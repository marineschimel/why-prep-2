open Owl
module M = Arm.Make (Arm.Defaults)
open Lib
open Base
include Model_typ

let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")
let data_dir = Cmdargs.(get_string "-data" |> force ~usage:"-data [dir in which data is]")

let lambda =
  Cmdargs.(get_float "-lambda" |> force ~usage:"-lambda [dir in which data is]")

let scale_mov = Cmdargs.(get_float "-scale_mov" |> default 1.)
let scale_prep = Cmdargs.(get_float "-scale_prep" |> default 1.)

let results_dir =
  Cmdargs.(get_string "-rdir" |> force ~usage:"-rdir [where to save the key results]")

let nonlin = if Cmdargs.(check "-tanh") then "tanh" else "relu"
let save_all = Cmdargs.(check "-save_all")
let soc = Cmdargs.(check "-soc")
let skew = Cmdargs.(check "-skew")
let triang = Cmdargs.(check "-triang")
let lr = Cmdargs.(check "-lr")
let rad_c = Cmdargs.(get_float "-rad_c" |> default 0.5)
let rad_w = Cmdargs.(get_float "-rad_w" |> default 0.5)
let tau_mov = Cmdargs.(get_float "-tau_mov" |> force ~usage:"tau_mov")
let t_coeff = Cmdargs.(get_float "-t_coeff" |> default 1.)
let exponent = Cmdargs.(get_int "-exponent" |> force ~usage:"exponent")
let seed = Cmdargs.(get_int "-seed" |> default 1)
let in_dir s = Printf.sprintf "%s/%s" dir s
let in_data_dir s = Printf.sprintf "%s/%s" data_dir s
let n_targets = 1

let targets =
  C.broadcast' (fun () ->
      Array.init n_targets ~f:(fun i ->
          let radius = 0.12 in
          let theta = Float.(of_int i *. Const.pi *. 2. /. of_int n_targets) in
          let x = Maths.(radius *. cos theta) in
          let y = Maths.(radius *. sin theta) in
          (* let x = Stats.uniform_rvs ~a:(Maths.neg radius) ~b:radius in
      let y = Stats.uniform_rvs ~a:(Maths.neg radius) ~b:radius in *)
          let t1, t2 = Misc.pos_to_angles x (y +. 0.199112) in
          Mat.of_array [| t1; t2; 0.; 0. |] 1 (-1)))

let hand_targets =
  C.broadcast' (fun () ->
      Array.init n_targets ~f:(fun i ->
          let radius = 0.12 in
          let theta = Float.(of_int i *. Const.pi *. 2. /. of_int n_targets) in
          let x = Maths.(radius *. cos theta) in
          let y = Maths.(radius *. sin theta) in
          Mat.of_array [| x; y |] 1 (-1)))

let _ =
  C.root_perform (fun () ->
      Mat.save_txt ~out:(in_dir "targets") (Mat.concatenate ~axis:0 targets);
      Mat.save_txt ~out:(in_dir "hand_targets") (Mat.concatenate ~axis:0 hand_targets))

let beta = AD.F 1E-2
let exponent = AD.F (Float.of_int exponent)
let phi_t t = AD.Maths.(t ** exponent)
let phi_x x = x
let d_phi_x x = AD.Maths.(F 0.5 * (F 1. + (F 0. * x)))
let d2_phi_x x = AD.Maths.(diagm (F 0. * x))
let link_f x = phi_x x
let target i = targets.(i)
let dt = 2E-3
let lambda_prep = lambda *. scale_prep
let lambda_mov = lambda *. scale_mov
let n_out = 2
let _n = 204
let m = 200
let tau = 150E-3
let n_output = 2

let _ =
  Mat.save_txt
    ~out:(in_dir "prms")
    (Mat.of_array
       [| tau; lambda_prep; lambda_mov; dt; AD.unpack_flt beta; rad_c |]
       1
       (-1))

let theta0 = Mat.of_arrays [| [| 0.174533; 2.50532; 0.; 0. |] |] |> AD.pack_arr
let t_preps = [| 0.3 |]
let w_soc seed = Mat.(load_txt (Printf.sprintf "%s/w_rec_%i" data_dir seed))
let w_rdn seed = Mat.gaussian ~sigma:Float.(0.9 /. sqrt (of_int m)) m m

let w_skew seed =
  let m = Mat.gaussian ~sigma:Float.(5. /. sqrt (of_int m)) m m in
  Mat.((m - transpose m) /$ 2.)

let w_nonnorm seed =
  let m = Mat.gaussian ~sigma:Float.(rad_w /. sqrt (of_int m)) m m in
  Mat.triu ~k:1 m

let c seed = AD.Mat.gaussian ~sigma:Float.(rad_c / sqrt (of_int m)) 2 m
let seeds = C.broadcast' (fun () -> [| 0 |])
let nets = C.broadcast' (fun () -> [| "soc"; "skew"; "nonnorm"; "random" |])
let lambdas = C.broadcast' (fun () -> [| 1., 1.; 1., 1000.; 1000., 1. |])
let x0 = C.broadcast' (fun () -> AD.Maths.(F 0.5 * AD.Mat.uniform ~a:5. ~b:15. m 1))

let eigenvalues m =
  let v = Linalg.D.eigvals m in
  let re = Dense.Matrix.Z.re v in
  let im = Dense.Matrix.Z.im v in
  Mat.(concat_horizontal (transpose re) (transpose im))

let b = Mat.eye m

let x0 =
  AD.Maths.concatenate [| AD.Maths.transpose theta0; x0 |] ~axis:0 |> AD.Maths.transpose

let _ = Stdio.printf "N targets : %i %!" n_targets

let tasks =
  Array.init n_targets ~f:(fun i ->
      let n_target = i in
      ( n_target
      , Model.
          { t_prep = t_preps.(0)
          ; x0
          ; t_movs = [| 0.4 |]
          ; dt
          ; t_hold = Some 0.2
          ; t_pauses = None
          ; scale_lambda = None
          ; target = AD.pack_arr (target n_target)
          ; theta0
          ; tau = 150E-3
          } ))

let task0 = tasks.(0)
let _ = Stdio.printf "Size of tasks : %i %!" (Array.length tasks)
let save_prms suffix prms = Misc.save_bin (Printf.sprintf "%s/prms_%s" dir suffix) prms
let save_task suffix task = Misc.save_bin (Printf.sprintf "%s/prms_%s" dir suffix) task
let epsilon = 1E-1

module U = Priors.Gaussian

module D0 = Dynamics.Arm_Plus (struct
  let phi_x x = phi_x x
  let d_phi_x x = d_phi_x x
  let phi_u x = x
  let d_phi_u _ = AD.Mat.eye m
end)

module L0 = Likelihoods.Ramping (struct
  let label = "output"
  let phi_x x = phi_x x
  let phi_t t = phi_t t
end)

module R = Readout

let parameters =
  let open Owl_parameters in
  C.broadcast' (fun () ->
      Array.init
        Int.(Array.length nets * Array.length lambdas * Array.length seeds)
        ~f:(fun i ->
          let n = Int.(Array.length nets * Array.length lambdas) in
          let n_lambdas = Array.length lambdas in
          let n_seed = i / n in
          let n_rest = Int.rem i n in
          let n_net = Int.(n_rest / Array.length lambdas) in
          let n_lambda = Int.rem n_rest n_lambdas in
          let c = c seed in
          let x0_base = AD.Maths.get_slice [ []; [ 4; -1 ] ] x0 |> AD.Maths.transpose in
          let u0 = phi_x x0_base in
          let norm_u0 = AD.Maths.(l2norm_sqr' u0) in
          let _ =
            Stdio.printf
              "%i %i %i %i%!"
              (AD.Mat.row_num c)
              (AD.Mat.col_num c)
              (AD.Mat.row_num u0)
              (AD.Mat.col_num u0)
          in
          let c = AD.Maths.(c - (c *@ u0 *@ transpose u0 / norm_u0)) in
          let net_name = nets.(n_net) in
          let w =
            if net_name == "soc"
            then w_soc seed
            else if net_name == "skew"
            then w_skew seed
            else if net_name == "nonnorm"
            then w_nonnorm seed
            else w_rdn seed
          in
          let baseline_input =
            AD.Maths.(neg ((AD.pack_arr w *@ link_f x0_base) - x0_base))
            |> AD.Maths.transpose
          in
          let scale_prep, scale_mov = lambdas.(n_lambda) in
          let lambda_prep = Float.(scale_prep *. lambda) in
          let lambda_mov = Float.(scale_mov *. lambda) in
          let likelihood =
            Likelihoods.Ramping_P.
              { c = (pinned : setter) c
              ; c_mask = None
              ; qs_coeff = (pinned : setter) (AD.F 1.)
              ; t_coeff = (pinned : setter) (AD.F t_coeff)
              ; g_coeff = (pinned : setter) (AD.F 1.)
              ; tau_mov = (pinned : setter) (AD.F Float.(0.001 *. tau_mov))
              }
          in
          let dynamics =
            Dynamics.Arm_Plus_P.
              { a = (pinned : setter) (AD.pack_arr Mat.(transpose w))
              ; b = (pinned : setter) (AD.pack_arr b)
              ; bias = (pinned : setter) baseline_input
              }
          in
          let prior =
            Priors.Gaussian_P.
              { lambda_prep = (pinned : setter) (AD.F lambda_prep)
              ; lambda_mov = (pinned : setter) (AD.F lambda_mov)
              }
          in
          let readout = R.Readout_P.{ c = (pinned : setter) c } in
          let generative = Model.Generative_P.{ prior; dynamics; likelihood } in
          Model.Full_P.{ generative; readout }, net_name, scale_prep, scale_mov, seed))

module I = Model.ILQR (U) (D0) (L0)

let sample_from_q ~mu ~q_cov_chol ~n_samples =
  (*cov = List of blocks of the covariance*)
  let _ =
    Stdio.printf
      "%i %i %i %i %!"
      (AD.shape mu).(0)
      (AD.shape mu).(1)
      (AD.Mat.row_num q_cov_chol.(0))
      (AD.Mat.row_num q_cov_chol.(1))
  in
  let m = (AD.Arr.shape mu).(1) in
  let mu = AD.Maths.reshape mu [| 1; (AD.Arr.shape mu).(0); m |] in
  let mu = AD.Maths.get_slice [ []; [ 0; -2 ] ] mu in
  let z =
    Array.map q_cov_chol ~f:(fun ell_t ->
        (* let z_t =
          AD.Linalg.linsolve ~trans:true ~typ:`l ell_t (AD.Mat.gaussian m n_samples)
        in *)
        let eta = AD.Mat.gaussian m n_samples in
        let z_t = AD.Maths.(ell_t *@ eta) in
        let z_t = AD.Maths.transpose z_t in
        AD.Maths.reshape z_t [| 1; n_samples; m |])
    |> AD.Maths.concatenate ~axis:0
  in
  let z = AD.Maths.transpose ~axis:[| 1; 0; 2 |] z in
  AD.Maths.(mu + z)

let q_cov_chol ?(epsilon = 0.) inv_q_covs =
  (*inv_q_covs is a list with each element being the hessian of that block of inputs*)
  inv_q_covs
  |> Array.of_list
  |> Array.map ~f:(fun x ->
         assert (Owl.Linalg.D.is_posdef (AD.unpack_arr x));
         let n = AD.Mat.row_num x in
         let x = AD.Maths.(x + (F 1E-5 * AD.Mat.eye n)) in
         let y = AD.Maths.(AD.Linalg.inv x + (F epsilon * AD.Mat.eye n)) in
         AD.Linalg.chol y)

(* inv_q_cov = U^T U = LL^T *)
(* AD.Linalg.chol x |> AD.Maths.transpose) *)

let log_lik_term ~prms ~readout ~task u =
  let module I = Dynamics.Integrate (D0) in
  let integrate = I.integrate ~readout:(Owl_parameters.extract prms.Full_P.readout.c) in
  let u = AD.Maths.reshape u [| (AD.shape u).(1); (AD.shape u).(2) |] in
  let xs = integrate ~prms:prms.generative.dynamics ~task ~n:204 ~u in
  let xs = AD.unpack_arr xs in
  let n_prep = Int.(of_float (task.t_prep /. dt)) in
  let le, lm =
    Analysis_funs.cost_x
      ~f:(fun k x ->
        let c =
          L0.neg_logp_t
            ~readout
            ~prms:prms.generative.likelihood
            ~task
            ~k
            ~z_t:(AD.pack_arr xs)
        in
        AD.unpack_flt c)
      ~n_prep
      xs
  in
  AD.pack_flt Float.(le +. lm)

let log_p ~prms ~task us =
  let _ =
    Stdio.printf "%i %i %i %!" (Arr.shape us).(1) (Arr.shape us).(2) (Arr.shape us).(0)
  in
  let n_prep = Int.(of_float (task.t_prep /. dt)) in
  let _, _, l =
    let f k u =
      Priors.Gaussian_Prior.neg_logp_t
        ~prms:prms.Full_P.generative.prior
        ~task
        ~k
        ~x:(AD.pack_arr u)
        ~u:(AD.pack_arr u)
    in
    Analysis_funs.cost_u
      ~f:(fun k x -> AD.unpack_flt (f k x))
      ~n_prep
      (Mat.reshape us [| (Arr.shape us).(1); (Arr.shape us).(2) |])
  in
  l |> AD.pack_flt

let logq us mu_u q_cov_chol =
  let us = AD.pack_arr us in
  let q_cov0 = q_cov_chol.(0) in
  let m_ = AD.Mat.col_num q_cov0 in
  let t_ = Array.length q_cov_chol in
  let m = Float.of_int m_ in
  let t = Float.of_int t_ in
  let cst = Float.(m * t * log Const.pi2) in
  let log_det_term =
    Array.fold q_cov_chol ~init:(AD.F 0.) ~f:(fun accu x ->
        AD.Maths.(accu + (F 2. * sum' (log (diag x)))))
  in
  let u_s = AD.shape us in
  assert (Array.length u_s = 3);
  (*us is of size 1 x M x T*)
  let n_samples = u_s.(0) in
  let mu_u = AD.Maths.get_slice [ [ 0; -2 ]; [] ] mu_u in
  let du =
    AD.Maths.(us - AD.Arr.reshape mu_u [| 1; (AD.shape mu_u).(0); (AD.shape mu_u).(1) |])
  in
  (* quadratic term: 
            assuming vec is stacking columns, du = vec(dU) and dU = is T x N
              du^t ((S^t S)⊗(T^t T))^{-1} du
          =  du^t ((S^{-1} S^{-t})⊗(T^{-1} T^{-t})) du 
          =  du^t (S^{-1}⊗T^{-1}) (S^{-t}⊗T^{-t}) du
          =  || (S^{-t}⊗T^{-t}) du ||^2 
          =  || vec(T^{-t} dU S^{-1}) ||^2 *)
  let quadratic_term =
    Array.foldi q_cov_chol ~init:(AD.F 0.) ~f:(fun t accu ell_t ->
        let du =
          AD.Maths.get_slice [ []; [ t ] ] du |> fun z -> AD.Maths.reshape z [| -1; m_ |]
          (*now it's n_samples x 1 x m*)
        in
        (*ell_t is upper trig*)
        let du_ell = AD.Linalg.linsolve ~typ:`u ell_t (AD.Maths.transpose du) in
        AD.Maths.(accu + l2norm_sqr' du_ell))
  in
  AD.Maths.(
    F (-0.5) * ((F Float.(of_int n_samples) * (F cst + log_det_term)) + quadratic_term))

let eval_ml ?(epsilon = 0.) ~mu_u ~q_cov_chol ~n_samples ~task ~prms =
  let u_samples = sample_from_q ~mu:mu_u ~q_cov_chol ~n_samples |> AD.unpack_arr in
  let u_samples = Arr.split ~axis:0 (Array.init n_samples ~f:(fun _ -> 1)) u_samples in
  let v =
    Array.foldi u_samples ~init:[] ~f:(fun i acc u ->
        if Int.(i % C.n_nodes = C.rank)
        then (
          let llik =
            log_lik_term
              ~prms
              ~readout:(Owl_parameters.extract prms.Full_P.readout.c)
              ~task
              (AD.pack_arr u)
          in
          let lp = log_p ~prms ~task u in
          let lq = logq u mu_u q_cov_chol in
          AD.Maths.(llik + lp - lq) :: acc)
        else acc)
    |> C.gather
  in
  C.broadcast' (fun () ->
      let v =
        v
        |> Array.to_list
        |> List.concat
        |> Array.of_list
        |> fun z -> AD.Maths.of_arrays [| z |]
      in
      AD.Maths.(log_sum_exp' v - log (AD.F (Float.of_int n_samples))))

let get_mll xs us quus n_target n_prep task prms =
  let q_cov_chol = q_cov_chol quus in
  let mll = eval_ml ~epsilon:0. ~mu_u:us ~q_cov_chol ~n_samples:1000 ~task ~prms in
  AD.unpack_flt mll

let approx_posteriors =
  let x0 = x0 in
  let v =
    Array.foldi
      parameters
      ~init:[]
      ~f:(fun i acc (prms, net_name, scale_prep, scale_mov, seed) ->
        if Int.(i % C.n_nodes = C.rank)
        then (
          let _, t = task0 in
          let n_prep = Float.to_int (t.t_prep /. dt) in
          let t_prep_int = Float.to_int (1000. *. t.t_prep) in
          let xs, us, l, quus, _ =
            I.solve ~u_init:Mat.(gaussian ~sigma:0. 2001 m) ~n:(m + 4) ~m ~x0 ~prms t
          in
          (prms, net_name, scale_prep, scale_mov, seed, xs, us, quus) :: acc)
        else acc)
    |> C.gather
  in
  C.broadcast' (fun () ->
      let _ = Stdio.printf "Broadcasting approx posteriors %!" in
      v |> Array.to_list |> List.concat)

(* should we broadcast here? *)

let _ = Stdio.printf "posteriors broadcasted %!"

let get_mll =
  let rec compute_mll ap =
    match ap with
    | hd :: tl ->
      let prms, net_name, scale_prep, scale_mov, seed, xs, us, quus = hd in
      let n_target, t = task0 in
      let n_prep = Float.to_int (t.t_prep /. dt) in
      let t_prep_int = Float.to_int (1000. *. t.t_prep) in
      let mll = get_mll xs us quus n_target n_prep t prms in
      C.root_perform (fun () ->
          Mat.save_txt
            ~out:
              (in_dir
                 (Printf.sprintf "mll_%s_%.1f_%.1f_%i" net_name scale_prep scale_mov seed))
            (Mat.of_arrays [| [| mll |] |]);
          save_prms
            (Printf.sprintf "%s_%.1f_%.1f_%i" net_name scale_prep scale_mov seed)
            prms;
          save_task (Printf.sprintf "%i_%i" n_target Int.(2 * n_prep)) t);
      compute_mll tl
    | [] -> Stdio.printf "done %!"
  in
  compute_mll approx_posteriors
