open Owl
module AD = Algodiff.D
module M = Arm.Make (Arm.Defaults)
open Lib
open Base

(* Setting up the parameters/directories
*)
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
let sim_soc = Cmdargs.(check "-sim_soc")
let skew = Cmdargs.(check "-skew")
let triang = Cmdargs.(check "-triang")
let lr = Cmdargs.(check "-lr")
let rad_c = Cmdargs.(get_float "-rad_c" |> default 0.5)
let rad_w = Cmdargs.(get_float "-rad_w" |> default 0.5)
let tau_mov = Cmdargs.(get_float "-tau_mov" |> default 600.)
let t_coeff = Cmdargs.(get_float "-t_coeff" |> default 1.)
let exponent = Cmdargs.(get_int "-exponent" |> force ~usage:"exponent")
let n_lr = Cmdargs.(get_int "-n_lr" |> default 1)
let sa_tgt = Cmdargs.(get_float "-sa_tgt" |> default 1.)
let seed = Cmdargs.(get_int "-seed" |> default 1)
let in_dir s = Printf.sprintf "%s/%s" dir s
let in_data_dir s = Printf.sprintf "%s/%s" data_dir s
let n_targets = 8

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

(* let phi_x x = x
let d_phi_x x = AD.Maths.(F 1. + (F 0. * x))
let d2_phi_x x = AD.Maths.(diagm (F 0. * x)) *)
let _ = if nonlin == "tanh" then Stdio.printf "tanh nonlinearity"
let phi_x x = if nonlin == "tanh" then AD.Maths.(F 5. * tanh x) else AD.Maths.relu x

let d_phi_x x =
  if nonlin == "tanh"
  then AD.Maths.(F 5. * (F 1. - (F 2. * sqr (tanh x))))
  else AD.Maths.(F 0.5 * (F 1. + signum x))

let d2_phi_x x =
  if nonlin == "tanh"
  then AD.Maths.(F (-10.) * tanh x * d_phi_x x)
  else AD.Maths.(diagm (F 0. * x))

(* let phi_x x = AD.Maths.(F 5. * tanh x)
let d_phi_x x = AD.Maths.(F 5. * (F 1. - F 2. * (sqr (tanh x))))
let d2_phi_x x = AD.Maths.(F (-10.) * tanh x * (d_phi_x x)) *)
(* let phi_x x = AD.Maths.(beta * AD.requad (x / beta))
let d_phi_x x = AD.Maths.(AD.d_requad (x / beta))
let d2_phi_x x = AD.Maths.(F 1. / beta * AD.d2_requad (x / beta)) *)
let link_f x = phi_x x
let target i = targets.(i)
let dt = 1E-3
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
let t_preps = [| 0.; 0.5 |]
(* 0.025; 0.05; 0.07; 0.1; 0.15; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9; 1.0 |] *)

let u, v =
  let m = Mat.gaussian ~sigma:Float.(rad_w /. sqrt (of_int m)) m m in
  let q, r, _ = Linalg.D.qr m in
  ( Mat.get_slice [ []; [ 0; n_lr - 1 ] ] m
  , Mat.get_slice [ []; [ n_lr; (2 * n_lr) - 1 ] ] m )

let _ =
  Stdio.printf
    "%i %i %i %i %!"
    (Mat.row_num u)
    (Mat.col_num u)
    (Mat.row_num v)
    (Mat.col_num v)

let eigenvalues m =
  let v = Linalg.D.eigvals m in
  let re = Dense.Matrix.Z.re v in
  let im = Dense.Matrix.Z.im v in
  Mat.(concat_horizontal (transpose re) (transpose im))

let get_sa m =
  let v = Linalg.D.eigvals m in
  let v = Dense.Matrix.Z.re v in
  Mat.max' v

let w =
  C.broadcast' (fun () ->
      if soc
      then Mat.(load_txt (Printf.sprintf "%s/w_rec_%i" data_dir seed))
      else if sim_soc
      then (
        let w = Mat.(load_txt (Printf.sprintf "%s/w_rec_%i" data_dir seed)) in
        let transformed_w, s = Misc.transform w in
        let _ =
          Mat.save_txt ~out:(Printf.sprintf "%s/sim_norm" data_dir) (Mat.l2norm ~axis:1 s)
        in
        transformed_w)
      else if lr
      then (
        let m = Mat.(u *@ transpose v) in
        let max_eig = get_sa m in
        if Float.(sa_tgt < 1.) then Mat.(m /$ Float.(max_eig *. sa_tgt)) else m)
      else (
        let m = Mat.gaussian ~sigma:Float.(rad_w /. sqrt (of_int m)) m m in
        if skew
        then Mat.((m - transpose m) /$ 2.)
        else if triang
        then Mat.triu ~k:1 m
        else m))

(* let _ =
  C.root_perform (fun () ->
      Mat.save_txt ~out:(Printf.sprintf "%s/w" dir) (Mat.gaussian ~sigma:0.05 m m))


let w = C.broadcast' (fun () -> Mat.(load_txt (Printf.sprintf "%s/w" dir))) *)

let c =
  C.broadcast' (fun () -> AD.Mat.gaussian ~sigma:Float.(rad_c / sqrt (of_int m)) 2 m)

let x0 = C.broadcast' (fun () -> AD.Maths.(F 0.5 * AD.Mat.uniform ~a:5. ~b:15. m 1))

let _ =
  C.root_perform (fun () ->
      Mat.save_txt ~out:(in_dir "eigs") (eigenvalues w);
      Mat.save_txt
        ~out:(in_dir "eigs_wrec")
        (eigenvalues Mat.(load_txt (Printf.sprintf "%s/w_rec_%i" data_dir seed)));
      Mat.save_txt ~out:(in_dir "w") w)
(* let c = C.broadcast' (fun () -> AD.pack_arr Mat.(load_txt (Printf.sprintf "%s/c" dir))) *)

(* let x0 =
  C.broadcast' (fun () ->
      AD.Maths.transpose (AD.pack_arr Mat.(load_txt (Printf.sprintf "%s/x0" dir)))) *)

let u0 = phi_x x0
let norm_u0 = AD.Maths.(l2norm_sqr' u0)
let c = AD.Maths.(c - (c *@ u0 *@ transpose u0 / norm_u0))

let _ =
  C.root_perform (fun () ->
      Mat.save_txt
        ~out:(Printf.sprintf "%s/nullspace" dir)
        (AD.unpack_arr AD.Maths.(c *@ u0)))

let _ =
  C.root_perform (fun () ->
      Mat.save_txt ~out:(Printf.sprintf "%s/c" dir) (AD.unpack_arr c))

(* 
let _ =
  C.root_perform (fun () ->
      Mat.save_txt ~out:(Printf.sprintf "%s/x0" dir) (AD.unpack_arr x0)) *)

let b = Mat.eye m

(*structure : get inputs from task, then from the same SOC just generating inputs, 
then from yet another SOC
Maybe do the same thing using 200 - 100 - 50 
what about one population receiving modulated inputs about the target? 
*)
(* x(t+1)- x(t) = Wx(t) + baseline => 0 when x = -W^(-1)*baseline  *)
(* let baseline_input = AD.Mat.ones 1 m *)

(* let x0 =
  let winv = AD.Linalg.linsolve (AD.pack_arr w) (AD.Mat.eye m) in
  AD.Maths.(neg (winv *@ transpose baseline_input)) |> AD.Maths.transpose *)

(* let x0 =
  C.broadcast' (fun () ->
      AD.Maths.transpose (AD.pack_arr (Mat.load_txt (Printf.sprintf "%s/x0" dir)))) *)

let baseline_input =
  C.broadcast' (fun () ->
      AD.Maths.(neg ((AD.pack_arr w *@ link_f x0) - x0)) |> AD.Maths.transpose)

let x0 =
  AD.Maths.concatenate [| AD.Maths.transpose theta0; x0 |] ~axis:0 |> AD.Maths.transpose

let _ = Stdio.printf "N targets : %i %!" n_targets

let tasks =
  Array.init
    (n_targets * Array.length t_preps)
    ~f:(fun i ->
      let n_time = i / n_targets in
      let n_target = Int.rem i n_targets in
      ( n_target
      , Model.
          { t_prep = t_preps.(n_time)
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

let _ = Stdio.printf "Size of tasks : %i %!" (Array.length tasks)
let save_prms suffix prms = Misc.save_bin (Printf.sprintf "%s/prms_%s" dir suffix) prms
let save_task suffix task = Misc.save_bin (Printf.sprintf "%s/prms_%s" dir suffix) task
let epsilon = 1E-1

(* let phi_x x = let x = AD.unpack_arr x in let y = Mat.map (fun x -> if Float.(x > epsilon) then x else if Float.(x > neg epsilon) 
    then Float.((x +. epsilon)**2. /. ((4. *. epsilon))) 
else 0.) x in AD.pack_arr y 

let d_phi_x x = let x = AD.unpack_arr x in let y = Mat.map (fun x -> if Float.(x > epsilon) then 1. 
else if Float.(x > neg epsilon)  then ((x +. epsilon)/.(2. *. epsilon)) else 0.) x in AD.pack_arr y

let d2_phi_x x = let x = AD.unpack_arr x in let y = Mat.map (fun x -> if Float.(x > epsilon) then 0. 
else if Float.(x > neg epsilon) then Float.(1./(2. * epsilon)) else 0.) x in AD.pack_arr y *)

module U = Priors.Gaussian
(* 
_Phi (struct
let phi_u x = phi_x x
  let d_phi_u x =  d_phi_x x
  let d2_phi_u x = d2_phi_x x
end) *)

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

let t_tot = 0.6

module R = Readout

let dt_scaling = Float.(dt /. 1E-3 *. 1000.)

let prms =
  let open Owl_parameters in
  C.broadcast' (fun () ->
      let likelihood =
        Likelihoods.Ramping_P.
          { c = (pinned : setter) c
          ; c_mask = None
          ; qs_coeff = (pinned : setter) (AD.F Float.(dt_scaling))
          ; t_coeff = (pinned : setter) (AD.F Float.(t_coeff *. dt_scaling))
          ; g_coeff = (pinned : setter) (AD.F Float.(dt_scaling))
          ; tau_mov = (pinned : setter) (AD.F Float.(t_tot))
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
          { lambda_prep = (pinned : setter) (AD.F Float.(lambda_prep *. dt_scaling))
          ; lambda_mov = (pinned : setter) (AD.F Float.(lambda_mov *. dt_scaling))
          }
      in
      let readout = R.Readout_P.{ c = (pinned : setter) c } in
      let generative = Model.Generative_P.{ prior; dynamics; likelihood } in
      Model.Full_P.{ generative; readout })

module I = Model.ILQR (U) (D0) (L0)

let summary_tasks =
  Array.init (Array.length t_preps) ~f:(fun _ ->
      Array.init n_targets ~f:(fun _ -> Mat.zeros 1 1, false))

let get_idx t =
  let _ = Stdio.printf "ts are %f %f %!" t t_preps.(0) in
  let idx, _ = Array.findi_exn t_preps ~f:(fun _ tp -> Float.(t = tp)) in
  idx

let save_results suffix xs us quus n_target n_prep task =
  let file s = Printf.sprintf "%s/%s_%s" dir s suffix in
  let _ = Stdio.printf "nprep is %i %!" n_prep in
  let xs = AD.unpack_arr xs in
  let us = AD.unpack_arr us in
  let _ = Stdio.printf "%i %i %!" (Mat.row_num xs) (Mat.col_num xs) in
  let torque_err, target_err =
    Analysis_funs.cost_x
      ~f:(fun k x ->
        let c =
          L0.neg_logp_t
            ~prms:prms.generative.likelihood
            ~task
            ~k
            ~z_t:(AD.pack_arr x)
            ~readout:(Owl_parameters.extract prms.readout.c)
        in
        AD.unpack_flt c)
      ~n_prep
      xs
  in
  let t_to_target = Analysis_funs.time_to_end xs targets.(n_target) n_prep in
  let thetas, xs, us =
    Mat.get_slice [ []; [ 0; 3 ] ] xs, Mat.get_slice [ []; [ 4; -1 ] ] xs, us
  in
  let x0 = Mat.get_slice [ []; [ 4; -1 ] ] (AD.unpack_arr x0) in
  let rates = AD.unpack_arr (link_f (AD.pack_arr xs)) in
  let ue_prep, ue_mov, ue_tot =
    Analysis_funs.cost_u ~f:(fun _ x -> Mat.l2norm_sqr' x) ~n_prep us
  in
  let input_cost_prep, input_cost_mov, input_cost_tot =
    let f k u =
      U.neg_logp_t
        ~prms:prms.generative.prior
        ~task
        ~k
        ~x:(AD.pack_arr u)
        ~u:(AD.pack_arr u)
    in
    Analysis_funs.cost_u ~f:(fun k x -> AD.unpack_flt (f k x)) ~n_prep us
  in
  let t_prep = Float.of_int n_prep *. dt in
  let loss = torque_err +. target_err +. input_cost_tot in
  let prep_idx = ue_prep /. ue_mov in
  let ratio_u_cost = input_cost_prep /. input_cost_mov in
  let tr_quus =
    List.map ~f:(fun x -> AD.unpack_arr x |> Mat.diag) quus
    |> Array.of_list
    |> fun z -> Mat.concatenate z
  in
  let summary =
    ( Mat.of_array
        [| t_prep; prep_idx; loss; input_cost_tot; torque_err; target_err; ratio_u_cost |]
        1
        (-1)
    , true )
  in
  Owl.Mat.save_txt
    ~out:(file "quus")
    (Mat.of_array [| input_cost_prep; input_cost_mov; input_cost_tot |] 1 (-1));
  Owl.Mat.save_txt
    ~out:(file "summary")
    (Mat.of_array
       [| t_prep; prep_idx; loss; input_cost_tot; torque_err; target_err; ratio_u_cost |]
       1
       (-1));
  Owl.Mat.save_txt
    ~out:(file "u_cost")
    (Mat.of_array [| input_cost_prep; input_cost_mov; input_cost_tot |] 1 (-1));
  Owl.Mat.save_txt
    ~out:(file "task_cost")
    (Mat.of_array [| torque_err; target_err |] 1 (-1));
  Owl.Mat.save_txt
    ~out:(file "u_energy")
    (Mat.of_array [| ue_prep; ue_mov; ue_tot |] 1 (-1));
  Owl.Mat.save_txt
    ~out:(file "t_to_tgt")
    (Mat.of_array [| Float.of_int t_to_target |] 1 (-1));
  let hands =
    let h =
      AD.Mat.map_by_row
        (fun x -> Arm.unpack_state_diff (M.hand_of (Arm.pack_state x)))
        (AD.pack_arr thetas)
    in
    AD.unpack_arr h
  in
  Owl.Mat.save_txt ~out:(file "hands") hands;
  Owl.Mat.save_txt ~out:(file "thetas") thetas;
  Owl.Mat.save_txt ~out:(file "xs") xs;
  Owl.Mat.save_txt ~out:(file "us") us;
  Owl.Mat.save_txt ~out:(file "rates") rates;
  Owl.Mat.save_txt ~out:(file "eff_us") (AD.unpack_arr (link_f (AD.pack_arr us)));
  Owl.Mat.save_txt
    ~out:(file "torques")
    Mat.((rates - AD.unpack_arr (link_f (AD.pack_arr x0))) *@ transpose (AD.unpack_arr c))

let () =
  let x0 = x0 in
  let _ = save_prms "" prms in
  Array.iteri tasks ~f:(fun i (n_target, t) ->
      if Int.(i % C.n_nodes = C.rank)
      then (
        try
          let n_prep = Float.to_int (t.t_prep /. dt) in
          let t_prep_int = Float.to_int (1000. *. t.t_prep) in
          let xs, us, l, quus, _ =
            I.solve ~u_init:Mat.(gaussian ~sigma:0. 2001 m) ~n:(m + 4) ~m ~x0 ~prms t
          in
          let () =
            save_results
              (Printf.sprintf "%i_%i" n_target t_prep_int)
              xs
              us
              quus
              n_target
              n_prep
              t
          in
          if save_all
          then (
            let _ = Stdio.printf "success %i %i" n_target t_prep_int in
            Mat.save_txt
              ~out:(in_dir (Printf.sprintf "loss_%i_%i" n_target t_prep_int))
              (Mat.of_array [| AD.unpack_flt l |] 1 (-1));
            save_task (Printf.sprintf "%i_%i" n_target t_prep_int) t)
        with
        | _ -> Stdio.printf "fail %i" n_target))
