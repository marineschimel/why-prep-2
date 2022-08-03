open Owl
module AD = Algodiff.D
module M = Arm.Make (Arm.Defaults)
open Lib
open Base
open Accessor.O

(* Setting up the parameters/directories
*)
let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")
let data_dir = Cmdargs.(get_string "-data" |> force ~usage:"-data [dir in which data is]")

let results_dir =
  Cmdargs.(get_string "-rdir" |> force ~usage:"-rdir [where to save the key results]")

let save_all = Cmdargs.(check "-save_all")

let lambda =
  Cmdargs.(get_float "-lambda" |> force ~usage:"-lambda [dir in which data is]")

let perturb = Cmdargs.(check "-perturb")
let annealing = Cmdargs.(check "-annealing")
let rad_c = Cmdargs.(get_float "-rad_c" |> default 0.05)
let rad_w = Cmdargs.(get_float "-rad_w" |> default 0.5)
let amp = Cmdargs.(get_float "-amp" |> default 1.)
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

let _ =
  C.root_perform (fun () ->
      Mat.save_txt ~out:(in_dir "targets") (Mat.concatenate ~axis:0 targets))

let beta = AD.F 1E-2
let phi_x x = AD.Maths.(F amp * tanh x)
let d_phi_x x = AD.Maths.(F amp * (F 1. - (F 2. * sqr (tanh x))))
let d2_phi_x x = AD.Maths.(F amp * F (-2.) * tanh x * d_phi_x x)
let link_f x = phi_x x
let target i = targets.(i)
let dt = 2E-3
let lambda_prep = lambda
let lambda_mov = lambda
let n_out = 2
let _n = 204
let m = 200
let tau = 150E-3
let n_output = 2

let _ =
  Mat.save_txt
    ~out:(in_dir "prms")
    (Mat.of_array
       [| tau; lambda_prep; lambda_mov; dt; AD.unpack_flt beta; rad_c; rad_w; amp |]
       1
       (-1))

let theta0 = Mat.of_arrays [| [| 0.174533; 2.50532; 0.; 0. |] |] |> AD.pack_arr
let t_preps = [| 0.; 0.05; 0.1; 0.15; 0.2; 0.3; 0.45; 0.5; 0.6 |]

let eigenvalues m =
  let v = Linalg.D.eigvals m in
  let re = Dense.Matrix.Z.re v in
  let im = Dense.Matrix.Z.im v in
  Mat.(concat_horizontal (transpose re) (transpose im))

let w = C.broadcast' (fun () -> Mat.gaussian ~sigma:Float.(rad_w /. sqrt (of_int m)) m m)
let eigs = C.broadcast' (fun () -> eigenvalues w)

let _ =
  C.root_perform (fun () ->
      Mat.save_txt ~out:(Printf.sprintf "%s/w" dir) w;
      Mat.save_txt ~out:(Printf.sprintf "%s/eigs" dir) eigs)

let c =
  C.broadcast' (fun () -> AD.Mat.gaussian ~sigma:Float.(rad_c / sqrt (of_int m)) 2 m)

let x0 = C.broadcast' (fun () -> AD.Mat.uniform ~a:(-2.) ~b:2. m 1)
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

let b = Mat.eye m

(*structure : get inputs from task, then from the same SOC just generating inputs, 
then from yet another SOC
Maybe do the same thing using 200 - 100 - 50 
what about one population receiving modulated inputs about the target? 
*)

let baseline_input =
  C.broadcast' (fun () ->
      AD.Maths.(neg ((AD.pack_arr w *@ link_f x0) - x0)) |> AD.Maths.transpose)

let x0 =
  let x0 =
    if perturb then AD.Maths.(x0 + AD.Arr.gaussian ~sigma:0.01 (AD.shape x0)) else x0
  in
  AD.Maths.concatenate [| AD.Maths.transpose theta0; x0 |] ~axis:0 |> AD.Maths.transpose

let tasks =
  Array.init
    (n_targets * Array.length t_preps)
    ~f:(fun i ->
      let n_time = i / n_targets in
      let n_target = Int.rem i n_targets in
      Model.
        { t_prep = t_preps.(n_time)
        ; x0
        ; t_movs = [| 0.3 |]
        ; dt
        ; t_hold = Some 0.2
        ; t_pauses = None
        ; scale_lambda = None
        ; target = AD.pack_arr (target n_target)
        ; theta0
        ; tau = 150E-3
        })

let save_prms suffix prms = Misc.save_bin (Printf.sprintf "%s/prms_%s" dir suffix) prms
let save_task suffix task = Misc.save_bin (Printf.sprintf "%s/prms_%s" dir suffix) task
let epsilon = 1E-1

module U = Priors.Gaussian

module D0 = Dynamics.Arm_Gated_Plus (struct
  let phi_x x = phi_x x
  let d_phi_x x = d_phi_x x
  let phi_u x = x
  let d_phi_u _ = AD.Mat.eye m
end)

(*  *)
module L0 = Likelihoods.End_Phi (struct
  let label = "output"
  let phi_x x = phi_x x
  let d_phi_x x = d_phi_x x
  let d2_phi_x x = d2_phi_x x
  let speed_end_penalty = 0.1
end)

module R = Readout

let prms =
  let open Owl_parameters in
  C.broadcast' (fun () ->
      let likelihood =
        Likelihoods.End_Phi_P.
          { c = (pinned : setter) c
          ; c_mask = None
          ; qs_coeff = (pinned : setter) (AD.F 1.)
          ; t_coeff = (pinned : setter) (AD.F 0.5)
          ; g_coeff = (pinned : setter) (AD.F 5.)
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
      Model.Full_P.{ generative; readout })

module I = Model.ILQR (U) (D0) (L0)

let summary_tasks =
  Array.init (Array.length t_preps) ~f:(fun _ ->
      Array.init n_targets ~f:(fun _ -> Mat.zeros 1 1, false))

let get_idx t =
  let _ = Stdio.printf "ts are %f %f %!" t t_preps.(0) in
  let idx, _ = Array.findi_exn t_preps ~f:(fun _ tp -> Float.(t = tp)) in
  idx

let save_results suffix xs us n_target n_prep task =
  let file s = Printf.sprintf "%s/%s_%s" dir s suffix in
  let _ = Stdio.printf "nprep is %i %!" n_prep in
  let xs = AD.unpack_arr xs in
  let us = AD.unpack_arr us in
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
  let summary =
    ( Mat.of_array
        [| t_prep; prep_idx; loss; input_cost_tot; torque_err; target_err; ratio_u_cost |]
        1
        (-1)
    , true )
  in
  if save_all
  then (
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
    Owl.Mat.save_txt ~out:(file "thetas") thetas;
    Owl.Mat.save_txt ~out:(file "xs") xs;
    Owl.Mat.save_txt ~out:(file "us") us;
    Owl.Mat.save_txt ~out:(file "rates") rates;
    Owl.Mat.save_txt ~out:(file "eff_us") (AD.unpack_arr (link_f (AD.pack_arr us)));
    Owl.Mat.save_txt
      ~out:(file "torques")
      Mat.(
        (rates - AD.unpack_arr (link_f (AD.pack_arr x0))) *@ transpose (AD.unpack_arr c));
    get_idx t_prep, n_target, summary)
  else get_idx t_prep, n_target, summary

let init_us =
  Array.init (Array.length tasks) ~f:(fun _ -> Mat.(gaussian ~sigma:0. 2001 m))

let prep_initial_conditions =
  if annealing
  then
    let open Model in
    let open Full_P.A in
    let open Generative_P.A in
    let open Owl_parameters in
    let open Dynamics_typ in
    let init_w = Owl_parameters.extract prms.generative.dynamics.a in
    let rescaled_prms =
      Accessor.map
        (generative @> dynamics @> Dynamics.Arm_Plus_P.A.a)
        ~f:(fun _ -> (pinned : setter) AD.Maths.(init_w / F 5.))
        prms
    in
    Array.foldi tasks ~init:[] ~f:(fun i accu t ->
        if Int.(i % C.n_nodes = C.rank)
        then (
          let n_target = Int.rem i n_targets in
          let n_prep = Float.to_int (t.t_prep /. dt) in
          let t_prep_int = Float.to_int (1000. *. t.t_prep) in
          let _ = Stdio.printf "trying weak rad" in
          let xs, us, l, success =
            I.solve
              ~u_init:Mat.(gaussian ~sigma:0. 2001 m)
              ~n:(m + 4)
              ~m
              ~x0
              ~prms:rescaled_prms
              t
          in
          if success then (i, AD.unpack_arr us) :: accu else accu)
        else accu)
    |> C.gather
    |> fun v ->
    C.root_perform (fun () ->
        let us = v |> Array.to_list |> List.concat |> Array.of_list in
        Array.iter us ~f:(fun (i, us) -> init_us.(i) <- us))

let full_tasks =
  C.broadcast' (fun () ->
      Array.init (Array.length tasks) ~f:(fun i -> tasks.(i), init_us.(i)))

let _ = Stdio.printf "ran the u initialization %!"

let summaries, try_tasks =
  let x0 = x0 in
  let _ = save_prms "" prms in
  Array.foldi full_tasks ~init:([], []) ~f:(fun i accu t ->
      let accu1, accu2 = accu in
      if Int.(i % C.n_nodes = C.rank)
      then (
        let t, u0 = t in
        let n_target = Int.rem i n_targets in
        let n_prep = Float.to_int (t.t_prep /. dt) in
        let t_prep_int = Float.to_int (1000. *. t.t_prep) in
        let xs, us, l, success = I.solve ~u_init:u0 ~n:(m + 4) ~m ~x0 ~prms t in
        let idx, n_target, summary =
          save_results
            (Printf.sprintf "%i_%i" n_target t_prep_int)
            xs
            us
            n_target
            n_prep
            t
        in
        if save_all
        then (
          Mat.save_txt
            ~out:(in_dir (Printf.sprintf "loss_%i_%i" n_target t_prep_int))
            (Mat.of_array [| AD.unpack_flt l |] 1 (-1));
          save_task (Printf.sprintf "%i_%i" n_target t_prep_int) t);
        if success then (idx, n_target, summary) :: accu1, accu2 else accu1, t :: accu2)
      else accu)
  |> C.gather
  |> fun v ->
  C.broadcast' (fun () ->
      let v1 = Array.map ~f:fst v in
      let v2 = Array.map ~f:snd v in
      ( v1 |> Array.to_list |> List.concat |> Array.of_list
      , v2 |> Array.to_list |> List.concat |> Array.of_list ))

let save_summaries =
  C.root_perform (fun () ->
      Array.iter summaries ~f:(fun (i, n, s) -> summary_tasks.(i).(n) <- s))

let _ = Stdio.printf "saved the summaries so far %!"

let final_save =
  C.root_perform (fun () ->
      let eigs = eigenvalues w in
      let norm_eigs = Mat.l2norm ~axis:1 eigs in
      let max_eig = Mat.max' eigs in
      let prefix =
        if perturb
        then Printf.sprintf "%s/rad_%.3f_amp_%.3f_seed_%i_y" results_dir rad_w amp seed
        else Printf.sprintf "%s/rad_%.3f_amp_%.3f_seed_%i_n" results_dir rad_w amp seed
      in
      let _ =
        Misc.save_bin (Printf.sprintf "%s_summaries" prefix) summary_tasks;
        Mat.save_txt ~out:(Printf.sprintf "%s_eigs" prefix) eigs
      in
      let mean_across_tgts =
        Array.map summary_tasks ~f:(fun x ->
            let arr_mat = Array.map ~f:fst x in
            let mat = Mat.concatenate ~axis:0 arr_mat in
            let m = Mat.mean ~axis:0 mat in
            let s = Mat.var ~axis:0 mat in
            Mat.(m @|| s))
        |> Mat.concatenate ~axis:0
      in
      let n_preps = Array.length t_preps in
      let full_summary =
        let vals_0 = summary_tasks.(0) |> Array.map ~f:fst |> Mat.concatenate ~axis:0 in
        let val_f =
          summary_tasks.(n_preps - 1) |> Array.map ~f:fst |> Mat.concatenate ~axis:0
        in
        let loss_ratio =
          Mat.(Mat.get_slice [ []; [ 2 ] ] val_f / Mat.get_slice [ []; [ 2 ] ] vals_0)
          |> Mat.mean'
        in
        let prep_idx_f = Mat.get_slice [ []; [ 1 ] ] val_f |> Mat.mean' in
        let norm_w = Mat.l2norm_sqr' w in
        Mat.of_array
          [| Mat.max' norm_eigs; max_eig; norm_w; loss_ratio; prep_idx_f |]
          1
          (-1)
      in
      Mat.save_txt ~out:(Printf.sprintf "%s_mean_across" prefix) mean_across_tgts;
      Mat.save_txt ~out:(Printf.sprintf "%s_full_summary" prefix) full_summary)

(* let full_summary =   *)

let _ = C.barrier ()
