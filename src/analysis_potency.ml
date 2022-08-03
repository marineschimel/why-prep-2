open Owl
module AD = Algodiff.D
module M = Arm.Make (Arm.Defaults)
open Lib
open Base

(* Setting up the parameters/directories
*)
let dir =
  "/home/mmcs3/rds/rds-t2-cs156-T7o4pEA8QoU/mmcs3/linear_lambda_0.000001/ramping_soc/seed_0_mixed"

let lambda = 0.000001
let scale_mov = 1.
let scale_prep = 1.
let results_dir = dir
let nonlin = "relu"
let save_all = Cmdargs.(check "-save_all")
let soc = true
let lr = Cmdargs.(check "-lr")
let rad_c = 0.05
let tau_mov = 200.
let t_coeff = 1.5
let exponent = 2
let seed = 1
let in_dir s = Printf.sprintf "%s/%s" dir s
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
let phi_x x = x
let d_phi_x x = AD.Maths.((F 0. * x) + F 1.)
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

let t_preps =
  [| 0.; 0.025; 0.05; 0.07; 0.1; 0.15; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9; 1.0 |]

let w = Mat.load_txt (in_dir "w")


let c = Mat.load_txt (in_dir "w")


let x0 = C.broadcast' (fun () -> AD.pack_arr (Mat.(load_txt (in_dir "rates_0_300")) |> fun m -> Mat.row m 0 |> Mat.transpose))

let baseline_input =
  C.broadcast' (fun () ->
      AD.Maths.(neg ((AD.pack_arr w *@ link_f x0) - x0)) |> AD.Maths.transpose)
(* 
let x0 =
  AD.Maths.concatenate [| AD.Maths.transpose theta0; x0 |] ~axis:0 |> AD.Maths.transpose *)

let eigenvalues m =
  let v = Linalg.D.eigvals m in
  let re = Dense.Matrix.Z.re v in
  let im = Dense.Matrix.Z.im v in
  Mat.(concat_horizontal (transpose re) (transpose im))

let _ =
  C.root_perform (fun () ->
      Mat.save_txt ~out:(in_dir "eigs") (eigenvalues w);
      Mat.save_txt ~out:(in_dir "w") w)
(* let c = C.broadcast' (fun () -> AD.pack_arr Mat.(load_txt (Printf.sprintf "%s/c" dir))) *)

(* let x0 =
  C.broadcast' (fun () ->
      AD.Maths.transpose (AD.pack_arr Mat.(load_txt (Printf.sprintf "%s/x0" dir)))) *)

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


let x0 = AD.Maths.transpose x0
let _ = Stdio.printf "N targets : %i %!" n_targets

let tasks =
  Array.init n_targets ~f:(fun i ->
      let n_time = i / n_targets in
      let n_target = Int.rem i n_targets in
      Model.
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
        })

let _ = Stdio.printf "Size of tasks : %i %!" (Array.length tasks)
let save_prms suffix prms = Misc.save_bin (Printf.sprintf "%s/prms_%s" dir suffix) prms
let save_task suffix task = Misc.save_bin (Printf.sprintf "%s/prms_%s" dir suffix) task
let epsilon = 1E-1

module U = Priors.Gaussian

module D0 = Dynamics.Nonlinear (struct
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

let prms =
  let open Owl_parameters in
  C.broadcast' (fun () ->
      let likelihood =
        Likelihoods.Ramping_P.
          { c = (pinned : setter) (AD.pack_arr c)
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
          ; bias = (pinned : setter) (AD.Maths.transpose baseline_input)
          }
      in
      let prior =
        Priors.Gaussian_P.
          { lambda_prep = (pinned : setter) (AD.F lambda_prep)
          ; lambda_mov = (pinned : setter) (AD.F lambda_mov)
          }
      in
      let readout = R.Readout_P.{ c = (pinned : setter) (AD.pack_arr c) } in
      let generative = Model.Generative_P.{ prior; dynamics; likelihood } in
      Model.Full_P.{ generative; readout })

module I = Model.ILQR (U) (D0) (L0)

let n = 200
let in_dir s = Printf.sprintf "%s/%s" dir s

let eigenvalues m =
  let v = Linalg.D.eigvals m in
  let re = Dense.Matrix.Z.re v in
  let im = Dense.Matrix.Z.im v in
  Mat.(concat_horizontal (transpose re) (transpose im))

let w = Mat.load_txt (in_dir "w")
let eigs_w = eigenvalues w
let _ = Mat.save_txt ~out:(in_dir "eigs_w") eigs_w
let c = Mat.load_txt (in_dir "c")

(*normalize c*)
let c = Mat.(c /$ Owl.Linalg.D.norm c)
let a = Mat.(w - eye n)
let disc_a = Mat.(zeros n n - eye n)
let obs_gramian a c = Linalg.D.lyapunov Mat.(transpose a) Mat.(neg (transpose c) *@ c)
let ctr_gramian a b = Linalg.D.lyapunov a Mat.(neg b *@ transpose b)
let obs_id = obs_gramian a (Mat.eye 200)
let obs_disc = obs_gramian disc_a Mat.(neg (transpose c) *@ c)
let ctr_disc = ctr_gramian disc_a Mat.(eye n)
let tr_p_disc = Mat.trace ctr_disc
let eig_disc_o, eig_disc_q = eigenvalues obs_disc, eigenvalues ctr_disc

let _ =
  Mat.save_txt (in_dir "obs_eig_disc") eig_disc_o;
  Mat.save_txt (in_dir "ctr_eig_disc") eig_disc_q

let obs_c = obs_gramian a c
let ctr = ctr_gramian a (Mat.eye 200)
let ctr = Mat.(ctr /$ Float.(tr_p_disc /. 200.))
(* obs_gramian a c *)

let obs_c_modes, obs_c_eigs =
  let obs_c_modes, _, _ = Linalg.D.svd obs_c in
  obs_c_modes, eigenvalues obs_c

let obs_id_modes, obs_id_eigs =
  let obs_id_modes, _, _ = Linalg.D.svd obs_id in
  obs_id_modes, eigenvalues obs_id

let ctr_modes, ctr_eigs =
  let ctr_modes, _, _ = Linalg.D.svd ctr in
  Mat.(get_slice [ []; [] ] ctr_modes), eigenvalues ctr

let traj i = Mat.load_txt (in_dir (Printf.sprintf "rates_%i_300" i))
let n_reaches = 8

module Int = Dynamics.Integrate (D0)
let rates =
  Array.init n_targets ~f:(fun i ->
      let us =  Mat.load_txt (in_dir (Printf.sprintf "us_%i_300" i)) in 
        Int.integrate
          ~readout:(AD.pack_arr c)
          ~prms:prms.Model.Full_P.generative.dynamics
          ~task:tasks.(i)
          ~n
          ~u:(AD.pack_arr us) |> AD.unpack_arr)

let mean_rates = 
  Array.map rates ~f:(fun t -> 
Arr.reshape t [| 1; -1; 200 |])
  |> fun z ->
  Arr.concatenate ~axis:0 z
  |> fun e -> Arr.mean ~axis:0 e |> fun e -> Mat.reshape e [| -1; 200 |]


let ctrl_num_unshuffled =
  Array.map rates ~f:(fun r ->
    let t = Mat.(r - mean_rates) in 
      Mat.(t *@ ctr_modes) |> fun z -> Arr.reshape z [| 1; -1; 200 |])
  |> fun z ->
  Arr.concatenate ~axis:0 z
  |> Mat.sqr
  |> fun e ->
  Arr.mean ~axis:0 e
  |> Mat.sqrt
  |> fun e -> Mat.reshape e [| -1; 200 |] |> Mat.sum ~axis:0

let _ = Mat.save_txt (in_dir "unshuffled_num_all") ctrl_num_unshuffled

let shuffle us idces =
  Mat.get_fancy [ R [ 0; -1 ]; L idces ] us


let ctrl_num_shuffled idces =
  let new_rates =   Array.init n_targets ~f:(fun i ->
    let us =  Mat.load_txt (in_dir (Printf.sprintf "us_%i_300" i)) in 
    let us = shuffle us idces in
      Int.integrate
        ~readout:(AD.pack_arr c)
        ~prms:prms.Model.Full_P.generative.dynamics
        ~task:tasks.(i)
        ~n
        ~u:(AD.pack_arr us) |> AD.unpack_arr)
    in let mean_rates = 
      Array.map new_rates ~f:(fun t -> 
    Arr.reshape t [| 1; -1; 200 |])
      |> fun z ->
      Arr.concatenate ~axis:0 z
      |> fun e -> Arr.mean ~axis:0 e |> fun e -> Mat.reshape e [| -1; 200 |] in  Array.map new_rates ~f:(fun r ->
        let t = Mat.(r - mean_rates) in 
          Mat.(t *@ ctr_modes) |> fun z -> Arr.reshape z [| 1; -1; 200 |])
      |> fun z ->
      Arr.concatenate ~axis:0 z
      |> Mat.sqr
      |> fun e ->
      Arr.mean ~axis:0 e
      |> Mat.sqrt
      |> fun e -> Mat.reshape e [| -1; 200 |] |> Mat.sum ~axis:0

let _ =
  Array.init 200 ~f:(fun tree ->
    let _ = Stdio.printf "%i %!" tree in 
     let shuffled_idces = List.init 200 ~f:(fun i -> i) |> List.permute in 
    ctrl_num_shuffled shuffled_idces)
  |> Mat.concatenate ~axis:0
  |> fun z -> Mat.save_txt (in_dir "shuffled_num_all") z
