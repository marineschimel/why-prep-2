open Owl
module AD = Algodiff.D
module M = Arm.Make (Arm.Defaults)
open Lib
open Base

(* Setting up the parameters/directories
*)
let dir i =
  Printf.sprintf "/home/mmcs3/rds/rds-t2-cs156-T7o4pEA8QoU/mmcs3/linear_lambda_0.000001/ramping_soc/seed_1_%i" i


let in_dir i =
    Printf.sprintf "%s/%s" (dir i)
  
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


let theta0 = Mat.of_arrays [| [| 0.174533; 2.50532; 0.; 0. |] |] |> AD.pack_arr

let t_preps =
  [| 0.; 0.025; 0.05; 0.07; 0.1; 0.15; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9; 1.0 |]

let w = Mat.load_txt (in_dir 0 "w")


let c i = Mat.load_txt (in_dir i "c")


let x0 = C.broadcast' (fun () -> AD.pack_arr (Mat.(load_txt (in_dir 0 "rates_0_300" )) |> fun m -> Mat.row m 0 |> Mat.transpose))

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


let b = Mat.eye m

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
let save_prms suffix prms = Misc.save_bin (Printf.sprintf "%s/prms_%s" (dir 0) suffix) prms
let save_task suffix task = Misc.save_bin (Printf.sprintf "%s/prms_%s" (dir 0) suffix) task
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
          {  c = (pinned : setter) (AD.pack_arr (c 0))
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
      let readout = R.Readout_P.{ c = (pinned : setter) (AD.pack_arr (c 0)) } in
      let generative = Model.Generative_P.{ prior; dynamics; likelihood } in
      Model.Full_P.{ generative; readout })

module I = Model.ILQR (U) (D0) (L0)

let n = 200

let eigenvalues m =
  let v = Linalg.D.eigvals m in
  let re = Dense.Matrix.Z.re v in
  let im = Dense.Matrix.Z.im v in
  Mat.(concat_horizontal (transpose re) (transpose im))

let w = Mat.load_txt (in_dir 0 "w")
let eigs_w = eigenvalues w
let _ = Mat.save_txt ~out:(in_dir 0 "eigs_w") eigs_w
let  idx_c i = let idx_c = Mat.load_txt (in_dir i " idx_c") in Mat.( idx_c /$ Owl.Linalg.D.norm  idx_c)

let a = let w = Mat.load_txt (in_dir 0 "w") in Mat.(w - eye n)
let disc_a = Mat.(zeros n n - eye n)
let obs_gramian a  idx_c = Linalg.D.lyapunov Mat.(transpose a) Mat.(neg (transpose  idx_c) *@  idx_c)
let ctr_gramian a b = Linalg.D.lyapunov a Mat.(neg b *@ transpose b)
let obs_id = obs_gramian a (Mat.eye 200)
let ctr_disc = ctr_gramian disc_a Mat.(eye n)
let tr_p_disc = Mat.trace ctr_disc


let ctr = ctr_gramian a (Mat.eye 200)
let ctr = Mat.(ctr /$ Float.(tr_p_disc /. 200.))
(* obs_gramian a  idx_c *)



let obs_id_modes, obs_id_eigs =
  let obs_id_modes, _, _ = Linalg.D.svd obs_id in
  obs_id_modes, eigenvalues obs_id

let ctr_modes, ctr_eigs =
  let ctr_modes, _, _ = Linalg.D.svd ctr in
  Mat.(get_slice [ []; [] ] ctr_modes), eigenvalues ctr

let traj idx_c i = Mat.load_txt (in_dir idx_c (Printf.sprintf "rates_%i_300" i))
let n_reaches = 8

module Int = Dynamics.Integrate (D0)
let rates idx_c=
  Array.init n_targets ~f:(fun i ->
      let us =  Mat.load_txt (in_dir idx_c (Printf.sprintf "us_%i_300" i)) in 
        Int.integrate
          ~readout:(AD.pack_arr (c idx_c))
          ~prms:prms.Model.Full_P.generative.dynamics
          ~task:tasks.(i)
          ~n
          ~u:(AD.pack_arr us) |> AD.unpack_arr)

let mean_rates idx_c= 
  Array.map (rates idx_c) ~f:(fun t -> 
Arr.reshape t [| 1; -1; 200 |])
  |> fun z ->
  Arr.concatenate ~axis:0 z
  |> fun e -> Arr.mean ~axis:0 e |> fun e -> Mat.reshape e [| -1; 200 |]


let ctrl_num_unshuffled idx_c =
  Array.map (rates idx_c) ~f:(fun r ->
    let t = Mat.(r - (mean_rates idx_c)) in 
      Mat.(t *@ ctr_modes) |> fun z -> Arr.reshape z [| 1; -1; 200 |])
  |> fun z ->
  Arr.concatenate ~axis:0 z
  |> Mat.sqr
  |> fun e ->
  Arr.mean ~axis:0 e
  |> Mat.sqrt
  |> fun e -> Mat.reshape e [| -1; 200 |] |> Mat.sum ~axis:0

let _ = Array.init 10 ~f:(fun i -> Mat.save_txt (in_dir i "unshuffled_num_all") (ctrl_num_unshuffled i)) 

let shuffle us idces =
  Mat.get_fancy [ R [ 0; -1 ]; L idces ] us


let ctrl_num_shuffled idces idx_c =
  let new_rates =   Array.init n_targets ~f:(fun i ->
    let us =  Mat.load_txt (in_dir idx_c (Printf.sprintf "us_%i_300" i)) in 
    let us = shuffle us idces in
      Int.integrate
        ~readout:(AD.pack_arr (c 0))
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
  Array.init 10 (fun idx_c -> Array.init 200 ~f:(fun tree ->
     let shuffled_idces = List.init 200 ~f:(fun i -> i) |> List.permute in 
    ctrl_num_shuffled shuffled_idces idx_c)
  |> Mat.concatenate ~axis:0
  |> fun z -> Mat.save_txt (in_dir idx_c "shuffled_num_all") z)
