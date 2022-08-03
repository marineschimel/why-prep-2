open Owl
module AD = Algodiff.D
module M = Arm.Make (Arm.Defaults)
open Lib
open Base
open Accessor.O

let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")
let data_dir = Cmdargs.(get_string "-data" |> force ~usage:"-data [dir in which data is]")
let in_dir s = Printf.sprintf "%s/%s" dir s
let in_data_dir s = Printf.sprintf "%s/%s" data_dir s
let n_targets = 32
let n_times = 4
let beta = AD.F 1E-2
let phi_x x = AD.Maths.relu x
let d_phi_x x = AD.Maths.(F 0.5 * (F 1. + signum x))
let d2_phi_x _x = AD.F 0.
(* let phi_x x = AD.Maths.(beta * AD.requad (x / beta))
let d_phi_x x = AD.Maths.(AD.d_requad (x / beta))
let d2_phi_x x = AD.Maths.(F 1. / beta * AD.d2_requad (x / beta)) *)

(* let targets = C.broadcast' (fun () -> Mat.load_txt (Printf.sprintf "%s/target_thetas" dir))  *)
let targets =
  C.broadcast' (fun () ->
      Array.init n_targets ~f:(fun _ ->
          let radius = Stats.uniform_rvs ~a:0.08 ~b:0.14 in
          let theta = Stats.uniform_rvs ~a:0. ~b:(2. *. Const.pi) in
          let x = Maths.(radius *. cos theta) in
          let y = Maths.(radius *. sin theta) in
          (* let x = Stats.uniform_rvs ~a:(Maths.neg radius) ~b:radius in
          let y = Stats.uniform_rvs ~a:(Maths.neg radius) ~b:radius in *)
          let t1, t2 = Misc.pos_to_angles x (y +. 0.199112) in
          [| t1; t2; 0.; 0. |])
      |> Mat.of_arrays)

let _ = C.root_perform (fun () -> Mat.save_txt ~out:(in_dir "target_thetas") targets)
let target i = Mat.row targets i
let t_prep = 0.3
let dt = 10E-3
let lambda_prep = 5E-3
let lambda_mov = 5E-3
let n_out = 2
let n = 204
let m = 200
let tau = 150E-3
let n_output = 2
let w = C.broadcast' (fun () -> Mat.(load_txt (Printf.sprintf "%s/w_rec" data_dir)))

(* let w = C.broadcast' (fun () -> Mat.gaussian ~sigma:0.001 m m) *)
let theta0 = Mat.of_arrays [| [| 0.174533; 2.50532; 0.; 0. |] |] |> AD.pack_arr
let t_preps = [| 0. |]
let c = C.broadcast' (fun () -> AD.Mat.gaussian ~sigma:0.003 2 m)
let x0 = C.broadcast' (fun () -> AD.Mat.uniform ~a:1. ~b:2. m 1)

let _ =
  C.root_perform (fun () ->
      Mat.save_txt ~out:(Printf.sprintf "%s/x0" dir) (AD.unpack_arr x0))

let u0 = phi_x x0
let norm_u0 = AD.Maths.(l2norm_sqr' u0)
let c = AD.Maths.(c - (c *@ u0 *@ transpose u0 / norm_u0))
let norm_c0 = AD.Maths.l2norm' c

let _ =
  C.root_perform (fun () ->
      Mat.save_txt ~out:(Printf.sprintf "%s/c0" dir) (AD.unpack_arr c))

let baseline_input =
  C.broadcast' (fun () ->
      AD.Maths.(neg ((AD.pack_arr w *@ phi_x x0) - x0)) |> AD.Maths.transpose)

let x0 =
  AD.Maths.concatenate [| AD.Maths.transpose theta0; x0 |] ~axis:0 |> AD.Maths.transpose

let tasks =
  C.broadcast' (fun () ->
      Array.init n_targets ~f:(fun i ->
          (* let n_time = i / n_targets in *)
          let n_target = Int.rem i n_targets in
          Model.
            { t_prep = Stats.uniform_rvs ~a:0. ~b:0.
            ; t_mov = Stats.uniform_rvs ~a:0.3 ~b:0.4
            ; dt
            ; x0
            ; t_hold = Some 0.1
            ; scale_lambda = None
            ; target = AD.pack_arr (target n_target)
            ; t_pauses = None
            ; theta0
            ; tau = 150E-3
            }))

let _ = C.print (Printf.sprintf "array len : %i %!" (Array.length tasks))

(* let save_prms suffix prms = Misc.save_bin (Printf.sprintf "%s/prms_%s" dir suffix) prms
let save_task suffix task = Misc.save_bin (Printf.sprintf "%s/prms_%s" dir suffix) task *)

module U = Priors.Gaussian
module D = Dynamics.Arm_Linear

(* module D = Dynamics.Arm_Plus (struct
  let phi_x x = phi_x x
  let d_phi_x x = d_phi_x x
  let phi_u x = x
  let d_phi_u _ = AD.Mat.eye m
end) *)
module L = Likelihoods.End (struct
  let label = "output"
end)
(* module L = Likelihoods.End_Phi (struct
  let label = "output"
  let phi_x x = phi_x x
  let d_phi_x x = d_phi_x x
  let d2_phi_x x = d2_phi_x x
  let speed_end_penalty = 0.1
end) *)

module R = Readout

let prms =
  let open Owl_parameters in
  C.broadcast' (fun () ->
      let likelihood =
        Likelihoods.End_P.
          { c = (pinned : setter) c
          ; c_mask = None
          ; qs_coeff = (pinned : setter) (AD.F 1.)
          ; t_coeff = (pinned : setter) (AD.F 0.5)
          ; g_coeff = (pinned : setter) (AD.F 5.)
          }
      in
      let dynamics =
        Dynamics.Arm_Linear_P.
          { a = (pinned : setter) (AD.pack_arr Mat.(transpose w))
          ; b = (pinned : setter) (AD.Mat.eye m)
          }
      in
      let prior =
        Priors.Gaussian_P.
          { lambda_prep = (pinned : setter) (AD.F lambda_prep)
          ; lambda_mov = (pinned : setter) (AD.F lambda_mov)
          }
      in
      let readout = R.Readout_P.{ c = (learned : setter) c } in
      let generative = Model.Generative_P.{ prior; dynamics; likelihood } in
      Model.Full_P.{ generative; readout })

module I = Model.ILQR (U) (D) (L)

let project_c c =
  let open Owl_parameters in
  let c = extract c in
  let norm_c = AD.Maths.l2norm' c in
  (learned : setter)
    AD.Maths.((c - (c *@ u0 *@ transpose u0 / norm_u0)) / norm_c * norm_c0)

let loss ~u_init ~prms t =
  let prms =
    C.broadcast
      (Accessor.map (Model.Full_P.A.readout @> Readout.Readout_P.A.c) prms ~f:project_c)
  in
  let _, us, l =
    match u_init with
    | Some u_init ->
      let _ = C.print "Some u_init  : " in
      I.solve ~u_init ~opt:true ~n ~m ~x0 ~prms t
    | None -> I.solve ~opt:true ~n ~m ~x0 ~prms t
  in
  let c = Owl_parameters.extract prms.readout.c in
  AD.Maths.((l + (F 0.01 * l2norm_sqr' c)) / F (Float.of_int n_targets)), AD.unpack_arr us

let save_results suffix prms tasks =
  Array.iteri tasks ~f:(fun i t ->
      if Int.(i % C.n_nodes = C.rank)
      then (
        let file s = Printf.sprintf "%s.%s_%i" suffix s i in
        let xs, us, _ = I.solve ~x0 ~n ~m ~prms t in
        let thetas, xs, us =
          ( Mat.get_slice [ []; [ 0; 3 ] ] (AD.unpack_arr xs)
          , Mat.get_slice [ []; [ 4; -1 ] ] (AD.unpack_arr xs)
          , AD.unpack_arr us )
        in
        Owl.Mat.save_txt ~out:(file "thetas") thetas;
        Owl.Mat.save_txt ~out:(file "xs") xs;
        Owl.Mat.save_txt ~out:(file "us") us))

let final_prms =
  let in_each_iteration ~prms k =
    if Int.(k % 30 = 0) then save_results (in_dir "train") prms tasks;
    if Int.(k % 30 = 0)
    then
      C.root_perform (fun () ->
          Mat.save_txt
            ~out:(in_dir "nullspace")
            (AD.unpack_arr AD.Maths.(x0 *@ transpose c)))
  in
  I.train
    ~max_iter:2000
    ~loss
    ~recycle_u:true
    ~eta:(`of_iter (fun k -> Float.(0.008 / (1. + sqrt (of_int k / 10.)))))
    ~init_prms:prms
    ~in_each_iteration
    ~save_progress_to:(1, 2, in_dir "progress")
    tasks

(* 
let _ =
  Mat.save_txt ~out:"c_test" (AD.unpack_arr (Owl_parameters.extract final_prms.readout.c)) *)
