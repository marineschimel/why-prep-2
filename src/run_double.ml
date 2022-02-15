open Owl
module AD = Algodiff.D
module M = Arm.Make (Arm.Defaults)
open Lib
open Base

let _ = Backtrace.Exn.set_recording true

(* Setting up the parameters/directories
*)
let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")
let subdir = Cmdargs.(get_string "-subdir" |> force ~usage:"-d [dir to save in]")
let data_dir = Cmdargs.(get_string "-data" |> force ~usage:"-data [dir in which data is]")
let in_dir s = Printf.sprintf "%s/%s" dir s
let in_data_dir s = Printf.sprintf "%s/%s" data_dir s
let n_targets = 8

let lambda =
  Cmdargs.(get_float "-lambda" |> force ~usage:"-lambda [dir in which data is]")


let rad = Cmdargs.(get_float "-rad" |> default 0.12)

let targets =
  C.broadcast' (fun () ->
      Array.init n_targets ~f:(fun i ->
          let radius = rad in
          let theta = Float.(of_int i *. Const.pi *. 2. /. of_int n_targets) in
          let x = Maths.(radius *. cos theta) in
          let y = Maths.(radius *. sin theta) in
          (* let x = Stats.uniform_rvs ~a:(Maths.neg radius) ~b:radius in
      let y = Stats.uniform_rvs ~a:(Maths.neg radius) ~b:radius in *)
          let t1, t2 = Misc.pos_to_angles x (y +. 0.199112) in
          Mat.of_array [| t1; t2; 0.; 0. |] 1 (-1)))


(* let double_targets =
  Array.init n_targets ~f:(fun n ->
      let ti = targets.(n) in
      let tj = Mat.of_arrays [| [| 0.174533; 2.50532; 0.; 0. |] |] in
      Mat.(ti @= tj)) *)

let double_targets =
  Array.init (n_targets ** 2) ~f:(fun n ->
      let i = Int.(n / n_targets) in
      let j = n - (i * n_targets) in
      let ti = targets.(i) in
      let tj = targets.(j) in
      Mat.(ti @= tj))


(* Array.map targets ~f:(fun ti ->
      Array.map targets ~f:(fun t -> Mat.(ti @= t))
      |> fun a -> Array.sub a ~pos:1 ~len:(n_targets - 1)) *)

let targets = C.broadcast targets
let n_targets = n_targets - 1

let _ =
  C.root_perform (fun () ->
      Mat.save_txt ~out:(in_dir "targets") (Mat.concatenate ~axis:0 targets))


let beta = AD.F 1E-2
let phi_x x = AD.Maths.relu x
let d_phi_x x = AD.Maths.(F 0.5 * (F 1. + signum x))
let d2_phi_x x = AD.Maths.(diagm (F 0. * x))
(* let phi_x x = AD.Maths.relu x
let d_phi_x x = AD.Maths.(F 0.5 * (F 1. + signum x))
let d2_phi_x x = AD.Maths.(diagm (F 0. * x)) *)

(* let phi_x x = AD.Maths.(beta * AD.requad (x / beta))
let d_phi_x x = AD.Maths.(AD.d_requad (x / beta))
let d2_phi_x x = AD.Maths.(F 1. / beta * AD.d2_requad (x / beta)) *)
let link_f x = phi_x x
let double_target i = double_targets.(i)
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
    (Mat.of_array [| tau; lambda_prep; lambda_mov; dt; AD.unpack_flt beta |] 1 (-1))


let theta0 = Mat.of_arrays [| [| 0.174533; 2.50532; 0.; 0. |] |] |> AD.pack_arr
let t_preps = [| 0.5 |]
let w = C.broadcast' (fun () -> Mat.(load_txt (Printf.sprintf "%s/w_rec" data_dir)))
let c = C.broadcast' (fun () -> AD.pack_arr Mat.(load_txt (Printf.sprintf "%s/c" dir)))

(* let c = C.broadcast' (fun () -> AD.Mat.gaussian ~sigma:0.003 2 m) *)
(* C.broadcast' (fun () ->
      AD.pack_arr Mat.(load_txt (Printf.sprintf "%s/opt_readout" dir))) *)

let b = Mat.eye m

(*structure : get inputs from task, then from the same SOC just generating inputs, 
then from yet another SOC
Maybe do the same thing using 200 - 100 - 50 
what about one population receiving modulated inputs about the target? 
*)
(* x(t+1)- x(t) = Wx(t) + baseline => 0 when x = -W^(-1)*baseline  *)

let x0 =
  C.broadcast' (fun () ->
      AD.Maths.transpose
        (AD.pack_arr
           (Mat.get_slice [ [ 0 ] ] Mat.(load_txt (Printf.sprintf "%s/xs_1_500" dir)))))


let baseline_input =
  C.broadcast' (fun () ->
      AD.Maths.(neg ((AD.pack_arr w *@ link_f x0) - x0)) |> AD.Maths.transpose)


let u0 = phi_x x0
let norm_u0 = AD.Maths.(l2norm_sqr' u0)
let c = AD.Maths.(c - (c *@ u0 *@ transpose u0 / norm_u0))

let x0 =
  AD.Maths.concatenate [| AD.Maths.transpose theta0; x0 |] ~axis:0 |> AD.Maths.transpose


let tasks =
  Array.init
    (n_targets * n_targets * Array.length t_preps)
    ~f:(fun i ->
      let k = i / (n_targets * Array.length t_preps) in
      let next_i = Int.rem i (n_targets * Array.length t_preps) in
      let j = next_i / Array.length t_preps in
      let n_time = Int.rem next_i (Array.length t_preps) in
      let double_target =
        let ti = targets.(k) in
        let tj = targets.(j) in
        Mat.(ti @= tj)
      in
      Model.
        { t_prep = t_preps.(n_time)
        ; x0
        ; t_movs = [| 0.3; 0.35 |]
        ; dt
        ; t_hold = None
        ; t_pauses = Some [| 0.1; 0.2 |]
        ; scale_lambda = None
        ; target = AD.pack_arr double_target
        ; theta0
        ; tau = 150E-3
        })


let save_prms suffix prms =
  Misc.save_bin (Printf.sprintf "%s/%s/prms_%s" dir subdir suffix) prms


let save_task suffix task =
  Misc.save_bin (Printf.sprintf "%s/%s/task_%s" dir subdir suffix) task


let epsilon = 1E-1

module U = Priors.Gaussian

module D0 = Dynamics.Arm_Plus (struct
  let phi_x x = phi_x x
  let d_phi_x x = d_phi_x x
  let phi_u x = x
  let d_phi_u _ = AD.Mat.eye m
end)

(*  *)
module L0 = Likelihoods.Successive (struct
  let label = "output"
  let phi_x x = phi_x x
  let d_phi_x x = d_phi_x x
  let d2_phi_x x = d2_phi_x x
end)

module R = Readout

let prms =
  let open Owl_parameters in
  C.broadcast' (fun () ->
      let likelihood =
        Likelihoods.Successive_P.
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

let save_results suffix xs us n_prep task =
  let file s = Printf.sprintf "%s/%s/%s_%s" dir subdir s suffix in
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
  let thetas, xs, us =
    Mat.get_slice [ []; [ 0; 3 ] ] xs, Mat.get_slice [ []; [ 4; -1 ] ] xs, us
  in
  let x0 = Mat.get_slice [ []; [ 4; -1 ] ] (AD.unpack_arr x0) in
  let rates = AD.unpack_arr (link_f (AD.pack_arr xs)) in
  Owl.Mat.save_txt
    ~out:(file "task_cost")
    (Mat.of_array [| torque_err; target_err |] 1 (-1));
  Owl.Mat.save_txt ~out:(file "thetas") thetas;
  Owl.Mat.save_txt ~out:(file "xs") xs;
  Owl.Mat.save_txt ~out:(file "us") us;
  Owl.Mat.save_txt ~out:(file "rates") rates;
  Owl.Mat.save_txt ~out:(file "eff_us") (AD.unpack_arr (link_f (AD.pack_arr us)));
  Owl.Mat.save_txt
    ~out:(file "torques")
    Mat.((rates - AD.unpack_arr (link_f (AD.pack_arr x0))) *@ transpose (AD.unpack_arr c))


let _ =
  let x0 = x0 in
  let _ = save_prms "" prms in
  try
    Array.iteri tasks ~f:(fun i t ->
        if Int.(i % C.n_nodes = C.rank)
        then (
          let k = i / (n_targets * Array.length t_preps) in
          let next_i = Int.rem i (n_targets * Array.length t_preps) in
          let j = next_i / Array.length t_preps in
          let t_prep = Float.to_int (1000. *. t.t_prep) in
          let xs, us, l =
            I.solve ~u_init:Mat.(gaussian ~sigma:0. 2001 m) ~n:(m + 4) ~m ~x0 ~prms t
          in
          save_results (Printf.sprintf "%i_%i_%i" k j t_prep) xs us t_prep t;
          Mat.save_txt
            ~out:(in_dir (Printf.sprintf "loss_%i_%i_%i" k j t_prep))
            (Mat.of_array [| AD.unpack_flt l |] 1 (-1))))
  with
  | _ -> ()


let _ = C.barrier ()
