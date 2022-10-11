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
let t1 = 0.3
let t2 = 0.4
let exponent = Cmdargs.(get_int "-exponent" |> force ~usage:"exponent")
let t_coeff = Cmdargs.(get_float "-t_coeff" |> default 1.0)
let g_coeff = Cmdargs.(get_float "-g_coeff" |> default 1.0)
let in_dir s = Printf.sprintf "%s/%s" dir s
let seed = Cmdargs.(get_int "-seed" |> force ~usage:"-seed")
let in_data_dir s = Printf.sprintf "%s/%s" data_dir s
let n_targets = 3
let pause = Cmdargs.(get_float "-pause" |> default 0.5)
let pause_coeff = Cmdargs.(get_float "-pause_coeff" |> default 1.)
let save_all = Cmdargs.(check "-save_all") 
let lambda =
  Cmdargs.(get_float "-lambda" |> force ~usage:"-lambda [dir in which data is]")

let exponent = AD.F (Float.of_int exponent)
let phi_t t = AD.Maths.(t ** exponent)
let scale_mov = Cmdargs.(get_float "-scale_mov" |> default 1.)
let rad = Cmdargs.(get_float "-rad" |> default 0.12)
let t_tot = 0.6

let targets =
  C.broadcast' (fun () ->
      Array.init n_targets ~f:(fun i ->
          let radius = rad in
          let theta = Float.(of_int i *. Const.pi *. 2. /. of_int n_targets) in
          let x = Maths.(radius *. cos theta) in
          let y = Maths.(radius *. sin theta) in
          let t1, t2 = Misc.pos_to_angles x (y +. 0.199112) in
          Mat.of_array [| t1; t2; 0.; 0. |] 1 (-1)))

let double_targets =
  Array.init (n_targets ** 2) ~f:(fun n ->
      let i = Int.(n / n_targets) in
      let j = n - (i * n_targets) in
      let ti = targets.(i) in
      let tj = targets.(j) in
      Mat.(ti @= tj))


let targets = C.broadcast targets
let n_targets = n_targets - 1

let _ =
  C.root_perform (fun () ->
      Mat.save_txt ~out:(in_dir "targets") (Mat.concatenate ~axis:0 targets))

let beta = AD.F 1E-2
let phi_x x = AD.Maths.relu x
let d_phi_x x = AD.Maths.(F 0.5 * (F 1. + signum x))
let d2_phi_x x = AD.Maths.(diagm (F 0. * x))
let link_f x = phi_x x
let double_target i = double_targets.(i)
let dt = 2E-3

(*angle is in hand space (arctan ((y2-y1)/(x2-x1))) and
peak_speed is peak radial speed and length = *)

let lambda_prep = lambda
let lambda_mov = lambda *. scale_mov

let _ =
  Mat.save_txt
    ~out:(in_dir (Printf.sprintf "%s_lambda" subdir))
    (Mat.of_arrays [| [| lambda |] |])

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

let w =
  C.broadcast' (fun () -> Mat.(load_txt (Printf.sprintf "%s/w_rec_%i" data_dir seed)))

let c = C.broadcast' (fun () -> AD.pack_arr Mat.(load_txt (Printf.sprintf "%s/c" dir)))


let b = Mat.eye m

(*structure : get inputs from task, then from the same SOC just generating inputs, 
then from yet another SOC
Maybe do the same thing using 200 - 100 - 50 
what about one population receiving modulated inputs about the target? 
*)
(* x(t+1)- x(t) = Wx(t) + baseline => 0 when x = -W^(-1)*baseline  *)

let x0 = C.broadcast' (fun () -> AD.Maths.(F 0.5 * AD.Mat.uniform ~a:5. ~b:15. m 1))

let x0 =
  C.broadcast' (fun () ->
      let m = Mat.load_txt (in_dir "rates_0_0") in
      AD.pack_arr (Mat.get_slice [ [ 0 ] ] m |> fun z -> Arr.reshape z [| -1; 1 |]))

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
      ( k
      , j
      , Model.
          { t_prep = t_preps.(n_time)
          ; x0
          ; t_movs = [| t1; t2 |]
          ; dt
          ; t_hold = None
          ; t_pauses = Some [| pause; 0.2 |]
          ; scale_lambda = None
          ; target = AD.pack_arr double_target
          ; theta0
          ; tau = 150E-3
          } ))

let save_prms suffix prms =
  Misc.save_bin (Printf.sprintf "%s/%s_prms_%s" dir subdir suffix) prms

let save_task suffix task =
  Misc.save_bin (Printf.sprintf "%s/%s_task_%s" dir subdir suffix) task

let epsilon = 1E-1

module U = Priors.Gaussian

module D0 = Dynamics.Arm_Plus (struct
  let phi_x x = phi_x x
  let d_phi_x x = d_phi_x x
  let phi_u x = x
  let d_phi_u _ = AD.Mat.eye m
end)

(*  *)
module L0 = Likelihoods.Successive_Ramping (struct
  let label = "output"
  let phi_x x = phi_x x
  let phi_t t = phi_t t
end)

module R = Readout

let dt_scaling = Float.(dt /. 1E-3 *. 1000.)

let prms =
  let open Owl_parameters in
  C.broadcast' (fun () ->
      let likelihood =
        Likelihoods.Successive_Ramping_P.
          { c = (pinned : setter) c
          ; c_mask = None
          ; qs_coeff = (pinned : setter) (AD.F Float.(dt_scaling))
          ; t_coeff = (pinned : setter) (AD.F Float.(t_coeff *. dt_scaling))
          ; g_coeff = (pinned : setter) (AD.F Float.(g_coeff *. dt_scaling))
          ; tau_mov_1 = (pinned : setter) (AD.F Float.(t_tot))
          ; tau_mov_2 = (pinned : setter) (AD.F Float.(t_tot))
          ; pause_coeff = (pinned : setter) (AD.F Float.(pause_coeff))
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

let save_results suffix xs us n_prep task =
  let file s = Printf.sprintf "%s/%s_%s_%s" dir subdir s suffix in
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
  let _, _, input_cost_tot =
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
  let rates = AD.unpack_arr (link_f (AD.pack_arr xs)) in
  let _ = 
  Owl.Mat.save_txt ~out:(file "u_cost") (Mat.of_array [| input_cost_tot |] 1 (-1));
  Owl.Mat.save_txt
    ~out:(file "task_cost")
    (Mat.of_array [| torque_err; target_err |] 1 (-1));
  in if save_all then 
  let hands =
    let h =
      AD.Mat.map_by_row
        (fun x -> Arm.unpack_state_diff (M.hand_of (Arm.pack_state x)))
        (AD.pack_arr thetas)
    in
    AD.unpack_arr h
  in
  Owl.Mat.save_txt ~out:(file "max_rates") (Mat.of_array [|(Mat.max' rates)|] 1 (-1))

let attempt_1 =
  let x0 = x0 in
  let _ = save_prms "" prms in
  Array.foldi tasks ~init:[] ~f:(fun i acc (k, j, t) ->
      if Int.(i % C.n_nodes = C.rank)
      then (
        try
          let t_prep = Float.to_int (1000. *. t.t_prep) in
          let xs, us, l, _, success =
            I.solve ~u_init:Mat.(gaussian ~sigma:0. 2001 m) ~n:(m + 4) ~m ~x0 ~prms t
          in
          if success then 
          save_results (Printf.sprintf "%i_%i_%i" k j t_prep) xs us t_prep t;
          Mat.save_txt
            ~out:(in_dir (Printf.sprintf "loss_%i_%i_%i" k j t_prep))
            (Mat.of_array [| AD.unpack_flt l |] 1 (-1));
          (true, (k, j, t)) :: acc
        with
        | _ -> (false, (k, j, t)) :: acc)
      else acc)
  |> C.gather
  |> fun a -> C.broadcast' (fun () -> a |> Array.to_list |> List.concat |> Array.of_list)

let _ =
  Array.iteri attempt_1 ~f:(fun i (b, (k, j, t)) ->
      if Int.(i % C.n_nodes = C.rank)
      then
        if not b
        then (
          let t_prep = Float.to_int (1000. *. t.t_prep) in
          let xs, us, l, _, success =
            I.solve ~u_init:Mat.(gaussian ~sigma:0.0001 2001 m) ~n:(m + 4) ~m ~x0 ~prms t
          in
          if success then 
          save_results (Printf.sprintf "%i_%i_%i" k j t_prep) xs us Int.(t_prep/2) t;
          Mat.save_txt
            ~out:(in_dir (Printf.sprintf "%s_loss_%i_%i_%i" subdir k j t_prep))
            (Mat.of_array [| AD.unpack_flt l |] 1 (-1))))

let _ = C.barrier ()