open Owl
module AD = Algodiff.D
module M = Arm.Make (Arm.Defaults)
open Lib
open Base

(* Setting up the parameters/directories
*)
let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")
let data_dir = Cmdargs.(get_string "-data" |> force ~usage:"-data [dir in which data is]")
let nu = Cmdargs.(get_int "-nu" |> force ~usage:"-nu [value of the scale coeff]")
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
          let t1, t2 = Misc.pos_to_angles x (y +. 0.199112) in
          Mat.of_array [| t1; t2; 0.; 0. |] 1 (-1)))

let _ =
  C.root_perform (fun () ->
      Mat.save_txt ~out:(in_dir "targets") (Mat.concatenate ~axis:0 targets))

let beta = AD.F 1E-2
(* let phi_x x = x
let d_phi_x x = AD.Maths.((F 0. * x) + F 1.) *)

let phi_x x = AD.Maths.relu x
let d_phi_x x = AD.Maths.(F 0.5 * (F 1. + signum x))
let d2_phi_x _x = AD.F 0.
let link_f x = phi_x x
let target i = targets.(i)
let dt = 2E-3
let lambda_prep = 1E-6
let lambda_mov = 1E-6
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
let t_preps = [| 0.; 0.05; 0.1; 0.15; 0.3; 0.6; 0.8 |]
let w = C.broadcast' (fun () -> Mat.(load_txt (Printf.sprintf "%s/w_rec" data_dir)))
let c = C.broadcast' (fun () -> AD.Mat.gaussian ~sigma:0.05 2 m)
let x0 = C.broadcast' (fun () -> AD.Mat.uniform ~a:2. ~b:10. m 1)
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

let _ =
  C.root_perform (fun () ->
      Mat.save_txt ~out:(Printf.sprintf "%s/x0" dir) (AD.unpack_arr x0))

let b = Mat.eye m

(*structure : get inputs from task, then from the same SOC just generating inputs, 
then from yet another SOC
Maybe do the same thing using 200 - 100 - 50 
what about one population receiving modulated inputs about the target? 
*)
(* x(t+1)- x(t) = Wx(t) + baseline => 0 when x = -W^(-1)*baseline  *)
let baseline_input =
  C.broadcast' (fun () ->
      AD.Maths.(neg ((AD.pack_arr w *@ link_f x0) - x0)) |> AD.Maths.transpose)

let x0 =
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
        ; t_mov = 0.4
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

(* let phi_x x = let x = AD.unpack_arr x in let y = Mat.map (fun x -> if Float.(x > epsilon) then x else if Float.(x > neg epsilon) 
    then Float.((x +. epsilon)**2. /. ((4. *. epsilon))) 
else 0.) x in AD.pack_arr y 

let d_phi_x x = let x = AD.unpack_arr x in let y = Mat.map (fun x -> if Float.(x > epsilon) then 1. 
else if Float.(x > neg epsilon)  then ((x +. epsilon)/.(2. *. epsilon)) else 0.) x in AD.pack_arr y

let d2_phi_x x = let x = AD.unpack_arr x in let y = Mat.map (fun x -> if Float.(x > epsilon) then 0. 
else if Float.(x > neg epsilon) then Float.(1./(2. * epsilon)) else 0.) x in AD.pack_arr y *)

module U = Priors.Student

module D0 = Dynamics.Arm_Plus (struct
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
        (* Priors.Sparse_P.
          { lambda_prep = (pinned : setter) (AD.F lambda_mov)
          ; lambda_mov = (pinned : setter) (AD.F lambda_mov)
          } *)
        Priors.Student_P.
          { nu = (pinned : setter) (AD.F (Float.of_int nu))
          ; spatial_stds = (pinned : setter) (AD.Mat.create 1 m (1. /. lambda_mov))
          }
      in
      let readout = R.Readout_P.{ c = (pinned : setter) c } in
      let generative = Model.Generative_P.{ prior; dynamics; likelihood } in
      Model.Full_P.{ generative; readout })

let save_results suffix xs us n_target n_prep task =
  let file s = Printf.sprintf "%s/%s_%s" dir s suffix in
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
  let t_to_target = Analysis_funs.time_to_end xs targets.(n_target) in
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
    Mat.((rates - AD.unpack_arr (link_f (AD.pack_arr x0))) *@ transpose (AD.unpack_arr c))

module I = Model.ILQR (U) (D0) (L0)

let _ =
  let x0 = x0 in
  let _ = save_prms "" prms in
  (* try *)
  Array.iteri tasks ~f:(fun i t ->
      if Int.(i % C.n_nodes = C.rank)
      then (
        let n_target = Int.rem i n_targets in
        let t_prep = Float.to_int (1000. *. t.t_prep) in
        let xs, us, l =
          I.solve ~u_init:Mat.(gaussian ~sigma:0. 2001 m) ~n:(m + 4) ~m ~x0 ~prms t
        in
        save_results (Printf.sprintf "%i_%i" n_target t_prep) xs us n_target t_prep t;
        Mat.save_txt
          ~out:(in_dir (Printf.sprintf "loss_%i_%i" n_target t_prep))
          (Mat.of_array [| AD.unpack_flt l |] 1 (-1));
        save_task (Printf.sprintf "%i_%i" n_target t_prep) t))

(* with
  | e -> Stdio.printf "%s" (Exn.to_string e) *)

let _ = C.barrier ()
