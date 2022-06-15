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

let lambda =
  Cmdargs.(get_float "-lambda" |> force ~usage:"-lambda [dir in which data is]")


let rad = Cmdargs.(get_float "-rad" |> default 0.5)
let in_dir s = Printf.sprintf "%s/%s" dir s
let in_data_dir s = Printf.sprintf "%s/%s" data_dir s
let n_targets = 1520

let targets =
  C.broadcast' (fun () ->
      Array.init n_targets ~f:(fun i ->
          let radius = Stats.uniform_rvs ~a:0.1 ~b:0.16 in
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

(* let phi_x x = x
let d_phi_x x = AD.Maths.(F 1. + (F 0. * x))
let d2_phi_x x = AD.Maths.(diagm (F 0. * x)) *)
let phi_x x = x
let d_phi_x x = AD.Maths.((F 1. + F 0. * x))
let d2_phi_x x = AD.Maths.(diagm (F 0. * x))

(* let phi_x x = AD.Maths.(beta * AD.requad (x / beta))
let d_phi_x x = AD.Maths.(AD.d_requad (x / beta))
let d2_phi_x x = AD.Maths.(F 1. / beta * AD.d2_requad (x / beta)) *)
let link_f x = phi_x x
let target i = targets.(i)
let dt = 2E-3
let lambda_prep = lambda
let lambda_mov = lambda *. 1000. 
let n_out = 2
let n = 54
let m = 50
let tau = 150E-3
let n_output = 2

let _ =
  Mat.save_txt
    ~out:(in_dir "prms")
    (Mat.of_array [| tau; lambda_prep; lambda_mov; dt; AD.unpack_flt beta; rad |] 1 (-1))


let theta0 = Mat.of_arrays [| [| 0.174533; 2.50532; 0.; 0. |] |] |> AD.pack_arr

(* let t_preps = [| 0.; 0.05; 0.1; 0.15; 0.2; 0.3; 0.45; 0.5; 0.6; 0.8; 1. |] *)
let t_preps = [| 0.4 |]
let w = C.broadcast' (fun () -> Mat.(load_txt (Printf.sprintf "%s/mini_w_rec" data_dir)))

(* let _ =
  C.root_perform (fun () ->
      Mat.save_txt ~out:(Printf.sprintf "%s/w" dir) (Mat.gaussian ~sigma:0.05 m m))


let w = C.broadcast' (fun () -> Mat.(load_txt (Printf.sprintf "%s/w" dir))) *)
let n_ = Int.(n-4)
let c = C.broadcast' (fun () -> AD.Mat.gaussian ~sigma:Float.(rad / sqrt (of_int n_)) 2 n_)

(* let c = C.broadcast' (fun () -> AD.pack_arr Mat.(load_txt (Printf.sprintf "%s/c" dir))) *)

let _ =
  C.root_perform (fun () ->
      Mat.save_txt ~out:(Printf.sprintf "%s/c" dir) (AD.unpack_arr c))


let b = Mat.eye n_

let x0 = C.broadcast' (fun () -> 
  AD.Maths.(F 0.5 * AD.Mat.uniform ~a:5. ~b:20. n_ 1))
(*structure : get inputs from task, then from the same SOC just generating inputs, 
then from yet another SOC
Maybe do the same thing using 200 - 100 - 50 
what about one population receiving modulated inputs about the target? 
*)
(* x(t+1)- x(t) = Wx(t) + baseline => 0 when x = -W^(-1)*baseline  *)
(* let baseline_input = AD.Mat.ones 1 m *)

let tasks =
  Array.init n_targets ~f:(fun i ->
      (* let x0 = AD.Maths.(F 0.5 * AD.Mat.uniform ~a:5. ~b:20. m 1) in *)
      let n_target = Int.rem i n_targets in
      let t_hold = Stats.uniform_rvs ~a:0.15 ~b:0.3 in
      let t_mov = Float.(0.3 +. 0.2 -. t_hold) in
      Model.
        { t_prep = 0.4
        ; x0 =
            AD.Maths.concatenate [| AD.Maths.transpose theta0; x0 |] ~axis:0
            |> AD.Maths.transpose
        ; t_movs = [| t_mov |]
        ; dt
        ; t_hold = Some t_hold
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

module U = Priors.Gaussian
(* 
_Phi (struct
let phi_u x = phi_x x
  let d_phi_u x =  d_phi_x x
  let d2_phi_u x = d2_phi_x x
end) *)

module D0 = Dynamics.Arm_Discrete (struct
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
          ; bias = (pinned : setter) (AD.Mat.zeros 1 m)
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

let save_results suffix xs us =
  let file s = Printf.sprintf "%s/%s_%s" dir s suffix in
  let xs = AD.unpack_arr xs in
  let rates = AD.unpack_arr (link_f (AD.pack_arr xs)) in
  Owl.Mat.save_txt ~out:(file "rates") rates;
  Owl.Mat.save_txt ~out:(file "us") (AD.unpack_arr us)


(* let torque_err, target_err =
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
    (Mat.of_array [| Float.of_int t_to_target |] 1 (-1)); *)
(* Owl.Mat.save_txt ~out:(file "thetas") thetas;
  Owl.Mat.save_txt ~out:(file "xs") xs;
  Owl.Mat.save_txt ~out:(file "us") us; *)

(* Owl.Mat.save_txt ~out:(file "eff_us") (AD.unpack_arr (link_f (AD.pack_arr us))) *)

let replace_baseline prms x0 =
  let open Model in
  let open Full_P.A in
  let open Generative_P.A in
  let open Owl_parameters in
  let open Dynamics_typ in 
  let b = Owl_parameters.extract prms.Full_P.generative.Generative_P.dynamics.Arm_Plus_P.b in
  let a = Owl_parameters.extract prms.Full_P.generative.Generative_P.dynamics.Arm_Plus_P.a in
  let a_discrete =
      let a = let n = AD.Mat.row_num a in
      AD.Maths.(a - AD.Mat.eye n) 
    in
    Linalg.D.(expm Mat.(AD.unpack_arr a *$ Float.(dt /. tau))) |> AD.pack_arr in 
  let bbinv = AD.Linalg.linsolve AD.Maths.((transpose b)*@(b)) (AD.Mat.eye n_) in 
  let x0 = AD.Maths.transpose x0 in 
  let u_base = AD.Maths.(x0 + neg ((link_f x0)*@a_discrete)) |> fun z -> AD.Maths.(b*@bbinv*@(transpose z)) |> AD.Maths.transpose in
  let _ = AD.Mat.print AD.Maths.(u_base*@b - x0 + ((link_f x0)*@a_discrete)) in 
  let new_prms =
    Accessor.map
      (generative @> dynamics @> Dynamics.Arm_Plus_P.A.bias)
      ~f:(fun _ -> (pinned : setter) u_base)
  in
  new_prms prms


let _ =
  let _ = save_prms "" prms in
  Array.mapi tasks ~f:(fun i t ->
      if Int.(i % C.n_nodes = C.rank)
      then (
        (* try *)
          let xs, us, _l, _ =
            I.solve
              ~u_init:Mat.(gaussian ~sigma:0. 2001 m)
              ~n:n
              ~m
              ~x0:t.x0
              ~prms:
                (replace_baseline
                   prms
                   AD.Maths.(get_slice [ [ 4; -1 ] ] (transpose t.x0)))
              t
          in
          save_results (Printf.sprintf "%i" i) xs us))
        (* with
        | _ -> ())) *)


(* save_results (Printf.sprintf "%i" i) xs us n_target t_prep t;
        Mat.save_txt
          ~out:(in_dir (Printf.sprintf "loss_%i" i))
          (Mat.of_array [| AD.unpack_flt l |] 1 (-1));
        save_task (Printf.sprintf "%i" i) t)) *)

let _ = C.barrier ()

(* let xs0, us0,_ =
  I0.solve
    ~u_init:(Mat.((sqr (Mat.gaussian ~sigma:0.01 901 m))))
    ~arm:true
    ~n:_n
    ~m
    ~prms
    task_noprep *)

(* let _ =
  Mat.save_txt ~out:(in_dir "eye_us_0") (AD.unpack_arr us0);
  Mat.save_txt ~out:(in_dir "eye_xs_0") (AD.unpack_arr xs0)


let xs_300, us_300,_ =
  I0.solve
    ~u_init:(Mat.( (sqr (Mat.gaussian ~sigma:0.01 901 m))))
    ~arm:true
    ~n:_n
    ~m
    ~prms
    task_prep


let _ =
  Mat.save_txt ~out:(in_dir "eye_us_300") (AD.unpack_arr us_300);
  Mat.save_txt ~out:(in_dir "eye_xs_300") (AD.unpack_arr xs_300) *)
