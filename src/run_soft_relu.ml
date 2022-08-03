open Owl
module AD = Algodiff.D
module M = Arm.Make (Arm.Defaults)
open Lib
open Base

(* Setting up the parameters/directories
*)
let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")
let data_dir = Cmdargs.(get_string "-data" |> force ~usage:"-data [dir in which data is]")
let in_dir s = Printf.sprintf "%s/%s" dir s
let in_data_dir s = Printf.sprintf "%s/%s" data_dir s

let targets =
  C.broadcast' (fun () -> Mat.load_txt (Printf.sprintf "%s/target_thetas" data_dir))

let target i = Mat.row targets i
let dt = 1E-3
let lambda_prep = 1E-6
let lambda_mov = 1E-6
let n_out = 2
let _n = 204
let m = 200
let tau = 150E-3
let n_output = 2
let n_targets = 8
let theta0 = Mat.of_arrays [| [| 0.174533; 2.50532; 0.; 0. |] |] |> AD.pack_arr
let t_preps = [| 0.; 0.05; 0.1; 0.15; 0.2; 0.3; 0.4; 0.5; 0.6; 0.8; 1. |]
let target_num = 0
let w = C.broadcast' (fun () -> Mat.(load_txt (Printf.sprintf "%s/w_rec" data_dir)))

let c =
  C.broadcast' (fun () -> AD.pack_arr Mat.(load_txt (Printf.sprintf "%s/c" data_dir)))

(* let c = C.broadcast' (fun () -> AD.pack_arr ((Mat.(load_txt (Printf.sprintf "%s/opt_c_soc" data_dir))))) *)
let b =
  Mat.concatenate
    ~axis:0
    [| Mat.concatenate ~axis:1 [| Mat.eye 120; Mat.zeros 120 80 |]
     ; Mat.concatenate ~axis:1 [| Mat.(neg (eye 80)); Mat.zeros 80 120 |]
    |]

(*structure : get inputs from task, then from the same SOC just generating inputs, 
then from yet another SOC
Maybe do the same thing using 200 - 100 - 50 
what about one population receiving modulated inputs about the target? 
*)

let tasks =
  Array.init
    (Array.length t_preps * n_targets)
    ~f:(fun i ->
      let n_time = i / n_targets in
      let n_target = Int.rem i n_targets in
      Model.
        { t_prep = t_preps.(n_time)
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

let beta = AD.F 1E-2
let phi_x x = AD.Maths.(beta * AD.requad (x / beta))
let d_phi_x x = AD.Maths.(AD.d_requad (x / beta))
let d2_phi_x x = AD.Maths.(F 1. / beta * AD.d2_requad (x / beta))
let link_f x = phi_x x

let save_results suffix xs us =
  let file s = Printf.sprintf "%s/%s_%s" dir s suffix in
  let thetas, xs, us =
    ( Mat.get_slice [ []; [ 0; 3 ] ] (AD.unpack_arr xs)
    , Mat.get_slice [ []; [ 4; -1 ] ] (AD.unpack_arr xs)
    , AD.unpack_arr us )
  in
  let rates = AD.unpack_arr (link_f (AD.pack_arr xs)) in
  Owl.Mat.save_txt ~out:(file "thetas") thetas;
  Owl.Mat.save_txt ~out:(file "xs") xs;
  Owl.Mat.save_txt ~out:(file "us") us;
  Owl.Mat.save_txt ~out:(file "rates") rates;
  Owl.Mat.save_txt ~out:(file "eff_us") (AD.unpack_arr (link_f (AD.pack_arr us)));
  Owl.Mat.save_txt ~out:(file "torques") Mat.(rates *@ transpose (AD.unpack_arr c))

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
  let phi_u x = AD.Maths.relu x
  let d_phi_u u = AD.Maths.(diagm (F 0.5 * (signum u + F 1.)))
end)

(*  *)
module L0 = Likelihoods.End_Phi (struct
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
          { a = (pinned : setter) (AD.pack_arr Mat.(transpose w - eye m))
          ; b = (pinned : setter) (AD.pack_arr b)
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

let _ =
  let x0 =
    AD.Maths.concatenate ~axis:1 [| theta0; AD.Maths.(F (-1.) * AD.Mat.ones 1 m) |]
  in
  let _ = save_prms "" prms in
  Array.mapi tasks ~f:(fun i t ->
      if Int.(i % C.n_nodes = C.rank)
      then (
        let n_target = Int.rem i n_targets in
        let t_prep = Float.to_int (1000. *. t.t_prep) in
        let xs, us, _ =
          I.solve ~u_init:Mat.(gaussian ~sigma:0.01 2001 m) ~n:(m + 4) ~m ~x0 ~prms t
        in
        save_results (Printf.sprintf "%i_%i" n_target t_prep) xs us;
        save_task (Printf.sprintf "%i_%i" n_target t_prep) t))

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
