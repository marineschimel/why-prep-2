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
let targets = Mat.load_txt (Printf.sprintf "%s/target_thetas" data_dir)
let target i = Mat.row targets i
let dt = 1E-3
let lambda_prep = 1E-4
let lambda_mov = 1E-4
let n_out = 2
let _n = 204
let m = 200
let tau = 150E-3
let n_output = 2
let n_targets = 8
let theta0 = Mat.of_arrays [| [| 0.174533; 2.50532; 0.; 0. |] |] |> AD.pack_arr
let t_prep = 0.
let target_num = 0
let c = Mat.gaussian ~sigma:0.5 n_out m
let b = Mat.gaussian ~sigma:0.5 m m
let _ = Mat.save_txt ~out:(in_dir "c") c
let _ = Mat.save_txt ~out:(in_dir "b") b
let b = Mat.eye 200

(*structure : get inputs from task, then from the same SOC just generating inputs, 
then from yet another SOC
Maybe do the same thing using 200 - 100 - 50 
what about one population receiving modulated inputs about the target? 
*)

let task_noprep =
  Model.
    { t_prep = 0.
    ; t_mov = 0.4
    ; dt
    ; t_hold = Some 0.2
    ; scale_lambda = None
    ; t_pauses = None
    ; target = AD.pack_arr (target target_num)
    ; theta0
    ; tau = 150E-3
    }


let task_prep =
  Model.
    { t_prep = 0.3
    ; t_mov = 0.4
    ; dt
    ; t_hold = Some 0.2
    ; scale_lambda = None
    ; t_pauses = None
    ; target = AD.pack_arr (target target_num)
    ; theta0
    ; tau = 150E-3
    }


let save_results suffix xs us =
  let file s = Printf.sprintf "%s/%s_%s" dir s suffix in
  let thetas, xs, us =
    ( Mat.get_slice [ []; [ 0; 3 ] ] (AD.unpack_arr xs)
    , Mat.get_slice [ []; [ 4; -1 ] ] (AD.unpack_arr xs)
    , AD.unpack_arr us )
  in
  Owl.Mat.save_txt ~out:(file "thetas") thetas;
  Owl.Mat.save_txt ~out:(file "xs") xs;
  Owl.Mat.save_txt ~out:(file "us") us


let save_prms suffix prms = Misc.save_bin (Printf.sprintf "%s/prms_%s" dir suffix) prms
let save_task suffix task = Misc.save_bin (Printf.sprintf "%s/prms_%s" dir suffix) task

module U = Priors.Gaussian

module D0 = Dynamics.Arm_Plus (struct
  let phi_x x = x
  let d_phi_x _ = AD.Mat.eye 200
  let phi_u x = AD.Maths.relu x
  let d_phi_u _ = AD.Mat.ones 1 200
end)

module L0 = Likelihoods.End (struct
  let label = "output"
end)

let prms =
  let open Owl_parameters in
  let likelihood =
    Likelihoods.End_P.
      { c =
          (pinned : setter) (AD.pack_arr (Mat.load_txt (Printf.sprintf "%s/c" data_dir)))
      ; c_mask = None
      ; qs_coeff = (pinned : setter) (AD.F 1.)
      ; t_coeff = (pinned : setter) (AD.F 0.5)
      ; g_coeff = (pinned : setter) (AD.F 1.)
      }
  in
  let dynamics =
    Dynamics.Arm_Plus_P.
      { a =
          (pinned : setter)
            (AD.pack_arr
               (Mat.transpose Mat.(load_txt (Printf.sprintf "%s/w_rec" data_dir) - eye m)))
      ; b = (pinned : setter) (AD.pack_arr b)
      ; c =
          (pinned : setter) (AD.pack_arr (Mat.load_txt (Printf.sprintf "%s/c" data_dir)))
      }
  in
  let prior =
    Priors.Gaussian_P.
      { lambda_prep = (pinned : setter) (AD.F lambda_prep)
      ; lambda_mov = (pinned : setter) (AD.F lambda_mov)
      }
  in
  Model.Generative_P.{ prior; dynamics; likelihood }


module I0 = Model.ILQR (U) (D0) (L0)

let xs0, us0 =
  I0.solve
    ~u_init:(Mat.sqr (Mat.gaussian ~sigma:0.01 901 m))
    ~arm:true
    ~n:_n
    ~m
    ~prms
    task_noprep


let _ =
  Mat.save_txt ~out:"log_barrier/eye_us_0" (AD.unpack_arr us0);
  Mat.save_txt ~out:"log_barrier/eye_xs_0" (AD.unpack_arr xs0)


let xs_300, us_300 =
  I0.solve
    ~u_init:(Mat.sqr (Mat.gaussian ~sigma:0.01 901 m))
    ~arm:true
    ~n:_n
    ~m
    ~prms
    task_prep


let _ =
  Mat.save_txt ~out:"log_barrier/eye_us_300" (AD.unpack_arr us_300);
  Mat.save_txt ~out:"log_barrier/eye_xs_300" (AD.unpack_arr xs_300)
