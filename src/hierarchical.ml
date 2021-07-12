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
let lambda_prep = 1E-6
let lambda_mov = 1E-6
let n_out = 2
let n = 204
let m = 200
let tau = 150E-3
let n_output = 2
let n_targets = 8
let theta0 = Mat.of_arrays [| [| 0.174533; 2.50532; 0.; 0. |] |] |> AD.pack_arr
let t_prep = 0.3
let target_num = 0

(* let c = Mat.gaussian ~sigma:0.1 n_out m
let _ = Mat.save_txt ~out:(in_data_dir "c") c *)

(*structure : get inputs from task, then from the same SOC just generating inputs, 
then from yet another SOC
Maybe do the same thing using 200 - 100 - 50 
what about one population receiving modulated inputs about the target? 
*)

let task =
  Model.
    { t_prep = 0.3
    ; t_mov = 0.4
    ; dt
    ; t_hold = Some 0.2
    ; scale_lambda = None
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
module D0 = Dynamics.Arm_Linear
module D1 = Dynamics.Linear

module L0 = Likelihoods.End (struct
  let label = "output"
end)

module L1 = Likelihoods.Match (struct
  let label = "output"
end)

let prms0 =
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
    Dynamics.Arm_Linear_P.
      { a =
          (pinned : setter)
            (AD.pack_arr
               (Mat.transpose Mat.(load_txt (Printf.sprintf "%s/w_rec" data_dir) - eye m)))
      ; b = (pinned : setter) (AD.Mat.eye m)
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


let prms1 =
  let open Owl_parameters in
  let likelihood = Likelihoods.Match_P.{ q_coeff = (pinned : setter) (AD.F 1.) } in
  let dynamics =
    Dynamics.Linear_P.
      { a =
          (pinned : setter)
            (AD.pack_arr
               (Mat.transpose Mat.(load_txt (Printf.sprintf "%s/w_rec" data_dir) - eye m)))
      ; b = (pinned : setter) (AD.Mat.eye m)
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
module I1 = Model.ILQR (U) (D1) (L1)

let xs0, us0 = I0.solve ~arm:true ~n ~m ~prms:prms0 task

let _ =
  Mat.save_txt ~out:(Printf.sprintf "us_0") (AD.unpack_arr us0);
  Mat.save_txt ~out:(Printf.sprintf "xs_0") (AD.unpack_arr xs0)


let n = 200

let _ =
  let rec hier k us =
    let task_next = Model.{ task with target = us } in
    let xs, next_us = I1.solve ~arm:false ~n ~m ~prms:prms1 task_next in
    let _ = Mat.save_txt ~out:(Printf.sprintf "us_%i" (k + 1)) (AD.unpack_arr next_us) in
    let _ = Mat.save_txt ~out:(Printf.sprintf "xs_%i" (k + 1)) (AD.unpack_arr xs) in
    if k = 10 then next_us else hier (k + 1) next_us
  in
  hier 0 us0
