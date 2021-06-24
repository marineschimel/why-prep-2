open Owl
module AD = Algodiff.D
module M = Arm.Make (Arm.Defaults)
open Lib

let _ = Printexc.record_backtrace true
let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")
let data_dir = Cmdargs.(get_string "-data" |> force ~usage:"-data [dir in which data is]")
let in_dir s = Printf.sprintf "%s/%s" dir s
let targets = Mat.load_txt (Printf.sprintf "%s/target_thetas" data_dir)
let target i = Mat.row targets i
let t_prep = 0.3
let dt = 1E-3
let lambda_prep = 1E-3
let lambda_mov = 1E-3
let n_out = 2
let n = 24
let m = 20
let n_output = 2
let theta0 = Mat.of_arrays [| [| 0.174533; 2.50532; 0.; 0. |] |] |> AD.pack_arr

let task =
  Model.
    { t_prep; t_mov = 0.4; dt; t_hold = None; target = AD.pack_arr (target 0); theta0 }


module U = Priors.Gaussian
module D = Dynamics.Arm_Linear

module L = Likelihoods.End (struct
  let label = "output"
end)

let prms =
  let open Owl_parameters in
  let likelihood =
    Likelihoods.End_P.
      { c = (pinned : setter) (AD.Mat.ones n_out n)
      ; c_mask =
          Some
            (AD.Maths.concatenate ~axis:1 [| AD.Mat.zeros n_out 4; AD.Mat.ones n_out m |])
      ; qs_coeff = (pinned : setter) (AD.F 1.)
      ; t_coeff = (pinned : setter) (AD.F 1.)
      }
  in
  let dynamics =
    Dynamics.Arm_Linear_P.
      { a = (pinned : setter) (AD.Mat.ones m m)
      ; b = (pinned : setter) (AD.Mat.eye m)
      ; c = (pinned : setter) (AD.Mat.ones n_out m)
      }
  in
  let prior =
    Priors.Gaussian_P.
      { lambda_prep = (pinned : setter) (AD.F lambda_prep)
      ; lambda_mov = (pinned : setter) (AD.F lambda_mov)
      }
  in
  Model.Generative_P.{ prior; dynamics; likelihood }


module I = Model.ILQR (U) (D) (L)

let _, us = I.solve ~n ~m ~prms task
