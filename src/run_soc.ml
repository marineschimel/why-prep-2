open Owl
module AD = Algodiff.D
module M = Arm.Make (Arm.Defaults)
open Lib
open Base

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
let n = 204
let m = 200
let tau = 150E-3
let n_output = 2
let n_targets = 8
let theta0 = Mat.of_arrays [| [| 0.174533; 2.50532; 0.; 0. |] |] |> AD.pack_arr
let t_preps = [| 0.; 0.01; 0.02; 0.05; 0.1; 0.15; 0.2; 0.3; 0.4; 0.5; 0.6; 0.8; 1. |]

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
        ; target = AD.pack_arr (target n_target)
        ; theta0
        ; tau = 150E-3
        })


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


module U = Priors.Gaussian
module D = Dynamics.Arm_Linear

module L = Likelihoods.End (struct
  let label = "output"
end)

let prms =
  let open Owl_parameters in
  let likelihood =
    Likelihoods.End_P.
      { c =
          (pinned : setter)
            (AD.Maths.concatenate
               ~axis:1
               [| AD.Mat.zeros n_out 4
                ; AD.pack_arr (Mat.load_txt (Printf.sprintf "%s/c" data_dir))
               |])
      ; c_mask =
          Some
            (AD.Maths.concatenate ~axis:1 [| AD.Mat.zeros n_out 4; AD.Mat.ones n_out m |])
      ; qs_coeff = (pinned : setter) (AD.F 1.)
      ; t_coeff = (pinned : setter) (AD.F 1.)
      }
  in
  let dynamics =
    Dynamics.Arm_Linear_P.
      { a =
          (pinned : setter)
            (AD.pack_arr (Mat.load_txt (Printf.sprintf "%s/w_rec" data_dir)))
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


module I = Model.ILQR (U) (D) (L)

let _ =
  Array.mapi tasks ~f:(fun i t ->
      if Int.(i % C.n_nodes = C.rank)
      then (
        let n_target = Int.rem i n_targets in
        let t_prep = Float.to_int (1000. *. t.t_prep) in
        let xs, us = I.solve ~n ~m ~prms t in
        save_results (Printf.sprintf "%i_%i" n_target t_prep) xs us))
