open Owl
open Lib

(*want to test for all noise levels the robustness to it
We choose a value epsilon and for each time step use gaussian noise of standard deviation  \epsilon u*)

let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir]")
let in_dir = Printf.sprintf "%s/%s" dir
let n = 204
let m = 200

module U = Priors.Gaussian
module D = Dynamics.Arm_Linear

module L = Likelihoods.End (struct
  let label = "output"
end)

let run_noisy ~i ~t_prep ~epsilon =
  let us =
    Mat.load_txt (in_dir (Printf.sprintf "us_%i_%i" i Float.(to_int (1000. *. t_prep))))
  in
  let n_steps = Mat.row_num us in
  let prms =
    Misc.read_bin
      (in_dir (Printf.sprintf "prms_%i_%i" i Float.(to_int (1000. *. t_prep))))
  in
  let task =
    Misc.read_bin
      (in_dir (Printf.sprintf "task_%i_%i" i Float.(to_int (1000. *. t_prep))))
  in
  let module I = Model.ILQR (U) (D) (L) in
  let pert_us = Mat.(us + epsilon $* gaussian n_steps n * us) in
  let xs, us = I.run ~ustars:pert_us ~n ~m ~prms ~task in
  let xs = AD.Maths.get_slice [ []; [ n; -1 ] ] new_taus in
  let thetas = Mat.get_slice [ []; [ 0; 3 ] ] (AD.unpack_arr xs) in
  Mat.save_txt
    ~out:(Printf.sprintf "noisy_thetas_%i_%i" i Float.(to_int (1000. *. t_prep)))
    thetas


let t_preps = [| 0.; 0.01; 0.02; 0.05; 0.2; 0.3; 0.5; 0.6; 0.8; 1. |]
