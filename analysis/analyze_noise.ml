open Owl
open Lib
open Base
module AD = Algodiff.D 

(*want to test for all noise levels the robustness to it
We choose a value epsilon and for each time step use gaussian noise of standard deviation  \epsilon u*)

let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir]")
let data_dir = Cmdargs.(get_string "-data" |> force ~usage:"-data [data_dir]")
let in_dir = Printf.sprintf "%s/%s" dir
let c = C.broadcast' (fun ()-> Mat.load_txt (Printf.sprintf "%s/c" data_dir))
let w = C.broadcast' (fun ()-> Mat.load_txt (Printf.sprintf "%s/w_rec" data_dir))
let a = Mat.((w - eye (Mat.row_num w))/$(150E-3)) |> Mat.transpose 
let targets = 
  C.broadcast' (fun () -> let x = Mat.load_txt (Printf.sprintf "%s/target_thetas" data_dir) in Mat.get_slice [[];[0;1]] x)
let n = 204
let m = 200
let noise_level = 0.00005

let epsilons = [|20.;15.;10.;8.;2.;1.;0.1|]
let t_preps = [| 0.; 0.01; 0.02; 0.05; 0.2; 0.3; 0.5; 0.6; 0.8; 1. |]

let task_prms = Array.init (8*(Array.length t_preps)) ~f:(fun i -> let t_idx = i/8 in let n_tgt = Int.rem i 8 in n_tgt,t_preps.(t_idx))
module U = Priors.Gaussian
module D = Dynamics.Arm_Linear

module L = Likelihoods.End (struct
  let label = "output"
end)

let get_accuracy n_target traj= 
let tgt = Mat.row targets n_target in 
let end_pt = Mat.get_slice [[-190;-1];[0;1]] traj
in Mat.((l2norm_sqr ~axis:0 Mat.(end_pt - tgt))/$190.)


let run_noisy ~i ~t_prep ~epsilon =
  let us =
    Mat.load_txt (in_dir (Printf.sprintf "us_%i_%i_%.6f" i Float.(to_int (1000. *. t_prep)) noise_level))
  in
  let n_steps = Mat.row_num us in
  let prms =
    Misc.read_bin
      (in_dir (Printf.sprintf "prms" ))
  in
  let task =
    Misc.read_bin
      (in_dir (Printf.sprintf "task_%i_%i_%.6f" i Float.(to_int (1000. *. t_prep)) noise_level))
  in
  let module I = Model.ILQR (U) (D) (L) in
  let pert_us = Mat.(us + (epsilon $* gaussian n_steps m * us)) in
  let xs, us = I.run ~ustars:pert_us ~n ~m ~prms task in
  let thetas = Mat.get_slice [ []; [ 0; 3 ] ] (AD.unpack_arr xs) in
  let xs =  Mat.get_slice [ []; [ 4;-1 ] ] (AD.unpack_arr xs) in 
  let us = AD.unpack_arr us in 
  let accuracy = get_accuracy i thetas in 
  Mat.save_txt ~append:true ~out:(in_dir (Printf.sprintf "added_noise/accuracy_%i_%i_%.6f_%.3f" i Float.(to_int (1000. *. t_prep)) noise_level epsilon)) accuracy;
  Mat.save_txt
    ~out:(in_dir (Printf.sprintf "added_noise/noisy_thetas_%i_%i_%.6f_%.3f" i Float.(to_int (1000. *. t_prep)) noise_level epsilon))
    thetas;
    Mat.save_txt
    ~out:(in_dir (Printf.sprintf "added_noise/noisy_xs_%i_%i_%.6f_%.3f" i Float.(to_int (1000. *. t_prep)) noise_level epsilon))
    xs
;  Mat.save_txt
~out:(in_dir (Printf.sprintf "added_noise/noisy_torques_%i_%i_%.6f_%.3f" i Float.(to_int (1000. *. t_prep)) noise_level epsilon))
Mat.(xs*@(transpose c));
Mat.save_txt
~out:(in_dir (Printf.sprintf "added_noise/rec_us_%i_%i_%.6f_%.3f" i Float.(to_int (1000. *. t_prep)) noise_level epsilon))
Mat.(xs*@(a));
Mat.save_txt
~out:(in_dir (Printf.sprintf "added_noise/noisy_us_%i_%i_%.6f_%.3f" i Float.(to_int (1000. *. t_prep)) noise_level epsilon))
us




let _ = Array.init 50 ~f:(fun _ -> Array.map epsilons ~f:(fun eps -> Array.iteri task_prms ~f:(fun j (i,x) -> 
  if Int.(j % C.n_nodes = C.rank) then run_noisy ~i ~t_prep:x ~epsilon:eps)))