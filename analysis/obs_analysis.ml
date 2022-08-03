open Owl
open Lib
open Defaults

(*let tests_reshape = let w = Mat.(load_txt "w") in let w_in = Mat.reshape  w [|size_net;size_net|]
in let flat_w = Mat.reshape w_in [|1;size_net*size_net|] in let reshaped_w  = Mat.reshape flat_w [|size_net;size_net|]
in Mat.print (Mat.(reshaped_w - w))*)

let nn_evol =
  Array.init 26 (fun i ->
      nonnormality (Mat.load_txt (Printf.sprintf "from_soc/w_new_%i" (succ i))))
  |> fun z -> Mat.of_array z 1 (-1)

let _ = Mat.save_txt ~out:"nn_evol" nn_evol
let w = Mat.(load_txt "w")

let w_1, w_2, w_72 =
  Mat.load_txt "soc_50/w", Mat.load_txt "soc_50/w_new_2", Mat.load_txt "soc_50/w_new_72"

let _ = Mat.save_txt ~out:"dw_small" Mat.(w_72 - w_1)

let _ =
  let module Z = Dense.Matrix.Z in
  let _, evals = Linalg.D.eig w_72 in
  let re, im = Z.re evals, Z.im evals in
  Mat.save_txt ~out:"eigs_w_small_72" Mat.(transpose (re @= im))

let _ = Printf.printf "%f %f %f" (nonnormality w) (nonnormality w_1) (nonnormality w_2)

(*let _ = Printf.printf "%f %f" (Mat.l2norm' w) (Mat.l2norm' w_3)*)

let _ = Printexc.record_backtrace true
let obs_gramian a c = Linalg.D.lyapunov Mat.(transpose a) Mat.(neg (transpose c) *@ c)
let ctrl_gramian a b = Linalg.D.lyapunov Mat.(transpose a) Mat.(neg (transpose b) *@ b)
let w = Mat.load_txt "soc_optc/size_200/w"

let traj, inputs =
  let x = Mat.load_txt "soc_optc/size_200/taus_ilqr" in
  Mat.get_slice [ []; [ 4; 203 ] ] x, Mat.get_slice [ []; [ 204; -1 ] ] x

let a = Mat.((w - eye n) /$ tau)

let proj_ac c =
  let og = obs_gramian a c in
  let modes, eigs, _ = Linalg.D.svd og in
  Mat.(traj *@ modes), Mat.(inputs *@ modes), Mat.sqr eigs |> Mat.transpose

let proj_nonnull c =
  let og = obs_gramian a c in
  let modes, _, _ = Linalg.D.svd og in
  Mat.(c *@ modes) |> Mat.transpose

let proj_nonnull_ctrl c =
  let cg = ctrl_gramian a (Mat.eye n) in
  let modes, _, _ = Linalg.D.svd cg in
  Mat.(c *@ modes) |> Mat.transpose

let proj_wc w c =
  let modes_w, eigs_w, _ = Linalg.D.svd w in
  let proj_modes_c = Mat.(c *@ modes_w) |> Mat.transpose in
  proj_modes_c, Mat.transpose eigs_w

let in_dir s = Printf.sprintf "projs_dprep/%s" s

let _ =
  Array.init 6 (fun i ->
      let c = Mat.load_txt (Printf.sprintf "soc_optc/size_200/c_new_%i" i) in
      let ac, inpts, eigs = proj_ac c in
      let read = proj_nonnull c in
      let read_ctrl = proj_nonnull_ctrl c in
      let proj_wc, eigs_w = proj_wc w c in
      ( Mat.save_txt ~out:(in_dir (Printf.sprintf "proj_wc_%i" i)) proj_wc
      , Mat.save_txt ~out:(in_dir (Printf.sprintf "eigs_w_%i" i)) eigs_w
      , Mat.save_txt ~out:(in_dir (Printf.sprintf "ac_%i" i)) ac
      , Mat.save_txt ~out:(in_dir (Printf.sprintf "inpts_%i" i)) inpts
      , Mat.save_txt ~out:(in_dir (Printf.sprintf "read_%i" i)) read
      , Mat.save_txt ~out:(in_dir (Printf.sprintf "eigs_obs_%i" i)) eigs
      , Mat.save_txt ~out:(in_dir (Printf.sprintf "read_ctrl_%i" i)) read_ctrl ))
