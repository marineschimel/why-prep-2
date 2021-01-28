open Owl
open Lib
open Defaults


let _ = Printexc.record_backtrace true
let obs_gramian a c =
  Linalg.D.lyapunov Mat.(transpose a) Mat.(neg (transpose c) *@ c)

let ctrl_gramian a b =  Linalg.D.lyapunov Mat.(transpose a) Mat.(neg (transpose b) *@ b)

let w = Mat.load_txt "soc_optc/size_200/w"

let traj,inputs = let x =  Mat.load_txt "soc_optc/size_200/taus_ilqr"
in Mat.get_slice [[];[4;203]] x, Mat.get_slice [[];[204;-1]] x

let a = Mat.((w - eye n)/$tau)

let proj_ac c = let og = obs_gramian a c in let modes,eigs,_ = Linalg.D.svd og in Mat.(traj*@modes),  Mat.(inputs*@modes), Mat.sqr eigs |> Mat.transpose

let proj_nonnull c = let og = obs_gramian a c in let modes,_,_ = Linalg.D.svd og in Mat.(c*@modes) |> Mat.transpose 


let proj_nonnull_ctrl c = let cg = ctrl_gramian a (Mat.eye n) in let modes,_,_ = Linalg.D.svd cg in Mat.(c*@modes) |> Mat.transpose

let proj_wc w c= let modes_w,eigs_w,_ = Linalg.D.svd w in let 
proj_modes_c = Mat.(c*@(modes_w)) |>Mat.transpose in proj_modes_c,Mat.transpose eigs_w
let in_dir s= Printf.sprintf "projs_dprep/%s" s
let _ = Array.init 6 (fun i -> let c = 
  Mat.load_txt (Printf.sprintf "soc_optc/size_200/c_new_%i" i)
in let ac,inpts,eigs = proj_ac c in let read = proj_nonnull c
in let read_ctrl = proj_nonnull_ctrl c in 
let proj_wc, eigs_w = proj_wc w c in 
Mat.save_txt ~out:(in_dir (Printf.sprintf "proj_wc_%i" i)) proj_wc, 
Mat.save_txt ~out:(in_dir (Printf.sprintf "eigs_w_%i" i)) eigs_w,
Mat.save_txt ~out:(in_dir (Printf.sprintf "ac_%i" i)) ac, 
Mat.save_txt ~out:(in_dir (Printf.sprintf "inpts_%i" i)) inpts,
Mat.save_txt ~out:(in_dir (Printf.sprintf "read_%i" i)) read,
Mat.save_txt ~out:(in_dir (Printf.sprintf "eigs_obs_%i" i)) eigs,
Mat.save_txt ~out:(in_dir (Printf.sprintf "read_ctrl_%i" i)) read_ctrl)





