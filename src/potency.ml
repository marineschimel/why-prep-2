open Owl
open Base

(*get obs and controllabiloty gramians*)
(*get the modes, as well as their eigenvalues*)
(*compute the mean energy across movements as a function of time *)
let n = 200

let dir = "/home/mmcs3/rds/rds-t2-cs156-T7o4pEA8QoU/mmcs3/linear_lambda_0.000001/ramping_soc/seed_0_mov"

let in_dir s = Printf.sprintf "%s/%s" dir s 

let eigenvalues m =
  let v = Linalg.D.eigvals m in
  let re = Dense.Matrix.Z.re v in
  let im = Dense.Matrix.Z.im v in
  Mat.(concat_horizontal (transpose re) (transpose im))


let w = Mat.load_txt (in_dir "w")

let eigs_w = eigenvalues w

let _ = Mat.save_txt ~out:(in_dir "eigs_w") eigs_w
let c =  Mat.load_txt (in_dir "c")

(*normalize c*)
let c = Mat.(c/$ Owl.Linalg.D.norm c)

let a = Mat.(w - eye n)

let disc_a = Mat.(zeros n n - eye n)

let obs_gramian a c =
  Linalg.D.lyapunov Mat.(transpose a) Mat.(neg (transpose c) *@ c)

let ctr_gramian a b = Linalg.D.lyapunov a Mat.(neg b *@ transpose b)

let obs_id = obs_gramian a (Mat.eye 200)

let obs_disc = obs_gramian disc_a Mat.(neg (transpose c) *@ c)

let ctr_disc = ctr_gramian disc_a Mat.(eye n)

let tr_p_disc = Mat.trace ctr_disc
let eig_disc_o, eig_disc_q = eigenvalues (obs_disc), eigenvalues (ctr_disc)

let _ = Mat.save_txt (in_dir "obs_eig_disc") eig_disc_o;
Mat.save_txt (in_dir "ctr_eig_disc") eig_disc_q 
(*these are all 0.5, so we multiply P by 2 to havethese be at 1*)


let obs_c = obs_gramian a c

let ctr = ctr_gramian a (Mat.eye 200)

let ctr = Mat.(ctr /$ Float.(tr_p_disc /.200.))
  (* obs_gramian a c *)


let obs_c_modes,obs_c_eigs = let obs_c_modes, _,_ = Linalg.D.svd obs_c
in obs_c_modes, (eigenvalues obs_c)

let obs_id_modes,obs_id_eigs = let obs_id_modes, _,_ = Linalg.D.svd obs_id
in obs_id_modes, (eigenvalues obs_id)

let ctr_modes,ctr_eigs = let ctr_modes, _,_ = Linalg.D.svd ctr
in Mat.(get_slice [[]; []] ctr_modes), (eigenvalues ctr)

let traj i = Mat.load_txt (in_dir (Printf.sprintf "rates_%i_300" i))

let n_reaches = 7

let mean_trajs = Array.init n_reaches ~f:(fun i -> let i = succ i in let t = traj i in Arr.reshape t [|1;-1;200|]) |> fun z -> Arr.concatenate ~axis:0 z |> fun e -> Arr.mean ~axis:0 e |> fun e -> Mat.reshape e [|-1;200|]

let obs_id_proj = 
  Array.init n_reaches ~f:(fun i -> let t = Mat.((traj i) - mean_trajs) in Mat.(t*@obs_id_modes) |> fun z -> Arr.reshape z [|1;-1;200|]) |> fun z -> Arr.concatenate ~axis:0 z |> fun e -> Arr.std ~axis:0 e |> fun e -> Mat.reshape e [|-1;200|]

let obs_c_proj = 
    Array.init n_reaches ~f:(fun i -> let t = Mat.((traj i) - mean_trajs) in Mat.(t*@obs_c_modes) |> fun z -> Arr.reshape z [|1;-1;200|]) |> fun z -> Arr.concatenate ~axis:0 z |> fun e -> Arr.std ~axis:0 e |> fun e -> Mat.reshape e [|-1;200|]

let ctrl_proj = 
  Array.init n_reaches ~f:(fun i ->  let t = Mat.((traj i) - mean_trajs) in Mat.(t*@ctr_modes) |> fun z -> Arr.reshape z [|1;-1;200|]) |> fun z -> Arr.concatenate ~axis:0 z |> fun e -> Arr.std ~axis:0 e |> fun e -> Mat.reshape e [|-1;200|]

let _ = Mat.save_txt (in_dir "obs_id_proj") obs_id_proj;
Mat.save_txt (in_dir "obs_c_proj") obs_c_proj;
Mat.save_txt (in_dir "ctrl_proj") ctrl_proj;
Mat.save_txt (in_dir "obs_id_eigs") obs_id_eigs;
Mat.save_txt (in_dir "obs_c_eigs") obs_c_eigs;
Mat.save_txt (in_dir "ctr_eigs") ctr_eigs