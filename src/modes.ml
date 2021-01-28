open Owl
open Lib
open Defaults

let _ = Printexc.record_backtrace true
let t_prep = Cmdargs.(get_int "-t" |> force ~usage:"-prep_time")

  let obs_gramian a c = Linalg.D.lyapunov Mat.(transpose a) Mat.(neg (transpose c) *@ c)
  let ctrl_gramian a b = Linalg.D.lyapunov a Mat.(neg b *@ transpose b)
  let seed = 1
  
  let n_reaches = 4
  let in_dir s = Printf.sprintf "reach_1/%s" s
  let size_prep = t_prep - 50

  let size_mov = 400
  let duration = t_prep + 600
  let saving_dir s = Printf.sprintf "modes/%s" s

  let softmax x =
    let m = Mat.(std ~axis:0 x +$ 0.5) in
    Mat.(x / m)


let preprocessed_data =
  let load i =
    let m = Mat.load_txt (Printf.(sprintf "reach_1/reach_%i/traj_%i" i t_prep)) in
    let ma = Mat.get_fancy [ R [ 0; duration-1 ]; R [ 4; -1 ] ] m in
    Arr.reshape ma [| duration; size_net; 1 |]
  in
  let data = Array.init n_reaches (fun i -> softmax (load (succ i))) in
  let mean_data = Arr.concatenate ~axis:2 data |> Arr.mean ~axis:2 in
  let dat =
    let a = Arr.concatenate ~axis:0 (Array.map (fun x -> Arr.(x - mean_data)) data) in
    Arr.reshape a [| (Arr.shape a).(0); (Arr.shape a).(1) |]
  in
  mean_data, Mat.of_arrays (Arr.to_arrays dat)



let c = Mat.load_txt (in_dir "c")

let w = Mat.load_txt (in_dir "w")

let a = Mat.((w - eye 200)/$tau)

let x_prep, x_mov, mean_data, x_all,x_allall =
  let mu, dat,all =
    let load i =
      Mat.(
        get_slice [ []; [ 4; -1 ] ] (load_txt (in_dir Printf.(sprintf "reach_%i/traj_%i" ( i) t_prep))))
      (*in
        Mat.get_fancy [ R [ 0; duration - 1 ]; R [ 0; -1 ] ] m*)
    in
    let data =
      if n_reaches > 1
      then Array.init n_reaches (fun i -> load (succ i)) |> stack ~axis:0
      else load 1
    in
    let mean_data = Arr.mean ~axis:0 data in
    let dat = Arr.(data - mean_data) in
    let all = Array.init n_reaches (fun i -> load (succ i)) |> Mat.concatenate ~axis:0 in 
    Arr.squeeze ~axis:[| 0 |] mean_data, dat, all
  in
  let x_prep =
    if n_reaches > 1
    then (
      let m = Arr.get_slice [ []; [ 200; 200 + size_prep - 1 ]; [] ] dat in
      Arr.reshape m [| -1; Arr.(shape m).(2) |])
    else Mat.get_slice [ [ 200; 200 + size_prep - 1 ]; [] ] dat
  in
  let x_mov =
    if n_reaches > 1
    then (
      let m = Arr.get_slice [ []; [ t_prep - 50; t_prep + size_mov - 50 - 1 ]; [] ] dat in
      Arr.reshape m [| -1; Arr.(shape m).(2) |])
    else Mat.get_slice [ [ t_prep - 50; t_prep + size_mov - 50 - 1 ]; [] ] dat
  in
  x_prep, x_mov, mu, dat,all


let inputs = Mat.load_txt (in_dir Printf.(sprintf "reach_1/results_us_%i" t_prep))
let top_pc_mov, _, _ = Linalg.D.svd Mat.(transpose (x_mov - mean ~axis:0 x_mov))
let top_pc_prep, _, _ = Linalg.D.svd Mat.(transpose (x_prep - mean ~axis:0 x_prep))
let top_pc, _, _ =  Linalg.D.svd Mat.(transpose (x_allall - mean ~axis:0 x_allall))

let _ = Printf.printf "%i %i %!" (Mat.row_num top_pc) (Mat.col_num top_pc)
let matmul x y =
  if n_reaches > 1
  then (
    let shp = Arr.shape x in
    let _ = Printf.printf "shape : %i %i %i %!" shp.(0) shp.(1) shp.(2) in 
    let x = Arr.reshape x [| -1; shp.(2) |] in
    let z = Mat.(x *@ y) in
    shp.(2) <- -1;
    Arr.reshape z shp)
  else Mat.(x *@ y)


let normalized x = Mat.get_slice [ []; [ 0; 10 ] ] Mat.(x / l2norm ~axis:0 x)

let () =
  let obs = obs_gramian a c in
  let obs_i = obs_gramian a (Mat.eye 200) in
  let con = ctrl_gramian a (Mat.eye 200) in
  let eigv_obsi, _, _ = Linalg.D.svd obs_i in
  let eigv_obs, _, _ = Linalg.D.svd obs in
  let _ =
    Mat.save_txt ~out:(saving_dir "obs_modes") eigv_obs;
    Mat.save_txt ~out:(saving_dir "pc_modes") top_pc
  in
  let eigv_con, _, _ = Linalg.D.svd con in
  let ac_obs_proj = matmul x_all eigv_obs in
  let ac_obsi_proj = matmul x_all eigv_obsi in
  let ac_ctrl_proj = matmul x_all eigv_con in
  let ac_pc_proj = matmul x_all top_pc in
  (* let inp_obs_proj = matmul inputs eigv_obs in
  let inp_ctrl_proj = matmul inputs eigv_con in
  let inp_pc_proj = matmul inputs top_pc in
  let inp_obsi_proj = matmul inputs eigv_obsi in *)
  let pc_obs, pc_ctrl=
  (* pc_prep_obs, pc_prep_ctrl, pc_mov_obs, pc_mov_ctrl *)
    ( Mat.(sqr (transpose (normalized top_pc) *@ normalized eigv_obs))
    , Mat.(sqr (transpose (normalized top_pc) *@ normalized eigv_con)))
    (* , Mat.(sqr (transpose (normalized top_pc_prep) *@ normalized eigv_obs))
    , Mat.(sqr (transpose (normalized top_pc_prep) *@ normalized eigv_con))
    , Mat.(sqr (transpose (normalized top_pc_mov) *@ normalized eigv_obs))
    , Mat.(sqr (transpose (normalized top_pc_mov) *@ normalized eigv_con)) ) *)
  in
  let _ =
    ( Mat.save_txt ~out:(saving_dir "pc_obs") pc_obs
    , Mat.save_txt ~out:(saving_dir "pc_ctrl") pc_ctrl)
  in
  (let get_dat i x =
     if n_reaches > 1 then Arr.(get_slice [ [ i ] ] x |> Arr.squeeze ~axis:[| 0 |]) else x
   in
   for i = 0 to n_reaches - 1 do
     Mat.save_txt
       ~out:(saving_dir (Printf.sprintf "ac_obs_%i" (succ i)))
       (get_dat i ac_obs_proj);
     Mat.save_txt
       ~out:(saving_dir (Printf.sprintf "ac_obsi_%i" (succ i)))
       (get_dat i ac_obsi_proj);
     Mat.save_txt
       ~out:(saving_dir (Printf.sprintf "ac_ctrl_%i" (succ i)))
       (get_dat i ac_ctrl_proj);
     Mat.save_txt
       ~out:(saving_dir (Printf.sprintf "ac_pc_%i_%i" i (t_prep)))
       (get_dat i ac_pc_proj);
     (* Mat.save_txt
       ~out:(saving_dir (Printf.sprintf "inp_obs_%i" (succ i)))
       (get_dat i inp_obs_proj);
     Mat.save_txt
       ~out:(saving_dir (Printf.sprintf "inp_obsi_%i" (succ i)))
       (get_dat i inp_obsi_proj);
     Mat.save_txt
       ~out:(saving_dir (Printf.sprintf "inp_ctrl_%i" (succ i)))
       (get_dat i inp_ctrl_proj);
     Mat.save_txt
       ~out:(saving_dir (Printf.sprintf "in_pc_%i" (succ i)))
       (get_dat i inp_pc_proj) *)
   done);
  let get_var x =
    if n_reaches > 1
    then Arr.var ~axis:0 x |> Arr.squeeze ~axis:[| 0 |]
    else Mat.var x ~axis:0
  in
  get_var ac_obs_proj |> Mat.save_txt ~out:(saving_dir (Printf.sprintf "ac_obs_proj_var"));
  get_var ac_obsi_proj
  |> Mat.save_txt ~out:(saving_dir (Printf.sprintf "ac_obsi_proj_var"));
  get_var ac_ctrl_proj
  |> Mat.save_txt ~out:(saving_dir (Printf.sprintf "ac_ctrl_proj_var"));
  get_var ac_pc_proj |> Mat.save_txt ~out:(saving_dir (Printf.sprintf "ac_pc_proj_var"));
  (* get_var inp_obs_proj
  |> Mat.save_txt ~out:(saving_dir (Printf.sprintf "inp_obs_proj_var"));
  get_var inp_obsi_proj
  |> Mat.save_txt ~out:(saving_dir (Printf.sprintf "inp_obsi_proj_var"));
  get_var inp_ctrl_proj
  |> Mat.save_txt ~out:(saving_dir (Printf.sprintf "inp_ctrl_proj_var"));
  get_var inp_pc_proj |> Mat.save_txt ~out:(saving_dir (Printf.sprintf "inp_pc_proj_var")); *)
  Mat.(c *@ eigv_con)
  |> Mat.transpose
  |> Mat.l2norm ~axis:1
  |> Mat.save_txt ~out:(saving_dir "test_read");
  Mat.(transpose (get_slice [ []; [ 0 ] ] top_pc_mov) *@ eigv_con)
  |> Mat.transpose
  |> Mat.l2norm ~axis:1
  |> Mat.save_txt ~out:(saving_dir "test_proxi");
  eigv_con |> Mat.save_txt ~out:(saving_dir "ctrl_eig")
