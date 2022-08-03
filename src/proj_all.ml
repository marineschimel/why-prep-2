open Owl
open Lib
open Defaults

let _ = Printexc.record_backtrace true
let times = [| 0; 100; 300; 500 |]
let n_reaches = 4
let in_dir s = Printf.sprintf "reach_1/projs/%s" s
let duration t = t + 400
let saving_dir s = Printf.sprintf "modes/%s" s
let softmax x = x
(* let m = Mat.(std ~axis:0 x +$ 0.01) in
    Mat.(x / m) *)

let data_mov =
  let load i t_prep =
    let m = Mat.load_txt Printf.(sprintf "reach_1/reach_%i/traj_%i" i t_prep) in
    let ma = Mat.get_fancy [ R [ t_prep; t_prep + 400 - 1 ]; R [ 4; -1 ] ] m in
    Arr.reshape ma [| 400; size_net; 1 |]
  in
  let get_data t_prep =
    let data = Array.init n_reaches (fun i -> softmax (load (succ i) t_prep)) in
    let mean_data = Arr.concatenate ~axis:2 data |> Arr.mean ~axis:2 in
    let dat =
      let a = Arr.concatenate ~axis:0 (Array.map (fun x -> Arr.(x - mean_data)) data) in
      Arr.reshape a [| (Arr.shape a).(0); (Arr.shape a).(1) |]
    in
    Mat.of_arrays (Arr.to_arrays dat)
  in
  Array.map (fun t_prep -> get_data t_prep) times |> Mat.concatenate ~axis:0

let data_prep =
  let load i t_prep =
    let m = Mat.load_txt Printf.(sprintf "reach_1/reach_%i/traj_%i" i t_prep) in
    let ma = Mat.get_fancy [ R [ 0; t_prep ]; R [ 4; -1 ] ] m in
    Arr.reshape ma [| t_prep + 1; size_net; 1 |]
  in
  let get_data t_prep =
    let data = Array.init n_reaches (fun i -> softmax (load (succ i) t_prep)) in
    let mean_data = Arr.concatenate ~axis:2 data |> Arr.mean ~axis:2 in
    let dat =
      let a = Arr.concatenate ~axis:0 (Array.map (fun x -> Arr.(x - mean_data)) data) in
      Arr.reshape a [| (Arr.shape a).(0); (Arr.shape a).(1) |]
    in
    Mat.of_arrays (Arr.to_arrays dat)
  in
  Array.map (fun t_prep -> get_data t_prep) times |> Mat.concatenate ~axis:0

let top_pc_mov, _, _ = Linalg.D.svd Mat.(transpose (data_mov - mean ~axis:0 data_mov))
let top_pc_prep, _, _ = Linalg.D.svd Mat.(transpose (data_prep - mean ~axis:0 data_prep))

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

let projs =
  Array.map
    (fun t ->
      let m =
        Mat.load_txt Printf.(sprintf "reach_1/reach_1/traj_%i" t)
        |> Mat.get_slice [ [ 0; duration t ]; [ 4; -1 ] ]
      in
      let top_pc_mov, top_pc_prep =
        Mat.get_slice [ []; [ 0; 1 ] ] top_pc_mov, Mat.get_slice [ []; [ 2 ] ] top_pc_prep
      in
      Mat.(m *@ (top_pc_mov @|| top_pc_prep))
      |> Mat.save_txt ~out:(in_dir (Printf.sprintf "proj_all_%i" t)))
    times
