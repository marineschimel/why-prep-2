open Owl
open Lib
open Defaults

let _ = Printexc.record_backtrace true
let times = [| 0; 100; 300; 500 |]
let n_reaches = 8
let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")
let t_prep = Cmdargs.(get_int "-prep" |> force ~usage:"-prep")
let in_dir s = Printf.sprintf "%s/%s" dir s
let softmax x = x

(* let m = Mat.(std ~axis:0 x +$ 0.01) in
    Mat.(x / m) *)
let preprocessed_data =
  let load i =
    let m = Mat.load_txt (in_dir Printf.(sprintf "rates_%i_%i" i t_prep)) in
    let ma = Mat.get_fancy [ R [ 0; duration - 1 ]; R [ 0; -1 ] ] m in
    Arr.reshape ma [| duration; -1; 1 |]
  in
  let data = Array.init n_reaches (fun i -> softmax (load i)) in
  let mean_data = Arr.concatenate ~axis:2 data |> Arr.mean ~axis:2 in
  let dat =
    let a = Arr.concatenate ~axis:0 (Array.map (fun x -> Arr.(x - mean_data)) data) in
    Arr.reshape a [| (Arr.shape a).(0); (Arr.shape a).(1) |]
  in
  mean_data, Mat.of_arrays (Arr.to_arrays dat)


let data_prep =
  let m = Arr.reshape (snd preprocessed_data) [| n_reaches; duration; -1 |] in
  let f i =
    Arr.reshape
      (Arr.get_slice [ [ i ]; [ n_prep - size_prep; n_prep - 1 ]; [] ] m)
      [| size_prep; -1 |]
  in
  Array.init n_reaches f |> Arr.concatenate ~axis:0


let data_mov =
  let m = Arr.reshape (snd preprocessed_data) [| n_reaches; duration; -1 |] in
  let f i =
    Arr.reshape
      (Arr.get_slice [ [ i ]; [ n_prep + 25; n_prep + size_mov - 1 + 25 ]; [] ] m)
      [| size_mov; -1 |]
  in
  Array.init n_reaches f |> Arr.concatenate ~axis:0


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
