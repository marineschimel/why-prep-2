open Owl
open Lib
open Defaults
open Base

(* /home/mmcs3/rds/hpc-work/_results/why_prep/results/uniform_1E-6/soc/seed_9 *)
let _ = Printexc.record_backtrace true
let times = [| 0; 50; 100; 300; 500; 800 |]

(* 500 ; 600 |] *)
let n_reaches = 8
let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")
let t_prep = Cmdargs.(get_int "-prep" |> force ~usage:"-prep")
let in_dir s = Printf.sprintf "%s/%s" dir s

let softmax x =
  let m = Mat.(max ~axis:0 x - min ~axis:0 x +$ (0.2 *. max' x)) in
  Mat.(x / m)

let n_prep = t_prep / 2
let size_prep = 40
let size_mov = 100

let duration t =
  let rates0 = Mat.load_txt (in_dir Printf.(sprintf "rates_0_%i" t)) in
  Mat.row_num rates0

let preprocessed_data t =
  let load i =
    let m = Mat.load_txt (in_dir Printf.(sprintf "rates_%i_%i" i t)) in
    let ma = Mat.get_fancy [ R [ 0; -1 ]; R [ 0; -1 ] ] m in
    Arr.reshape ma [| Mat.row_num m; -1; 1 |]
  in
  let data = Array.init n_reaches (fun i -> load i) in
  let concat_data = Arr.concatenate ~axis:0 data in
  let norm_factor =
    Mat.(max ~axis:0 concat_data - min ~axis:0 concat_data +$ (0.1 *. max' concat_data))
    |> fun z -> Arr.reshape z [| 1; -1; 1 |]
  in
  let norm_data = Array.map ~f:(fun x -> Arr.(x / norm_factor)) data in
  let mean_data = Arr.concatenate ~axis:2 norm_data |> Arr.mean ~axis:2 in
  let dat =
    let a =
      Arr.concatenate ~axis:0 (Array.map ~f:(fun x -> Arr.(x - mean_data)) norm_data)
    in
    Arr.reshape a [| (Arr.shape a).(0); (Arr.shape a).(1) |]
  in
  mean_data
  |> fun a ->
  ( Arr.reshape a [| (Arr.shape a).(0); (Arr.shape a).(1) |]
  , Mat.of_arrays (Arr.to_arrays dat) )

let data_prep =
  let m =
    Arr.reshape (snd (preprocessed_data t_prep)) [| n_reaches; duration t_prep; -1 |]
  in
  let f i =
    Arr.reshape
      (Arr.get_slice [ [ i ]; [ n_prep - size_prep; n_prep - 1 ]; [] ] m)
      [| size_prep; -1 |]
  in
  Array.init n_reaches f |> Arr.concatenate ~axis:0

let data_mov =
  let m =
    Arr.reshape (snd (preprocessed_data t_prep)) [| n_reaches; duration t_prep; -1 |]
  in
  let f i =
    Arr.reshape
      (Arr.get_slice [ [ i ]; [ n_prep + 25; n_prep + size_mov - 1 + 25 ]; [] ] m)
      [| size_mov; -1 |]
  in
  Array.init n_reaches f |> Arr.concatenate ~axis:0

let data t =
  let m = Arr.reshape (snd (preprocessed_data t)) [| n_reaches; duration t; -1 |] in
  let f i = Arr.reshape (Arr.get_slice [ [ i ]; []; [] ] m) [| -1; 200 |] in
  Array.init n_reaches f |> Arr.concatenate ~axis:0

let all_data = data t_prep
(* Array.map times (fun t -> data t) |> Arr.concatenate ~axis:0 *)

let top_pc, _, _ = Linalg.D.svd Mat.(transpose (all_data - mean ~axis:0 all_data))
let top_pc_mov, _, _ = Linalg.D.svd Mat.(transpose (data_mov - mean ~axis:0 data_mov))
let top_pc_prep, _, _ = Linalg.D.svd Mat.(transpose (data_prep - mean ~axis:0 data_prep))

let matmul x y =
  if n_reaches > 1
  then (
    let shp = Arr.shape x in
    let _ = Stdio.printf "shape : %i %i %i %!" shp.(0) shp.(1) shp.(2) in
    let x = Arr.reshape x [| -1; shp.(2) |] in
    let z = Mat.(x *@ y) in
    shp.(2) <- -1;
    Arr.reshape z shp)
  else Mat.(x *@ y)

let normalized x = Mat.get_slice [ []; [ 0; 10 ] ] Mat.(x / l2norm ~axis:0 x)

let projs reach =
  Array.map
    ~f:(fun t ->
      let m =
        let m = Arr.reshape (snd (preprocessed_data t)) [| n_reaches; duration t; -1 |] in
        Arr.reshape (Arr.get_slice [ [ reach ]; []; [] ] m) [| -1; 200 |]
      in
      let top_pc_mov, top_pc_prep =
        ( Mat.get_slice [ []; [ 0; 2 ] ] top_pc_mov
        , Mat.get_slice [ []; [ 0; 2 ] ] top_pc_prep )
      in
      Mat.(m *@ top_pc_mov)
      |> Mat.save_txt ~out:(in_dir (Printf.sprintf "pca/proj_reach_mov_%i_%i" reach t));
      Mat.(m *@ top_pc_prep)
      |> Mat.save_txt ~out:(in_dir (Printf.sprintf "pca/proj_reach_prep_%i_%i" reach t));
      let top_pc = Mat.get_slice [ []; [ 0; 5 ] ] top_pc in
      Mat.(m *@ top_pc)
      |> Mat.save_txt ~out:(in_dir (Printf.sprintf "pca/proj_reach_%i_%i" reach t)))
    times

let _ = Array.init n_reaches ~f:(fun i -> projs i)
