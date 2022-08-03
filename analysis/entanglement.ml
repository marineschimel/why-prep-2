open Owl
open Lib
open Defaults

let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")
let in_dir s = Printf.sprintf "%s/%s" dir s

(*Q(t) = max_t' ||(\dot{x}(t)-\dot{x}(t')||^2/(||(\dot{x}(t)-\dot{x}(t')||^2+\epsilon)*)

let epsilon =
  0.1
  *. (Mat.(get_slice [ []; [ 4; -1 ] ] (Mat.load_txt (in_dir "traj_400")))
     |> fun z -> Mat.(sum' (var ~axis:0 z)))

let dt = sampling_dt
let sampling_window = 1

let tangling t traj_1 traj_2 dt =
  let epsilon = 0.001 *. (traj_1 |> fun z -> Mat.(sum' (var ~axis:0 z))) in
  let get_row i y = Mat.get_slice [ [ i ]; [] ] y in
  let state_t = Mat.get_slice [ [ t ]; [] ] traj_1
  and dx_t =
    Mat.((get_slice [ [ t ]; [] ] traj_1 - get_slice [ [ pred t ]; [] ] traj_1) /$ dt)
  and dx_2 =
    Mat.((get_slice [ [ 1; -1 ]; [] ] traj_2 - get_slice [ [ 0; -2 ]; [] ] traj_2) /$ dt)
  in
  let entangled =
    Mat.mapi_rows
      (fun i r ->
        let dx' = get_row i dx_2 in
        let top = Mat.l2norm_sqr' Mat.(dx' - dx_t)
        and bottom = Mat.l2norm_sqr' Mat.(r - state_t) +. epsilon in
        top /. bottom)
      (Mat.get_slice [ [ 1; -1 ]; [] ] traj_2)
  in
  let max_idx = Utils.Array.max_i entangled in
  Mat.max' (Mat.of_array entangled 1 (-1)), max_idx

(* || Tests || *)

let activity ?reach:_reach time =
  let _reach =
    match _reach with
    | Some a -> a
    | None -> 1
  in
  let idces = List.init (400 / sampling_window) (fun i -> time + (i * sampling_window)) in
  let x =
    Mat.get_fancy
      [ L idces; R [ 4; -1 ] ]
      (Mat.load_txt (Printf.sprintf "%s/traj_%i" dir time))
  in
  Mat.(x /$ (max' x -. min' x)) |> fun z -> Mat.(z - mean ~axis:0 z)

let muscles ?reach:_reach time =
  let _reach =
    match _reach with
    | Some a -> a
    | None -> 1
  in
  let idces = List.init (400 / sampling_window) (fun i -> i * sampling_window) in
  let x =
    Mat.get_fancy
      [ L idces; R [ 0; 1 ] ]
      (Mat.load_txt (Printf.sprintf "%s/torques_%i" dir time))
  in
  Mat.(x /$ (Mat.max' x -. Mat.min' x)) |> fun z -> Mat.(z - mean ~axis:0 z)

let _ =
  let tang, i = tangling 1 (activity 600) (activity 600) dt in
  Printf.printf "%f %i \n %!" tang i

let pca_proj x dim =
  let cov = Mat.(transpose x *@ x) in
  let modes, _, _ = Linalg.D.svd cov in
  let top_modes = Mat.get_slice [ []; [ 0; dim ] ] modes in
  Mat.(x *@ top_modes)

let dim = 8

let get_all_x pc pm =
  Array.map
    (fun i ->
      let x =
        Mat.get_slice
          [ [ 400; -200 ]; [ 4; -1 ] ]
          (Mat.load_txt
             (Printf.sprintf "results/weighing_pm/w_%i_%i/reach_%i/traj_400" pc pm i))
      in
      let m = Mat.(max ~axis:0 x - min ~axis:0 x +$ 0.) in
      let x = Mat.(x / m) in
      let modes, _, _ = Linalg.D.svd (Mat.transpose Mat.(x - Mat.mean ~axis:1 x)) in
      let pcs = Mat.get_slice [ []; [ 0; dim - 1 ] ] modes in
      let x = Mat.(x *@ pcs) in
      let dx = Mat.(get_slice [ [ 1; -1 ]; [] ] x - get_slice [ [ 0; -2 ]; [] ] x) in
      Mat.(get_slice [ [ 1; -1 ]; [] ] x @|| dx))
    [| 1; 2; 3; 4; 5; 6; 7 |]
  |> Mat.concatenate ~axis:0

let get_torques_all pc pm =
  Array.map
    (fun i ->
      let x =
        Mat.get_slice
          [ [ 400; -200 ]; [] ]
          (Mat.load_txt
             (Printf.sprintf "results/weighing_pm/w_%i_%i/reach_%i/torques_400" pc pm i))
      in
      let m = Mat.(max ~axis:0 x - min ~axis:0 x) in
      let x = Mat.(x / m) in
      let dx = Mat.(get_slice [ [ 1; -1 ]; [] ] x - get_slice [ [ 0; -2 ]; [] ] x) in
      Mat.(get_slice [ [ 1; -1 ]; [] ] x @|| dx))
    [| 1; 2; 3; 4; 5; 6; 7 |]
  |> Mat.concatenate ~axis:0

let tangling_all y dim =
  let epsilon = 0.1 *. (y |> fun z -> Mat.(sum' (var ~axis:0 z))) in
  [| Mat.map_rows
       (fun x ->
         let x, dx =
           ( Mat.get_slice [ []; [ 0; dim - 1 ] ] Mat.(x - y) |> Mat.l2norm_sqr ~axis:1
           , Mat.get_slice [ []; [ dim; -1 ] ] Mat.((x - y) /$ dt)
             |> Mat.l2norm_sqr ~axis:1 )
         in
         Mat.max' Mat.(dx / (x +$ epsilon)))
       y
  |]
  |> Mat.of_arrays
  |> Mat.transpose

let tangling pc pm =
  let tangling_ac = tangling_all (get_all_x pc pm) dim in
  let tangling_to = tangling_all (get_torques_all pc pm) 2 in
  Mat.save_txt
    ~out:(Printf.sprintf "results/weighing_pm/entanglement/test_tangling_all_%i_%i" pc pm)
    Mat.(tangling_ac @|| tangling_to)

let _ =
  tangling 1 1;
  tangling 1 5;
  tangling 1 10;
  tangling 1 100;
  tangling 1 500;
  tangling 1 1000;
  tangling 5 1;
  tangling 10 1;
  tangling 100 1;
  tangling 500 1;
  tangling 1000 1

let _ =
  Array.map
    (fun (cp, cm) ->
      let x =
        Mat.load_txt
          (Printf.sprintf
             "results/weighing_pm/entanglement/test_tangling_all_%i_%i"
             (int_of_float cp)
             (int_of_float cm))
      in
      let y = Mat.(of_arrays [| [| cm /. cp |] |] @|| Mat.mean ~axis:0 x) in
      let _ = Printf.printf "%i %i %!" (Mat.row_num y) (Mat.col_num y) in
      y)
    [| 1000., 1.
     ; 500., 1.
     ; 100., 1.
     ; 10., 1.
     ; 5., 1.
     ; 1., 1.
     ; 1., 5.
     ; 1., 10.
     ; 1., 100.
     ; 1., 500.
     ; 1., 1000.
    |]
  |> Mat.concatenate ~axis:0
  |> Mat.save_txt ~out:"results/weighing_pm/entanglement/vary_prep"

let _ =
  let m_mov =
    Mat.col (Mat.load_txt "results/weighing_pm/entanglement/test_tangling_all_1000_1") 0
    |> Mat.mean'
  in
  let m_mixed =
    Mat.col (Mat.load_txt "results/weighing_pm/entanglement/test_tangling_all_1_1") 0
    |> Mat.mean'
  in
  let m_prep =
    Mat.col (Mat.load_txt "results/weighing_pm/entanglement/test_tangling_all_1_1000") 0
    |> Mat.mean'
  in
  Mat.save_txt
    ~out:"results/weighing_pm/entanglement/comp_entanglement"
    (Mat.transpose (Mat.of_arrays [| [| m_mov; m_mixed; m_prep |] |]))

(*do across conditions if I can!
Compute entanglement for projection of activity/compare w full state
get results for multiple netwtorks with inputs at different times and make a plot? 
For visualization make PCA plot across conditions
 *)
