open Owl
open Base
module AD = Algodiff.D
module M = Arm.Make (Arm.Defaults)

let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir]")
let n = 7
let in_dir = Printf.sprintf "%s/%s" dir

let us i t_prep =
  Mat.load_txt (in_dir (Printf.sprintf "us_%i_%i" i Float.(to_int (1000. *. t_prep))))


let us_mov us t_prep = Mat.get_slice [ [ Float.(to_int (1000. *. t_prep)); -1 ]; [] ] us

let us_prep us t_prep =
  Mat.get_slice [ [ 0; Float.(to_int (1000. *. t_prep)) - 1 ]; [] ] us


let get_entropy us =
  let _, sis, _ = Linalg.D.svd us in
  let pi = Mat.(sqr sis /$ sum' (sqr sis)) in
  let e = Mat.(pi *@ log (transpose pi)) in
  assert (Mat.row_num e = 1);
  Mat.sum' e


let energy us = Mat.l2norm_sqr' us
let abs us = Mat.mean' (Mat.abs us)

let get_mov_onset ~threshold ~thetas =
  let hands =
    AD.Mat.map_by_row
      (fun x -> Arm.unpack_state_diff (M.hand_of (Arm.pack_state x)))
      (AD.pack_arr thetas)
  in
  let vel = Mat.get_fancy [ R [ 0; -1 ]; L [ 1; 3 ] ] (AD.unpack_arr hands) in
  let norm_vel = Mat.l2norm ~axis:1 vel in
  let idces = Mat.filter (fun x -> Float.(x > threshold)) norm_vel in
  let flt_idces =
    Array.map ~f:(fun x -> Float.of_int x) idces |> fun z -> Mat.of_array z 1 (-1)
  in
  Mat.min' flt_idces


let gather_energies =
  Array.map
    ~f:(fun t ->
      let us = us n t in
      [| t; energy (us_prep us t) /. energy (us_mov us t) |])
    [| 0.; 0.01; 0.02; 0.05; 0.1; 0.15; 0.2; 0.3; 0.4; 0.5; 0.6 |]


let gather_entropies =
  Array.map
    ~f:(fun t ->
      let us = us n t in
      let ep, em, et =
        get_entropy (us_prep us t), get_entropy (us_mov us t), get_entropy us
      in
      [| t; ep; em; et |])
    [| 0.; 0.01; 0.02; 0.05; 0.1; 0.15; 0.2; 0.3; 0.4; 0.5; 0.6; 0.8; 1. |]


let _ =
  Mat.save_txt ~out:(in_dir "entropies_7") (Mat.of_arrays gather_entropies);
  Mat.save_txt ~out:(in_dir "prep_idx_7") (Mat.of_arrays gather_energies)
