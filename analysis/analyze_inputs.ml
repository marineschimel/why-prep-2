open Owl
open Base
module AD = Algodiff.D
module M = Arm.Make (Arm.Defaults)

let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir]")
let n = 7
let in_dir = Printf.sprintf "%s/%s" dir
let t_preps = [| 0.; 0.01; 0.02; 0.05; 0.2; 0.3; 0.5; 0.6; 0.8; 1. |]

let us i t_prep =
  Mat.load_txt (in_dir (Printf.sprintf "us_%i_%i" i Float.(to_int (1000. *. t_prep))))


let xs i t_prep =
  Mat.load_txt (in_dir (Printf.sprintf "xs_%i_%i" i Float.(to_int (1000. *. t_prep))))

let thetas i t_prep =
  Mat.load_txt (in_dir (Printf.sprintf "thetas_%i_%i" i Float.(to_int (1000. *. t_prep))))


let hands i t_prep =
  let thetas = 
  Mat.load_txt (in_dir (Printf.sprintf "thetas_%i_%i" i Float.(to_int (1000. *. t_prep))))
  in AD.unpack_arr (AD.Mat.map_by_row (fun x -> Arm.unpack_state_diff (M.hand_of (Arm.pack_state x))) (AD.pack_arr thetas))
  
let us_mov us t_prep =
  Mat.get_slice [ [ Float.(to_int (1000. *. t_prep)) + 1; -1 ]; [] ] us


let us_prep us t_prep = Mat.get_slice [ [ 0; Float.(to_int (1000. *. t_prep)) ]; [] ] us

let get_entropy us =
  let us = Mat.(us /$ Float.of_int (Mat.row_num us)) in
  let _, sis, _ = Linalg.D.svd us in
  let pi = Mat.(sis /$ sum' sis) in
  let e = Mat.(pi *@ neg (log (transpose pi))) in
  assert (Mat.col_num e = 1 && Mat.row_num e = 1);
  Mat.sum' e


let get_norm us = Mat.l2norm_sqr' us

let get_rank us =
  let _, sis, _ = Linalg.D.svd us in
  Mat.(sum' sis)


let energy us = Mat.l2norm_sqr' us
let abs us = Mat.mean' (Mat.abs us)

let get_mov_onset ~thetas =
  let hands =
    AD.Mat.map_by_row
      (fun x -> Arm.unpack_state_diff (M.hand_of (Arm.pack_state x)))
      (AD.pack_arr thetas)
  in
  let vel = Mat.get_fancy [ R [ 0; -1 ]; L [ 1; 3 ] ] (AD.unpack_arr hands) in
  let norm_vel = Mat.l2norm ~axis:1 vel in
  let threshold = 0.05 *. Mat.max' norm_vel in 
  let idces = Mat.filter (fun x -> Float.(x > threshold)) norm_vel in
  let flt_idces =
    Array.map ~f:(fun x -> Float.of_int x) idces |> fun z -> Mat.of_array z 1 (-1)
  in
  Mat.min' flt_idces


let gather_energies n =
  Array.map
    ~f:(fun t ->
      let us = us n t in
      [| t; energy (us_prep us t) /. energy (us_mov us t) |])
    t_preps


let gather_ranks n =
  Array.map
    ~f:(fun t ->
      let us = us n t in
      let rp, rm, rt = get_rank (us_prep us t), get_rank (us_mov us t), get_rank us in
      [| t; rp; rm; rt |])
    t_preps


let gather_entropies f n =
  Array.map
    ~f:(fun t ->
      let us = f n t in
      let ep, em, et =
        get_entropy (us_prep us t), get_entropy (us_mov us t), get_entropy us
      in
      [| t; ep; em; et |])
    t_preps


let gather_norms n =
  Array.map
    ~f:(fun t ->
      let us = us n t in
      let ep, em, et = get_norm (us_prep us t), get_norm (us_mov us t), get_norm us in
      [| t; ep; em; et |])
    t_preps


let get_move_onset n = 
  Array.map
  ~f:(fun t ->
    let thetas = thetas n t in
  let time_onset = get_mov_onset ~thetas
  in  [| t; time_onset |]) t_preps

let save_hands n =  Array.iter
~f:(fun t ->
  let hands = hands n t in 
Mat.save_txt ~out:(in_dir (Printf.sprintf "hands_%i_%i" n Float.(to_int (1000. *. t)))) hands) t_preps

let _ = Array.init 8 ~f:(fun n ->
  save_hands n;
  Mat.save_txt
    ~out:(in_dir (Printf.sprintf "analysis/prep_idx_%i" n))
    (Mat.of_arrays (gather_energies n));
  Mat.save_txt
    ~out:(in_dir (Printf.sprintf "analysis/norms_%i" n))
    (Mat.of_arrays (gather_norms n));
  Mat.save_txt
    ~out:(in_dir (Printf.sprintf "analysis/onset_%i" n))
    (Mat.of_arrays (get_move_onset n)))
(* let _ =
  Array.init 8 ~f:(fun n ->
      Mat.save_txt
        ~out:(in_dir (Printf.sprintf "entropies_%i" n))
        (Mat.of_arrays (gather_entropies us n));
      Mat.save_txt
        ~out:(in_dir (Printf.sprintf "entropies_x_%i" n))
        (Mat.of_arrays (gather_entropies xs n));
      Mat.save_txt
        ~out:(in_dir (Printf.sprintf "prep_idx_%i" n))
        (Mat.of_arrays (gather_energies n));
      Mat.save_txt
        ~out:(in_dir (Printf.sprintf "ranks_%i" n))
        (Mat.of_arrays (gather_ranks n));
      Mat.save_txt
        ~out:(in_dir (Printf.sprintf "norms_%i" n))
        (Mat.of_arrays (gather_norms n))) *)
