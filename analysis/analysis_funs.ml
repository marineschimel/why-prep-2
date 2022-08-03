open Owl
module AD = Algodiff.D
module M = Arm.Make (Arm.Defaults)

let _ = Printexc.record_backtrace true
let m = size_net

let max_speed x =
  let speed = Mat.get_fancy [ R [ 0; -1 ]; L [ 1; 3 ] ] x |> Mat.l2norm ~axis:1 in
  Mat.max' speed

let time_to_end x =
  let target_hands = Mat.load_txt "data/target_hand" in
  let target_pos = Mat.get_fancy [ R [ 0; -1 ]; L [ 0; 2 ] ] target_hands in
  let idces =
    Mat.filter_rows
      (fun x ->
        let x_pos = Mat.get_fancy [ R [ 0; -1 ]; L [ 0; 2 ] ] x in
        Mat.(l2norm' (x_pos - target_pos)) < 0.01)
      x
  in
  idces.(0)

let get_final_pos x =
  let m = Mat.get_fancy [ R [ -20; -1 ]; R [ 0; 1 ] ] x in
  Mat.mean ~axis:0 m

let energy ~beta ~alpha ~t_prep x =
  let r = Mat.(eye m *$ (r_coeff *. 0.5)) in
  let cost_function t u =
    Mat.(sum' (u *@ r *@ transpose u))
    *. alpha
    *. (1. +. ((beta -. 1.) *. Maths.sigmoid ((t -. t_prep) /. 2E-3)))
  in
  let m =
    Mat.mapi_rows
      (fun i u ->
        let t = sampling_dt *. float_of_int i in
        cost_function t u)
      x
  in
  Mat.of_array m (-1) 1

let energies ~t_prep e =
  let i_trans = int_of_float (t_prep /. sampling_dt) in
  let _ = Printf.printf "%i %i %!" (Mat.row_num e) (Mat.col_num e) in
  let preparation = Mat.get_fancy [ R [ 0; i_trans ]; R [] ] e in
  let movement = Mat.get_fancy [ R [ i_trans + 1; -1 ]; R [] ] e in
  Mat.mean' preparation, Mat.mean' movement

let compute_energies_beta ~beta ~alpha ~t_prep x =
  let ep, em = energies ~t_prep (energy ~beta ~alpha ~t_prep x) in
  Mat.(
    save_txt
      ~append:true
      ~out:(Printf.sprintf "results/energy/mean_energy_beta")
      (of_array [| beta; ep; em |] 1 (-1)))

let inputs ~t_prep x =
  let i_prep, i_mov =
    ( (Mat.of_array
         (Mat.mapi_rows
            (fun i u ->
              let t = sampling_dt *. float_of_int i in
              if t < t_prep then Mat.abs u |> Mat.sum' else 0.)
            x))
        1
        (-1)
      |> Mat.sum'
    , (Mat.of_array
         (Mat.mapi_rows
            (fun i u ->
              let t = sampling_dt *. float_of_int i in
              if t > t_prep then Mat.abs u |> Mat.sum' else 0.)
            x))
        1
        (-1)
      |> Mat.sum' )
  in
  i_prep, i_mov

let cost_inputs ~t_prep x =
  let r = Mat.(eye size_inputs *$ (r_coeff *. 0.5)) in
  let i_prep, i_mov =
    ( (Mat.of_array
         (Mat.mapi_rows
            (fun i u ->
              let t = sampling_dt *. float_of_int i in
              if t < t_prep then Mat.(u *@ r * u) |> Mat.sum' else 0.)
            x))
        1
        (-1)
      |> Mat.sum'
    , (Mat.of_array
         (Mat.mapi_rows
            (fun i u ->
              let t = sampling_dt *. float_of_int i in
              if t > t_prep then Mat.(u *@ r * u) |> Mat.sum' else 0.)
            x))
        1
        (-1)
      |> Mat.sum' )
  in
  i_prep, i_mov

let accuracy_theta x =
  let tp = Mat.(get_fancy [ I 0; L [ 0; 1 ] ] (load_txt "data/target_thetas")) in
  let fp =
    let m = Mat.get_slice [ [ -200 ]; [ 0; 1 ] ] x in
    Mat.mean ~axis:0 m
  in
  Mat.(fp - tp) |> Mat.l2norm_sqr'
