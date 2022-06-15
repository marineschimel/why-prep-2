open Owl
module AD = Algodiff.D

let _ = Printexc.record_backtrace true
let m = 200
let sampling_dt = 2E-3

let max_speed x =
  let speed = Mat.get_fancy [ R [ 0; -1 ]; L [ 1; 3 ] ] x |> Mat.l2norm ~axis:1 in
  Mat.max' speed


let time_to_end x target_hands n_prep =
  let target_pos = Mat.get_fancy [ R [ 0; -1 ]; L [ 0; 2 ] ] target_hands in
  let idces =
    Mat.filter_rows
      (fun x ->
        let x_pos = Mat.get_fancy [ R [ 0; -1 ]; L [ 0; 2 ] ] x in
        Mat.(l2norm' (x_pos - target_pos)) < 0.01)
      x
  in
  try
    let i = idces.(0) - n_prep in
    i
  with
  | _ -> -100


let get_final_pos x =
  let m = Mat.get_fancy [ R [ -20; -1 ]; R [ 0; 1 ] ] x in
  Mat.mean ~axis:0 m


let cost_u ~f ~n_prep x =
  let arr = Mat.mapi_rows f x in
  let mat = Mat.of_array arr (-1) 1 in
  if n_prep = 0
  then 0., Mat.sum' mat, Mat.sum' mat
  else
    ( Mat.sum' (Mat.get_slice [ [ 0; n_prep - 1 ] ] mat)
    , Mat.sum' (Mat.get_slice [ [ n_prep; -1 ] ] mat)
    , Mat.sum' mat )


let cost_x ~f ~n_prep x =
  let arr = Mat.mapi_rows f x in
  let mat = Mat.of_array arr (-1) 1 in
  if n_prep = 0
  then 0., Mat.sum' mat
  else
    ( Mat.sum' (Mat.get_slice [ [ 0; n_prep - 1 ] ] mat)
    , Mat.sum' (Mat.get_slice [ [ n_prep; -1 ] ] mat) )


let energies ~t_prep ~dt e =
  let i_trans = Int.of_float (t_prep /. dt) in
  let preparation = Mat.get_fancy [ R [ 0; i_trans ]; R [] ] e in
  let movement = Mat.get_fancy [ R [ i_trans + 1; -1 ]; R [] ] e in
  Mat.mean' preparation, Mat.mean' movement


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


let accuracy_theta x =
  let tp = Mat.(get_fancy [ I 0; L [ 0; 1 ] ] (load_txt "data/target_thetas")) in
  let fp =
    let m = Mat.get_slice [ [ -200 ]; [ 0; 1 ] ] x in
    Mat.mean ~axis:0 m
  in
  Mat.(fp - tp) |> Mat.l2norm_sqr'
