open Owl
module AD = Algodiff.D
open Lib
open Defaults

let _ = Printexc.record_backtrace true
let m = size_net

let max_speed x =
  let speed = Mat.get_fancy [ R [ 0; -1 ]; L [ 1; 3 ] ] x |> Mat.l2norm ~axis:1 in
  Mat.max' speed

let get_final_pos x =
  let m = Mat.get_slice [ [ -200 ]; [ 0; 1 ] ] x in
  Mat.mean ~axis:0 m

let accuracy_theta x =
  let tp = Mat.(get_fancy [ I 0; L [ 0; 1 ] ] (load_txt "data/target_thetas")) in
  let fp = get_final_pos x in
  let _ = Mat.print fp
in let _ = Mat.print tp in 
  Mat.(fp - tp) |> Mat.l2norm_sqr'


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
  let r = Mat.(eye m *$ (r_coeff *. 0.5)) in
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

  let dir = Cmdargs.(get_string "-d" |> force ~usage:"-dp [dir where param search is saved]")

 let ts =  [|0.;0.1;0.5;1.;2.5;5.|]
 let rs =   [|1.;1.2;1.5;2.;3.;5.;7.;10.;100.;1000.|] 

let torques_cost t_prep x = Mat.(l2norm_sqr' (get_slice [[0;int_of_float (t_prep*.1000.)];[]] x))
let load_traj r t  = Mat.load_txt (Printf.sprintf "%s/r_%.3f/t_%.3f/traj_300" dir r t)

let load_torques r t  = Mat.load_txt (Printf.sprintf "%s/r_%.3f/t_%.3f/torques_300" dir r t)

let _ = Printf.printf "checkpoint 1 %!"
let load_us r t = let _ = "true tr" in Mat.load_txt (Printf.sprintf "%s/r_%.3f/t_%.3f/results_us_300" dir r t)

let load_linesearch r t = Mat.load_txt (Printf.sprintf "%s/r_%.3f/t_%.3f/linesearch" dir r t)

let load_loss r t = Mat.load_txt (Printf.sprintf "%s/r_%.3f/t_%.3f/loss_time" dir r t)

let summary = 
(Array.map (fun r -> (Array.map (fun t -> 
  let loss_time = load_loss r t in 
 let n = Mat.row_num loss_time in 
 let _ = Printf.printf "%f %f %!" r t in 
  let inpts = Mat.get loss_time (n-1) 3 
in let tot_loss =  Mat.get loss_time (n-1) 2
in let torques =  Mat.get loss_time  (n-1) 4
in let acc = Mat.get loss_time  (n-1) 5 in 
(* let us = load_us r t in *)
(* let torques = load_torques r t in  *)
(* in let accuracy = accuracy_theta (Mat.get_slice [[];[0;3]] traj) in  *)
(* let ip,im = inputs ~t_prep:0.3 us in  *)
(* let torques = torques_cost 0.3 torques in  *)
Mat.of_arrays [|[|0.1*.(r)/.(float size_net);t;acc; torques; inpts; tot_loss;torques|]|]
)ts) |> fun y -> Mat.concatenate ~axis:0 y ) rs) |> fun z -> Mat.concatenate ~axis:0 z

let _ = Mat.save_txt ~out:(Printf.sprintf "%s/summary" dir) summary

(*add something st if it's not defined then "linearsearch error = true"*)