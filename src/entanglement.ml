open Owl
open Lib
open Defaults

(*Q(t) = max_t' ||(\dot{x}(t)-\dot{x}(t')||^2/(||(\dot{x}(t)-\dot{x}(t')||^2+\epsilon)*)
let epsilon = 1E-6
let dt = sampling_dt

let tangling t traj_1 traj_2 =
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

let activity ?reach:_reach net time =
  let reach =
    match _reach with
    | Some a -> a
    | None   -> 1
  in
  Mat.get_slice
    [ [ time; -1 ]; [ 4; -1 ] ]
    (Mat.load_txt (Printf.sprintf "results/%s/reach_%i/traj_%i" net reach time))


let muscles ?reach:_reach net time =
  let reach =
    match _reach with
    | Some a -> a
    | None   -> 1
  in
  Mat.load_txt (Printf.sprintf "results/%s/reach_%i/torques_%i" net reach time)


let _ =
  let tang, i = tangling 1 (activity "soc" 300) (activity "soc" 300) in
  Printf.printf "%f %i \n %!" tang i


let save_tangling net t_prep =
  let nr = Mat.row_num (activity net t_prep) - 1 in
  Mat.save_txt
    ~out:(Printf.sprintf "entanglement/%s/traj_%i" net t_prep)
    Mat.(
      init nr 1 (fun i ->
          let tang, _ = tangling (succ i) (activity net t_prep) (activity net t_prep) in
          tang)
      @|| init nr 1 (fun i ->
              let tang, _ =
                tangling (succ i) (activity net t_prep) (activity ~reach:2 net t_prep)
              in
              tang));
  Mat.save_txt
    ~out:(Printf.sprintf "entanglement/%s/torques_%i" net t_prep)
    (Mat.init
       (Mat.row_num (muscles net t_prep) - 1)
       1
       (fun i ->
         let tang, _ = tangling (succ i) (muscles net t_prep) (muscles net t_prep) in
         tang))


let _ = save_tangling "soc" 300

(*do across conditions if I can!
Compute entanglement for projection of activity/compare w full state
get results for multiple netwtorks with inputs at different times and make a plot? 
For visualization make PCA plot across conditions
 *)
