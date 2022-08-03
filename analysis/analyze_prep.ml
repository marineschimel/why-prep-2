open Owl

let _ = Printexc.record_backtrace true

(*First, analysis of the scaled trajectories*)
let scaled_traj t_prep =
  Mat.load_txt (Printf.sprintf "scaling/scaled_traj_%i" (int_of_float (1000. *. t_prep)))
  |> fun x ->
  Mat.get_slice [ [ 1; -1 ]; [ 1; -1 ] ] x
  |> fun m -> Arr.reshape m [| (Arr.shape m).(0); (Arr.shape m).(1); 1 |]

let tpreps = [| 0.05; 0.1; 0.3; 0.5; 0.6; 0.8 |]

let scaled_trajs =
  Array.map (fun t -> scaled_traj t) tpreps |> fun z -> Arr.concatenate ~axis:2 z

let std_st = Arr.(var ~axis:2 scaled_trajs)
let mean_st = Arr.l2norm_sqr ~axis:2 scaled_trajs
let scaled_diff = Arr.((scaled_trajs - mean_st) / std_st)
let n, m, i = (Arr.shape std_st).(0), (Arr.shape std_st).(1), (Arr.shape std_st).(2)
let _ = Printf.printf "%i %i %i %!" n m i
let err_neurons = Arr.mean ~axis:0 Arr.(std_st / mean_st)

let _ =
  Mat.save_txt ~out:"scaling/err_neurons" (Mat.sqrt (Arr.reshape err_neurons [| m; 1 |]))

let corr =
  let get_ac ~tprep =
    let y =
      Mat.load_txt
        (Printf.sprintf "scaling/scaled_traj_%i" (int_of_float (tprep *. 1000.)))
    in
    Mat.reshape y [| 1; Mat.row_num y * Mat.col_num y |]
  in
  let ar =
    Array.map
      (fun tprep ->
        let ac = get_ac ~tprep in
        let x = Mat.(ac - mean ~axis:1 ac) in
        Mat.to_array Mat.(x / l2norm ~axis:1 x))
      tpreps
  in
  Mat.of_arrays ar

let _ =
  Mat.save_txt
    ~out:(Printf.sprintf "scaling/correlation_scaled")
    Mat.(corr *@ transpose corr)

(*let _ = (Arr.split ~axis:0 (Array.init (Array.length tpreps) (fun _ -> 1)) (Arr.transpose ~axis:[|2;0;1|] scaled_diff)) |> 
(Array.map (fun z -> Arr.reshape z [|n;m|]) ) |> Mat.concatenate ~axis:0 |> (Mat.sum ~axis:0) |> fun x -> Mat.save_txt ~out:"err_neurons" x*)
