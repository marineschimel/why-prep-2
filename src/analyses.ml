open Owl

let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")
let in_dir s = Printf.sprintf "%s/%s" dir s

let ls_double_reaches =
  [| 0, 3
   ; 0, 4
   ; 0, 5
   ; 1, 2
   ; 1, 5
   ; 2, 3
   ; 2, 1
   ; 2, 4
   ; 2, 5
   ; 3, 2
   ; 4, 2
   ; 4, 3
   ; 4, 6
   ; 5, 1
   ; 5, 2
   ; 5, 3
   ; 5, 4
   ; 6, 0
   ; 6, 1
   ; 6, 2
  |]

let average_torques =
  Array.map
    (fun (i, j) ->
      let x = Mat.load_txt (in_dir (Printf.sprintf "torques_%i_%i_500" i j)) in
      Arr.reshape x [| 1; -1; Mat.col_num x |])
    ls_double_reaches
  |> Arr.concatenate ~axis:0
  |> Arr.mean ~axis:0
  |> fun z -> Arr.reshape z [| -1; 2 |]

let average_vel =
  Array.map
    (fun (i, j) ->
      let x = Mat.load_txt (in_dir (Printf.sprintf "thetas_%i_%i_500" i j)) in
      let y = Mat.get_slice [ []; [ 2; 3 ] ] x in
      Mat.l2norm ~axis:1 y)
    ls_double_reaches
  |> Mat.concatenate ~axis:1
  |> Mat.mean ~axis:1

let _ = Mat.save_txt ~out:(in_dir "analysis/mean_torques") average_torques
let _ = Mat.save_txt ~out:(in_dir "analysis/mean_vel") average_vel
