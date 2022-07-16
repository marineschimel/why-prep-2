open Owl
module AD = Algodiff.D
open Lib
open Base

let _ = Printexc.record_backtrace true
let dir = Cmdargs.(get_string "-d" |> default "/home/mmcs3/rds/hpc-work/_results/why_prep/results/uniform_1E-6")
let in_dir s = Printf.sprintf "%s/%s" dir s
let net_typ = "rdn/rdn"
let rads = [|0.2;0.5;0.95; 1.5|]
let seeds = [|6|]

let filename net rad seed = Printf.sprintf "%s/%s/rad_%.3f_seed_%i_full_summary" dir net rad seed

let make_summary_file = 
  Array.map ~f:(fun s -> Array.map ~f:(fun r -> let f = filename net_typ r s in 
  let summary = Mat.load_txt f in 
  let prep_idx = Mat.get summary 0 4 in 
  let norm = Mat.get summary 0 2
in [|r; norm; prep_idx|]
  ) rads |> fun x -> Mat.of_arrays x
  |> fun z -> Arr.reshape z [|(Mat.row_num z); (Mat.col_num z);1|]) seeds 
  |> fun z -> Arr.concatenate ~axis:2 z
 |> Arr.mean ~axis:2 |> fun z -> Arr.reshape z [|(Arr.shape z).(0); -1|] 

let _ = Mat.save_txt ~out:(in_dir (Printf.sprintf "%s/across_summary" net_typ)) make_summary_file


