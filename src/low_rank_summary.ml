open Owl
module AD = Algodiff.D
open Lib
open Base

let _ = Printexc.record_backtrace true
let dir = Cmdargs.(get_string "-d" |> default "/home/mmcs3/rds/rds-t2-cs133-hh9aMiOkJqI/mmcs3/low_rank_monkeys_lambda_1E-7")
let in_dir s = Printf.sprintf "%s/%s" dir s
let rank =  Cmdargs.(get_int "-rank" |> default 50)
let seed = 1
let n_preps = [|0;50;100;150;200;300|]
let get_loss tgt prep = Mat.load_txt (in_dir (Printf.sprintf "b_%i/seed_%i/loss_%i_%i" rank seed tgt prep)) |> Mat.sum'
let n_reaches = [|0;1|]
let get_losses = Array.map n_reaches ~f:(fun i -> 
  Array.map n_preps ~f:(fun n -> [|Float.of_int n; get_loss i n|]
  ) |> Mat.of_arrays |> fun z -> Arr.reshape z [|(-1);2;1|]) 
  |> fun z -> Arr.concatenate ~axis:2 z |>
   fun z -> Arr.mean ~axis:2 z |> fun z -> Arr.reshape z [|-1;2|] 

let _ = Mat.save_txt ~out:(in_dir (Printf.sprintf "b_%i/seed_%i/losses_all" rank seed)) get_losses