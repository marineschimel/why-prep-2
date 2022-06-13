open Owl
open Base

let dir = "/home/mmcs3/rds/rds-t2-cs133-hh9aMiOkJqI/mmcs3/random_monkeys_lambda_1E-7"
let vprep i = Mat.load_txt (Printf.sprintf "%s/seed_%i/projs_300_350_short/vprep" dir i)

let vpreps = (Array.map [|1;2;3;4|] ~f:(fun i -> let v = vprep i in v)) |> fun z -> Mat.concatenate ~axis:1 z |> fun m -> let me = Mat.mean ~axis:1 m in let v = Mat.std ~axis:1 m in Mat.(me@||v) |> fun z -> Mat.save_txt ~out:(Printf.sprintf "%s/vpreps_prep" dir) z