open Owl


let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir]")
let epsilon = Cmdargs.(get_float "-eps" |> force ~usage:"-eps [dir]")
let in_dir = Printf.sprintf "%s/%s" dir
let noise_level = 0.00005
let ts = [| 0.; 0.01; 0.02; 0.05; 0.2; 0.3; 0.5; 0.6; 0.8; 1. |]
let accuracy i t_prep = 
let acc = Mat.load_txt (in_dir (Printf.sprintf "added_noise/accuracy_%i_%i_%.6f_%.3f" i Float.(to_int (1000. *. t_prep)) noise_level epsilon)) in Mat.mean' (Mat.sqr acc)


let accuracy_all t_prep = (Array.init 7 (fun i -> accuracy i t_prep)) |> fun a -> Mat.of_arrays [|a|] |> fun m -> Mat.sum' m

let _ = 
(Array.map (fun t -> [|t; accuracy_all t|]) ts) |> fun z -> 
  Mat.of_arrays z |> fun m -> Mat.save_txt ~out:(in_dir (Printf.sprintf "plotting/accuracies_%.6f_%.3f" noise_level epsilon)) m
