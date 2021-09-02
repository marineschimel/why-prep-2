open Owl
open Base
open Lib

let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")
let data_dir = Cmdargs.(get_string "-data" |> force ~usage:"-data [dir in which data is]")
let in_dir s = Printf.sprintf "%s/%s" dir s
let in_data_dir s = Printf.sprintf "%s/%s" data_dir s
let c = Mat.gaussian ~sigma:0.2 40 50
let _ = Mat.save_txt ~out:(in_dir "c") c

let gather_hands = Array.init 144 ~f:(fun i -> let x = 
  Mat.load_txt (in_data_dir (Printf.sprintf "subsamp_hands_%i_300" i))
in let x = Mat.get_fancy [R [0;-1];L [0;2]] x 
in Arr.reshape x [|1;-1;2|]) |> Arr.concatenate ~axis:0


let gather_ac = Array.init 144 ~f:(fun i -> let x = 
  Mat.load_txt (in_data_dir (Printf.sprintf "subsamp_xs_%i_300" i))
in Arr.reshape x [|1;-1;40|]) |> Arr.concatenate ~axis:0


let gather_neural = Array.init 144 ~f:(fun i -> let x = 
  Mat.load_txt (in_data_dir (Printf.sprintf "subsamp_xs_%i_300" i))
in let x = Mat.(x*@c) in Arr.reshape x [|1;-1;50|]) |> Arr.concatenate ~axis:0

let gather_us = Array.init 144 ~f:(fun i -> let x = 
  Mat.load_txt (in_data_dir (Printf.sprintf "subsamp_us_%i_300" i))
in Arr.reshape x [|1;-1;40|]) |> Arr.concatenate ~axis:0


let _ = Arr.save_npy ~out:(in_dir "true_hands.npy") gather_hands;
Arr.save_npy ~out:(in_dir "true_xs.npy") gather_ac;
Arr.save_npy ~out:(in_dir "true_neural.npy") gather_neural;
Arr.save_npy ~out:(in_dir "true_us.npy") gather_us;
Misc.save_bin (in_dir "true_hands") gather_hands;
Misc.save_bin (in_dir "true_xs") gather_ac;
Misc.save_bin (in_dir "true_us") gather_us;
Misc.save_bin (in_dir "true_neural") gather_neural


