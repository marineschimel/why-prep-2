open Owl
module AD = Algodiff.D
module M = Arm.Make (Arm.Defaults)
open Lib
open Base

let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")
let in_dir s = Printf.sprintf "%s/%s" dir s
let t_preps = [|0; 100; 200; 300; 500|]
let hands t = let h = AD.Mat.map_by_row (fun x -> Arm.unpack_state_diff (M.hand_of (Arm.pack_state x))) (AD.pack_arr t) in AD.unpack_arr h

(* let _ = Array.map t_preps ~f:(fun t -> Array.init 8 ~f:(fun i ->let thetas = Mat.load_txt (in_dir (Printf.sprintf "thetas_%i_%i" i t)) in Mat.save_txt ~out:(in_dir (Printf.sprintf "hands_%i_%i" i t)) (hands thetas))) *)

let _ = let tgts = Mat.load_txt (in_dir (Printf.sprintf "targets")) in Mat.save_txt ~out:(in_dir (Printf.sprintf "target_hands")) (hands tgts)