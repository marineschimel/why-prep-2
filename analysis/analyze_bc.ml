open Owl
open Lib
open Defaults

let in_dir s = Printf.sprintf "results_bc_r2/bc_opt_prep/%s" s

let w = Mat.load_txt "data/w_rec"

let ctrl_gramian b=
  Linalg.D.lyapunov
    Mat.(transpose Mat.((w - eye 200) /$ tau))
    Mat.(neg (transpose b) *@ b)

    let obs_gramian c=
  Linalg.D.lyapunov
    Mat.(transpose Mat.((w - eye 200) /$ tau))
    Mat.(neg (transpose c) *@ c)
let proj i = let b = Mat.load_txt (in_dir (Printf.sprintf "b_new_%i" i)) in 
let c = Mat.load_txt (in_dir (Printf.sprintf "c_new_%i" i)) 
in let _proj = Mat.(c*@(transpose b)) in let cg = ctrl_gramian b in let co = obs_gramian c
in Mat.(trace cg), Mat.trace co


(* let evol_projs = Array.init 15 (fun i -> (proj (succ i)) |> Mat.sqr |> Mat.sum ~axis:1 |> Mat.transpose) |> Mat.concatenate ~axis:0 *)


let evol_traces = Array.init 16 (fun i -> let cg, co = proj (succ i) in [|cg;co|]) |> Mat.of_arrays
let _ = Mat.save_txt ~out:(in_dir "evol_traces") evol_traces