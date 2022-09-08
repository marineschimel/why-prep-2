open Owl
open Base

let n_targets = 8

let t_prep = 300
let n_prep = Int.(t_prep/2)

(*get obs and controllabiloty gramians*)
(*get the modes, as well as their eigenvalues*)
(*compute the mean energy across movements as a function of time *)
let n = 200

let gen_traj f1 f2 = 
let test t = let x = Mat.of_array [|Float.(Maths.cos (t /. f1) *. exp (neg t/.f2)); Float.(Maths.sin (t /. f1) *. exp (neg t/.f2))|] 1 2 in Mat.(x@||zeros 1 Int.(n-2))
in 
let traj = Array.init 250 ~f:(fun i -> test (Float.of_int i)) |> Mat.concatenate ~axis:0
in 
let traj_rev =  Array.init n_prep ~f:(fun i -> test (Float.of_int Int.(250 - i)))|> Mat.concatenate ~axis:0 in 
Mat.(traj_rev @= traj)


let movs = Array.init n_targets ~f:(fun i ->
  let f1 =Float.(5.*. (of_int Int.(i + 1)))
in let f2 = 150. in let t = gen_traj f1 f2
in Mat.save_txt ~out:(Printf.sprintf "ctrls/rates_%i_%i" i t_prep) t)