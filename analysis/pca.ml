open Owl

let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [directory we're getting files from]")

let t_prep = Cmdargs.(get_int "-t" |> force ~usage:"-t [prep times considered]")

let n_reaches = 4

let dim = 8

let pca_decomp x = let cov = Mat.((transpose x)*@x) in let modes,_,_ = Linalg.D.svd cov
in Mat.get_slice [[];[0;dim]] modes 

let pca_proj x = Mat.(x*@(pca_decomp x))
let full_modes = Array.init n_reaches (fun i -> Mat.(get_slice [[];[4;-1]] (load_txt (Printf.sprintf "%s/reach_%i/traj_%i" dir (succ i) t_prep)))) |> Mat.concatenate ~axis:0 |> fun x -> pca_decomp x
let _ = Array.init n_reaches (fun i -> let x = Mat.(get_slice [[];[4;-1]] (load_txt (Printf.sprintf "%s/reach_%i/traj_%i" dir (succ i) t_prep))) in Mat.save_txt ~out:(Printf.sprintf "%s/individual_proj_%i_%i" dir t_prep i) (pca_proj x); Mat.save_txt ~out:(Printf.sprintf "%s/full_proj_%i_%i" dir t_prep i) (Mat.(x*@(full_modes))))
