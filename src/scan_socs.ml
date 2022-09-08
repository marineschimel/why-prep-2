open Owl
module AD = Algodiff.D
module M = Arm.Make (Arm.Defaults)
open Lib
open Printf
module P = Stdlib

(* Setting up the parameters/directories
*)
let dir = Cmdargs.(get_string "-data" |> force ~usage:"-d [dir to save in]")
let in_dir s = Printf.sprintf "%s/%s" dir s

let data_dir = Cmdargs.(get_string "-data" |> force ~usage:"-data [dir in which data is]")
let n = 200
let sa_tgt = Cmdargs.(get_float "-sa" |> default 0.8)
let eta = Cmdargs.(get_float "-eta" |> default 1.)
let n_e = Maths.round (0.8 *. float n) |> int_of_float
let n_i = n - n_e
let p_con = Cmdargs.(get_float "-p_con" |> default 0.2) (* exc. connection density *)
let rhs = Mat.(neg (eye n))
let dc_eval = -10.
let e = [ []; [ 0; pred n_e ] ]
let i = [ []; [ n_e; -1 ] ]
let ee = [ [ 0; pred n_e ]; [ 0; pred n_e ] ]
let ei = [ [ 0; pred n_e ]; [ n_e; -1 ] ]
let ie = [ [ n_e; -1 ]; [ 0; pred n_e ] ]
let ii = [ [ n_e; -1 ]; [ n_e; -1 ] ]
let ninh = Mat.init_2d n 1 (fun i _ -> float (if i < n_e then n_i else pred n_i))

let eigenvalues m =
  let v = Linalg.D.eigvals m in
  let re = Dense.Matrix.Z.re v in
  let im = Dense.Matrix.Z.im v in
  Mat.(concat_horizontal (transpose re) (transpose im))

let spectral_abscissa m = Linalg.D.eig m |> snd |> Dense.Matrix.Z.re |> Mat.max'

let save ?step label w =
  let open Printf in 
  let w_eigs = eigenvalues w in
  Mat.(save_txt ~out:(in_dir Printf.(sprintf "socs/w_%s" label)) w);
  Mat.save_txt w_eigs ~out:(in_dir (sprintf "socs/w_%s_eig" label))

let normalize w =
  let z =
    Mat.((dc_eval $- sum ~axis:1 (get_slice e w) + sum ~axis:1 (get_slice i w)) / ninh)
  in
  Mat.(set_slice i w (neg (relu (neg (z + get_slice i w)))));
  for i = 0 to n - 1 do
    Mat.set w i i 0.
  done
    
let _ =
Array.mapi (fun j init_radius ->
  let open Stdlib in 
  if Base.Int.(j % C.n_nodes = C.rank) then
    let w =
      let w = Mat.zeros n n in
      let ids_e = Array.init n_e (fun i -> i) in
      let k_e = Maths.round (p_con *. float n_e) |> int_of_float in
      let rec draw_ids i =
        let ids = Stats.choose ids_e k_e in
        if Array.exists (( = ) i) ids then draw_ids i else ids
      in
      let w0 = init_radius /. sqrt (float n *. p_con *. (1. -. p_con)) in
      for i = 0 to pred n do
        draw_ids i
        |> Array.iter (fun j ->
               Mat.set w i j (exp (log w0 +. (0.5 *. Stats.gaussian_rvs ~mu:0. ~sigma:1.))));
        for j = n_e to pred n do
          if i <> j then Mat.set w i j (-.Random.float 1.)
        done
      done;
      normalize w;
      w
  in 
let rec iterate k sas =
  if k mod 10 = 0
  then (
    save ~step:(k / 10) (Printf.sprintf "rad_%.1f_sa_%.2f" init_radius sa_tgt) w);
  let sa = spectral_abscissa w in
  printf "\riteration %5i | sa = %.5f%!" k sa;
  let shift = max 1. (1.2 *. sa) in
  let w_s = Mat.(add_diag w Maths.(neg shift)) in
  let p = Linalg.D.lyapunov w_s rhs in
  let q = Linalg.D.lyapunov (Mat.transpose w_s) rhs in
  let grad = Mat.(q *@ p) in
  let eta = eta /. Mat.trace grad in
  let grad = Mat.(grad - diagm (diag grad)) in
  Mat.(set_slice i w (neg (relu (neg (get_slice i w - (eta $* get_slice i grad))))));
  normalize w;
  if sa > sa_tgt then iterate (k + 1) (sa :: sas) else ()
in 
  iterate 0 []) [|0.3; 0.5; 0.7; 0.8; 1.0; 1.1; 1.3; 1.4; 1.6; 1.7; 1.8; 1.9; 2.1|]
  
  (* 2.2; 2.7; 3.2; 3.7; 1.0; 1.5; 2.0; 2.5; 3.0; 3.5; 4.0; 4.5; 5.0; 5.5; 6.0; 6.5; 10.0; 1.3; 3.1; 3.4; 3.6; 3.8; 4.1; 4.2;4.3; 4.6; 4.9; 5.2; 5.3; 5.4; 5.7; 5.9; 6.1; 6.3; 6.7; 6.9|] *)
   (* [|3.1; 3.4; 3.6; 3.8; 4.1; 4.2;4.3; 4.6; 4.9; 5.2; 5.3; 5.4; 5.7; 5.9; 6.1; 6.3; 6.7; 6.9|] 
      [|1.0; 1.5; 2.0; 2.5; 3.0; 3.5; 4.0; 4.5; 5.0; 5.5; 6.0; 6.5; 10.0|]*)
(* 
  [|0.2; 0.4; 0.6; 0.9; 1.3; 2.2; 2.7; 3.2; 3.7|] *)

(* 
  open Printf
open Owl
module P = Stdlib

let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")
let in_dir s = Printf.sprintf "%s/%s" dir s
let save_all = Cmdargs.(check "-save_all")
let eta = Cmdargs.(get_float "-eta" |> default 1.)
let n = Cmdargs.(get_int "-n" |> default 200)
let sa_tgt = Cmdargs.(get_float "-sa" |> default 0.8)
let n_e = Maths.round (0.8 *. float n) |> int_of_float
let n_i = n - n_e
let p_con = Cmdargs.(get_float "-p_con" |> default 0.2) (* exc. connection density *)

let radius = Cmdargs.(get_float "-radius" |> default 1.5)
let rhs = Mat.(neg (eye n))
let dc_eval = -10.

(* useful shortcut for slicing *)
let e = [ []; [ 0; pred n_e ] ]
let i = [ []; [ n_e; -1 ] ]
let ee = [ [ 0; pred n_e ]; [ 0; pred n_e ] ]
let ei = [ [ 0; pred n_e ]; [ n_e; -1 ] ]
let ie = [ [ n_e; -1 ]; [ 0; pred n_e ] ]
let ii = [ [ n_e; -1 ]; [ n_e; -1 ] ]
let ninh = Mat.init_2d n 1 (fun i _ -> float (if i < n_e then n_i else pred n_i))

(* normalise I connectivity to preserve DC mode *)
let normalize w =
  let z =
    Mat.((dc_eval $- sum ~axis:1 (get_slice e w) + sum ~axis:1 (get_slice i w)) / ninh)
  in
  Mat.(set_slice i w (neg (relu (neg (z + get_slice i w)))));
  for i = 0 to n - 1 do
    Mat.set w i i 0.
  done

(* initial W *)
let w =
  let w = Mat.zeros n n in
  let ids_e = Array.init n_e (fun i -> i) in
  let k_e = Maths.round (p_con *. float n_e) |> int_of_float in
  let rec draw_ids i =
    let ids = Stats.choose ids_e k_e in
    if Array.exists (( = ) i) ids then draw_ids i else ids
  in
  let w0 = radius /. sqrt (float n *. p_con *. (1. -. p_con)) in
  for i = 0 to pred n do
    draw_ids i
    |> Array.iter (fun j ->
           Mat.set w i j (exp (log w0 +. (0.5 *. Stats.gaussian_rvs ~mu:0. ~sigma:1.))));
    for j = n_e to pred n do
      if i <> j then Mat.set w i j (-.Random.float 1.)
    done
  done;
  normalize w;
  w

let eigenvalues m =
  let v = Linalg.D.eigvals m in
  let re = Dense.Matrix.Z.re v in
  let im = Dense.Matrix.Z.im v in
  Mat.(concat_horizontal (transpose re) (transpose im))

let spectral_abscissa m = Linalg.D.eig m |> snd |> Dense.Matrix.Z.re |> Mat.max'

let save ?step label =
  let w_eigs = eigenvalues w in
  Mat.(save_txt w ~out:(in_dir (sprintf "w_%s" label)));
  Mat.save_txt w_eigs ~out:(in_dir (sprintf "w_%s_eig" label));
  if save_all
  then (
    match step with
    | Some s -> Mat.save_txt w_eigs ~out:(in_dir (sprintf "movie/eig%i" s))
    | None -> ())

let save_training sas =
  [| sas |> List.rev |> Array.of_list |]
  |> Mat.of_arrays
  |> Mat.transpose
  |> Mat.save_txt ~out:(in_dir "training_info")

let rec iterate k sas =
  if k mod 10 = 0
  then (
    save ~step:(k / 10) (Printf.sprintf "rec_11");
    if Cmdargs.(check "-save_training") then save_training sas);
  let sa = spectral_abscissa w in
  printf "\riteration %5i | sa = %.5f%!" k sa;
  let shift = max 1. (1.2 *. sa) in
  let w_s = Mat.(add_diag w Maths.(neg shift)) in
  let p = Linalg.D.lyapunov w_s rhs in
  let q = Linalg.D.lyapunov (Mat.transpose w_s) rhs in
  let grad = Mat.(q *@ p) in
  let eta = eta /. Mat.trace grad in
  let grad = Mat.(grad - diagm (diag grad)) in
  Mat.(set_slice i w (neg (relu (neg (get_slice i w - (eta $* get_slice i grad))))));
  normalize w;
  if sa > sa_tgt then iterate (k + 1) (sa :: sas) else ()

let () =
  if save_all then save "init";
  iterate 0 [] *)
