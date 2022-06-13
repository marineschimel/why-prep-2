open Owl
module AD = Algodiff.D
open Printf

let n_mov = 6

let get_mskew ~dx ~x =
  let max_steps = 50000 in
  let n_lat = Mat.col_num x.(0) in
  let x = x |> Array.map (Mat.get_slice [ [ 0; -2 ] ]) in
  let open Algodiff.D in
  let dx = Array.map pack_arr dx in
  let x = Array.map pack_arr x in
  let cost w =
    let rec accu i c =
      if i < n_mov
      then accu (succ i) Maths.(c + l2norm_sqr' (dx.(i) - (x.(i) *@ w)))
      else Maths.sqrt c
    in
    accu 0 (F 0.)
  in
  let dcost = grad cost in
  let w = Mat.gaussian n_lat n_lat in
  let alpha = F 0.1 in
  let rec optimise w step old_c e =
    let w = Maths.(F 0.5 * (w - transpose w)) in
    if e > 1E-9 && step < max_steps
    then (
      let c = cost w |> unpack_flt in
      let e = Stdlib.((old_c -. c) /. old_c) in
      if step mod 10 = 0 then printf "\rstep %i | cost: %f | epsilon: %f%!" step c e;
      let dw = dcost w in
      let w = Maths.(w - (alpha * dw)) in
      optimise w (succ step) c e)
    else w
  in
  optimise w 0 1E9 1. |> unpack_arr


let get_modes_jpca mskew =
  let module Z = Dense.Matrix.Z in
  let u, s = Linalg.D.eig mskew in
  let er = Z.re s in
  let ei = Z.im s in
  printf "real eigenvalues \n";
  Mat.print er;
  printf "imag eigenvalues \n";
  Mat.print ei;
  let v1 = Z.col u 0 in
  let v2 = Z.col u 1 in
  let j1 = Z.(re Mat.(v1 + v2)) in
  let j1 = Z.(j1 /$ Mat.l2norm' j1) in
  let j2 = Z.(im (v1 - v2)) in
  let j2 = Mat.(j2 /$ l2norm' j2) in
  Z.(j1 @|| j2)
