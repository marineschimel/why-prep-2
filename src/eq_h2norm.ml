open Owl
module AD = Algodiff.D

(* optimize \alpha to have H2norm equal to H_tgt
H2norm given by Tr(CPC^T) where AP + PA^T + BB^T = 0 *)

let dir = Cmdargs.(get_string "-d" |> default "/home/mmcs3/rds/hpc-work/_results/why_prep/results/uniform_1E-7/skew")
let data_dir = Cmdargs.(get_string "-data" |> default "/home/mmcs3/rds/hpc-work/_results/why_prep/data")
let in_dir s = Printf.sprintf "%s/%s" dir s
let n = 200
let c0 = AD.pack_arr (Mat.load_txt (in_dir "c"))
let b = AD.Mat.eye n
let w_tgt = AD.pack_arr (Mat.load_txt (Printf.sprintf "%s/w_rec_6" data_dir))
let c_tgt =  AD.pack_arr (Mat.load_txt (Printf.sprintf "%s/c" data_dir))
let w0 = let m = AD.Mat.gaussian ~sigma:(0.6/.sqrt(200.)) n n
in AD.Maths.(m - transpose m)

let h2norm w c = let a = AD.Maths.((w - AD.Mat.eye n)) 
in let p = AD.Linalg.lyapunov (AD.Maths.transpose a) (AD.Maths.(neg (AD.Mat.eye n)))
in AD.Maths.trace (AD.Maths.(c*@p*@(transpose c)))
let h2_tgt = h2norm w_tgt c_tgt
let mu = AD.F 0.


module Prms = struct
  type 'a t = { alpha : 'a; mu : 'a } [@@deriving prms]
end

open Prms
module P = Owl_opt_lbfgs.Make (Prms)

let cost prms = 
  let alpha = prms.alpha in 
  let mu = prms.mu in 
  let w = AD.Maths.(alpha * w0 + mu*AD.Mat.eye n) in 
  let h2norm = h2norm w c0 in AD.Maths.(sqr (h2_tgt - h2norm))



let learn alpha0 =
  let prms0 = { alpha = AD.F alpha0; mu = AD.F 0.3 } in
  let stop =
    let fv_prev = ref 1E9 in
    fun s ->
      let k = P.iter s in
      let fv = P.(fv s) in
      let pct_change = abs_float ((fv -. !fv_prev) /. !fv_prev) in
      fv_prev := fv;
      Printf.printf "\r iter %i | fv %f | pct change %f %!" k fv pct_change;
      pct_change < 1E-3
  in
  let s0 = P.init ~f:cost ~prms0 () in
  let sf = P.min ~stop s0 in
  let prms = P.prms sf in
  let alpha = prms.alpha |> AD.unpack_flt
in let mu = prms.mu |> AD.unpack_flt 
in let _ = Printf.printf "alpha = %f || mu %f || val : %f %!" alpha mu (P.fv sf) in  alpha, mu


let _ = 
  let alpha, mu = learn 0. in
  let new_w =
  AD.Maths.(w0*AD.F alpha + AD.F mu * AD.Mat.eye n) |> AD.unpack_arr 
  in 
  let eigenvalues m =
    let v = Linalg.D.eigvals m in
    let re = Dense.Matrix.Z.re v in
    let im = Dense.Matrix.Z.im v in
    Mat.(concat_horizontal (transpose re) (transpose im))
  in Mat.save_txt ~out:(in_dir "test_w") new_w; 
  Mat.save_txt ~out:(in_dir "test_w_eigs") (eigenvalues new_w)
(* bisection : f a is always supposed to have negative values and f b always postivie *)
(* let n_iter = 20

let opt_alpha = 
  let h2_tgt = AD.unpack_flt h2_tgt in 
  let w0 = AD.unpack_arr w0 in 
  let _ = Printf.printf "H_tgt = %f \n %!" h2_tgt in 
  let f alpha = 
  let w = Mat.(alpha $* w0) in 
  let h2norm = (h2norm (AD.pack_arr w) c0) |> AD.unpack_flt in Float.((h2norm -. h2_tgt))
  in 
let rec bisection a b k = 
  let c = Float.((a +. b)/.2.) in 
  if k > n_iter then f a, f b else
  let _ = Printf.printf "|| Diff = %f %f %i \n ||%!" (f a -. f 0.) (f b -. f 0.) k in 
  if f c < 0. then bisection c b (succ k) else bisection c b (succ k)
in bisection 0.1 20. 0

let fa, fb = opt_alpha
let _ = Printf.printf "%f %f %!" fa fb *)
