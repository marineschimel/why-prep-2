open Owl
module AD = Algodiff.D
open Lib
open Defaults

let _ = Printexc.record_backtrace true

(*we take the nonlin and the w from the Defaults file *)

let f x = Defaults.g x
let w_rec = Defaults.w
let __w = Defaults.__w
let tau = AD.F Defaults.tau
let dt = AD.F Defaults.sampling_dt

(*find candidate fixed points for the networks by running the dynamics many times*)

(*then check whether they are fixed points or very slow points?

Use the routine from Sussillo : 
Add noise to the fixed point candidates ('noise_var')
    Optimize to find the closest fixed points / slow points (many hps, 
      see optimize_fps)
    Exclude any fixed points whose fixed point loss is above threshold ('fp_tol')
    Exclude any non-unique fixed points according to a tolerance ('unique_tol')
    Exclude any far-away "outlier" fixed points ('outlier_tol')*)
let run_net n_steps x0 =
  let _ = Printf.printf "Running the dynamics %! \n" in
  let rec run_dyn k x =
    if k = n_steps
    then x
    else (
      let dx = AD.Maths.(((__w *@ f x) - x) / tau * dt) in
      run_dyn (k + 1) AD.Maths.(x + dx))
  in
  run_dyn 0 x0

let init_points = Mat.gaussian 200 10 (*result from running the dynamics many times *)

let h0s = Mat.map_cols (fun x -> run_net 1 (AD.pack_arr x)) init_points
let _ = Printf.printf "size h0s %i %!" (Array.length h0s)

module Prms = struct
  type 'a t = { h : 'a } [@@deriving prms]
end

open Prms
module P = Owl_opt_lbfgs.Make (Prms)

let cost prms = AD.Maths.(l2norm_sqr' (prms.h - (__w *@ f prms.h)))
let sig_noise = 0.001
let max_iter = 1000
let tolerance = 1E-6
let fixed_points = []

let learn h0 =
  let _ = AD.Mat.print h0 in
  let noise = AD.Mat.(gaussian ~sigma:sig_noise 200 1) in
  let prms0 = { h = AD.Maths.(h0 + noise) } in
  let stop =
    let fv_prev = ref 1E9 in
    fun s ->
      let k = P.iter s in
      let fv = P.(fv s) in
      let pct_change = abs_float ((fv -. !fv_prev) /. !fv_prev) in
      fv_prev := fv;
      if k mod 10 = 0
      then Printf.printf "iter %i | fv %f | pct change %f %!" k fv pct_change;
      k > max_iter || fv < tolerance
  in
  let f prms = cost prms in
  let s0 = P.init ~f ~prms0 () in
  let sf = P.min ~stop s0 in
  let prms = P.prms sf in
  let f_cost = cost prms in
  if AD.unpack_flt f_cost < tolerance then Some (AD.unpack_arr prms.h) else None

let final_fixed_points =
  let x, _ =
    Array.fold_left
      (fun (acc, k) x ->
        let _ = Printf.printf "num %i %!" k in
        match learn x with
        | Some a -> acc @ [ a ], k + 1
        | None -> acc, k + 1)
      ([], 0)
      h0s
  in
  Printf.printf "%i %!" (List.length x);
  Mat.concatenate ~axis:1 (Array.of_list x)

let _ = Mat.save_txt ~out:"ffpoints" final_fixed_points
