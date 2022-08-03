open Owl
open Lib
open Defaults
module A = Analysis

let _ = Printexc.record_backtrace true
let _ = Mat.save_txt ~out:"data/w_disc" (Mat.zeros 200 200)

let in_dir f s =
  Printf.sprintf "results_c_alpha/skew/alpha_%i/%s" (int_of_float (f *. 1000.)) s
(* Printf.sprintf "results_small/big_soc/alpha_%i/%s" (int_of_float (f *. 1000.)) s *)

module Z = Dense.Matrix.Z

let _ =
  Array.map
    (fun psi ->
      let losses = Mat.load_txt (in_dir psi "loss_time") in
      let loss_0 = Mat.get losses 0 2 in
      let loss_end = Mat.get losses (Mat.row_num losses - 1) 2 in
      let times = Mat.get_slice [ []; [ 0 ] ] losses in
      let norm_losses = Mat.(get_slice [ []; [ 2 ] ] losses /$ loss_0) in
      let nonnormality = Mat.get losses 0 3 in
      let w =
        let x = Mat.gaussian size_net size_net ~sigma:(psi /. sqrt (float size_net)) in
        Mat.(x - transpose x)
      in
      let a = Mat.((w - eye size_net) /$ tau) in
      let _, eigs = Linalg.D.eig a in
      let im, re = Z.im eigs, Z.re eigs in
      let max_im = Mat.max' im in
      let max_re = Mat.max' re in
      let c = Mat.load_txt (in_dir psi "c") in
      (* let c = Mat.load_txt (in_dir psi "new_c") in  *)
      let q = Linalg.D.lyapunov Mat.(transpose a) Mat.(neg (transpose c) *@ c) in
      let frob_norm = Mat.l2norm' w in
      let h2norm = Mat.(trace q) in
      let input_ratio =
        let inputs = Mat.load_txt (in_dir psi "results_us_300") in
        let ip, im = A.inputs ~t_prep:0.3 inputs in
        ip /. im
      in
      Mat.save_txt ~out:(in_dir psi "norm_losses") Mat.(times @|| norm_losses);
      [| psi
       ; (loss_0 /. loss_end) -. 1.
       ; nonnormality
       ; h2norm
       ; input_ratio
       ; frob_norm
       ; max_im
       ; max_re
       ; float (Owl.Linalg.D.rank w)
      |])
    [| 0.01; 0.03; 0.05; 0.07; 0.12; 0.15; 0.2 |]
  |> Mat.of_arrays
  |> Mat.save_txt ~out:"results_c_alpha/skew/prep_amount"

let _ =
  Array.map
    (fun psi ->
      let w =
        let x = Mat.gaussian size_net size_net ~sigma:(psi /. sqrt (float size_net)) in
        Mat.(x - transpose x)
      in
      let a = Mat.((w - eye size_net) /$ tau) in
      let _, eigs = Linalg.D.eig a in
      let im, re = Z.im eigs, Z.re eigs in
      let max_im = Mat.max' im in
      let max_re = Mat.max' re in
      let schur, _, _ = Linalg.D.schur a in
      Mat.save_txt ~out:(in_dir psi "eigs") Mat.(transpose re @|| transpose im);
      Mat.save_txt ~out:(in_dir psi "schur") schur;
      [| psi; max_im; max_re; nonnormality w |])
    [| 0.01; 0.03; 0.05; 0.07; 0.12; 0.15; 0.2 |]
  |> Mat.of_arrays
  |> Mat.save_txt ~out:"results_c_alpha/skew/data"

(* let ctrl =
  let schur, _, _ =
    let inv_ctrl = Linalg.D.linsolve ctrl_gramian (Mat.eye n) in
    let _ = Mat.save_txt ~out:"skew_term_wrec" Mat.(skew_term w_rec *@ inv_ctrl) in
    Linalg.D.schur Mat.(skew_term w_rec *@ inv_ctrl)
  in
  Printf.printf "max nonnormality %f %!" (nonnormality schur) *)

(* 
      let inv_ctrl = Linalg.D.linsolve ctrl_gramian (Mat.eye n) in
      let _ = Mat.save_txt ~out:"skew_term_wrec" Mat.(skew_term w_rec *@ inv_ctrl) in
      Linalg.D.schur Mat.(skew_term w_rec *@ inv_ctrl)
    in
    Printf.printf "max nonnormality %f %!" (nonnormality schur) *)
