open Owl
module AD = Algodiff.D
open Lib
open Defaults
open Typ
open Misc

(*for each reach we give a target-specific input at the beginning on prep, then we train
the RNN weights for the input (not at first but after) + the connectivity w with a fixed
readout  )*)

let _ = Printexc.record_backtrace true
let __c = Defaults.__c

module Prms = struct
  type 'a t =
    { x0 : 'a
    ; w : 'a
    }
  [@@deriving prms]
end

open Prms
module P = Owl_opt_lbfgs.Make (Prms)

let n_reaches = 7
let rdm = AD.Mat.gaussian ~sigma:(0.9 /. Maths.sqrt (float n)) n n

let trajectory ~n_steps ~dt ~w ~x0 =
  let a = AD.Maths.((w - AD.Mat.eye n) / F tau) in
  let rec accu k x xs =
    if k = n_steps
    then xs |> List.rev |> Array.of_list
    else (
      let x =
        let dx = AD.Maths.(a *@ x * F dt) in
        AD.Maths.(x + dx)
      in
      let xs = AD.Maths.transpose x :: xs in
      accu (k + 1) x xs)
  in
  accu 0 x0 []


let target_torques =
  Array.init n_reaches (fun i ->
      Mat.(load_txt (Printf.sprintf "results/soc/reach_%i/torques_0" (succ i))))


let packed_torques = Array.map (fun x -> AD.pack_arr x) target_torques

(*Misc.read_bin (in_dir "target_torques.bin")*)
let n_steps = Mat.row_num target_torques.(0)

(* let n_steps, n_torques = Mat.shape target_torques.(0) *)
(*let n_reaches = Array.length target_torques*)

let traj w x0 =
  AD.Maths.concatenate
    ~axis:0
    (Array.init n_reaches (fun i ->
         AD.Maths.concatenate
           ~axis:0
           (trajectory
              ~w
              ~n_steps
              ~dt:sampling_dt
              ~x0:(AD.Maths.get_slice [ []; [ i ] ] x0))))


let cost target_torques =
  let norm = Mat.of_arrays [| [| 1.; 3. |] |] |> Mat.transpose |> AD.pack_arr in
  fun prms ->
    let w = prms.w in
    let x0 = prms.x0 in
    let xs = traj w x0 in
    let ys = AD.Maths.(xs *@ transpose __c) in
    let t = AD.Maths.concatenate ~axis:0 target_torques in
    let dy = AD.Maths.(transpose (ys - t)) in
    let loss =
      AD.Maths.(
        (l2norm_sqr' (dy / norm) / F (float n_steps))
        + (F 1. / F (float n *. float n_reaches) * l2norm_sqr' prms.x0)
        + (F 1. / F (float n *. float n) * l2norm_sqr' prms.w))
    in
    loss


let in_dir s = in_dir Printf.(sprintf "end_to_end/%s" s)

let save_prms prms =
  let w = AD.unpack_arr prms.w in
  Mat.save_txt w ~out:"results/end_to_end/w";
  let x0 = AD.unpack_arr prms.x0 in
  Mat.save_txt x0 ~out:"results/end_to_end/x0"


let save_traj prms =
  let _ =
    Array.init n_reaches (fun i ->
        let xi = AD.unpack_arr (traj prms.w prms.x0) in
        let x = Mat.get_slice [ [ n_steps * i; (n_steps * i) + n_steps - 1 ]; [] ] xi in
        let torques = Mat.(x *@ transpose c) in
        let thetas = M.theta_trajectory ~dt:sampling_dt torques in
        let hands = Array.map M.hand_of thetas in
        Mat.save_txt x ~out:(in_dir Printf.(sprintf "x_%i" (succ i)));
        Mat.save_txt ~out:(in_dir "x0") (AD.unpack_arr prms.x0);
        Mat.save_txt torques ~out:(in_dir Printf.(sprintf "torques_%i" (succ i)));
        Mat.save_txt
          (Arm.unpack_sequence thetas)
          ~out:(in_dir Printf.(sprintf "thetas_%i" 1));
        Mat.save_txt
          (Arm.unpack_sequence hands)
          ~out:(in_dir Printf.(sprintf "hands_%i" (succ i))))
  in
  ()


let learn ~torques =
  let prms0 =
    { x0 = AD.Mat.gaussian ~sigma:0.001 n n_reaches
    ; w = AD.Mat.gaussian ~sigma:0.001 n n
    }
  in
  let stop =
    let fv_prev = ref 1E9 in
    fun s ->
      let k = P.iter s in
      let fv = P.(fv s) in
      let pct_change = abs_float (fv -. !fv_prev) /. !fv_prev in
      fv_prev := fv;
      Printf.printf "\riter %i | fv %f | pct change %f %!" k fv pct_change;
      if k mod 10 = 0
      then (
        let prms = P.prms s in
        save_traj prms;
        save_prms prms);
      fv < 0.0016
  in
  let f =
    let cost = cost torques in
    fun prms -> cost prms
  in
  let s0 = P.init ~f ~prms0 () in
  let sf = P.min ~stop s0 in
  let prms = P.prms sf in
  save_traj prms;
  save_prms prms;
  prms


let _ = learn ~torques:packed_torques
