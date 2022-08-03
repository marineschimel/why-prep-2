open Owl
module AD = Algodiff.D
open Lib

let _ = Printexc.record_backtrace true
let t_prep = 500
let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")
let in_dir seed s = Printf.sprintf "%s/seed_%i/%s" dir seed s
let n_dim = 8
let n_reaches = 8
let n_var = 8
let dt = 2E-3
let reach_0 = Mat.load_txt (in_dir 1 Printf.(sprintf "rates_%i_%i" 0 t_prep))
let size_prep = 200
let size_mov = 200
let duration = Mat.row_num reach_0
let n_prep = float_of_int t_prep /. 1000. /. dt |> int_of_float

let softmax x =
  let m = Mat.(max ~axis:0 x - min ~axis:0 x +$ (0.1 *. max' x)) in
  Mat.(x / m)

let preprocessed_data seed =
  let load i =
    let m = Mat.load_txt (in_dir seed Printf.(sprintf "rates_%i_%i" i t_prep)) in
    let ma = Mat.get_fancy [ R [ 0; duration - 1 ]; R [ 0; -1 ] ] m in
    Arr.reshape ma [| duration; -1; 1 |]
  in
  let data = Array.init n_reaches (fun i -> load i) in
  let concat_data = Arr.concatenate ~axis:0 data in
  let norm_factor =
    Mat.(max ~axis:0 concat_data - min ~axis:0 concat_data +$ (0.1 *. max' concat_data))
    |> fun z -> Arr.reshape z [| 1; -1; 1 |]
  in
  let norm_data = Array.map (fun x -> Arr.(x / norm_factor)) data in
  let mean_data = Arr.concatenate ~axis:2 norm_data |> Arr.mean ~axis:2 in
  let dat =
    let a =
      Arr.concatenate ~axis:0 (Array.map (fun x -> Arr.(x - mean_data)) norm_data)
    in
    Arr.reshape a [| (Arr.shape a).(0); (Arr.shape a).(1) |]
  in
  mean_data, Mat.of_arrays (Arr.to_arrays dat)

let x_prep seed =
  let m = Arr.reshape (snd (preprocessed_data seed)) [| n_reaches; duration; -1 |] in
  let f i =
    Arr.reshape
      (Arr.get_slice [ [ i ]; [ n_prep - size_prep; n_prep - 1 ]; [] ] m)
      [| size_prep; -1 |]
  in
  Array.init n_reaches f |> Arr.concatenate ~axis:0

let x_mov seed =
  let m = Arr.reshape (snd (preprocessed_data seed)) [| n_reaches; duration; -1 |] in
  let f i =
    Arr.reshape
      (Arr.get_slice [ [ i ]; [ n_prep + 25; n_prep + size_mov - 1 + 25 ]; [] ] m)
      [| size_mov; -1 |]
  in
  Array.init n_reaches f |> Arr.concatenate ~axis:0

let n = 200

let cov_p seed =
  let r = Mat.(x_prep seed - mean ~axis:0 (x_prep seed)) in
  Mat.(transpose r *@ r)

let cov_m seed =
  let r = Mat.(x_mov seed - mean ~axis:0 (x_mov seed)) in
  Mat.(transpose r *@ r)

module Prms = struct
  type 'a t = { w : 'a } [@@deriving prms]
end

open Prms
module P = Owl_opt_lbfgs.Make (Prms)

let cost seed prms =
  let cov_p = cov_p seed in
  let cov_m = cov_m seed in
  let __cp, __cm = AD.pack_arr cov_p, AD.pack_arr cov_m in
  let _, sprep, _ = Linalg.D.svd cov_p in
  let _, smov, _ = Linalg.D.svd cov_m in
  let sp = Mat.sum' (Mat.get_slice [ []; [ 0; n_dim - 1 ] ] sprep)
  and sm = Mat.sum' (Mat.get_slice [ []; [ 0; n_dim - 1 ] ] smov) in
  let orth_modes, _ = AD.Linalg.qr prms.w in
  let new_wp = AD.Maths.get_slice [ []; [ 0; n_dim - 1 ] ] orth_modes
  and new_wm = AD.Maths.get_slice [ []; [ n_dim; -1 ] ] orth_modes in
  let obj =
    AD.Maths.(
      (trace (transpose new_wp *@ __cp *@ new_wp) / F sp)
      + (trace (transpose new_wm *@ __cm *@ new_wm) / F sm))
  in
  AD.Maths.(neg obj), new_wp, new_wm

let learn seed =
  let prms0 = { w = AD.Mat.gaussian n (2 * n_dim) } in
  let stop =
    let fv_prev = ref 1E9 in
    fun s ->
      let k = P.iter s in
      let fv = P.(fv s) in
      let pct_change = abs_float ((fv -. !fv_prev) /. !fv_prev) in
      fv_prev := fv;
      Printf.printf "\r iter %i | fv %f | pct change %f %!" k fv pct_change;
      pct_change < 1E-4
  in
  let f prms =
    let f, _, _ = cost seed prms in
    f
  in
  let s0 = P.init ~f ~prms0 () in
  let sf = P.min ~stop s0 in
  let prms = P.prms sf in
  let _, wp, wm = cost seed prms in
  prms, wp, wm

let ws seed =
  let _, new_wp, new_wm = learn seed in
  new_wp |> AD.unpack_arr, new_wm |> AD.unpack_arr

let captured_variance seed =
  let wp, wm = ws seed in
  let cov_p = cov_p seed in
  let cov_m = cov_m seed in
  let _mp, sp, _ = Linalg.D.svd (Mat.transpose (x_prep seed)) in
  let _mv, sm, _ = Linalg.D.svd (Mat.transpose (x_mov seed)) in
  let sp, sm = Mat.sqr sp, Mat.sqr sm in
  let _, _ = Mat.(cumsum ~axis:1 sp /$ sum' sp), Mat.(cumsum ~axis:1 sm /$ sum' sm) in
  let top_mp, top_mv =
    Mat.get_slice [ []; [ 0; n_var - 1 ] ] wp, Mat.get_slice [ []; [ 0; n_var - 1 ] ] wm
  in
  let m = Mat.(transpose top_mp *@ cov_p *@ top_mp) in
  let x = Mat.(x_mov seed *@ top_mp *@ transpose top_mp) in
  let y = Mat.(x_prep seed *@ top_mv *@ transpose top_mv) in
  let z = Mat.(transpose top_mv *@ cov_p *@ top_mv) in
  let zz = Mat.(transpose top_mp *@ cov_m *@ top_mp) in
  let xcent = Mat.(x - mean ~axis:0 x) in
  let ycent = Mat.(y - mean ~axis:0 y) in
  let _, ss, _ = Linalg.D.svd xcent in
  let _, tt, _ = Linalg.D.svd ycent in
  let _ =
    Printf.printf
      "pct_var %f %f %f %!"
      (Mat.trace z /. Mat.sum' (Mat.get_slice [ []; [ 0; n_var - 1 ] ] sp))
      (Mat.trace zz /. Mat.sum' (Mat.get_slice [ []; [ 0; n_var - 1 ] ] sm))
      (Mat.trace m /. Mat.sum' (Mat.get_slice [ []; [ 0; n_var - 1 ] ] sp))
  in
  ( Mat.save_txt
      ~out:(in_dir seed (Printf.sprintf "analysis/pct_var_%i" seed))
      Mat.(
        (cumsum ~axis:1 (Mat.sqr tt) /$ sum' sp)
        @= (cumsum ~axis:1 (Mat.sqr ss) /$ sum' sm))
  , top_mp
  , top_mv )

let proj x modes =
  (*x is T*N and modes is N*n_dim*)
  let proj = Mat.(x *@ modes) in
  proj

let reconstructed proj modes = Mat.(proj *@ transpose modes)

let occupancy seed =
  let load i =
    Mat.get_slice
      [ [ i * duration; ((i + 1) * duration) - 1 ]; [] ]
      (snd (preprocessed_data seed))
    (* let m = Mat.load_txt (in_dir Printf.(sprintf "rates_%i_%i" i t_prep)) in
    let ma = Mat.get_fancy [ R [ 0; duration - 1 ]; R [ 0; -1 ] ] m in
    ma *)
  in
  let wp, wm = ws seed in
  let projections x =
    let rp = reconstructed (proj x wp) wp in
    let rm = reconstructed (proj x wm) wm in
    Mat.(rp + rm), rp, rm, x
  in
  let rps =
    Arr.concatenate
      ~axis:2
      (Array.init n_reaches (fun i ->
           let _, rp, _, _ = projections (load i) in
           let m = Arr.reshape rp [| (Arr.shape rp).(0); (Arr.shape rp).(1); 1 |] in
           m))
  and rms =
    Arr.concatenate
      ~axis:2
      (Array.init n_reaches (fun i ->
           let _, _, rm, _ = projections (load i) in
           let m = Arr.reshape rm [| (Arr.shape rm).(0); (Arr.shape rm).(1); 1 |] in
           m))
  in
  let vprep, vmov =
    let x = Arr.var ~axis:2 rps |> Arr.sum ~axis:1 |> Arr.to_array in
    ( Mat.of_array x (-1) 1
    , let x = Arr.var ~axis:2 rms |> Arr.sum ~axis:1 |> Arr.to_array in
      Mat.of_array x (-1) 1 )
  in
  Mat.(vprep /$ max' vprep), Mat.(vmov /$ max' vmov)

let _ =
  let a = Array.init 4 (fun seed -> occupancy (succ seed)) in
  let vprep_mean = a |> Array.map fst |> Mat.concatenate ~axis:1 |> Mat.mean ~axis:1 in
  let vmov_mean = a |> Array.map snd |> Mat.concatenate ~axis:1 |> Mat.mean ~axis:1 in
  let vprep_var = a |> Array.map fst |> Mat.concatenate ~axis:1 |> Mat.var ~axis:1 in
  let vmov_var = a |> Array.map snd |> Mat.concatenate ~axis:1 |> Mat.var ~axis:1 in
  Mat.save_txt ~out:(Printf.sprintf "%s/analysis/vmov" dir) Mat.(vmov_mean @|| vmov_var);
  Mat.save_txt
    ~out:(Printf.sprintf "%s/analysis/vprep" dir)
    Mat.(vprep_mean @|| vprep_var)
