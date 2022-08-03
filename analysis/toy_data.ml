open Owl
open Lib
open Defaults

let dims = 16
let n_reaches = 8
let t_prep = 500
let size_prep = 400
let size_mov = 400
let duration = 1300
let n_dim = 8

let init_prep =
  Array.init n_reaches (fun i ->
      let o = Mat.gaussian ~mu:(float_of_int i) 500 1 in
      let m =
        let m0 = Mat.zeros 500 n_reaches in
        let _ = Mat.set_slice [ []; [ i ] ] m0 o in
        m0
      and z = Mat.zeros 500 n_reaches in
      Arr.reshape Mat.(m @|| z) [| 500; dims; 1 |])
  |> Arr.concatenate ~axis:2

let init_mov =
  Array.init n_reaches (fun i ->
      let o = Mat.gaussian ~mu:(float_of_int i) 500 1 in
      let m =
        let m0 = Mat.zeros 500 8 in
        Mat.set_slice [ []; [ i ] ] m0 o;
        m0
      and z = Mat.zeros 500 8 in
      Arr.reshape Mat.(z @|| m) [| 500; dims; 1 |])
  |> Arr.concatenate ~axis:2

let data = Arr.(rotate (init_prep @= init_mov) 180)

let softmax x =
  let m = Mat.(std ~axis:0 x +$ 0.5) in
  Mat.(x / m)

let x_prep =
  let a =
    let x =
      Array.init n_reaches (fun i ->
          Mat.get_slice [ [ 50; 50 + size_prep - 1 ]; []; [ i ] ] data)
    in
    Arr.concatenate ~axis:0 x
  in
  let shp = Arr.shape a in
  Arr.reshape a [| shp.(0); shp.(1) |] |> Arr.to_arrays |> Mat.of_arrays

let x_mov =
  let a =
    let x =
      Array.init n_reaches (fun i ->
          Mat.get_slice [ [ t_prep; t_prep + size_mov ]; []; [ i ] ] data)
    in
    Arr.concatenate ~axis:0 x
  in
  let shp = Arr.shape a in
  Arr.reshape a [| shp.(0); shp.(1) |] |> Arr.to_arrays |> Mat.of_arrays

let _ = Mat.save_txt ~out:"x" Mat.(x_prep @= x_mov)
let n = Mat.col_num x_mov

let cov_p =
  let r = Mat.(x_prep - mean ~axis:0 x_prep) in
  Mat.(transpose r *@ r)

let cov_m =
  let r = Mat.(x_mov - mean ~axis:0 x_mov) in
  Mat.(transpose r *@ r)

let __cp, __cm = AD.pack_arr cov_p, AD.pack_arr cov_m

(*Elsayed et al analysis : are the subspaces orthogonal?*)

let corr_p =
  let r = Mat.(x_prep - mean ~axis:0 x_prep) in
  let rr = Mat.(r / l2norm ~axis:0 r) in
  Mat.(transpose rr *@ rr)

let corr_m =
  let r = Mat.(x_mov - mean ~axis:0 x_mov) in
  let rr = Mat.(r / l2norm ~axis:0 r) in
  Mat.(transpose rr *@ rr)

let _ = Mat.save_txt ~out:"test_cp" corr_p, Mat.save_txt ~out:"test_cm" corr_m

let pairwise_corr =
  let cp = Mat.reshape corr_p [| n * n; 1 |]
  and cm = Mat.reshape corr_m [| n * n; 1 |] in
  Mat.(cp @|| cm) |> Mat.save_txt ~out:"test_pairwise"

module Prms = struct
  type 'a t = { w : 'a } [@@deriving prms]
end

open Prms
module P = Owl_opt_lbfgs.Make (Prms)

let cost prms =
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
  AD.Maths.(neg obj)

let learn =
  let prms0 = { w = AD.Mat.gaussian n (2 * n_dim) } in
  let stop =
    let fv_prev = ref 1E9 in
    fun s ->
      let k = P.iter s in
      let fv = P.(fv s) in
      let pct_change = abs_float ((fv -. !fv_prev) /. !fv_prev) in
      fv_prev := fv;
      Printf.printf "\r iter %i | fv %f | pct change %f %!" k fv pct_change;
      pct_change < 1E-8
  in
  let f prms = cost prms in
  let s0 = P.init ~f ~prms0 () in
  let sf = P.min ~stop s0 in
  let prms = P.prms sf in
  prms

let wp, wm =
  let prms = learn in
  let orth_modes, _ = AD.Linalg.qr prms.w in
  let new_wp = AD.Maths.get_slice [ []; [ 0; n_dim - 1 ] ] orth_modes
  and new_wm = AD.Maths.get_slice [ []; [ n_dim; -1 ] ] orth_modes in
  new_wp |> AD.unpack_arr, new_wm |> AD.unpack_arr

let _ =
  Mat.print (Mat.l2norm ~axis:0 wm);
  Mat.print (Mat.l2norm ~axis:0 wp)

let proj x modes =
  (*x is T*N and modes is N*n_dim*)
  let proj = Mat.(x *@ modes) in
  proj

(*proj is T*n_dim*)
let reconstructed proj modes = Mat.(proj *@ transpose modes)

(*T*N*)
let proj2d i =
  let x =
    let a =
      Arr.reshape (Arr.get_slice [ []; []; [ i - 1 ] ] data) [| 1000; dims |]
      |> Arr.to_arrays
    in
    Mat.of_arrays a
  in
  let proj_mov = reconstructed (proj x wm) wm in
  let proj_prep = reconstructed (proj x wp) wp in
  Mat.save_txt ~out:(Printf.sprintf "test_proj_prep_%i" i) proj_prep;
  Mat.save_txt ~out:(Printf.sprintf "test_proj_mov_%i" i) proj_mov

let captured_variance, top_mp, top_mv =
  let mp, sp, _ = Linalg.D.svd cov_p in
  let mv, sm, _ = Linalg.D.svd cov_m in
  let pct_varp, pct_varm =
    Mat.(cumsum ~axis:1 sp /$ sum' sp), Mat.(cumsum ~axis:1 sm /$ sum' sm)
  in
  let top_mp, top_mv =
    Mat.get_slice [ []; [ 0; 4 - 1 ] ] mp, Mat.get_slice [ []; [ 0; 4 - 1 ] ] mv
  in
  Mat.save_txt ~out:"pct_var_own" Mat.(pct_varp @= pct_varm);
  let x = Mat.(x_mov *@ top_mp *@ transpose top_mp) in
  let y = Mat.(x_prep *@ top_mv *@ transpose top_mv) in
  let z = Mat.(transpose top_mv *@ cov_p *@ top_mv) in
  let zz = Mat.(transpose top_mp *@ cov_m *@ top_mp) in
  let xcent = Mat.(x - mean ~axis:0 x) in
  let ycent = Mat.(y - mean ~axis:0 y) in
  let _, ss, _ = Linalg.D.svd xcent in
  let _, tt, _ = Linalg.D.svd ycent in
  let _ =
    Printf.printf
      "pct_var %f %f %!"
      (Mat.trace z /. Mat.sum' sp)
      (Mat.trace zz /. Mat.sum' sm)
  in
  ( Mat.save_txt
      ~out:"test_pct_var"
      Mat.(
        (cumsum ~axis:1 (Mat.sqr tt) /$ sum' sp)
        @= (cumsum ~axis:1 (Mat.sqr ss) /$ sum' sm))
  , top_mp
  , top_mv )

let occupancy =
  let rs i =
    let x =
      let a =
        Arr.reshape (Arr.get_slice [ []; []; [ i - 1 ] ] data) [| 1000; dims |]
        |> Arr.to_arrays
      in
      Mat.of_arrays a
    in
    reconstructed (proj x wp) wp, reconstructed (proj x wm) wm
  in
  let rps =
    Arr.concatenate
      ~axis:2
      (Array.init n_reaches (fun i ->
           let rp, _ = rs i in
           let m = Arr.reshape rp [| (Arr.shape rp).(0); (Arr.shape rp).(1); 1 |] in
           m))
  and rms =
    Arr.concatenate
      ~axis:2
      (Array.init n_reaches (fun i ->
           let _, rm = rs i in
           let m = Arr.reshape rm [| (Arr.shape rm).(0); (Arr.shape rm).(1); 1 |] in
           m))
  in
  let vprep, vmov =
    let x = Arr.var ~axis:2 rps |> Arr.sum ~axis:1 |> Arr.to_array in
    ( Mat.of_array x (-1) 1
    , let x = Arr.var ~axis:2 rms |> Arr.sum ~axis:1 |> Arr.to_array in
      Mat.of_array x (-1) 1 )
  in
  ( Mat.save_txt ~out:"vprep" Mat.(vprep /$ max' vprep)
  , Mat.save_txt ~out:"vmov" Mat.(vmov /$ max' vmov) )

let _ =
  Array.init n_reaches (fun i ->
      proj2d (succ i);
      0.)
