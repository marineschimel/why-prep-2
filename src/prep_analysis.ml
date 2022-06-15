open Owl
module AD = Algodiff.D
open Lib

let _ = Printexc.record_backtrace true
let t_prep = Cmdargs.(get_int "-prep" |> force ~usage:"-prep")
let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")
let in_dir s = Printf.sprintf "%s/%s" dir s
let n_dim = 8
let reaches = [|1;2;3;4;6|]
let n_reaches = (Array.length reaches)
let n_var = 8
let dt = 2E-3
let reach_0 = Mat.load_txt (in_dir Printf.(sprintf "rates_%i_%i" 1 t_prep))
let size_prep = 100
let size_mov = 100
let duration = Mat.row_num reach_0
let n_prep = float_of_int t_prep /. 1000. /. dt |> int_of_float

let softmax x =
  let m = Mat.(max ~axis:0 x - min ~axis:0 x +$ (0.1 *. max' x)) in
  Mat.(x / m)


let preprocessed_data =
  let load i =
    let m = Mat.load_txt (in_dir Printf.(sprintf "rates_%i_%i" i t_prep)) in
    let ma = Mat.get_fancy [ R [ 0; duration - 1 ]; R [ 0; -1 ] ] m in
    Arr.reshape ma [| duration; -1; 1 |]
  in
  let data = Array.map (fun i -> load i) reaches  in
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


let x_prep =
  let m = Arr.reshape (snd preprocessed_data) [| n_reaches; duration; -1 |] in
  let f i _ =
    let _ = Stdio.printf "%i %i %i %i %i %i %!" n_prep size_prep (Array.length reaches) (Arr.shape m).(0) (Arr.shape m).(1) (Arr.shape m).(2) in
    Arr.reshape
      (Arr.get_slice [ [ i ]; [ n_prep - size_prep; n_prep - 1 ]; [] ] m)
      [| size_prep; -1 |]
  in
  Array.mapi f reaches |> Arr.concatenate ~axis:0


let x_mov =
  let m = Arr.reshape (snd preprocessed_data) [| n_reaches; duration; -1 |] in
  let f i _ =
    Arr.reshape
      (Arr.get_slice [ [ i ]; [ n_prep + 25; n_prep + size_mov - 1 + 25 ]; [] ] m)
      [| size_mov; -1 |]
  in
  Array.mapi f reaches |> Arr.concatenate ~axis:0


let n = Mat.col_num x_mov

let cov_p =
  let r = Mat.(x_prep - mean ~axis:0 x_prep) in
  Mat.(transpose r *@ r)


let cov_m =
  let r = Mat.(x_mov - mean ~axis:0 x_mov) in
  Mat.(transpose r *@ r)


let __cp, __cm = AD.pack_arr cov_p, AD.pack_arr cov_m

let capt_var dim =
  let cp, varp, _ = Linalg.D.svd cov_p in
  let cm, varm, _ = Linalg.D.svd cov_m in
  let cm = Mat.get_slice [ []; [ 0; dim - 1 ] ] cm in
  let cp = Mat.get_slice [ []; [ 0; dim - 1 ] ] cp in
  let _ = Printf.printf "%i %i %!" (Mat.row_num cp) (Mat.col_num cp) in
  let var_p_mov = Mat.(trace Mat.(transpose cm *@ cov_p *@ cm)) /. Mat.sum' varp in
  let var_m_prep = Mat.(trace Mat.(transpose cp *@ cov_m *@ cp)) /. Mat.sum' varm in
  let var_p_prep = Mat.(trace Mat.(transpose cp *@ cov_p *@ cp)) /. Mat.sum' varp in
  let var_m_mov = Mat.(trace Mat.(transpose cm *@ cov_m *@ cm)) /. Mat.sum' varm in
  [| var_p_mov; var_m_prep; var_p_prep; var_m_mov |]


(*this saves : var of prep in move subspace, var move in prep subspace, var prep in prep subspace, var move in move subspace
  for various dimensions*)
let _ =
  Array.map (fun i -> capt_var i) reaches 
  |> Mat.of_arrays
  |> Mat.save_txt ~out:(in_dir "analysis/capt_var")


(*below, code to check how orthogonal the activity at the end of prep is to the rest of the movement*)

(*Elsayed et al analysis : are the subspaces orthogonal?*)

let corr_p =
  let r = Mat.(x_prep - mean ~axis:0 x_prep) in
  let rr = Mat.(r / l2norm ~axis:0 r) in
  Mat.(transpose rr *@ rr)


let corr_m =
  let r = Mat.(x_mov - mean ~axis:0 x_mov) in
  let rr = Mat.(r / l2norm ~axis:0 r) in
  Mat.(transpose rr *@ rr)


let _ =
  ( Mat.save_txt ~out:(in_dir "analysis/cp") corr_p
  , Mat.save_txt ~out:(in_dir "analysis/cm") corr_m )


let pairwise_corr =
  let cp = Mat.reshape corr_p [| n * n; 1 |]
  and cm = Mat.reshape corr_m [| n * n; 1 |] in
  Mat.(cp @|| cm) |> Mat.save_txt ~out:(in_dir "analysis/pairwise")


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
  AD.Maths.(neg obj), new_wp, new_wm


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
      pct_change < 1E-4
  in
  let f prms =
    let f, _, _ = cost prms in
    f
  in
  let s0 = P.init ~f ~prms0 () in
  let sf = P.min ~stop s0 in
  let prms = P.prms sf in
  let _, wp, wm = cost prms in
  prms, wp, wm


let wp, wm =
  let _, new_wp, new_wm = learn in
  new_wp |> AD.unpack_arr, new_wm |> AD.unpack_arr

let _ =
  Mat.save_txt ~out:(in_dir (Printf.sprintf "analysis/wp" )) wp;
  Mat.save_txt ~out:(in_dir (Printf.sprintf "analysis/wm" )) wm


let proj x modes =
  (*x is T*N and modes is N*n_dim*)
  let proj = Mat.(x *@ modes) in
  proj


(*proj is T*n_dim*)
let reconstructed proj modes = Mat.(proj *@ transpose modes)

(*T*N*)
let proj2d =
  let m = Arr.reshape (snd preprocessed_data) [| n_reaches; duration; -1 |] in
  fun i ->
    let x =
      let a = Arr.reshape (Arr.get_slice [ [ i ]; []; [] ] m) [| duration; -1 |] in
      a
    in
    let proj_prep = proj x wp in
    let proj_mov = proj x wm in 
    Mat.save_txt ~out:(in_dir (Printf.sprintf "analysis/proj_prep_%i" i)) proj_prep;
    Mat.save_txt ~out:(in_dir (Printf.sprintf "analysis/proj_mov_%i" i)) proj_mov


let _ = Array.mapi (fun i _ -> proj2d i) reaches

let captured_variance, top_mp, top_mv =
  let _mp, sp, _ = Linalg.D.svd (Mat.transpose x_prep) in
  let _mv, sm, _ = Linalg.D.svd (Mat.transpose x_mov) in
  let _ = Mat.save_txt ~out:"sm" Mat.(transpose (sqr sm)) in
  let sp, sm = Mat.sqr sp, Mat.sqr sm in
  let pct_varp, pct_varm =
    Mat.(cumsum ~axis:1 sp /$ sum' sp), Mat.(cumsum ~axis:1 sm /$ sum' sm)
  in
  let top_mp, top_mv =
    Mat.get_slice [ []; [ 0; n_var - 1 ] ] wp, Mat.get_slice [ []; [ 0; n_var - 1 ] ] wm
  in
  let m = Mat.(transpose top_mp *@ cov_p *@ top_mp) in
  Mat.save_txt ~out:(in_dir "analysis/pct_var_own") Mat.(pct_varp @= pct_varm);
  let x = Mat.(x_mov *@ top_mp *@ transpose top_mp) in
  let y = Mat.(x_prep *@ top_mv *@ transpose top_mv) in
  let z = Mat.(transpose top_mv *@ cov_p *@ top_mv) in
  let zz = Mat.(transpose top_mp *@ cov_m *@ top_mp) in
  let _ = Mat.save_txt ~out:(in_dir "analysis/top_mp") top_mp in
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
      ~out:(in_dir "analysis/pct_var")
      Mat.(
        (cumsum ~axis:1 (Mat.sqr tt) /$ sum' sp)
        @= (cumsum ~axis:1 (Mat.sqr ss) /$ sum' sm))
  , top_mp
  , top_mv )


let occupancy =
  let load i =
    Mat.get_slice
      [ [ i * duration; ((i + 1) * duration) - 1 ]; [] ]
      (snd preprocessed_data)
    (* let m = Mat.load_txt (in_dir Printf.(sprintf "rates_%i_%i" i t_prep)) in
    let ma = Mat.get_fancy [ R [ 0; duration - 1 ]; R [ 0; -1 ] ] m in
    ma *)
  in
  let projections x =
    let rp = reconstructed (proj x wp) wp in
    let rm = reconstructed (proj x wm) wm in
    Mat.(rp + rm), rp, rm, x
  in
  let rps =
    Arr.concatenate
      ~axis:2
      (Array.mapi (fun i _ ->
           let _, rp, _, _ = projections (load i) in
           let m = Arr.reshape rp [| (Arr.shape rp).(0); (Arr.shape rp).(1); 1 |] in
           m) reaches)
  and rms =
    Arr.concatenate
      ~axis:2
      (Array.mapi (fun i _ ->
           let _, _, rm, _ = projections (load i) in
           let m = Arr.reshape rm [| (Arr.shape rm).(0); (Arr.shape rm).(1); 1 |] in
           m) reaches)
  in
  let vprep, vmov =
    let x = Arr.var ~axis:2 rps |> Arr.sum ~axis:1 |> Arr.to_array in
    ( Mat.of_array x (-1) 1
    , let x = Arr.var ~axis:2 rms |> Arr.sum ~axis:1 |> Arr.to_array in
      Mat.of_array x (-1) 1 )
  in
  ( Mat.save_txt ~out:(in_dir "analysis/vprep") Mat.(vprep /$ max' vprep)
  , Mat.save_txt ~out:(in_dir "analysis/vmov") Mat.(vmov /$ max' vmov) )
