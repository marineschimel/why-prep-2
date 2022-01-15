open Owl
module AD = Algodiff.D
open Lib

let _ = Printexc.record_backtrace true
let t_prep = Cmdargs.(get_int "-prep" |> force ~usage:"-prep")
let t_prep_2 = Cmdargs.(get_int "-prep_2" |> default (t_prep + 900))
let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")
let double_dir = Cmdargs.(get_string "-dd" |> force ~usage:"-d [dir to save in]")
let compound_dir = Cmdargs.(get_string "-dc")

let compound =
  match compound_dir with
  | None -> false
  | Some _ -> true


let in_dir s = Printf.sprintf "%s/%s" dir s
let in_double_dir s = Printf.sprintf "%s/%s" double_dir s
let n_dim = 8
let n_reaches = 7
let n_var = 8
let dt = 2E-3
let reach_0 = Mat.load_txt (in_dir Printf.(sprintf "rates_%i_%i" 0 t_prep))

let double_reach_0 =
  Mat.load_txt (in_double_dir Printf.(sprintf "rates_%i_%i_%i" 0 2 t_prep))


let size_prep = 200
let size_prep_2 = 100
let size_mov = 200
let duration = Mat.row_num reach_0
let double_duration = Mat.row_num double_reach_0
let n_prep = float_of_int t_prep /. 1000. /. dt |> int_of_float
let n_prep_2 = float_of_int t_prep_2 /. 1000. /. dt |> int_of_float

let softmax x =
  let m = Mat.(max ~axis:0 x - min ~axis:0 x +$ (0.1 *. max' x)) in
  Mat.(x / m)


let ls_double_reaches =
  [| 0, 1
   ; 0, 2
   ; 0, 3
   ; 0, 4
   ; 0, 5
   ; 0, 6
   ; 1, 2
   ; 1, 5
   ; 2, 3
   ; 2, 1
   ; 2, 0
   ; 2, 4
   ; 2, 5
   ; 3, 1
   ; 3, 2
   ; 3, 3
   ; 3, 4
   ; 3, 5
   ; 3, 6
   ; 4, 2
   ; 4, 3
   ; 4, 5
   ; 4, 6
   ; 5, 0
   ; 5, 1
   ; 5, 2
   ; 5, 3
   ; 5, 4
   ; 5, 6
   ; 6, 0
   ; 6, 1
   ; 6, 2
   ; 6, 3
   ; 6, 5
  |]


let preprocessed_data =
  let load i =
    let m = Mat.load_txt (in_dir Printf.(sprintf "rates_%i_%i" i t_prep)) in
    let ma = Mat.get_fancy [ R [ 0; duration - 1 ]; R [ 0; -1 ] ] m in
    Arr.reshape ma [| duration; -1; 1 |]
  in
  let data = Array.init n_reaches (fun i -> load i) in
  let concat_data = Arr.concatenate ~axis:0 data in
  let norm_factor =
    Mat.(max ~axis:0 concat_data - min ~axis:0 concat_data +$ (0.2 *. max' concat_data))
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
  norm_factor, Mat.of_arrays (Arr.to_arrays dat)


let double_data =
  let load i j =
    let m = Mat.load_txt (in_double_dir Printf.(sprintf "rates_%i_%i_%i" i j t_prep)) in
    Arr.reshape m [| Mat.row_num m; -1; 1 |]
  in
  let data = Array.map (fun (i, j) -> load i j) ls_double_reaches in
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
  Mat.of_arrays (Arr.to_arrays dat)


let x_prep_from_single =
  let m = Arr.reshape (snd preprocessed_data) [| n_reaches; duration; -1 |] in
  let f i =
    Arr.reshape
      (Arr.get_slice [ [ i ]; [ n_prep - size_prep; n_prep - 1 ]; [] ] m)
      [| size_prep; -1 |]
  in
  Array.init n_reaches f |> Arr.concatenate ~axis:0


let x_prep_from_double =
  let m =
    Arr.reshape double_data [| Array.length ls_double_reaches; double_duration; -1 |]
  in
  let f i =
    let p1 =
      Arr.reshape
        (Arr.get_slice [ [ i ]; [ n_prep - size_prep; n_prep - 1 ]; [] ] m)
        [| size_prep; -1 |]
    in
    let p2 =
      Arr.reshape
        (Arr.get_slice [ [ i ]; [ n_prep_2 - size_prep_2; n_prep_2 - 1 ]; [] ] m)
        [| size_prep_2; -1 |]
    in
    Arr.concatenate ~axis:0 [| p1; p2 |]
  in
  Array.init n_reaches f |> Arr.concatenate ~axis:0


let x_mov_from_single =
  let m = Arr.reshape (snd preprocessed_data) [| n_reaches; duration; -1 |] in
  let f i =
    Arr.reshape
      (Arr.get_slice [ [ i ]; [ n_prep + 25; n_prep + 25 + size_mov - 1 ]; [] ] m)
      [| size_prep; -1 |]
  in
  Array.init n_reaches f |> Arr.concatenate ~axis:0


let x_mov_from_double =
  let m =
    Arr.reshape double_data [| Array.length ls_double_reaches; double_duration; -1 |]
  in
  let f i =
    let p1 =
      Arr.reshape
        (Arr.get_slice [ [ i ]; [ n_prep + 25; n_prep + 25 + size_mov - 1 ]; [] ] m)
        [| size_mov; -1 |]
    in
    let p2 =
      Arr.reshape
        (Arr.get_slice [ [ i ]; [ n_prep_2 + 25; n_prep_2 + 25 + size_mov - 1 ]; [] ] m)
        [| size_mov; -1 |]
    in
    Arr.concatenate ~axis:0 [| p1; p2 |]
  in
  Array.init (Array.length ls_double_reaches) f |> Arr.concatenate ~axis:0


let double_data =
  if compound
  then (
    let load i j =
      let m =
        Mat.load_txt
          Printf.(sprintf "%s/rates_%i_%i_%i" (Option.get compound_dir) i j t_prep)
      in
      Arr.reshape m [| Mat.row_num m; -1; 1 |]
    in
    let data = Array.map (fun (i, j) -> load i j) ls_double_reaches in
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
    Mat.of_arrays (Arr.to_arrays dat))
  else double_data


let n = Mat.col_num x_mov_from_single
let x_prep = Arr.(x_prep_from_single @= x_prep_from_double)
let x_mov = Arr.(x_mov_from_single @= x_mov_from_double)

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


let saving_dir = if compound then Option.get compound_dir else double_dir

(*this saves : var of prep in move subspace, var move in prep subspace, var prep in prep subspace, var move in move subspace
  for various dimensions*)
let _ =
  Array.init 9 (fun i -> capt_var i)
  |> Mat.of_arrays
  |> Mat.save_txt ~out:(Printf.sprintf "%s/analysis/capt_var" saving_dir)


(*below, code to check how orthogonal the activity at the end of prep is to the rest of the movement*)

(*Elsayed et al analysis : are the subspaces orthogonal?*)

let corr_p =
  let r = Mat.(x_prep_from_single - mean ~axis:0 x_prep_from_single) in
  let rr = Mat.(r / l2norm ~axis:0 r) in
  Mat.(transpose rr *@ rr)


let corr_m =
  let r = Mat.(x_mov_from_single - mean ~axis:0 x_mov_from_single) in
  let rr = Mat.(r / l2norm ~axis:0 r) in
  Mat.(transpose rr *@ rr)


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


let proj x modes =
  (*x is T*N and modes is N*n_dim*)
  let proj = Mat.(x *@ modes) in
  proj


(*proj is T*n_dim*)
let reconstructed proj modes = Mat.(proj *@ transpose modes)
let in_saving_dir = Printf.sprintf "%s/%s" saving_dir

(*T*N*)
let proj2d =
  let m = Arr.reshape (snd preprocessed_data) [| n_reaches; duration; -1 |] in
  fun i ->
    let x =
      let a = Arr.reshape (Arr.get_slice [ [ i ]; []; [] ] m) [| duration; -1 |] in
      a
    in
    let _ = Mat.save_txt ~out:(in_dir (Printf.sprintf "analysis/check_%i" i)) x in
    let proj_mov = proj x wm in
    let proj_prep = proj x wp in
    Mat.save_txt ~out:(in_saving_dir (Printf.sprintf "analysis/proj_prep_%i" i)) proj_prep;
    Mat.save_txt ~out:(in_saving_dir (Printf.sprintf "analysis/proj_mov_%i" i)) proj_mov


let _ = Array.init n_reaches (fun i -> proj2d i)

let occupancy_double =
  let double_duration =
    let reach_0 = Mat.load_txt (in_saving_dir Printf.(sprintf "rates_%i_%i" 0 t_prep)) in
    Mat.row_num reach_0
  in
  let load k =
    Mat.get_slice
      [ [ k * double_duration; ((k + 1) * double_duration) - 1 ]; [] ]
      double_data
    (* let m = Mat.load_txt Printf.(sprintf "%s/rates_%i_%i_%i" double_dir i j t_prep) in
    m *)
  in
  let projections x =
    let rp = Mat.(x *@ wp) in
    let rm = Mat.(x *@ wm) in
    Mat.(rp + rm), rp, rm, x
  in
  let rps =
    Arr.concatenate
      ~axis:2
      (Array.mapi
         (fun k _ ->
           let _, rp, _, _ = projections (load k) in
           let m = Arr.reshape rp [| (Arr.shape rp).(0); (Arr.shape rp).(1); 1 |] in
           m)
         ls_double_reaches)
  and rms =
    Arr.concatenate
      ~axis:2
      (Array.mapi
         (fun k _ ->
           let _, _, rm, _ = projections (load k) in
           let m = Arr.reshape rm [| (Arr.shape rm).(0); (Arr.shape rm).(1); 1 |] in
           m)
         ls_double_reaches)
  in
  let vprep, vmov =
    let x = Arr.var ~axis:2 rps |> Arr.sum ~axis:1 |> Arr.to_array in
    ( Mat.of_array x (-1) 1
    , let x = Arr.var ~axis:2 rms |> Arr.sum ~axis:1 |> Arr.to_array in
      Mat.of_array x (-1) 1 )
  in
  ( Mat.save_txt ~out:(in_double_dir "analysis/vprep") Mat.(vprep /$ max' vprep)
  , Mat.save_txt ~out:(in_double_dir "analysis/vmov") Mat.(vmov /$ max' vmov) )
