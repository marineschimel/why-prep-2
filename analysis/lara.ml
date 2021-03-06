open Owl
module AD = Algodiff.D
open Lib
open Defaults
module J = Jpca

let _ = Printexc.record_backtrace true

let t_prep = 600

let n_dim = 4

let n_reaches = 7

let n_var = 4

let in_dir s = Printf.sprintf "results/test_ort/random/%s" s

let size_prep = 400

let size_mov = 400

let duration = 1200

let n_shift = float_of_int t_prep /. 1000. /. sampling_dt |> int_of_float

let softmax x =
  let m = Mat.(max ~axis:0 x - min ~axis:1 x +$ (0.1 *. max' x)) in
  Mat.(x / m)

let preprocessed_data =
  let load i =
    let m =
      Mat.load_txt (in_dir Printf.(sprintf "reach_%i/traj_%i" i t_prep))
    in
    let ma = Mat.get_fancy [ R [ 0; duration ]; R [ 4; -1 ] ] m in
    Arr.reshape ma [| duration + 1; size_net; 1 |]
  in
  let data = Array.init n_reaches (fun i -> softmax (load (succ i))) in
  let mean_data = Arr.concatenate ~axis:2 data |> Arr.mean ~axis:2 in
  let dat =
    let a =
      Arr.concatenate ~axis:0 (Array.map (fun x -> Arr.(x - mean_data)) data)
    in
    Arr.reshape a [| (Arr.shape a).(0); (Arr.shape a).(1) |]
  in
  (mean_data, Mat.of_arrays (Arr.to_arrays dat))

(*
let x_prep =
  Mat.concatenate ~axis:0
    (Array.init n_reaches (fun i ->
         Mat.concatenate ~axis:0
           (Array.init (n_reaches - 1) (fun j ->
                softmax
                  (Mat.load_txt
                     (in_dir
                        Printf.(
                          sprintf "reach_%i%i/traj_%i" (succ i) (succ j) t_prep)))
                |> Mat.get_fancy [ R [ 0; t_prep - 150 ]; R [ 4; -1 ] ]))))

let x_mov =
  Mat.concatenate ~axis:0
    (Array.init n_reaches (fun i ->
         Mat.concatenate ~axis:0
           (Array.init (n_reaches - 1) (fun j ->
                softmax
                  (Mat.load_txt
                     (in_dir
                        Printf.(
                          sprintf "reach_%i%i/traj_%i" (succ i) (succ j) t_prep)))
                |> Mat.get_fancy
                     [ R [ t_prep - 50; t_prep + 250 ]; R [ 4; -1 ] ]))))
*)

let x_prep =
  let m =
    Arr.reshape (snd preprocessed_data) [| n_reaches; duration + 1; size_net |]
  in
  let f i =
    Arr.reshape
      (Arr.get_slice [ [ i ]; [ t_prep - size_prep; t_prep - 1 ]; [] ] m)
      [| size_prep; size_net |]
  in
  Array.init n_reaches f |> Arr.concatenate ~axis:0

let x_mov =
  let m =
    Arr.reshape (snd preprocessed_data) [| n_reaches; duration + 1; size_net |]
  in
  let f i =
    Arr.reshape
      (Arr.get_slice
         [ [ i ]; [ t_prep + 50; t_prep + size_mov - 1 + 50 ]; [] ]
         m)
      [| size_mov; size_net |]
  in
  Array.init n_reaches f |> Arr.concatenate ~axis:0

let n = Mat.col_num x_mov

let cov_p =
  let r = Mat.(x_prep - mean ~axis:0 x_prep) in
  Mat.(transpose r *@ r)

let cov_m =
  let r = Mat.(x_mov - mean ~axis:0 x_mov) in
  Mat.(transpose r *@ r)

let __cp, __cm = (AD.pack_arr cov_p, AD.pack_arr cov_m)

let capt_var dim =
  let cp, varp, _ = Linalg.D.svd cov_p in
  let cm, varm, _ = Linalg.D.svd cov_m in
  let cm = Mat.get_slice [ []; [ 0; dim - 1 ] ] cm in
  let cp = Mat.get_slice [ []; [ 0; dim - 1 ] ] cp in
  let _ = Printf.printf "%i %i %!" (Mat.row_num cp) (Mat.col_num cp) in
  let var_p_mov =
    Mat.(trace Mat.(transpose cm *@ cov_p *@ cm)) /. Mat.sum' varp
  in
  let var_m_prep =
    Mat.(trace Mat.(transpose cp *@ cov_m *@ cp)) /. Mat.sum' varm
  in
  let var_p_prep =
    Mat.(trace Mat.(transpose cp *@ cov_p *@ cp)) /. Mat.sum' varp
  in
  let var_m_mov =
    Mat.(trace Mat.(transpose cm *@ cov_m *@ cm)) /. Mat.sum' varm
  in
  [| var_p_mov; var_m_prep; var_p_prep; var_m_mov |]

(*this saves : var of prep in move subspace, var move in prep subspace, var prep in prep subspace, var move in move subspace
  for various dimensions*)
let _ =
  Array.init 9 (fun i -> capt_var (succ i))
  |> Mat.of_arrays
  |> Mat.save_txt ~out:(in_dir "analysis/capt_var")

(*below, code to check how orthogonal the activity at the end of prep is to the rest of the movement*)

(* let comp_ac_end =
  let x =
    Mat.load_txt (in_dir Printf.(sprintf "traj_%i" t_prep)) |> fun z ->
    Mat.get_slice [ []; [ 4; -1 ] ] z |> fun z -> Mat.(z / l2norm ~axis:1 z)
  in
  let ac_end_prep =
    Mat.get_slice [ [ t_prep; t_prep + 50 ]; [] ] x |> Mat.mean ~axis:0
  in
  let ac_move = Mat.get_slice [ [ 0; -1 ]; [] ] x in
  let comp = Mat.(ac_end_prep *@ transpose ac_move) |> Mat.transpose in
  Mat.save_txt ~out:"comp" comp *)

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
  ( Mat.save_txt ~out:(in_dir "analysis/cp") corr_p,
    Mat.save_txt ~out:(in_dir "analysis/cm") corr_m )

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
  (AD.Maths.(neg obj), new_wp, new_wm)

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
  let f _ prms =
    let f, _, _ = cost prms in
    f
  in
  let s0 = P.init ~f ~prms0 () in
  let sf = P.min ~stop s0 in
  let prms = P.prms sf in
  let _, wp, wm = cost prms in
  (prms, wp, wm)

let wp, wm =
  let _, new_wp, new_wm = learn in
  (new_wp |> AD.unpack_arr, new_wm |> AD.unpack_arr)

let proj x modes =
  (*x is T*N and modes is N*n_dim*)
  let proj = Mat.(x *@ modes) in
  proj

(*proj is T*n_dim*)
let reconstructed proj modes = Mat.(proj *@ transpose modes)

(*T*N*)
let proj2d =
  let m =
    Arr.reshape (snd preprocessed_data) [| n_reaches; duration + 1; size_net |]
  in
  fun i ->
    let x =
      let a =
        Arr.reshape
          (Arr.get_slice [ [ i - 1 ]; []; [] ] m)
          [| duration + 1; size_net |]
        |> Arr.to_arrays
      in
      Mat.of_arrays a
    in
    let proj_mov = proj x wm in
    let proj_prep = proj x wp in
    Mat.save_txt
      ~out:(in_dir (Printf.sprintf "analysis/proj_prep_%i" i))
      proj_prep;
    Mat.save_txt
      ~out:(in_dir (Printf.sprintf "analysis/proj_mov_%i" i))
      proj_mov

let modes_jpca x =
  let diff x =
    Mat.(get_slice [ [ 0; -2 ]; [] ] x - get_slice [ [ 1; -1 ]; [] ] x)
  in
  let dx = Array.map diff x in
  let wskew = J.get_mskew ~dx ~x in
  let modes = J.get_modes_jpca wskew in
  Mat.save_txt ~out:(in_dir "analysis/modes_jpca") modes

let captured_variance, top_mp, top_mv =
  let _mp, sp, _ = Linalg.D.svd (Mat.transpose x_prep) in
  let _mv, sm, _ = Linalg.D.svd (Mat.transpose x_mov) in
  let _ = Mat.save_txt ~out:"sm" Mat.(transpose (sqr sm)) in
  let sp, sm = (Mat.sqr sp, Mat.sqr sm) in
  let pct_varp, pct_varm =
    (Mat.(cumsum ~axis:1 sp /$ sum' sp), Mat.(cumsum ~axis:1 sm /$ sum' sm))
  in
  let top_mp, top_mv =
    ( Mat.get_slice [ []; [ 0; n_var - 1 ] ] wp,
      Mat.get_slice [ []; [ 0; n_var - 1 ] ] wm )
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
    Printf.printf "pct_var %f %f %f %!"
      (Mat.trace z /. Mat.sum' (Mat.get_slice [ []; [ 0; n_var - 1 ] ] sp))
      (Mat.trace zz /. Mat.sum' (Mat.get_slice [ []; [ 0; n_var - 1 ] ] sm))
      (Mat.trace m /. Mat.sum' (Mat.get_slice [ []; [ 0; n_var - 1 ] ] sp))
  in
  ( Mat.save_txt
      ~out:(in_dir "analysis/pct_var")
      Mat.(
        (cumsum ~axis:1 (Mat.sqr tt) /$ sum' sp)
        @= (cumsum ~axis:1 (Mat.sqr ss) /$ sum' sm)),
    top_mp,
    top_mv )

let occupancy =
  let load i =
    let m =
      Mat.load_txt (in_dir Printf.(sprintf "reach_%i/traj_%i" i t_prep))
    in
    let ma = Mat.get_fancy [ R [ 0; duration - 1 ]; R [ 4; -1 ] ] m in
    ma
  in
  let projections x =
    let rp = reconstructed (proj x wp) wp in
    let rm = reconstructed (proj x wm) wm in
    (Mat.(rp + rm), rp, rm, x)
  in
  let rps =
    Arr.concatenate ~axis:2
      (Array.init n_reaches (fun i ->
           let _, rp, _, _ = projections (load (succ i)) in
           let m =
             Arr.reshape rp [| (Arr.shape rp).(0); (Arr.shape rp).(1); 1 |]
           in
           m))
  and rms =
    Arr.concatenate ~axis:2
      (Array.init n_reaches (fun i ->
           let _, _, rm, _ = projections (load (succ i)) in
           let m =
             Arr.reshape rm [| (Arr.shape rm).(0); (Arr.shape rm).(1); 1 |]
           in
           m))
  in
  let vprep, vmov =
    let x = Arr.var ~axis:2 rps |> Arr.sum ~axis:1 |> Arr.to_array in
    ( Mat.of_array x (-1) 1,
      let x = Arr.var ~axis:2 rms |> Arr.sum ~axis:1 |> Arr.to_array in
      Mat.of_array x (-1) 1 )
  in
  ( Mat.save_txt ~out:(in_dir "analysis/vprep") Mat.(vprep /$ max' vprep),
    Mat.save_txt ~out:(in_dir "analysis/vmov") Mat.(vmov /$ max' vmov) )

let occupancy_double =
  let load i j =
    let m =
      Mat.load_txt
        Printf.(sprintf "results/double/soc/compound_%i%i/traj_600" i j)
    in
    let ma = Mat.get_fancy [ R [ 0; -1 ]; R [ 4; -1 ] ] m in
    ma
  in
  let projections x =
    let rp = Mat.(x *@ wp) in
    let rm = Mat.(x *@ wm) in
    (Mat.(rp + rm), rp, rm, x)
  in
  let rps =
    Arr.concatenate ~axis:2
      (Array.map
         (fun (i, j) ->
           let _, rp, _, _ = projections (load i j) in
           let m =
             Arr.reshape rp [| (Arr.shape rp).(0); (Arr.shape rp).(1); 1 |]
           in
           m)
         [|
           (0, 1);
           (0, 3);
           (0, 4);
           (0, 5);
           (0, 6);
           (1, 2);
           (1, 3);
           (1, 5);
           (2, 3);
           (3, 2);
           (3, 4);
           (4, 3);
           (5, 6);
           (6, 0);
           (6, 5);
         |])
  and rms =
    Arr.concatenate ~axis:2
      (Array.map
         (fun (i, j) ->
           let _, _, rm, _ = projections (load i j) in
           let m =
             Arr.reshape rm [| (Arr.shape rm).(0); (Arr.shape rm).(1); 1 |]
           in
           m)
         [|
           (0, 1);
           (0, 3);
           (0, 4);
           (0, 5);
           (0, 6);
           (1, 2);
           (1, 3);
           (1, 5);
           (2, 3);
           (3, 2);
           (3, 4);
           (4, 3);
           (5, 6);
           (6, 0);
           (6, 5);
         |])
  in
  let vprep, vmov =
    let x = Arr.var ~axis:2 rps |> Arr.sum ~axis:1 |> Arr.to_array in
    ( Mat.of_array x (-1) 1,
      let x = Arr.var ~axis:2 rms |> Arr.sum ~axis:1 |> Arr.to_array in
      Mat.of_array x (-1) 1 )
  in
  ( Mat.save_txt ~out:"results/double/analysis/vprep" Mat.(vprep /$ max' vprep),
    Mat.save_txt ~out:"results/double/analysis/vmov" Mat.(vmov /$ max' vmov) )

let _ =
  Array.init n_reaches (fun i ->
      proj2d (succ i);
      0.)

let proj_pc =
  let mp =
    let x, _, _ = Linalg.D.svd (Mat.transpose x_prep) in
    Mat.get_slice [ []; [ 0; 5 ] ] x
  in
  let mv =
    let x, _, _ = Linalg.D.svd (Mat.transpose x_mov) in
    Mat.get_slice [ []; [ 0; 5 ] ] x
  in
  let m =
    Arr.reshape (snd preprocessed_data) [| n_reaches; duration + 1; size_net |]
  in
  fun i ->
    let x =
      let a =
        Arr.reshape
          (Arr.get_slice [ [ i - 1 ]; []; [] ] m)
          [| duration + 1; size_net |]
        |> Arr.to_arrays
      in
      Mat.of_arrays a
    in
    let proj_mov = reconstructed (proj x mv) mv in
    let proj_prep = reconstructed (proj x mp) mp in
    let proj_mov_ort = reconstructed (proj x wm) wm in
    let proj_prep_ort = reconstructed (proj x wp) wp in
    Mat.save_txt
      ~out:(in_dir (Printf.sprintf "analysis/proj_pc_prep_%i" i))
      proj_prep;
    Mat.save_txt
      ~out:(in_dir (Printf.sprintf "analysis/proj_pc_all_%i" i))
      Mat.(proj_prep @|| proj_mov);
    Mat.save_txt
      ~out:(in_dir (Printf.sprintf "analysis/proj_pc_mov_%i" i))
      proj_mov;
    Mat.save_txt
      ~out:(in_dir (Printf.sprintf "analysis/proj_ort_prep_%i" i))
      proj_prep_ort;
    Mat.save_txt
      ~out:(in_dir (Printf.sprintf "analysis/proj_ort_mov_%i" i))
      proj_mov_ort

let _ =
  Array.init n_reaches (fun i ->
      proj_pc (succ i);
      0.)

let _ =
  Linalg.D.rank cov_m |> Printf.printf "analysis/rank cov_m %i %!";
  Linalg.D.rank cov_p |> Printf.printf "analysis/rank cov_p %i %!"

let x =
  Array.init n_reaches (fun i ->
      let y =
        Mat.load_txt
          (in_dir Printf.(sprintf "reach_%i/traj_%i" (succ i) t_prep))
        |> Mat.get_fancy [ R [ 0; -1 ]; R [ 4; -1 ] ]
      in
      proj y wp
      |> Mat.get_slice [ []; [ 0; 2 ] ]
      |> Mat.save_txt ~out:(in_dir (Printf.sprintf "analysis/proj_%i" i)))

let potency_in_nullspace x =
  let x = Mat.(transpose (x /$ (l2norm' x +. 1E-8))) in
  let w = Mat.load_txt (in_dir "w") in
  let c = Mat.load_txt (in_dir "c") in
  let a = Mat.((w - eye 200) /$ Defaults.tau) in
  let q = obs_gramian a c in
  let null_proj =
    Mat.((eye 200 - (transpose c *@ (inv (c *@ transpose c) *@ c))) *@ x)
  in
  Mat.(transpose null_proj *@ q *@ null_proj)

let test =
  Mat.map_rows potency_in_nullspace
    Mat.(
      get_slice [ []; [ 4; -1 ] ]
        (load_txt Printf.(sprintf "double/soc/compound_01/traj_600")))
  |> Mat.concatenate ~axis:0 |> Mat.save_txt ~out:"test"

(*
let ()  =   
let dat = (snd preprocessed_data) in let sub_x =  let x, _, _ = Linalg.D.svd (Mat.transpose dat) in
Mat.get_slice [ []; [ 0; 5 ] ] x in 
let m = Arr.reshape (dat) [| n_reaches; duration + 1; size_net |] in
let x i = 
  let a =
    Arr.reshape (Arr.get_slice [ [ i - 1 ]; []; [] ] m) [| duration + 1; size_net |]
    |> Arr.to_arrays
  in
  Mat.of_arrays a
in 
let red_x i =  proj
(x i) sub_x in let _ = Mat.save_txt ~out:"red_x" (red_x 0) in 
  modes_jpca
    (Array.init n_reaches (fun i -> 
        red_x i
           ))

let jpca_proj x =

  let modes = Mat.load_txt (in_dir "analysis/modes_jpca") in
  let _ = Printf.printf "%i %i %i %i %!" (Mat.row_num modes) (Mat.col_num modes) (Mat.row_num x) (Mat.col_num x) in 
  Mat.((x)*@modes)

let _ =let dat = (snd preprocessed_data) in let x, _, _ = Linalg.D.svd (Mat.transpose dat) in
let m =  Mat.get_slice [ []; [ 0; 5 ] ] x in 
let _ = Printf.printf "%i %i %i %i %!" (Mat.row_num m) (Mat.col_num m) (Mat.row_num dat) (Mat.col_num dat) in 
let m = Arr.reshape Mat.(dat*@m) [| n_reaches; duration+1; 6 |]  in 
  Array.init n_reaches (fun i -> let x = Arr.get_slice [[i];[];[]] m in let x = Arr.reshape x [|duration+1;6|] in 
      Mat.save_txt
        ~out:(in_dir (Printf.sprintf "proj_jpca_%i" (succ i)))
        (jpca_proj (x));
      0.)  *)
