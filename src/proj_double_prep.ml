open Owl
module AD = Algodiff.D
open Lib

let _ = Printexc.record_backtrace true
let t_prep = Cmdargs.(get_int "-prep" |> force ~usage:"-prep")
let t_prep_2 = Cmdargs.(get_int "-prep_2" |> default (t_prep + 800))
let n_prep_2 = t_prep_2 / 2
let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")
let dir = "/home/mmcs3/rds/rds-t2-cs156-T7o4pEA8QoU/mmcs3/random_monkeys_random_lambda_2E-6/ramping/seed_5_200_2"
let in_dir seed s = Printf.sprintf "%s/%s" dir s
  (*Printf.sprintf "%s/seed_%i/%s" dir seed s*)
let in_double_dir _seed s = Printf.sprintf "/home/mmcs3/rds/rds-t2-cs156-T7o4pEA8QoU/mmcs3/random_monkeys_random_lambda_2E-6/ramping/seed_5_200_2/double_ramping_1002002_pause_0.05_0._0.000002_300400/%s" s
let superscript = Cmdargs.(get_string "-superscript" |> default "")
let proj_dir = Printf.sprintf "single_projs_nopause%s" superscript
let compound = Cmdargs.check "-compound"
let in_compound_dir _seed s = Printf.sprintf "/home/mmcs3/rds/rds-t2-cs133-hh9aMiOkJqI/mmcs3/random_monkeys_lambda_1E-6/double_ramping_short/%s" s
  (* Printf.sprintf "%s/300_350%s" dir superscript *)
let ws_dir = "analysis"
(* let ws_dir = "double_projs_long" *)
let saving_dir = if compound then in_compound_dir else in_double_dir
let n_dim = 8
let n_var = 8
let dt = 2E-3
let n_neurons = 200

let seed = Cmdargs.(get_int "-seed" |> force ~usage:"-seed")

let reach_0 = Mat.load_txt (in_dir 1 Printf.(sprintf "rates_%i_%i" 1 t_prep))

let double_reach_0 =
  Mat.load_txt (in_double_dir seed Printf.(sprintf "rates_%i_%i_%i" 0 1 t_prep))

(* |> Mat.get_slice [ [ 0; size_double ] ] *)

let size_prep = 200
let size_prep_2 = 100
let size_mov = 200
let duration = Mat.row_num reach_0
let double_duration = Mat.row_num double_reach_0
let size_double = Int.(double_duration - 1)
let n_prep = float_of_int t_prep /. 1000. /. dt |> int_of_float
let n_prep_2 = float_of_int t_prep_2 /. 1000. /. dt |> int_of_float

let softmax x =
  let m = Mat.(max ~axis:0 x - min ~axis:0 x +$ (0.1 *. max' x)) in
  Mat.(x / m)


(* let ls_double_reaches = [| 0, 0; 0, 1; 0, 2; 0, 3; 0, 4; 0, 5; 0, 6 |] *)
let ls_double_reaches =
  [| 0, 1
   ; 0, 2
   ; 0, 4
   ; 0, 5
   ; 0, 6
   ; 1, 0
   ; 1, 2
   ; 1, 3
   ; 1, 4
   ; 1, 5
   ; 1, 6
   ; 2, 1
   ; 2, 3
   ; 2, 4
   ; 2, 5
   ; 2, 6
   ; 2, 0
   ; 3, 0
   ; 3, 1
   ; 3, 2
   ; 3, 4
   ; 3, 5
   ; 3, 6
   ; 4, 2
   ; 4, 1
   ; 4, 3
   ; 4, 5
   ; 4, 6
   ; 5, 1
   ; 5, 0
   ; 5, 2
   ; 5, 4
   ; 5, 3
   ; 5, 6
   ; 6, 0
   ; 6, 1
   ; 6, 2
   ; 6, 3
   ; 6, 4
   ; 6, 5 (* ; 6, 4 *)
  |]

let single_reaches = [|1;2;3;4;6|]

let n_reaches = (Array.length single_reaches)

let n_double_reaches = Array.length ls_double_reaches

let preprocessed_data seed =
  let load i =
    let m = Mat.load_txt (in_dir seed Printf.(sprintf "rates_%i_%i" i t_prep)) in
    let ma = Mat.get_fancy [ R [ 0; duration - 1 ]; R [ 0; -1 ] ] m in
    Arr.reshape ma [| duration; -1; 1 |]
  in
  let data = Array.map (fun i -> load i) single_reaches in
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


let get_double_data seed =
  let load i j =
    let m =
      Mat.load_txt (in_double_dir seed Printf.(sprintf "rates_%i_%i_%i" i j t_prep))
      |> fun z -> Mat.get_slice [ [ 0; size_double ] ] z
    in
    Arr.reshape m [| Mat.row_num m; -1; 1 |]
  in
  let data =
    let open Base in
    Array.fold ls_double_reaches ~init:[] ~f:(fun acc x ->
        let i, j = x in
        try
          let m = load i j in
          let new_acc = m :: acc in
          new_acc
        with
        | _ -> acc)
    |> Array.of_list
  in
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


let x_prep_from_single seed =
  let m = Arr.reshape (snd (preprocessed_data seed)) [| n_reaches; duration; -1 |] in
  let f i =
    Arr.reshape
      (Arr.get_slice [ [ i ]; [ n_prep - size_prep; n_prep - 1 ]; [] ] m)
      [| size_prep; -1 |]
  in
  Array.init n_reaches f |> Arr.concatenate ~axis:0


let x_mov_from_single seed =
  let m = Arr.reshape (snd (preprocessed_data seed)) [| n_reaches; duration; -1 |] in
  let f i =
    Arr.reshape
      (Arr.get_slice [ [ i ]; [ n_prep + 25; n_prep + size_mov - 1 + 25 ]; [] ] m)
      [| size_mov; -1 |]
  in
  Array.init n_reaches f |> Arr.concatenate ~axis:0


let x_prep_from_double seed =
  let dat = get_double_data seed in
  let m = Arr.reshape dat [| -1; double_duration;n_neurons |] in
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
  Array.init (Arr.shape m).(0) f |> Arr.concatenate ~axis:0


let x_mov_from_double seed =
  let m = Arr.reshape (get_double_data seed) [| -1; double_duration; n_neurons |] in
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
  Array.init (Arr.shape m).(0) f |> Arr.concatenate ~axis:0


let data_to_proj seed =
  let load i j =
    let m =
      Mat.load_txt (saving_dir seed Printf.(sprintf "rates_%i_%i_%i" i j t_prep))
      |> Mat.get_slice [ [ 0; -1 ] ]
    in
    Arr.reshape m [| Mat.row_num m; -1 |]
  in
  let data =
    let open Base in
    Array.fold ls_double_reaches ~init:[] ~f:(fun acc x ->
        let i, j = x in
        try
          let m = load i j in
          let new_acc = m :: acc in
          new_acc
        with
        | _ -> acc)
    |> Array.of_list
  in
  let norm_data =  let concat_data = Arr.concatenate ~axis:0 data in  
  let norm_factor =
    Mat.(max ~axis:0 concat_data - min ~axis:0 concat_data +$ (0.1 *. max' concat_data)) 
    |> fun z -> Arr.reshape z [| 1; -1 |]
  in 
  let _ = Stdio.printf "shape : %i %i %i %i %!" (Arr.shape norm_factor).(0)  (Arr.shape norm_factor).(1) (Mat.row_num data.(0)) (Mat.col_num data.(0))
in 
Array.map (fun x -> Arr.(x / norm_factor) |> fun z ->  Arr.reshape z [| (Arr.shape z).(0); -1 ; 1 |]) data in 
  let mean_data = Arr.concatenate ~axis:2 norm_data |> Arr.mean ~axis:2 in
  let dat =
    let a =
      Arr.concatenate ~axis:0 (Array.map (fun x -> Arr.(x - mean_data)) norm_data)
    in
    Arr.reshape a [| (Arr.shape a).(0); (Arr.shape a).(1) |]
  in
  Mat.of_arrays (Arr.to_arrays dat)


let n = Mat.col_num (x_mov_from_single 1)
let x_prep seed = Arr.(x_prep_from_single seed @= x_prep_from_double seed)
let x_mov seed = Arr.(x_mov_from_single seed @= x_mov_from_double seed)

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
      pct_change < 1E-3
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


let proj x modes =
  (*x is T*N and modes is N*n_dim*)
  let proj = Mat.(x *@ modes) in
  proj


let reconstructed proj modes = Mat.(proj *@ transpose modes)

let preprocessed_data seed =
  let load i =
    let m = Mat.load_txt (in_dir seed Printf.(sprintf "rates_%i_%i" i t_prep)) in
    let ma = Mat.get_fancy [ R [ 0; duration - 1 ]; R [ 0; -1 ] ] m in
    Arr.reshape ma [| duration; -1; 1 |]
  in
  let data = Array.map (fun i -> load i) single_reaches in
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



let wp, wm =
  try
    ( Mat.load_txt (in_dir seed (Printf.sprintf "%s/wp" ws_dir))
    , Mat.load_txt (in_dir seed (Printf.sprintf "%s/wm" ws_dir)) )
  with
  | _ ->
    let wp, wm = ws seed in
    Mat.save_txt ~out:(in_dir seed (Printf.sprintf "%s/wp" ws_dir)) wp;
    Mat.save_txt ~out:(in_dir seed (Printf.sprintf "%s/wm" ws_dir)) wm;
    wp, wm


let x_sub_prep_from_double =
  let open Base in
  let m = Arr.reshape (get_double_data seed) [| -1; double_duration; n_neurons |] in
  Array.init 6 ~f:(fun reach ->
      let ls =
        Array.foldi ls_double_reaches ~init:[] ~f:(fun j accu (n, i) ->
            try
              let x =
                Arr.reshape
                  (Arr.get_slice [ [ j ]; [ n_prep - size_prep; n_prep - 1 ]; [] ] m)
                  [| size_prep; -1 |]
                (* in
              if n = reach && (Int.(n + i = 6) || (n = 3 && i = 2) || (n = 5 && i = 2)) *)
                (* then x :: accu *)
              in
              x :: accu
            with
            | _ -> accu)
        (* else accu) *)
      in
      Mat.concatenate ~axis:0 (Array.of_list ls))
  |> Mat.concatenate ~axis:0


(*size n_reaches * N*)

let occupancy seed =
  let double_reach_0 =
    Mat.load_txt (saving_dir seed Printf.(sprintf "rates_%i_%i_%i" 0 1 t_prep))
    |> Mat.get_slice [ [ 0; -1 ] ]
  in
  let dat = data_to_proj seed in
  let n = (Arr.shape dat).(0) in
  let double_duration = Mat.row_num double_reach_0 in
  let load i =
    Mat.get_slice [ [ i * double_duration; ((i + 1) * double_duration) - 1 ]; [] ] dat
  in
  let projections x =
    let rp = reconstructed (proj x wp) wp in
    rp
  in
  let rps =
    Arr.concatenate
      ~axis:2
      (Array.init (n / double_duration) (fun i ->
           let rp = projections (load i) in
           let m = Arr.reshape rp [| (Arr.shape rp).(0); (Arr.shape rp).(1); 1 |] in
           m))
  in
  let vprep =
    let x = Arr.var ~axis:2 rps |> Arr.sum ~axis:1 |> Arr.to_array in
    Mat.of_array x (-1) 1
  in
  let projections x =
    let rm = reconstructed (proj x wm) wm in
    rm
  in
  let rms =
    Arr.concatenate
      ~axis:2
      (Array.init (n / double_duration) (fun i ->
           let rp = projections (load i) in
           let m = Arr.reshape rp [| (Arr.shape rp).(0); (Arr.shape rp).(1); 1 |] in
           m))
  in
  let vmov =
    let x = Arr.var ~axis:2 rms |> Arr.sum ~axis:1 |> Arr.to_array in
    Mat.of_array x (-1) 1
  in
  Mat.(vmov /$ max' vmov), Mat.(vprep /$ max' vprep)


let _ =
  let vmov, vprep = occupancy seed in
  Mat.save_txt ~out:(in_dir seed (Printf.sprintf "%s/vprep" proj_dir)) vprep;
  Mat.save_txt ~out:(in_dir seed (Printf.sprintf "%s/vmov" proj_dir)) vmov


let proj2d =
  let double_reach_0 =
    Mat.load_txt (saving_dir seed Printf.(sprintf "rates_%i_%i_%i" 0 1 t_prep))
    |> Mat.get_slice [ [ 0; -1 ] ]
  in
  let double_duration = Mat.row_num double_reach_0 in
  let m = Arr.reshape (data_to_proj seed) [| -1; double_duration; n_neurons |] in
  fun i ->
    let x =
      let a = Arr.reshape (Arr.get_slice [ [ i ]; []; [] ] m) [| double_duration; -1 |] in
      a
    in
    let proj_prep = proj x wp in
    let proj_mov = proj x wm in
    Mat.save_txt
      ~out:(Printf.sprintf "%s/%s/proj_prep_%i" dir proj_dir i)
      proj_prep;
    Mat.save_txt
      ~out:(Printf.sprintf "%s/%s/proj_mov_%i" dir proj_dir i)
      proj_mov


let _ =
  Array.init 6 (fun i ->
      proj2d i;
      i)


let first_reaches = Mat.(x_prep_from_single 1 @= x_sub_prep_from_double)

let b =
  let proj_wp = Mat.(first_reaches *@ wp) in
  let _, _, v = Linalg.D.svd proj_wp in
  v


let wp = Mat.(wp *@ b)

let single_projs reach =
  let m = Arr.reshape (preprocessed_data seed) [| n_reaches; duration; -1 |] in
  let x = Arr.reshape (Arr.get_slice [ [ reach ]; [ n_prep ]; [] ] m) [| 1; -1 |] in
  proj x wp


let double_projs reach =
  let open Base in
  let m = Arr.reshape (get_double_data seed) [| -1; double_duration; n_neurons |] in
  let ls =
    Array.foldi ls_double_reaches ~init:[] ~f:(fun j accu (n, _) ->
        try
          let x = Arr.reshape (Arr.get_slice [ [ j ]; [ n_prep ]; [] ] m) [| 1; -1 |] in
          if n = reach then proj x wp :: accu else accu
        with
        | _ -> accu)
  in
  Mat.concatenate ~axis:0 (Array.of_list ls)


(*rotate W : project the prep of single reaches and 7 double reaches onto it, then apply PCA to the N_d x (14xT) matrix and we get a matrix B by which we rotate Wp *)
(* let _ =
  Array.iter
    (fun i ->
      let m_doub = double_projs i in
      let m_sing = single_projs i in
      let conc = Mat.(transpose m_sing @|| transpose m_doub) in
      Mat.save_txt ~out:(in_dir seed (Printf.sprintf "%s/single_%i" proj_dir i)) m_sing;
      Mat.save_txt ~out:(in_dir seed (Printf.sprintf "%s/doub_%i" proj_dir i)) m_doub;
      Mat.save_txt ~out:(in_dir seed (Printf.sprintf "%s/conc_%i" proj_dir i)) conc)
    (Array.init (n_reaches) (fun i -> i))  *)
