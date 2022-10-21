open Owl
module AD = Algodiff.D
open Lib
open Base

let _ = Printexc.record_backtrace true
let t_prep = Cmdargs.(get_int "-prep" |> force ~usage:"-prep")
let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")
let in_dir s = Printf.sprintf "%s/%s" dir s
let n_dim = 6
let reaches = [| 0; 1; 3; 4; 5; 6; 7 |]
let n_reaches = Array.length reaches
let n_var = n_dim
let dt = 2E-3
let reach_0 = Mat.load_txt (in_dir Printf.(sprintf "rates_%i_%i" 1 t_prep))
let size_prep = 150
let size_mov = 150
let duration = Mat.row_num reach_0
let n_prep = float_of_int t_prep /. 1000. /. dt |> int_of_float

let softmax x =
  let m = Mat.(max ~axis:0 x - min ~axis:0 x +$ (0.1 *. max' x)) in
  Mat.(x / m)

let array_rates =
  C.broadcast' (fun () ->
      Array.init 8 (fun i ->
          let m = Mat.load_txt (in_dir Printf.(sprintf "rates_%i_%i" i t_prep)) in
          let ma = Mat.get_fancy [ R [ 0; duration - 1 ]; R [ 0; -1 ] ] m in
          Arr.reshape ma [| duration; -1; 1 |]))

module Prms = struct
  type 'a t = { w : 'a } [@@deriving prms]
end

open Prms
module P = Owl_opt_lbfgs.Make (Prms)

let occupancies =
  Array.foldi
    (Array.init 100 ~f:(fun _ -> 1))
    ~init:[]
    ~f:(fun i accu _ ->
      if Int.(i % C.n_nodes = C.rank)
      then (
        let preprocessed_data =
          let load i = array_rates.(i) in
          let data = Array.map ~f:(fun i -> load i) reaches in
          let concat_data = Arr.concatenate ~axis:0 data in
          let norm_factor =
            Mat.(
              max ~axis:0 concat_data
              - min ~axis:0 concat_data
              +$ (0.1 *. max' concat_data))
            |> fun z -> Arr.reshape z [| 1; -1; 1 |]
          in
          let slice_neurons = List.init 200 (fun _ -> Random.int 200) in
          let norm_data = Array.map ~f:(fun x -> Arr.(x / norm_factor)) data in
          let norm_data =
            Array.map
              ~f:(fun x -> Arr.get_fancy [ R [ 0; -1 ]; L slice_neurons; R [ 0; -1 ] ] x)
              norm_data
          in
          let mean_data = Arr.concatenate ~axis:2 norm_data |> Arr.mean ~axis:2 in
          let dat =
            let a =
              Arr.concatenate
                ~axis:0
                (Array.map ~f:(fun x -> Arr.(x - mean_data)) norm_data)
            in
            Arr.reshape a [| (Arr.shape a).(0); (Arr.shape a).(1) |]
          in
          mean_data, Mat.of_arrays (Arr.to_arrays dat)
        in
        let x_prep =
          let m = Arr.reshape (snd preprocessed_data) [| n_reaches; duration; -1 |] in
          let f i _ =
            Arr.reshape
              (Arr.get_slice [ [ i ]; [ n_prep - size_prep; n_prep - 1 ]; [] ] m)
              [| size_prep; -1 |]
          in
          Array.mapi ~f reaches |> Arr.concatenate ~axis:0
        in
        let x_mov =
          let m = Arr.reshape (snd preprocessed_data) [| n_reaches; duration; -1 |] in
          let f i _ =
            Arr.reshape
              (Arr.get_slice [ [ i ]; [ n_prep + 25; n_prep + size_mov - 1 + 25 ]; [] ] m)
              [| size_mov; -1 |]
          in
          Array.mapi ~f reaches |> Arr.concatenate ~axis:0
        in
        let n = Mat.col_num x_mov in
        let cov_p =
          let r = Mat.(x_prep - mean ~axis:0 x_prep) in
          Mat.(transpose r *@ r)
        in
        let cov_m =
          let r = Mat.(x_mov - mean ~axis:0 x_mov) in
          Mat.(transpose r *@ r)
        in
        let __cp, __cm = AD.pack_arr cov_p, AD.pack_arr cov_m in
        let capt_var dim =
          let cp, varp, _ = Linalg.D.svd cov_p in
          let cm, varm, _ = Linalg.D.svd cov_m in
          let cm = Mat.get_slice [ []; [ 0; dim ] ] cm in
          let cp = Mat.get_slice [ []; [ 0; dim ] ] cp in
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
        in
        (*this saves : var of prep in move subspace, var move in prep subspace, var prep in prep subspace, var move in move subspace
  for various dimensions*)
        let _ =
          Array.init n_dim ~f:(fun i -> capt_var i)
          |> Mat.of_arrays
          |> Mat.save_txt ~out:(in_dir (Printf.sprintf "analysis_%i/capt_var" n_dim))
        in
        let corr_p =
          let r = Mat.(x_prep - mean ~axis:0 x_prep) in
          let rr = Mat.(r / l2norm ~axis:0 r) in
          Mat.(transpose rr *@ rr)
        in
        let corr_m =
          let r = Mat.(x_mov - mean ~axis:0 x_mov) in
          let rr = Mat.(r / l2norm ~axis:0 r) in
          Mat.(transpose rr *@ rr)
        in
        let proj x modes =
          (*x is T*N and modes is N*n_dim*)
          let proj = Mat.(x *@ modes) in
          proj
        in
        (*proj is T*n_dim*)
        let reconstructed proj modes = Mat.(proj *@ transpose modes) in
        let cp, varp, _ = Linalg.D.svd cov_p in
        let cm, varm, _ = Linalg.D.svd cov_m in
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
        in
        let learn =
          let prms0 = { w = AD.Mat.gaussian n (2 * n_dim) } in
          let stop =
            let fv_prev = ref 1E9 in
            fun s ->
              let k = P.iter s in
              let fv = P.(fv s) in
              let pct_change = abs_float ((fv -. !fv_prev) /. !fv_prev) in
              fv_prev := fv;
              (* Stdio.printf " iter %i | fv %f | pct change %f %!" k fv pct_change; *)
              Float.(pct_change < 1E-4)
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
        in
        let wp, wm =
          let _, new_wp, new_wm = learn in
          new_wp |> AD.unpack_arr, new_wm |> AD.unpack_arr
        in
        let load i =
          Mat.get_slice
            [ [ i * duration; ((i + 1) * duration) - 1 ]; [] ]
            (snd preprocessed_data)
        in
        let projections x =
          let rp = reconstructed (proj x wp) wp in
          let rm = reconstructed (proj x wm) wm in
          Mat.(rp + rm), rp, rm, x
        in
        let rps =
          Arr.concatenate
            ~axis:2
            (Array.mapi
               ~f:(fun i _ ->
                 let _, rp, _, _ = projections (load i) in
                 let m = Arr.reshape rp [| (Arr.shape rp).(0); (Arr.shape rp).(1); 1 |] in
                 m)
               reaches)
        and rms =
          Arr.concatenate
            ~axis:2
            (Array.mapi
               ~f:(fun i _ ->
                 let _, _, rm, _ = projections (load i) in
                 let m = Arr.reshape rm [| (Arr.shape rm).(0); (Arr.shape rm).(1); 1 |] in
                 m)
               reaches)
        in
        let vprep =
          let x = Arr.var ~axis:2 rps |> Arr.sum ~axis:1 |> Arr.to_array in
          let v = Mat.of_array x (-1) 1 in
          Mat.(v /$ max' v)
        in
        let vmov =
          let x = Arr.var ~axis:2 rms |> Arr.sum ~axis:1 |> Arr.to_array in
          let v = Mat.of_array x (-1) 1 in
          Mat.(v /$ max' v)
        in
        (vprep, vmov) :: accu)
      else accu)
  |> C.gather
  |> fun v ->
  let vpreps, vmovs =
    C.broadcast' (fun () ->
        (* v is an array of lists *)
        let v = v |> Array.to_list |> List.concat |> Array.of_list in
        let vprep = Array.map v ~f:fst |> Mat.concatenate ~axis:1 in
        let vmov = Array.map v ~f:snd |> Mat.concatenate ~axis:1 in
        vprep, vmov)
  in
  C.root_perform (fun () ->
      let vprep_mean = Mat.mean ~axis:1 vpreps in
      let vmov_mean = Mat.mean ~axis:1 vmovs in
      let vprep_std = Mat.std ~axis:1 vpreps in
      let vmov_std = Mat.std ~axis:1 vmovs in
      Mat.save_txt
        ~out:(in_dir (Printf.sprintf "analysis_%i/vprep_mean" n_dim))
        vprep_mean;
      Mat.save_txt ~out:(in_dir (Printf.sprintf "analysis_%i/vmov_mean" n_dim)) vmov_mean;
      Mat.save_txt ~out:(in_dir (Printf.sprintf "analysis_%i/vprep_std" n_dim)) vprep_std;
      Mat.save_txt ~out:(in_dir (Printf.sprintf "analysis_%i/vmov_std" n_dim)) vmov_std)
(* ( Mat.save_txt
      ~out:(in_dir (Printf.sprintf "analysis_%i/vprep" n_dim))
      Mat.(vprep /$ max' vprep)
  , Mat.save_txt
      ~out:(in_dir (Printf.sprintf "analysis_%i/vmov" n_dim))
      Mat.(vmov /$ max' vmov) ) *)
