open Owl
open Base

let dir = "/home/mmcs3/rds/rds-t2-cs156-T7o4pEA8QoU/mmcs3/eval_rad_abg"
let in_dir = Printf.sprintf "%s/%s" dir
let rad_w = Cmdargs.(get_float "-rad" |> force ~usage:"-rad")
let sa_w = Cmdargs.(get_float "-sa" |> force ~usage:"-sa")
let net = Cmdargs.(get_string "-net" |> force ~usage:"-net")
let seed = Cmdargs.(get_int "-seed" |> force ~usage:"-seed")
let n_target = Cmdargs.(get_int "-n_target" |> force ~usage:"-n_target")
let w = Mat.load_txt (in_dir (Printf.sprintf "w_%s_%.1f_%.1f_%i" net rad_w sa_w seed))
let c = Mat.load_txt (in_dir (Printf.sprintf "c_%s_%.1f_%.1f_%i_0" net rad_w sa_w seed))
let c = Mat.(c /$ l2norm' c)
let n = 200
let tau = 150E-3
let dt = 2E-3

(*normalize c*)
let a = Mat.((w - eye 200) /$ tau)
let n_steps = 200

let c_bar =
  let q, r, _ = Linalg.D.qr ~thin:false Mat.(transpose c) in
  Mat.(transpose (get_slice [ []; [ 2; -1 ] ] q))

let c_ort =
  let q, r, _ = Linalg.D.qr ~thin:false Mat.(transpose c) in
  Mat.(transpose (get_slice [ []; [ 0; 2 ] ] q))

let prior_sigma = Mat.(eye 2 *$ 1E-3)
let lam0 = Mat.(transpose c_bar *@ c_bar)
let sigma lam = Mat.((c *@ lam *@ transpose c) + prior_sigma)
let new_lam lam = Mat.(a *@ lam *@ transpose a)

let all_sigmas =
  let rec sigs lam ls k =
    let new_l = new_lam lam in
    let new_sigma = sigma new_l in
    if Int.(k = n_steps) then List.rev ls else sigs new_l (new_sigma :: ls) Int.(k + 1)
  in
  sigs lam0 [] 0

let torques =
  Array.init 8 ~f:(fun i ->
      Mat.(
        load_txt
          (Printf.sprintf
             "/home/mmcs3/rds/rds-t2-cs156-T7o4pEA8QoU/mmcs3/final_results/ramping_soc/seed_0_mixed/torques_%i_0"
             i)))

let compute_ll i =
  let t =
    torques.(i)
    |> fun y ->
    Mat.get_slice [ [ 0; n_steps - 1 ] ] y
    |> fun z -> Mat.map_rows (fun x -> x) z |> Array.to_list
  in
  let logdet_term = List.mapi all_sigmas ~f:(fun i s -> Linalg.D.logdet s) in
  let quad_term =
    let l =
      List.map2 all_sigmas t ~f:(fun s m ->
          let inv_s = Linalg.D.linsolve s (Mat.transpose m) in
          Mat.(0.5 *. sum' (neg (m *@ inv_s))))
    in
    match l with
    | Ok l -> l
    | _ -> []
  in
  let qt = [| Array.of_list quad_term |] |> Mat.of_arrays |> Mat.sum' in
  let lt =
    [| Array.of_list logdet_term |]
    |> Mat.of_arrays
    |> Mat.sum'
    |> fun l ->
    l +. (0.5 *. 2. *. Float.of_int n_steps *. Float.log Float.(2. *. Const.pi))
  in
  let kappa = Float.(neg qt /. Float.(of_int Int.(200 * n_steps))) in
  Float.((qt /. kappa) - lt - (of_int Int.(n_steps * 200) *. log kappa)), qt, lt

let obs_gramian a c = Linalg.D.lyapunov Mat.(transpose a) Mat.(neg (transpose c) *@ c)
let ctr_gramian a b = Linalg.D.lyapunov a Mat.(neg b *@ transpose b)
let ctr = ctr_gramian a (Mat.eye 200)
let obs_c = obs_gramian a c

let get_alpha u =
  (*u is a vector of size 1x200*)
  Mat.(u *@ ctr *@ transpose u)

let get_beta u =
  (*u is a vector of size 1x200*)
  Mat.(u *@ obs_c *@ transpose u)

let tr_obs = Mat.trace obs_c
let tr_ctr = Mat.trace ctr

let all_alphas_c =
  Array.init 2 ~f:(fun i ->
      let u = Mat.get_slice [ [ i ] ] c_ort in
      get_alpha u)
  |> fun z -> Mat.concatenate z ~axis:0

let all_betas_c =
  Array.init 2 ~f:(fun i ->
      let u = Mat.get_slice [ [ i ] ] c_ort in
      get_beta u)
  |> fun z -> Mat.concatenate z ~axis:0

let all_alphas_cbar =
  Array.init 198 ~f:(fun i ->
      let u = Mat.get_slice [ [ i ] ] c_bar in
      get_alpha u)
  |> fun z -> Mat.concatenate z ~axis:0

let all_betas_cbar =
  Array.init 198 ~f:(fun i ->
      let u = Mat.get_slice [ [ i ] ] c_bar in
      get_beta u)
  |> fun z -> Mat.concatenate z ~axis:0

let _ = Stdio.printf "%f %f %f %!" (Mat.trace obs_c) (Mat.trace ctr) (Mat.l2norm' w)
(* let _ =
  Mat.save_txt
    ~out:
      (in_dir
         (Printf.sprintf "all_alphas_c_%s_%.1f_%.1f_%i_%i" net rad_w sa_w seed n_target))
    all_alphas_c;
  Mat.save_txt
    ~out:
      (in_dir
         (Printf.sprintf
            "all_alphas_cbar_%s_%.1f_%.1f_%i_%i"
            net
            rad_w
            sa_w
            seed
            n_target))
    all_alphas_cbar;
  Mat.save_txt
    ~out:
      (in_dir
         (Printf.sprintf "all_betas_c_%s_%.1f_%.1f_%i_%i" net rad_w sa_w seed n_target))
    all_betas_c;
  Mat.save_txt
    ~out:
      (in_dir
         (Printf.sprintf "all_betas_cbar_%s_%.1f_%.1f_%i_%i" net rad_w sa_w seed n_target))
    all_betas_cbar;
  Mat.save_txt
    ~out:
      (in_dir (Printf.sprintf "tr_ctr_%s_%.1f_%.1f_%i_%i" net rad_w sa_w seed n_target))
    (Mat.of_arrays [| [| tr_ctr |] |]);
  Mat.save_txt
    ~out:
      (in_dir (Printf.sprintf "tr_obs_%s_%.1f_%.1f_%i_%i" net rad_w sa_w seed n_target))
    (Mat.of_arrays [| [| tr_obs |] |])
 *)
(* open Owl
open Base

let dir = "/home/mmcs3/rds/rds-t2-cs156-T7o4pEA8QoU/mmcs3/hyperparams_continuous_abg_3"
let in_dir = Printf.sprintf "%s/%s" dir
let rad_w = Cmdargs.(get_float "-rad" |> force ~usage:"-rad")
let sa_w = Cmdargs.(get_float "-sa" |> force ~usage:"-sa")
let net = Cmdargs.(get_string "-net" |> force ~usage:"-net")
let seed = Cmdargs.(get_int "-seed" |> force ~usage:"-seed")
let n_target = Cmdargs.(get_int "-n_target" |> force ~usage:"-n_target")

let w =
  Mat.load_txt
    (in_dir (Printf.sprintf "w_%s_%.1f_%.1f_%i_%i" net rad_w sa_w seed n_target))

let c =
  Mat.load_txt
    (in_dir (Printf.sprintf "c_%s_%.1f_%.1f_%i_%i" net rad_w sa_w seed n_target))

let n = 200
let tau = 150E-3
let dt = 2E-3

(*normalize c*)
let a = Mat.(eye n + ((w - eye n) *$ (dt /. tau)))
let n_steps = 200

let c_bar =
  let q, r, _ = Linalg.D.qr ~thin:false Mat.(transpose c) in
  Mat.(transpose (get_slice [ []; [ 2; -1 ] ] q))

let prior_sigma = Mat.(eye 2 *$ 1E-3)
let lam0 = Mat.(transpose c_bar *@ c_bar)
let sigma lam = Mat.((c *@ lam *@ transpose c) + prior_sigma)
let new_lam lam = Mat.(a *@ lam *@ transpose a)

let all_sigmas =
  let rec sigs lam ls k =
    let new_l = new_lam lam in
    let new_sigma = sigma new_l in
    if Int.(k = n_steps) then List.rev ls else sigs new_l (new_sigma :: ls) Int.(k + 1)
  in
  sigs lam0 [] 0

let torques =
  Array.init 8 ~f:(fun i ->
      Mat.(
        load_txt
          (Printf.sprintf
             "/home/mmcs3/rds/rds-t2-cs156-T7o4pEA8QoU/mmcs3/final_results/ramping_soc/seed_0_mixed/torques_%i_0"
             i)))

let compute_ll i =
  let t =
    torques.(i)
    |> fun y ->
    Mat.get_slice [ [ 0; n_steps - 1 ] ] y
    |> fun z -> Mat.map_rows (fun x -> x) z |> Array.to_list
  in
  let logdet_term = List.mapi all_sigmas ~f:(fun i s -> Linalg.D.logdet s) in
  let quad_term =
    let l =
      List.map2 all_sigmas t ~f:(fun s m ->
          let inv_s = Linalg.D.linsolve s (Mat.transpose m) in
          Mat.(0.5 *. sum' (neg (m *@ inv_s))))
    in
    match l with
    | Ok l -> l
    | _ -> []
  in
  let qt = [| Array.of_list quad_term |] |> Mat.of_arrays |> Mat.sum' in
  let lt =
    [| Array.of_list logdet_term |]
    |> Mat.of_arrays
    |> Mat.sum'
    |> fun l ->
    l +. (0.5 *. 2. *. Float.of_int n_steps *. Float.log Float.(2. *. Const.pi))
  in
  let kappa = Float.(neg qt /. Float.(of_int Int.(200 * n_steps))) in
  Float.((qt /. kappa) - lt - (of_int Int.(n_steps * 200) *. log kappa)), qt, lt

let lll =
  Array.init 1 ~f:(fun i ->
      let ll, _, _ = compute_ll n_target in
      [| ll |])
  |> Mat.of_arrays
  |> Mat.mean'

let _ = Stdio.printf "%f" lll

let _ =
  Mat.save_txt
    ~out:
      (in_dir (Printf.sprintf "g_full_%s_%.1f_%.1f_%i_%i" net rad_w sa_w seed n_target))
    (Mat.of_arrays [| [| lll |] |])
(*Array.iter
    (Array.init 8 ~f:(fun i -> i))
    ~f:(fun i ->
      let ll, _, _ = compute_ll i in
      Mat.save_txt
        ~out:(in_dir (Printf.sprintf "g_full_%s_%.1f_%.1f_%i_%i" net rad_w sa_w seed i))
        (Mat.of_arrays [| [| ll |] |]))*) *)
