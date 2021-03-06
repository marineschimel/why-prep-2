open Owl
open Lib
open Defaults
module A = Analysis

let _ = Printexc.record_backtrace true

let tdir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")

let tprep_dir s = Printf.sprintf "%s/%s" tdir s

let windows = [| 0.; 0.001; 0.05; 0.15; 0.3 |]

let t_preps = [| 0.05; 0.1; 0.3; 0.5; 0.6; 0.8 |]

let get_ac ~win ~tprep ~seed =
  let traj =
    Mat.load_txt
      (tprep_dir
         (Printf.sprintf "seed_%i/traj_%i" seed (int_of_float (tprep *. 1000.))))
  in
  let nprep = int_of_float (1000. *. tprep) in
  let nwin = int_of_float (1000. *. win) in
  let nm = nprep + nwin in
  let y = Mat.get_slice [ [ nprep; nm ]; [ 4; -1 ] ] traj in
  Mat.reshape y [| 1; Mat.row_num y * Mat.col_num y |]

let dt = 1E-3

let init_theta = Mat.of_arrays [| [| 0.174533; 2.50532; 0.; 0. |] |]

let target_theta = Mat.of_arrays [| [| -0.287147; 2.39155; 0.; 0. |] |]

let _ =
  let res =
    Array.map
      (fun t_prep ->
        let traj =
          Mat.load_txt
            (tprep_dir
               (Printf.sprintf "seed_1/traj_%i"
                  (int_of_float (t_prep *. 1000.))))
        in
        let inputs =
          Mat.load_txt
            (tprep_dir
               (Printf.sprintf "seed_1/results_us_%i"
                  (int_of_float (t_prep *. 1000.))))
        in
        let r = Mat.get (Mat.load_txt (tprep_dir "seed_1/loss_time")) 0 1 in
        let input_energy = Mat.(l2norm_sqr' inputs) *. r in
        let torques =
          Mat.load_txt
            (tprep_dir
               (Printf.sprintf "seed_1/torques_%i"
                  (int_of_float (t_prep *. 1000.))))
        in
        let accuracy =
          let t1 =
            let dp =
              Mat.(
                get_slice [ []; [ 0; 1 ] ] traj
                - get_slice [ []; [ 0; 1 ] ] target_theta)
              |> Mat.l2norm_sqr ~axis:1
            in
            Mat.mapi
              (fun i x ->
                let t = float_of_int i *. dt in
                Maths.sigmoid ((t -. (t_prep +. 0.4)) /. 20E-3) *. x)
              dp
            |> Mat.sum'
          in
          let t2 =
            let dv =
              Mat.(
                get_slice [ []; [ 2; 3 ] ] traj
                - get_slice [ []; [ 2; 3 ] ] target_theta)
              |> Mat.l2norm_sqr ~axis:1
            in
            Mat.mapi
              (fun i x ->
                let t = float_of_int i *. dt in
                Maths.sigmoid ((t -. (t_prep +. 0.4)) /. 20E-3) *. x)
              dv
            |> Mat.sum'
          in
          t1 +. (0.1 *. t2)
        in
        let no_mov =
          let t1 =
            let dp =
              Mat.(
                get_slice [ []; [ 0; 1 ] ] traj
                - get_slice [ []; [ 0; 1 ] ] init_theta)
              |> Mat.l2norm_sqr ~axis:1
            in
            Mat.mapi
              (fun i x ->
                let t = float_of_int i *. dt in
                Maths.(sigmoid ((t_prep -. t) /. 2E-4)) *. x)
              dp
            |> Mat.sum'
          in
          let t2 =
            let dv =
              Mat.(get_slice [ []; [ 2; 3 ] ] traj) |> Mat.l2norm_sqr ~axis:1
            in
            Mat.mapi
              (fun i x ->
                let t = float_of_int i *. dt in
                Maths.(sigmoid ((t_prep -. t) /. 2E-4)) *. x)
              dv
            |> Mat.sum'
          in
          let t3 =
            let dto =
              Mat.(
                get_slice
                  [ [ 0; int_of_float (t_prep *. 1000.) ]; [ 0; 1 ] ]
                  torques)
              |> Mat.l2norm_sqr ~axis:1
            in
            Mat.mapi
              (fun i x ->
                let t = float_of_int i *. dt in
                Maths.(sigmoid ((t_prep -. t) /. 2E-4)) *. x)
              dto
            |> Mat.sum'
          in
          t1 +. t2 +. t3
        in
        [|
          t_prep;
          0.5 *. accuracy *. dt;
          input_energy /. 2. *. dt;
          0.5 *. no_mov *. dt;
        |])
      t_preps
  in
  Mat.save_txt ~out:(tprep_dir "summary") Mat.(of_arrays res)

let w = Mat.load_txt "reach_1/w"

let c = Mat.load_txt "reach_1/c"

let top_modes_obs =
  let obs_gramian a c =
    Linalg.D.lyapunov Mat.(transpose a) Mat.(neg (transpose c) *@ c)
  in
  let a = Mat.((w - eye size_net) /$ tau) in
  let obs_modes, obs_eigs, _ = Linalg.D.svd (obs_gramian a c) in
  let _ = Mat.save_txt ~out:(tprep_dir "obs_eigs") obs_eigs in
  Mat.get_slice [ []; [ 0; 10 ] ] obs_modes

let top_modes_ctrl =
  let ctrl_gramian =
    Linalg.D.lyapunov
      Mat.(transpose Mat.(w - eye size_net) /$ tau)
      Mat.(neg (transpose (eye size_net)) *@ eye n)
  in
  let ctrl_modes, _, _ = Linalg.D.svd ctrl_gramian in
  Mat.get_slice [ []; [ 0; 10 ] ] ctrl_modes

let corr ~win ~seed =
  let tps = Mat.of_array t_preps (-1) 1 in
  let t_preps =
    let idces = Mat.filter (fun x -> x >= win) tps |> Array.to_list in
    Mat.get_fancy [ L idces ] tps |> Mat.to_array
  in
  let ar =
    Array.map
      (fun tprep ->
        let ac = get_ac ~win ~tprep ~seed in
        let x = Mat.(ac - (0. $* mean ~axis:1 ac)) in
        Mat.to_array Mat.(x / l2norm ~axis:1 x))
      t_preps
  in
  Mat.of_arrays ar

(* let _ = Array.map (fun window -> let c = (corr ~win:window) in Mat.save_txt ~out:(tprep_dir (Printf.sprintf "correlation_prep_%i" (int_of_float (window*.1000.))))
Mat.(c*@(transpose c))) windows  *)

let fun_pc i1 i2 window =
  let c1 = corr ~win:window ~seed:i1 in
  let c2 = corr ~win:window ~seed:i2 in
  Mat.(c1 *@ transpose c2)

let _ = Mat.print (fun_pc 1 2 0.01)

let pair_corr i1 i2 =
  Array.map
    (fun window ->
      let c1 = corr ~win:window ~seed:i1 in
      let c2 = corr ~win:window ~seed:i2 in
      Mat.save_txt
        ~out:
          (tprep_dir
             (Printf.sprintf "correlation_mov_%i_%i%i"
                (int_of_float (window *. 1000.))
                i1 i2))
        Mat.(c1 *@ transpose c2))
    windows

let av_corrs i_min i_max wdw nt =
  let rec tot_corrs tot i j acc =
    if i = i_max then Mat.(tot /$ acc)
    else if j = i_max then
      tot_corrs Mat.(tot + fun_pc i j wdw) (i + 1) (i + 2) (acc +. 1.)
    else tot_corrs Mat.(tot + fun_pc i j wdw) i (j + 1) (acc +. 1.)
  in
  tot_corrs (Mat.zeros nt nt) i_min (i_min + 1) 0.

let _ =
  Array.map
    (fun win ->
      Mat.save_txt
        ~out:(tprep_dir (Printf.sprintf "av_correlation_times_win_%.3f" win))
        (av_corrs 1 5 win (Array.length t_preps)))
    [| 0.01 |]

(* 
let prod_subspace_obs x =
  let z = Mat.mean ~axis:0 x in 
Mat.(((z)/(l2norm ~axis:1 z ))*@top_modes_obs) *)

(* let prod_subspace_ctrl x =
  let z = Mat.mean ~axis:0 x in 
  Mat.(((z)/(l2norm ~axis:1 z ))*@top_modes_ctrl)

let prod_obs ~win = let y =  (Array.map (fun t -> let x = Mat.load_txt (tprep_dir (Printf.sprintf "traj_%i"(int_of_float (1000.*.t)))) |> fun z -> Mat.get_slice [[(int_of_float (t*.1000.));win+(int_of_float (t*.1000.))];[4;-1]] z in Mat.to_array (prod_subspace_obs x)) t_preps) |> Mat.of_arrays 
in Mat.save_txt ~out:(tprep_dir (Printf.sprintf "obs_occ_%i" win)) (Mat.abs y)

let prod_ctrl ~win = let y =  (Array.map (fun t -> let x = Mat.load_txt (tprep_dir (Printf.sprintf "traj_%i"(int_of_float (1000.*.t)))) |> fun z -> Mat.get_slice [[(int_of_float (t*.1000.));win+(int_of_float (t*.1000.))];[4;-1]] z in Mat.to_array (prod_subspace_ctrl x)) t_preps) |> Mat.of_arrays 
in Mat.save_txt ~out:(tprep_dir (Printf.sprintf "/ctrl_occ_%i" win)) (Mat.abs y)

let l2_inpt_ratio = [|(Array.map (fun t ->  let inpts = Mat.load_txt (tprep_dir (Printf.sprintf "results_us_%i" (int_of_float (1000.*.t)))) in let cp,cm = A.cost_inputs ~t_prep:t inpts in cp/.cm) t_preps )|] |> Mat.of_arrays |> fun f -> Mat.save_txt ~out:(tprep_dir "l2_inp_ratio") f
let _ = prod_obs ~win:0; prod_obs ~win:10; prod_obs ~win:50; prod_obs ~win:400; prod_ctrl ~win:0;prod_ctrl ~win:10; prod_ctrl ~win:50; prod_ctrl ~win:400 *)

let _ = Algodiff.D.Linalg.lq (AD.Mat.gaussian 10 10)
