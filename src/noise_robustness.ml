open Owl
module AD = Algodiff.D
module M = Arm.Make (Arm.Defaults)
open Lib
module T = Core.Time
open Defaults
open Printf
module A = Analysis_funs

let _ = Printexc.record_backtrace true

let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")

let in_dir s = Printf.sprintf "%s/%s" dir s

let rinv =
  Cmdargs.(get_int "-rinv" |> force ~usage:"-rinv [inverse of r_coeff]")

let evaluate ?target:_target ?qs_coeff:_t_coeff ?r_coeff:_r_coeff
    ?t_prep:_t_prep ?gamma:_gamma ?cost:_cost ?x0:_x0 ?weighing_pm:_wpm ~t_mov
    ~noise ~c ~w ~b subdir =
  let module PT = struct
    let qs_coeff = match _t_coeff with Some a -> a | None -> 1000.

    let cost = match _cost with Some a -> a | None -> "running"

    let r_coeff =
      match _r_coeff with
      | Some a -> Defaults.r_coeff *. a
      | None -> Defaults.r_coeff

    let t_prep = match _t_prep with Some a -> a | None -> 0.

    let gamma_exponent = match _gamma with Some a -> a | None -> 2.

    let x0 =
      match _x0 with
      | Some x0 -> x0
      | None ->
          AD.Maths.(
            concatenate ~axis:1
              [|
                initial_theta;
                AD.Mat.zeros 1 n
                (*; Mat.of_array [| -0.5 |] 1 (-1) |> AD.pack_arr*);
              |])

    let target_theta =
      match _target with
      | Some a -> a
      | None -> [| Mat.of_arrays [| [| -0.287147; 2.39155; 0.; 0. |] |] |]

    let weighing_pm =
      match _wpm with Some (wp, wm) -> (wp, wm) | None -> (1., 1.)

    let w = w

    let n = size_net + 4

    let m = size_inputs

    let aug = 0

    let t_mov = t_mov

    let duration = t_prep +. t_mov +. 0.2

    let saving_dir = Printf.sprintf "%s/%s/%s" dir subdir

    let c = c

    let __c = AD.pack_arr c

    let b = b
  end in
  let module P = Prep.Make (PT) (Costs.C_Running (PT)) in
  let module I = Ilqr.Default.Make (P) in
  let __w = AD.pack_arr PT.w in
  let angles_to_x x0 us =
    let traj = I.trajectory x0 us in
    AD.Mat.map_by_row
      (fun x -> Arm.unpack_state_diff (M.hand_of (Arm.pack_state x)))
      traj
  in
  let x0 = PT.x0 in
  let us0 =
    Mat.load_txt
      (Printf.sprintf "results/noise/soc_r_%i/results_us_%i" rinv
         (int_of_float (PT.t_prep *. 1000.)))
  in
  let us =
    List.init P.n_steps (fun i ->
        let u = Mat.get_slice [ [ i ]; [] ] us0 in
        AD.pack_arr Mat.(u + (noise $* u * gaussian 1 PT.m)))
  in
  I.trajectory x0 us |> AD.unpack_arr
  |> Mat.save_txt ~out:(PT.saving_dir "traj0");
  angles_to_x x0 us |> AD.unpack_arr |> Mat.save_txt ~out:(PT.saving_dir "x0");
  let t0 = T.now () in
  let stop =
    let cprev = ref 1E9 in
    fun k us ->
      let c = I.loss x0 us in
      let pct_change = abs_float (c -. !cprev) /. !cprev in
      let traj_ad = I.trajectory x0 us in
      let traj = traj_ad |> AD.unpack_arr in
      let inputs =
        us |> Array.of_list |> AD.Maths.concatenate ~axis:0 |> AD.unpack_arr
      in
      let recurrent =
        AD.Maths.(
          transpose
            ( __w
            *@ g (transpose (get_slice [ [ 0; -2 ]; [ 4; pred PT.n ] ] traj_ad))
            ))
        |> AD.unpack_arr
      in
      if k mod 1 = 0 then (
        Printf.printf "iter %i | cost %f | pct change %f\n%!" k c pct_change;
        cprev := c;
        Mat.save_txt
          ~out:
            (PT.saving_dir
               (sprintf "traj_%i" (int_of_float (PT.t_prep *. 1000.))))
          traj;
        Mat.save_txt
          ~out:
            (PT.saving_dir
               (sprintf "floored_traj_%i" (int_of_float (PT.t_prep *. 1000.))))
          (Mat.map (fun y -> if y > 0. then y else 0.) traj);
        Mat.save_txt
          ~out:
            (PT.saving_dir
               (sprintf "torques_%i" (int_of_float (PT.t_prep *. 1000.))))
          (Mat.concatenate ~axis:0
             (Mat.map_rows
                (fun x ->
                  Mat.(
                    transpose
                      (PT.c *@ AD.unpack_arr (g (AD.pack_arr (transpose x))))))
                (Mat.get_slice [ []; [ 4; PT.n - 1 ] ] traj)));
        AD.unpack_arr (angles_to_x x0 us)
        |> Mat.save_txt
             ~out:
               (PT.saving_dir
                  (sprintf "hands_%i" (int_of_float (PT.t_prep *. 1000.))));
        inputs
        |> Mat.save_txt
             ~out:
               (PT.saving_dir
                  (sprintf "results_us_%i" (int_of_float (PT.t_prep *. 1000.))));
        let t = T.now () in
        let dt = T.diff t t0 in
        Printf.printf "Time iter %i | %s | %!" k (T.Span.to_string dt) );
      recurrent
      |> Mat.save_txt
           ~out:
             (PT.saving_dir
                (sprintf "recurrent_us_%i" (int_of_float (PT.t_prep *. 1000.))));
      if pct_change > 0. then (
        Mat.(
          save_txt ~append:true
            ~out:(PT.saving_dir "loss_time")
            (Mat.of_array
               [| PT.t_prep; PT.r_coeff; I.loss x0 us; nonnormality w |]
               1 (-1)));
        let ip, im = A.inputs ~t_prep:PT.t_prep inputs in
        let cp, cm = A.cost_inputs ~t_prep:PT.t_prep inputs in
        Mat.save_txt ~append:true
          ~out:(PT.saving_dir "input_distrib")
          (Mat.of_arrays [| [| PT.t_prep; ip; im |] |]);
        Mat.save_txt ~append:true
          ~out:(PT.saving_dir "sqr_input_distrib")
          (Mat.of_arrays [| [| PT.t_prep; cp; cm |] |]);
        Mat.(
          save_txt ~append:true ~out:(PT.saving_dir "sum_us")
            (Mat.of_array
               [|
                 PT.t_prep;
                 PT.r_coeff;
                 inputs |> Mat.map (fun x -> Maths.abs x *. r_coeff) |> Mat.sum';
               |]
               1 (-1)));
        Mat.(
          save_txt ~append:true ~out:(PT.saving_dir "accuracy")
            (Mat.of_array
               [| PT.r_coeff; PT.t_prep; A.accuracy_theta traj |]
               1 (-1))) );
      pct_change > 0.
  in
  let final_us = I.learn ~stop x0 us in
  let energy =
    List.fold_left
      (fun acc x ->
        let x = AD.unpack_arr x in
        acc +. (Mat.(sqr x |> sum') *. sampling_dt))
      0. final_us
  in
  energy

let targets = Mat.load_txt "data/target_thetas"

let run_run =
  Array.map
    (fun (pc, mc) ->
      let w = Mat.load_txt "results/reach_1/w" in
      let c = Mat.load_txt "results/reach_1/c" in
      let x0 =
        AD.Maths.(concatenate ~axis:1 [| initial_theta; AD.Mat.zeros 1 n |])
      in
      Array.mapi
        (fun _ noise ->
          evaluate ~x0 ~t_mov:0.4 ~t_prep:0.3 ~gamma:2.
            ~target:[| Mat.row targets 0 |]
            ~noise ~c ~w ~b:Defaults.__b
            ~r_coeff:(1. /. float rinv)
            ~qs_coeff:1. ~weighing_pm:(pc, mc)
            (Printf.sprintf "noise/soc_r_%i/level_%.3f" rinv noise))
        [| 0.; 0.001; 0.01; 0.1; 1.; 2.; 10. |])
    [| (1., 1.) |]