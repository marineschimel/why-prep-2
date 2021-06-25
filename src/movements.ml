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

let n_reaches = 48

let pos_to_angles x y =
  let ct2 =
    Maths.((sqr x +. sqr y -. sqr M._L1 -. sqr M._L2) /. (2. *. M._L1 *. M._L2))
  in
  let st2 = Maths.(sqrt (1. -. sqr ct2)) in
  let t2 = Maths.acos ct2 in
  let alpha = (M._L1 +. (M._L2 *. ct2)) /. (M._L2 *. st2)
  and bet = M._L2 *. st2
  and g = Maths.(neg (M._L1 +. (M._L2 *. ct2))) in
  let st1 = (x -. (alpha *. y)) /. ((alpha *. g) -. bet) in
  let ct1 = Maths.(sqrt (1. -. sqr st1)) in
  let t1 = Maths.atan (st1 /. ct1) in
  (t1, t2)

let new_targets n radius =
  let module Task = struct
    let angle_i = 0.

    let angle_f = 360.

    let reach_angles ?(angle_i = angle_i) ?(angle_f = angle_f) n_angles =
      Mat.linspace angle_i angle_f n_angles |> Mat.to_array
  end in
  let reach_angles = Task.reach_angles n in
  Array.map
    (fun reach_angle ->
      let hand_trajectory, _ =
        M.straight_reach ~tau:0.140 ~dt:sampling_dt ~duration:0.5
          ~angle:(Arm.deg reach_angle) ~radius
      in
      hand_trajectory |> Arm.unpack_sequence |> Mat.get_fancy [ R [ -1 ]; R [] ]
      |> fun z ->
      pos_to_angles (Mat.get z 0 0) (Mat.get z 0 2) |> fun (x, y) ->
      [| x; y; 0.; 0. |])
    reach_angles

(* let targets =
  let tgts_1 = Mat.of_arrays (new_targets 48 0.08) in
  let tgts_2 = Mat.of_arrays (new_targets 48 0.12) in
  Mat.concatenate ~axis:0 [| tgts_1; tgts_2 |]

let _ = Mat.save_txt ~out:"data_arm/targets" targets *)

let targets = Mat.load_txt "data_arm/targets"

let evaluate ?target:_target ?qs_coeff:_t_coeff ?r_coeff:_r_coeff
    ?t_prep:_t_prep ?gamma:_gamma ?x0:_x0 ?weighing_pm:_wpm ~t_mov ~c ~w ~b
    ~reach_num =
  let module PT = struct
    let qs_coeff = match _t_coeff with Some a -> a | None -> 1000.

    let cost = "running"

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

    let t_mov = t_mov

    let duration = t_prep +. t_mov +. 0.2

    let saving_dir = Printf.sprintf "%s/%s" dir

    let c = c

    let __c = AD.pack_arr c

    let b = b

    let aug = 0
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
  let us = List.init P.n_steps (fun _ -> AD.Mat.zeros 1 size_inputs) in
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
          ~out:(PT.saving_dir (sprintf "x_%i" reach_num))
          (Mat.get_slice [ []; [ 4; -1 ] ] traj);
        Mat.save_txt
          ~out:(PT.saving_dir (sprintf "thetas_%i" reach_num))
          (Mat.get_slice [ []; [ 0; 3 ] ] traj);
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
        |> Mat.save_txt ~out:(PT.saving_dir (sprintf "hands_%i" reach_num));
        inputs
        |> Mat.save_txt ~out:(PT.saving_dir (sprintf "results_us_%i" reach_num));
        let t = T.now () in
        let dt = T.diff t t0 in
        Printf.printf "Time iter %i | %s | %!" k (T.Span.to_string dt) );
      recurrent
      |> Mat.save_txt ~out:(PT.saving_dir (sprintf "recurrent_us_%i" reach_num));
      pct_change < 1E-2
  in
  let final_us = I.learn ~stop x0 us in
  final_us

let _ =
  let w = Mat.load_txt (Printf.sprintf "data_arm/w_rec_100") in
  let c = Mat.gaussian 2 100 in
  let _ = Mat.save_txt ~out:"data_arm/c" c in
  (* let c = Mat.load_txt (Printf.sprintf "data_arm/c") in *)
  let x0 =
    AD.Maths.(concatenate ~axis:1 [| initial_theta; AD.Mat.zeros 1 n |])
  in
  let _ = Printf.printf "checkpoint %!" in
  Array.mapi
    (fun i tgt ->
      evaluate ~x0 ~t_mov:0.3 ~t_prep:0.4 ~gamma:2. ~target:[| tgt |] ~c ~w
        ~b:Defaults.__b ~r_coeff:0.1 ~qs_coeff:10. ~weighing_pm:(1., 10000.)
        ~reach_num:i)
    (Array.init 2 (fun i -> Mat.row targets i))
