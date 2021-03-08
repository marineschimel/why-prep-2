open Owl
open Lib
open Defaults
open Params
module AD = Algodiff.D
module M = Arm.Make (Arm.Defaults)

let _ = Printexc.record_backtrace true

let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")

let in_dir s = Printf.sprintf "%s/%s" dir s

let run ?target:_target ?qs_coeff:_t_coeff ?r_coeff:_r_coeff ?t_prep:_t_prep
    ?gamma:_gamma ?cost:_cost ~t_mov ~c ~w ~b subdir =
  let rec plan_traj step traj inputs x0 us =
    let module PT = struct
      let qs_coeff = match _t_coeff with Some a -> a | None -> 1000.

      let cost = match _cost with Some a -> a | None -> "running"

      let r_coeff =
        match _r_coeff with
        | Some a -> Defaults.r_coeff *. a
        | None -> Defaults.r_coeff

      let t_prep = match _t_prep with Some a -> a | None -> 0.

      let gamma_exponent = match _gamma with Some a -> a | None -> 2.

      let target_theta =
        match _target with
        | Some a -> a
        | None -> [| Mat.of_arrays [| [| -0.287147; 2.39155; 0.; 0. |] |] |]

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

      let step = step
    end in
    if step < step_go then
      let module P = Prep.Make (PT) (Mpc_cost.C_Early (PT)) in
      let module I = Ilqr.Default.Make (P) in
      let stop =
        let cprev = ref 1E9 in
        fun k us ->
          let c = I.loss x0 us in
          let pct_change = abs_float (c -. !cprev) /. !cprev in
          if step mod 10 = 0 then
            Printf.printf "loop %i | iter %i | cost %f | pct change %f\n%!" step
              k c pct_change;
          cprev := c;
          pct_change < 1E-3
      in
      let us = I.learn ~stop x0 us in
      let plan_x =
        I.trajectory x0 us |> AD.Maths.get_slice [ [ 0; n_planned - 1 ]; [] ]
      in
      let planned_traj =
        AD.unpack_arr (AD.Maths.get_slice [ []; [ 4; -1 ] ] plan_x)
      in
      let _ =
        if step = 0 then (
          Mat.save_txt ~out:(PT.saving_dir "planned") (AD.unpack_arr plan_x);

          Mat.save_txt
            ~out:(PT.saving_dir "plan_torques")
            Mat.(planned_traj *@ transpose c) )
        else
          Mat.save_txt ~append:true ~out:(PT.saving_dir "planned")
            (AD.unpack_arr plan_x);
        Mat.save_txt ~append:true
          ~out:(PT.saving_dir "plan_torques")
          Mat.(planned_traj *@ transpose c)
      in
      let next_x0 =
        I.trajectory x0 us |> AD.Maths.get_slice [ [ n_planned ]; [] ]
        (*traj inputs x0 us*)
      in
      let plan_u =
        us |> Array.of_list
        |> AD.Maths.concatenate ~axis:0
        |> AD.Maths.get_slice [ [ 0; n_planned - 1 ]; [] ]
      in
      let next_u0 =
        if Maths.abs (float (step_go - step)) < float n_planned then
          List.init n_horizon (fun _ -> AD.Mat.zeros 1 size_inputs)
        else
          Array.sub (us |> Array.of_list) n_planned (n_horizon - n_planned)
          |> Array.to_list
      in
      plan_traj (step + n_planned) (plan_x :: traj) (plan_u :: inputs) next_x0
        (next_u0 @ List.init n_planned (fun _ -> AD.Mat.zeros 1 size_inputs))
    else if step > last_step then (
      let module P = Prep.Make (PT) (Mpc_cost.C_Post (PT)) in
      let _ = Printf.printf "threshold passed %!" in
      let module I = Ilqr.Default.Make (P) in
      let stop =
        let cprev = ref 1E9 in
        fun k us ->
          let c = I.loss x0 us in
          let pct_change = abs_float (c -. !cprev) /. !cprev in
          if step mod 10 = 0 then
            Printf.printf "loop %i | iter %i | cost %f | pct change %f\n%!" step
              k c pct_change;
          cprev := c;
          pct_change < 1E-3
      in
      let us = I.learn ~stop x0 us in
      let traj =
        AD.Maths.concatenate ~axis:0
          [|
            List.rev traj |> Array.of_list |> AD.Maths.concatenate ~axis:0;
            I.trajectory x0 us;
          |]
      in
      let inputs =
        AD.Maths.concatenate ~axis:0
          [|
            List.rev inputs |> Array.of_list |> AD.Maths.concatenate ~axis:0;
            us |> Array.of_list |> AD.Maths.concatenate ~axis:0;
          |]
      in

      Mat.save_txt ~out:(PT.saving_dir "traj") (AD.unpack_arr traj);
      Mat.save_txt ~out:(PT.saving_dir "torques")
        Mat.(
          AD.unpack_arr (AD.Maths.get_slice [ []; [ 4; -1 ] ] traj)
          *@ transpose c);
      Mat.save_txt ~out:(PT.saving_dir "inputs") (AD.unpack_arr inputs)
      (* in
         (traj, inputs) *) )
    else
      let module P = Prep.Make (PT) (Mpc_cost.C_Post (PT)) in
      let _ = Printf.printf "threshold passed %!" in
      let module I = Ilqr.Default.Make (P) in
      let stop =
        let cprev = ref 1E9 in
        fun k us ->
          let c = I.loss x0 us in
          let pct_change = abs_float (c -. !cprev) /. !cprev in
          if step mod 10 = 0 then
            Printf.printf "loop %i | iter %i | cost %f | pct change %f\n%!" step
              k c pct_change;
          cprev := c;
          pct_change < 1E-3
      in
      let us = I.learn ~stop x0 us in
      let plan_x =
        I.trajectory x0 us |> AD.Maths.get_slice [ [ 0; n_planned - 1 ]; [] ]
      in
      let planned_traj =
        AD.unpack_arr (AD.Maths.get_slice [ []; [ 4; -1 ] ] plan_x)
      in
      let _ =
        Mat.save_txt ~append:true ~out:(PT.saving_dir "planned")
          (AD.unpack_arr plan_x);
        Mat.save_txt ~append:true
          ~out:(PT.saving_dir "plan_torques")
          Mat.(planned_traj *@ transpose c)
      in
      let next_x0 =
        I.trajectory x0 us |> AD.Maths.get_slice [ [ n_planned ]; [] ]
        (*traj inputs x0 us*)
      in
      let plan_u =
        us |> Array.of_list
        |> AD.Maths.concatenate ~axis:0
        |> AD.Maths.get_slice [ [ 0; n_planned - 1 ]; [] ]
      in
      let next_u0 =
        Array.sub (us |> Array.of_list) n_planned (n_horizon - n_planned)
        |> Array.to_list
      in
      plan_traj (step + n_planned) (plan_x :: traj) (plan_u :: inputs) next_x0
        (next_u0 @ List.init n_planned (fun _ -> AD.Mat.zeros 1 size_inputs))
  in
  let us0 = List.init n_horizon (fun _ -> AD.Mat.zeros 1 size_inputs) in
  let x0 =
    AD.Maths.(
      concatenate ~axis:1
        [|
          initial_theta;
          AD.Mat.zeros 1 n (*; Mat.of_array [| -0.5 |] 1 (-1) |> AD.pack_arr*);
        |])
  in
  plan_traj 0 [] [] x0 us0

let _ =
  run ~t_mov:0.3 ~t_prep:0.5 ~gamma:2.
    ~target:[| Mat.row targets 0 |]
    ~c:(Mat.load_txt "data/size_50/c_soc_50")
    ~w:(Mat.load_txt "data/size_50/w_rec_50")
    ~b:Defaults.__b ~r_coeff:0.001 ~qs_coeff:1. (Printf.sprintf "mpc")
