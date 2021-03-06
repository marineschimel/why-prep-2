open Owl
open Lib
open Defaults
module AD = Algodiff.D
module M = Arm.Make (Arm.Defaults)

let _ = Printexc.record_backtrace true

let step_exp = 100

let step_go = 300

let last_step = step_go + 400 + 100

let _T = AD.F 0.4

let step_moves = 400

let n_horizon = 500

let n_planned = 20

let m = 200

let t_mat = Mat.(eye 2 *$ (Defaults.t_coeff *. 0.5)) |> AD.pack_arr

let r = AD.Maths.(AD.Mat.eye 200 * F (r_coeff *. 0.5))

let q = AD.Maths.(AD.Mat.eye 2 * F (q_coeff *. 0.5))

let __dt = AD.F sampling_dt

let target = Mat.of_arrays [| [| -0.287147; 2.39155; 0.; 0. |] |]

let tgt_pos = Mat.get_slice [ []; [ 0; 1 ] ] target |> AD.pack_arr

let rec plan_traj step traj inputs x0 us =
  let rl ~x ~u ~k =
    let thetas = unpack_pos x in
    let vel = unpack_vel x in
    let _, x_state = unpack_full_state x 4 in
    let dx_p = AD.Maths.(tgt_pos - thetas) in
    let dx_vel = vel in
    let t = AD.Maths.(__dt * F (float_of_int k)) in
    let torques = AD.Maths.(__c *@ transpose x_state) in

    if step < step_exp then
      let ds = step_exp - step in
      let t_left = AD.Maths.(F (float ds) * __dt) in

      let tau = AD.Maths.(_T + t_left) in

      AD.Maths.(
        ( sum' (u *@ r * u)
        + (sum' (dx_p *@ q * dx_p) * sigmoid ((t - tau) / F 20E-3))
        + (F 0.1 * (sum' (dx_vel *@ q * dx_vel) * sigmoid ((t - tau) / F 20E-3)))
        + sum' (transpose torques *@ t_mat * torques)
          * sigmoid ((t_left - t) / F 2E-3) )
        * __dt)
    else if step < step_go then
      let t_left = AD.F 0.02 in
      let tau = AD.Maths.(_T + t_left) in
      AD.Maths.(
        ( sum' (u *@ r * u)
        + (sum' (dx_p *@ q * dx_p) * sigmoid ((t - tau) / F 20E-3))
        + (F 0.1 * (sum' (dx_vel *@ q * dx_vel) * sigmoid ((t - tau) / F 20E-3)))
        + F 10.
          * sum' (transpose torques *@ t_mat * transpose torques)
          * sigmoid ((t_left - t) / F 2E-3) )
        * __dt)
    else
      let end_mov = step_go + step_moves in
      if step < end_mov then
        let ds = end_mov - step in
        let tau = AD.Maths.(F (float ds) * __dt) in
        AD.Maths.(
          ( sum' (u *@ r * u)
          + (sum' (dx_p *@ q * dx_p) * sigmoid ((t - tau) / F 20E-3))
          + F 0.1
            * (sum' (dx_vel *@ q * dx_vel) * sigmoid ((t - tau) / F 20E-3)) )
          * __dt)
      else
        AD.Maths.(
          ( sum' (u *@ r * u)
          + sum' (dx_p *@ q * dx_p)
          + (F 0.1 * sum' (dx_vel *@ q * dx_vel)) )
          * __dt)
  in
  let module P = struct
    let n = 204

    let n_theta = 4

    let m = 200

    let fl_x =
      let f ~k:_k ~x:_x = AD.F 0. in
      Some f

    let fl_xx =
      let f ~k:_k ~x:_x = AD.F 0. in
      Some f

    let rl_u =
      (*Need to check this*)
      let rlu ~k:_k ~x:_x ~u =
        let r = Mat.(eye m *$ r_coeff) |> AD.pack_arr in
        AD.Maths.(u *@ r * __dt)
      in
      Some rlu

    let rl_uu =
      let rluu ~k:_k ~x:_x ~u:_u =
        let ma = Mat.(eye m *$ r_coeff) |> AD.pack_arr in
        AD.Maths.(ma * __dt)
      in
      Some rluu

    let rl_ux =
      let f ~k:_k ~x:_x ~u:_u = AD.F 0. in
      Some f

    let rl_x = None

    let rl_xx =
      if step < step_exp then
        let ds = step_exp - step in
        let t_left = AD.Maths.(F (float ds) * __dt) in

        let tau = AD.Maths.(_T + t_left) in
        let rlxx ~k ~x:_x ~u:_u =
          let t = AD.Maths.(__dt * F (float_of_int k)) in
          let mu =
            AD.Maths.(
              AD.Mat.eye 2 * AD.F q_coeff * sigmoid ((t - tau) / F 20E-3) * __dt)
          in
          let mv =
            AD.Maths.(
              __dt * AD.Mat.eye 2 * F q_coeff
              * (F 0.1 * sigmoid ((t - tau) / F 20E-3)))
          in
          let mx =
            AD.Maths.(
              ( sigmoid ((t_left - t) / F 2E-3)
                * F Defaults.t_coeff
                * (transpose __c *@ __c)
              + (F Defaults.a_coeff * AD.Mat.eye m) )
              * __dt)
          in
          let mf1 =
            AD.Maths.concatenate ~axis:1
              [| mu; AD.Mat.zeros (n_theta - 2) (m + 2) |]
          in
          let mf2 =
            AD.Maths.concatenate ~axis:1
              [|
                AD.Mat.zeros (n_theta - 2) 2; mv; AD.Mat.zeros (n_theta - 2) m;
              |]
          in
          let mf3 =
            AD.Maths.concatenate ~axis:1 [| AD.Mat.zeros m n_theta; mx |]
          in
          AD.Maths.concatenate ~axis:0 [| mf1; mf2; mf3 |]
        in
        Some rlxx
      else if step < step_go then
        let t_left = AD.F 0.02 in
        let tau = AD.Maths.(_T + t_left) in
        let rlxx ~k ~x:_x ~u:_u =
          let t = AD.Maths.(__dt * F (float_of_int k)) in
          let mu =
            AD.Maths.(
              AD.Mat.eye 2 * AD.F q_coeff * sigmoid ((t - tau) / F 20E-3) * __dt)
          in
          let mv =
            AD.Maths.(
              __dt * AD.Mat.eye 2 * F q_coeff
              * (F 0.1 * sigmoid ((t - tau) / F 20E-3)))
          in
          let mx =
            AD.Maths.(
              ( sigmoid ((t_left - t) / F 2E-3)
                * F Defaults.t_coeff
                * (transpose __c *@ __c)
              + (F Defaults.a_coeff * AD.Mat.eye m) )
              * __dt)
          in
          let mf1 =
            AD.Maths.concatenate ~axis:1
              [| mu; AD.Mat.zeros (n_theta - 2) (m + 2) |]
          in
          let mf2 =
            AD.Maths.concatenate ~axis:1
              [|
                AD.Mat.zeros (n_theta - 2) 2; mv; AD.Mat.zeros (n_theta - 2) m;
              |]
          in
          let mf3 =
            AD.Maths.concatenate ~axis:1 [| AD.Mat.zeros m n_theta; mx |]
          in
          AD.Maths.concatenate ~axis:0 [| mf1; mf2; mf3 |]
        in
        Some rlxx
      else
        let end_mov = step_moves + step_go in
        if step < end_mov then
          let ds = end_mov - step in
          let tau = AD.Maths.(F (float ds) * __dt) in
          let rlxx ~k ~x:_x ~u:_u =
            let t = AD.Maths.(__dt * F (float_of_int k)) in
            let mu =
              AD.Maths.(
                AD.Mat.eye 2 * AD.F q_coeff
                * sigmoid ((t - tau) / F 20E-3)
                * __dt)
            in
            let mv =
              AD.Maths.(
                __dt * AD.Mat.eye 2 * F q_coeff
                * (F 0.1 * sigmoid ((t - tau) / F 20E-3)))
            in
            let mx = AD.Mat.zeros m m in
            let mf1 =
              AD.Maths.concatenate ~axis:1
                [| mu; AD.Mat.zeros (n_theta - 2) (m + 2) |]
            in
            let mf2 =
              AD.Maths.concatenate ~axis:1
                [|
                  AD.Mat.zeros (n_theta - 2) 2; mv; AD.Mat.zeros (n_theta - 2) m;
                |]
            in
            let mf3 =
              AD.Maths.concatenate ~axis:1 [| AD.Mat.zeros m n_theta; mx |]
            in
            AD.Maths.concatenate ~axis:0 [| mf1; mf2; mf3 |]
          in
          Some rlxx
        else
          let rlxx ~k:_k ~x:_x ~u:_u =
            let mu = AD.Maths.(AD.Mat.eye 2 * AD.F q_coeff * __dt) in
            let mv = AD.Maths.(__dt * AD.Mat.eye 2 * F q_coeff * F 0.1) in
            let mf1 =
              AD.Maths.concatenate ~axis:1
                [| mu; AD.Mat.zeros (n_theta - 2) (m + 2) |]
            in
            let mf2 =
              AD.Maths.concatenate ~axis:1
                [|
                  AD.Mat.zeros (n_theta - 2) 2; mv; AD.Mat.zeros (n_theta - 2) m;
                |]
            in
            let mf3 =
              AD.Maths.concatenate ~axis:1
                [| AD.Mat.zeros m n_theta; AD.Mat.zeros m m |]
            in
            AD.Maths.concatenate ~axis:0 [| mf1; mf2; mf3 |]
          in
          Some rlxx

    let dyn ~k:_k ~x ~u =
      let thetas, xs = unpack_full_state x n_theta in
      let xs = AD.Maths.get_slice [ []; [ 0; m - 1 ] ] xs in
      let xst = AD.Maths.transpose xs in
      let s = Arm.pack_state thetas in
      let __w = Defaults.__w in
      let dx =
        AD.Maths.(
          ( (((Defaults.__w *@ g xst) - (AD.Mat.eye m *@ xst)) / F Defaults.tau)
          + (Defaults.__b *@ transpose u) )
          * __dt)
        |> AD.Maths.transpose
      in
      let tau = AD.Maths.(__c *@ g xst) |> AD.Maths.transpose in
      let dotdot = M.theta_dot_dot s tau in
      let new_thetas =
        AD.Maths.(
          of_arrays
            [|
              [|
                s.x1 + (__dt * s.x1_dot);
                s.x2 + (__dt * s.x2_dot);
                s.x1_dot + (__dt * AD.Maths.get_item dotdot 0 0);
                s.x2_dot + (__dt * AD.Maths.get_item dotdot 1 0);
              |];
            |])
      and new_x = AD.Maths.(xs + dx) in
      AD.Maths.(concatenate ~axis:1 [| new_thetas; new_x |])

    let minv ~x =
      let thetas, _ = unpack_full_state x n_theta in
      let st = Arm.pack_state thetas in
      AD.Linalg.inv (M.inertia st)

    let ms ~x =
      let open AD.Maths in
      let thetas, _ = unpack_full_state x n_theta in
      let st = Arm.pack_state thetas in
      let m12 =
        let cst1 = AD.F M._A2 * sin st.x2 * AD.F (-1.)
        and cst2 = AD.F M._A2 * cos st.x2
        and mat2 =
          [|
            [| AD.F 0.; neg st.x2_dot * ((F 2. * st.x1_dot) + st.x2_dot) |];
            [| AD.F 0.; sqr st.x1_dot |];
          |]
          |> of_arrays
        and mat1 =
          [|
            [| AD.F 0.; (F 2. * st.x1_dot) + st.x2_dot |]; [| F 0.; st.x1_dot |];
          |]
          |> of_arrays
        in
        (cst2 * mat2) + (cst1 * mat1)
      in
      let m22 =
        let cst1 = F M._A2 * sin st.x2 in
        let mat1 =
          [|
            [| neg st.x2_dot; neg (st.x2_dot + st.x1_dot) |];
            [| AD.F 2. * st.x1_dot; AD.F 0. |];
          |]
          |> of_arrays
        in
        (cst1 * mat1) + M._B
      in
      (m12, m22)

    let _lin_dyn_x ~k:_k ~x ~u:_u =
      let nminv = AD.Maths.(neg (minv ~x)) in
      let m21, m22 = ms ~x in
      let b1 =
        AD.Maths.concatenate ~axis:1
          [| AD.Mat.zeros 2 2; AD.Mat.eye 2; AD.Mat.zeros 2 m |]
      and b2 =
        AD.Maths.(
          concatenate ~axis:1
            [| nminv *@ m21; nminv *@ m22; neg AD.Maths.(nminv *@ __c) |])
      and b3 =
        AD.Maths.concatenate ~axis:1
          [|
            AD.Mat.zeros m 2;
            AD.Mat.zeros m 2;
            AD.Maths.((Defaults.__w - AD.Mat.eye m) / F Defaults.tau);
          |]
      in
      let mat = AD.Maths.(concatenate ~axis:0 [| b1; b2; b3 |]) in
      AD.Maths.((mat * __dt) + AD.Mat.eye n) |> AD.Maths.transpose

    let _lin_dyn_u ~k:_k ~x:_x ~u:_u =
      let m =
        AD.Maths.concatenate ~axis:0
          [| AD.Mat.zeros 2 m; AD.Mat.zeros 2 m; Defaults.__b |]
      in
      AD.Maths.(transpose (m * __dt))

    let running_loss ~k ~x ~u = rl ~u ~x ~k

    let final_loss ~k:_k ~x:_x = AD.F 0.

    let dyn_x = Some _lin_dyn_x

    let dyn_u = Some _lin_dyn_u
  end in
  let module I = Ilqr.Default.Make (P) in
  let stop =
    let cprev = ref 1E9 in
    fun k us ->
      let c = I.loss x0 us in
      let pct_change = abs_float (c -. !cprev) /. !cprev in
      if step mod 10 = 0 then
        Printf.printf "loop %i | iter %i | cost %f | pct change %f\n%!" step k c
          pct_change;
      cprev := c;
      pct_change < 1E-3
  in
  let us = I.learn ~stop x0 us in
  if step > last_step then
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
    let _ =
      Mat.save_txt ~out:"traj" (AD.unpack_arr traj);
      Mat.save_txt ~out:"inputs" (AD.unpack_arr inputs)
    in
    (traj, inputs)
  else
    let plan_x =
      I.trajectory x0 us |> AD.Maths.get_slice [ [ 0; n_planned - 1 ]; [] ]
    in
    let _ =
      if step mod 100 = 0 then Mat.save_txt ~out:"plan_x" (AD.unpack_arr plan_x)
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
      (next_u0 @ List.init n_planned (fun _ -> AD.Mat.zeros 1 m))

let _ =
  let xs_0 = AD.Mat.zeros 1 m in
  let us0 = List.init n_horizon (fun _ -> AD.Mat.zeros 1 m) in
  let x0 = AD.Maths.(concatenate ~axis:1 [| initial_theta; xs_0 |]) in
  plan_traj 0 [] [] x0 us0
