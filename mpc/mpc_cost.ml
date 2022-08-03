open Owl
module AD = Algodiff.D
open Lib
open Defaults
open Params
open Typ

module C_Early (P : Prms) = struct
  (*default setting of qcoeff = 1E-7*)
  let n_theta = 4
  let q = Mat.(eye 2 *$ (q_coeff *. 0.5)) |> AD.pack_arr
  let r = Mat.(eye size_net *$ (P.r_coeff *. 0.5)) |> AD.pack_arr
  let t_mat = Mat.(eye 2 *$ (P.qs_coeff *. 0.5)) |> AD.pack_arr
  let a_mat = Mat.(eye size_net *$ (Defaults.a_coeff *. 0.5)) |> AD.pack_arr
  let __dt = AD.F sampling_dt
  let target = P.target_theta.(0) |> AD.pack_arr
  let tgt_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] target
  let in_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] initial_theta
  let q_start = Mat.(eye 2 *$ (P.qs_coeff *. 0.5)) |> AD.pack_arr
  let q_start_coeff = P.qs_coeff
  let step = P.step

  let t_left =
    let ds = step_exp - step in
    if float ds > 0.
    then float ds *. sampling_dt
    else sampling_dt *. float (step_max - step) /. 2.

  let tau = AD.F (t_left +. t_mov)

  let cost ~u ~x ~k =
    let _ = if k mod 1000 = 0 then Printf.printf "T left %f step %i \n %!" t_left step in
    let thetas = unpack_pos x in
    let vel = unpack_vel x in
    let _, x_state = unpack_full_state x 4 in
    let dx_p = AD.Maths.(tgt_pos - thetas) in
    let dx_vel = vel in
    let _dx_start = AD.Maths.(in_pos - thetas) in
    let t = AD.Maths.(__dt * F (float_of_int k)) in
    let torques = AD.Maths.(P.__c *@ transpose x_state) in
    let start = AD.Maths.(sigmoid ((F t_left - t) / F 2E-4)) in
    let u = AD.Maths.(u *@ transpose P.b) in
    AD.Maths.(
      (sum' (u *@ r * u)
      + (sum' (dx_p *@ q * dx_p) * sigmoid ((t - tau) / F 20E-3))
      + (F 0.1 * (sum' (dx_vel *@ q * dx_vel) * sigmoid ((t - tau) / F 20E-3)))
      + ((F 1. * (sum' (dx_vel *@ q_start * dx_vel) * start))
        + (sum' (_dx_start *@ q_start * _dx_start) * start)
        + (sum' (transpose torques *@ t_mat * transpose torques)
          * sigmoid ((F t_left - t) / F 2E-3))))
      * __dt)

  let rl_u =
    (*Need to check this*)
    let rlu ~k:_ ~x:_x ~u =
      let r = Mat.(eye size_net *$ P.r_coeff) |> AD.pack_arr in
      let u = AD.Maths.(u *@ transpose P.b) in
      AD.Maths.(u *@ r *@ P.b * __dt)
    in
    Some rlu

  let rl_x =
    let rlx ~k ~x ~u:_u =
      let _, x_state = unpack_full_state x 4 in
      let thetas = unpack_pos x in
      let t = AD.Maths.(__dt * F (float_of_int k)) in
      let start = AD.Maths.(sigmoid ((F t_left - t) / F 2E-4)) in
      let dx_p = AD.Maths.(tgt_pos - thetas) in
      let __c = AD.pack_arr P.c in
      let dx_vel = unpack_vel x in
      let dx_start = AD.Maths.(in_pos - thetas) in
      let r_xstate =
        AD.Maths.(
          x_state
          *@ transpose __c
          *@ t_mat
          *@ __c
          * F 2.
          * sigmoid ((F t_left - t) / F 2E-3))
      in
      let r_xp1 = AD.Maths.(AD.F 2. * neg dx_p *@ q * sigmoid ((t - tau) / F 20E-3)) in
      let r_xp2 = AD.Maths.(AD.F 2. * neg dx_start *@ q_start * start) in
      let r_xv1 = AD.Maths.(dx_vel *@ q * sigmoid ((t - tau) / F 20E-3) * F 0.2) in
      let r_xv2 = AD.Maths.(dx_vel *@ q_start * start * F 2.) in
      AD.Maths.(concatenate ~axis:1 [| r_xp1 + r_xp2; r_xv1 + r_xv2; r_xstate |] * __dt)
    in
    Some rlx

  let rl_uu =
    let rluu ~k:_ ~x:_x ~u:_u =
      let ma = Mat.(eye size_net *$ P.r_coeff) |> AD.pack_arr in
      AD.Maths.(P.b *@ ma *@ P.b * __dt)
    in
    Some rluu

  let rl_ux =
    let f ~k:_k ~x:_x ~u:_u = AD.F 0. in
    Some f

  let rl_xx =
    let rlxx ~k ~x:_x ~u:_u =
      let t = AD.Maths.(__dt * F (float_of_int k)) in
      let start = AD.Maths.(sigmoid ((F t_left - t) / F 2E-4)) in
      let mu =
        AD.Maths.(
          (AD.Mat.eye 2 * AD.F Defaults.q_coeff * __dt * sigmoid ((t - tau) / F 20E-3))
          + (AD.F P.qs_coeff * AD.Mat.eye 2 * start * __dt))
      in
      let mv =
        AD.Maths.(
          (__dt
          * AD.Mat.eye 2
          * F Defaults.q_coeff
          * (F 0.1 * sigmoid ((t - tau) / F 20E-3)))
          + (__dt * AD.Mat.eye 2 * F P.qs_coeff * start * F 1.))
      in
      let mx =
        AD.Maths.(
          ((sigmoid ((F t_left - t) / F 2E-3) * (transpose P.__c *@ t_mat * F 2. *@ P.__c))
          + (F Defaults.a_coeff * AD.Mat.eye size_net))
          * __dt)
      in
      let mf1 =
        AD.Maths.concatenate ~axis:1 [| mu; AD.Mat.zeros (n_theta - 2) (size_net + 2) |]
      in
      let mf2 =
        AD.Maths.concatenate
          ~axis:1
          [| AD.Mat.zeros (n_theta - 2) 2; mv; AD.Mat.zeros (n_theta - 2) size_net |]
      in
      let mf3 = AD.Maths.concatenate ~axis:1 [| AD.Mat.zeros size_net n_theta; mx |] in
      AD.Maths.concatenate ~axis:0 [| mf1; mf2; mf3 |]
    in
    Some rlxx

  let final_cost ~x ~k:_k =
    let q = Owl.Mat.(eye 2 *$ Defaults.q_coeff) |> AD.pack_arr in
    let fl =
      let thetas = unpack_pos x in
      let thetas_dot = unpack_vel x in
      let dx_p = AD.Maths.(tgt_pos - thetas)
      and dx_vel = thetas_dot in
      AD.(Maths.((sum' (dx_p *@ q * dx_p) * F 0.) + (F 0. * sum' (dx_vel *@ q * dx_vel))))
    in
    fl

  let fl_x = None

  (*let f ~k:_k ~x:_x = AD.F 0. in
    Some f*)

  let fl_xx = None
end

module C_Post (P : Prms) = struct
  (*default setting of qcoeff = 1E-7*)
  let n_theta = 4
  let q = Mat.(eye 2 *$ (q_coeff *. 0.5)) |> AD.pack_arr
  let r = Mat.(eye size_net *$ (P.r_coeff *. 0.5)) |> AD.pack_arr
  let t_mat = Mat.(eye 2 *$ (P.qs_coeff *. 0.5)) |> AD.pack_arr
  let a_mat = Mat.(eye size_net *$ (Defaults.a_coeff *. 0.5)) |> AD.pack_arr
  let __dt = AD.F sampling_dt
  let target = P.target_theta.(0) |> AD.pack_arr
  let tgt_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] target
  let in_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] initial_theta
  let q_start = Mat.(eye 2 *$ (P.qs_coeff *. 0.5)) |> AD.pack_arr
  let q_start_coeff = P.qs_coeff
  let step = P.step
  let ds = end_mov - step
  let tau = AD.Maths.(F (float ds) * __dt)

  let cost ~u ~x ~k =
    let thetas = unpack_pos x in
    let vel = unpack_vel x in
    let dx_p = AD.Maths.(tgt_pos - thetas) in
    let dx_vel = vel in
    let t = AD.Maths.(__dt * F (float_of_int k)) in
    AD.Maths.(
      (sum' (u *@ r * u)
      + (sum' (dx_p *@ q * dx_p) * sigmoid ((t - tau) / F 20E-3))
      + (F 0.1 * (sum' (dx_vel *@ q * dx_vel) * sigmoid ((t - tau) / F 20E-3))))
      * __dt)

  let rl_u =
    (*Need to check this*)
    let _rlu ~k ~x:_x ~u =
      let _t = AD.Maths.(__dt * F (float_of_int k)) in
      let r = Mat.(eye size_net *$ P.r_coeff) |> AD.pack_arr in
      let u = AD.Maths.(u *@ transpose P.b) in
      AD.Maths.(u *@ r *@ P.b * __dt)
    in
    None

  let rl_x =
    let _rlx ~k ~x ~u:_u =
      let thetas = unpack_pos x in
      let t = AD.Maths.(__dt * F (float_of_int k)) in
      let dx_p = AD.Maths.(tgt_pos - thetas) in
      let __c = AD.pack_arr P.c in
      let dx_vel = unpack_vel x in
      let r_xp1 = AD.Maths.(AD.F 2. * neg dx_p *@ q * sigmoid ((t - tau) / F 20E-3)) in
      let r_xv1 = AD.Maths.(dx_vel *@ q * sigmoid ((t - tau) / F 20E-3) * F 0.2) in
      AD.Maths.(concatenate ~axis:1 [| r_xp1; r_xv1; AD.Mat.zeros 1 P.m |] * __dt)
    in
    None

  let rl_uu =
    let _rluu ~k ~x:_x ~u:_u =
      let _t = AD.Maths.(__dt * F (float_of_int k)) in
      let ma = Mat.(eye size_net *$ P.r_coeff) |> AD.pack_arr in
      AD.Maths.(P.b *@ ma *@ P.b * __dt)
    in
    None

  let rl_ux =
    let f ~k:_k ~x:_x ~u:_u = AD.F 0. in
    Some f

  let rl_xx =
    let _rlxx ~k ~x:_x ~u:_u =
      let t = AD.Maths.(__dt * F (float_of_int k)) in
      let start = AD.F 0. in
      let mu =
        AD.Maths.(
          (AD.Mat.eye 2 * AD.F Defaults.q_coeff * __dt * sigmoid ((t - tau) / F 20E-3))
          + (AD.F P.qs_coeff * AD.Mat.eye 2 * start * __dt))
      in
      let mv =
        AD.Maths.(
          (__dt
          * AD.Mat.eye 2
          * F Defaults.q_coeff
          * (F 0.1 * sigmoid ((t - tau) / F 20E-3)))
          + (__dt * AD.Mat.eye 2 * F P.qs_coeff * start * F 1.))
      in
      let mf1 =
        AD.Maths.concatenate ~axis:1 [| mu; AD.Mat.zeros (n_theta - 2) (size_net + 2) |]
      in
      let mf2 =
        AD.Maths.concatenate
          ~axis:1
          [| AD.Mat.zeros (n_theta - 2) 2; mv; AD.Mat.zeros (n_theta - 2) size_net |]
      in
      let mf3 =
        AD.Maths.concatenate ~axis:1 [| AD.Mat.zeros size_net (size_net + n_theta) |]
      in
      AD.Maths.concatenate ~axis:0 [| mf1; mf2; mf3 |]
    in
    None

  let final_cost ~x ~k:_k =
    let q = Owl.Mat.(eye 2 *$ Defaults.q_coeff) |> AD.pack_arr in
    let fl =
      let thetas = unpack_pos x in
      let thetas_dot = unpack_vel x in
      let dx_p = AD.Maths.(tgt_pos - thetas)
      and dx_vel = thetas_dot in
      AD.(Maths.((sum' (dx_p *@ q * dx_p) * F 0.) + (F 0. * sum' (dx_vel *@ q * dx_vel))))
    in
    fl

  let fl_x = None

  (*let f ~k:_k ~x:_x = AD.F 0. in
    Some f*)

  let fl_xx = None
end
