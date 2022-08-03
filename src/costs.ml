open Owl
module AD = Algodiff.D
open Lib
open Defaults
open Typ

let _ = Printexc.record_backtrace true

module C_Running_4D (P : Prms) = struct
  (*default setting of qcoeff = 1E-7*)
  let n_theta = 4
  let q = Mat.(eye 2 *$ (q_coeff *. 0.5)) |> AD.pack_arr
  let r = Mat.(eye size_net *$ (P.r_coeff *. 0.5)) |> AD.pack_arr
  let t_mat = Mat.(eye 2 *$ (P.qs_coeff *. 0.5)) |> AD.pack_arr
  let a_mat = Mat.(eye size_net *$ (Defaults.a_coeff *. 0.5)) |> AD.pack_arr
  let tau = AD.F (P.t_prep +. P.t_mov)
  let __dt = AD.F sampling_dt
  let target = P.target_theta.(0) |> AD.pack_arr
  let tgt_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] target
  let in_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] initial_theta
  let cost_function t = AD.Maths.(F 1. + (F 0. * t))
  let q_start = Mat.(eye 2 *$ (P.qs_coeff *. 0.5)) |> AD.pack_arr
  let q_start_coeff = P.qs_coeff

  (*AD.Maths.(
    F P.alpha
    * (F 1. + ((F P.beta - F 1.) * sigmoid ((t - F P.t_prep) / F 2E-3))))*)

  let power u g =
    AD.pack_flt
      (Mat.sum' (Mat.map (fun x -> Maths.pow x g) (AD.unpack_arr (AD.primal' u))))

  let cost ~u ~x ~k =
    let thetas = unpack_pos x in
    let vel = unpack_vel x in
    let _, x_state = unpack_full_state x 4 in
    let dx_p = AD.Maths.(tgt_pos - thetas) in
    let dx_vel = vel in
    let _dx_start = AD.Maths.(in_pos - thetas) in
    let t = AD.Maths.(__dt * F (float_of_int k)) in
    let torques = AD.Maths.(cmc *@ phi (P.__c *@ g (transpose x_state))) in
    let start = AD.Maths.(sigmoid ((F P.t_prep - t) / F 2E-4)) in
    (* let u = AD.Maths.(u*@(transpose __b)) in  *)
    AD.Maths.(
      ((sum' (u *@ r * u) * cost_function t)
      + (sum' (dx_p *@ q * dx_p) * sigmoid ((t - tau) / F 20E-3))
      + (F 0.1 * (sum' (dx_vel *@ q * dx_vel) * sigmoid ((t - tau) / F 20E-3)))
      + ((F 1. * (sum' (dx_vel *@ q_start * dx_vel) * start))
        + (sum' (_dx_start *@ q_start * _dx_start) * start)
        + (sum' (transpose torques *@ t_mat * transpose torques)
          * sigmoid ((F P.t_prep - t) / F 2E-3))))
      * __dt)

  let rl_u =
    (*Need to check this*)
    let rlu ~k ~x:_x ~u =
      let t = AD.Maths.(__dt * F (float_of_int k)) in
      let c = cost_function t in
      let r = Mat.(eye size_net *$ P.r_coeff) |> AD.pack_arr in
      let u = AD.Maths.(u *@ transpose P.b) in
      AD.Maths.(u *@ r *@ P.b * c * __dt)
    in
    Some rlu

  let rl_x =
    let rlx ~k ~x ~u:_u =
      let _, x_state = unpack_full_state x 4 in
      let thetas = unpack_pos x in
      let t = AD.Maths.(__dt * F (float_of_int k)) in
      let start = AD.Maths.(sigmoid ((F P.t_prep - t) / F 2E-4)) in
      let dx_p = AD.Maths.(tgt_pos - thetas) in
      let __c = AD.pack_arr P.c in
      let torques = AD.Maths.(cmc *@ phi (P.__c *@ g (transpose x_state))) in
      let dtorques = AD.Maths.(cmc *@ dphi (__c *@ transpose x_state) *@ __c) in
      let dx_vel = unpack_vel x in
      let dx_start = AD.Maths.(in_pos - thetas) in
      let r_xstate =
        AD.Maths.(
          transpose torques
          *@ t_mat
          *@ dtorques
          * F 2.
          * sigmoid ((F P.t_prep - t) / F 2E-3))
      in
      let r_xp1 = AD.Maths.(AD.F 2. * neg dx_p *@ q * sigmoid ((t - tau) / F 20E-3)) in
      let r_xp2 = AD.Maths.(AD.F 2. * neg dx_start *@ q_start * start) in
      let r_xv1 = AD.Maths.(dx_vel *@ q * sigmoid ((t - tau) / F 20E-3) * F 0.2) in
      let r_xv2 = AD.Maths.(dx_vel *@ q_start * start * F 2.) in
      AD.Maths.(concatenate ~axis:1 [| r_xp1 + r_xp2; r_xv1 + r_xv2; r_xstate |] * __dt)
    in
    Some rlx

  let rl_uu =
    let rluu ~k ~x:_x ~u:_u =
      let t = AD.Maths.(__dt * F (float_of_int k)) in
      let c = cost_function t in
      let ma = Mat.(eye size_net *$ P.r_coeff) |> AD.pack_arr in
      AD.Maths.(transpose P.b *@ ma *@ P.b * c * __dt)
    in
    Some rluu

  let rl_ux =
    let f ~k:_k ~x:_x ~u:_u = AD.F 0. in
    Some f

  let rl_xx =
    let rlxx ~k ~x:_x ~u:_u =
      let t = AD.Maths.(__dt * F (float_of_int k)) in
      let start = AD.Maths.(sigmoid ((F P.t_prep - t) / F 2E-4)) in
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
          ((sigmoid ((F P.t_prep - t) / F 2E-3)
           * (transpose P.__c *@ t_mat * F 2. *@ P.__c))
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

module C_Running (P : Prms) = struct
  (*default setting of qcoeff = 1E-7*)
  let n_theta = 4
  let q = Mat.(eye 2 *$ (q_coeff *. 0.5)) |> AD.pack_arr
  let r = Mat.(eye size_net *$ (P.r_coeff *. 0.5)) |> AD.pack_arr
  let t_mat = Mat.(eye 2 *$ (P.qs_coeff *. 0.5)) |> AD.pack_arr
  let a_mat = Mat.(eye size_net *$ (Defaults.a_coeff *. 0.5)) |> AD.pack_arr
  let tau = AD.F (P.t_prep +. P.t_mov)
  let __dt = AD.F sampling_dt
  let target = P.target_theta.(0) |> AD.pack_arr
  let tgt_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] target
  let in_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] initial_theta
  let cost_function t = AD.Maths.(F 1. + (F 0. * t))
  let q_start = Mat.(eye 2 *$ (P.qs_coeff *. 0.5)) |> AD.pack_arr
  let q_start_coeff = P.qs_coeff

  (*AD.Maths.(
    F P.alpha
    * (F 1. + ((F P.beta - F 1.) * sigmoid ((t - F P.t_prep) / F 2E-3))))*)

  let power u g =
    AD.pack_flt
      (Mat.sum' (Mat.map (fun x -> Maths.pow x g) (AD.unpack_arr (AD.primal' u))))

  let cost ~u ~x ~k =
    let wp, wm = P.weighing_pm in
    let cost_function t =
      AD.Maths.(
        (sigmoid ((F P.t_prep - t) / F 2E-4) * F wp)
        + (sigmoid ((t - F P.t_prep) / F 2E-4) * F wm))
    in
    let thetas = unpack_pos x in
    let vel = unpack_vel x in
    let _, x_state = unpack_full_state x 4 in
    let dx_p = AD.Maths.(tgt_pos - thetas) in
    let dx_vel = vel in
    let _dx_start = AD.Maths.(in_pos - thetas) in
    let t = AD.Maths.(__dt * F (float_of_int k)) in
    let torques = AD.Maths.(P.__c *@ transpose x_state) in
    let start = AD.Maths.(sigmoid ((F P.t_prep - t) / F 2E-4)) in
    let u = AD.Maths.(u *@ transpose P.b) in
    AD.Maths.(
      ((sum' (u *@ r * u) * cost_function t)
      + (sum' (dx_p *@ q * dx_p) * sigmoid ((t - tau) / F 20E-3))
      + (F 0.1 * (sum' (dx_vel *@ q * dx_vel) * sigmoid ((t - tau) / F 20E-3)))
      + ((F 1. * (sum' (dx_vel *@ q_start * dx_vel) * start))
        + (sum' (_dx_start *@ q_start * _dx_start) * start)
        + (sum' (transpose torques *@ t_mat * transpose torques)
          * sigmoid ((F P.t_prep - t) / F 2E-3))))
      * __dt)

  let rl_u =
    (*Need to check this*)
    let rlu ~k ~x:_x ~u =
      let t = AD.Maths.(__dt * F (float_of_int k)) in
      let c =
        let wp, wm = P.weighing_pm in
        AD.Maths.(
          (sigmoid ((F P.t_prep - t) / F 2E-4) * F wp)
          + (sigmoid ((t - F P.t_prep) / F 2E-4) * F wm))
      in
      let r = Mat.(eye size_net *$ P.r_coeff) |> AD.pack_arr in
      let u = AD.Maths.(u *@ transpose P.b) in
      AD.Maths.(u *@ r *@ P.b * c * __dt)
    in
    Some rlu

  let rl_x =
    let rlx ~k ~x ~u:_u =
      let _, x_state = unpack_full_state x 4 in
      let thetas = unpack_pos x in
      let t = AD.Maths.(__dt * F (float_of_int k)) in
      let start = AD.Maths.(sigmoid ((F P.t_prep - t) / F 2E-4)) in
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
          * sigmoid ((F P.t_prep - t) / F 2E-3))
      in
      let r_xp1 = AD.Maths.(AD.F 2. * neg dx_p *@ q * sigmoid ((t - tau) / F 20E-3)) in
      let r_xp2 = AD.Maths.(AD.F 2. * neg dx_start *@ q_start * start) in
      let r_xv1 = AD.Maths.(dx_vel *@ q * sigmoid ((t - tau) / F 20E-3) * F 0.2) in
      let r_xv2 = AD.Maths.(dx_vel *@ q_start * start * F 2.) in
      AD.Maths.(concatenate ~axis:1 [| r_xp1 + r_xp2; r_xv1 + r_xv2; r_xstate |] * __dt)
    in
    Some rlx

  let rl_uu =
    let rluu ~k ~x:_x ~u:_u =
      let t = AD.Maths.(__dt * F (float_of_int k)) in
      let c =
        let wp, wm = P.weighing_pm in
        AD.Maths.(
          (sigmoid ((F P.t_prep - t) / F 2E-4) * F wp)
          + (sigmoid ((t - F P.t_prep) / F 2E-4) * F wm))
      in
      let ma = Mat.(eye size_net *$ P.r_coeff) |> AD.pack_arr in
      AD.Maths.(P.b *@ ma *@ P.b * c * __dt)
    in
    Some rluu

  let rl_ux =
    let f ~k:_k ~x:_x ~u:_u = AD.F 0. in
    Some f

  let rl_xx =
    let rlxx ~k ~x:_x ~u:_u =
      let t = AD.Maths.(__dt * F (float_of_int k)) in
      let start = AD.Maths.(sigmoid ((F P.t_prep - t) / F 2E-4)) in
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
          ((sigmoid ((F P.t_prep - t) / F 2E-3)
           * (transpose P.__c *@ t_mat * F 2. *@ P.__c))
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

(*
module C_End (P : Prms) = struct
  let n_theta = 4
  let q = Mat.(eye 2 *$ Defaults.q_coeff *$ 0.5) |> AD.pack_arr
  let r = Mat.(eye P.m *$ (P.r_coeff *. 0.5)) |> AD.pack_arr
  let tau = AD.F (target_duration +. P.t_prep)
  let __dt = AD.F sampling_dt
  let target = P.target_theta.(0) |> AD.pack_arr
  let t_mat = Mat.(eye 2 *$ (P.t_coeff *. 0.5)) |> AD.pack_arr
  let a_mat = Mat.(eye 2 *$ (P.t_coeff *. 0.)) |> AD.pack_arr
  let tgt_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] target
  let tgt_vel = AD.Maths.get_slice [ []; [ 2; 3 ] ] target
  let q_start = AD.F q_start
  let in_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] initial_theta
  let cost_function t = AD.Maths.(F 0. * t)

  let cost ~u ~x ~k =
    let thetas = unpack_pos x in
    let dx_start = AD.Maths.(in_pos - thetas) in
    let _, x_state = unpack_full_state x n_theta in
    let t = AD.Maths.(__dt * F (float_of_int k)) in
    let torques = AD.Maths.(P.__c *@ transpose x_state) in
    AD.Maths.(
      (sum' (u *@ r * u)
      + (q_start * (sum' (dx_start *@ q * dx_start) * sigmoid ((F P.t_prep - t) / F 2E-4)))
      + (sum' (transpose torques *@ t_mat * transpose torques)
        * sigmoid ((F P.t_prep - t) / F 2E-3)))
      * __dt)


  let rl_u =
    let rlu ~k:_k ~x:_x ~u =
      let r = Mat.(eye P.m *$ P.r_coeff) |> AD.pack_arr in
      AD.Maths.(u *@ r * __dt)
    in
    Some rlu


  let rl_uu =
    let rluu ~k:_k ~x:_x ~u:_u =
      let ma = Mat.(eye P.m *$ P.r_coeff) |> AD.pack_arr in
      AD.Maths.(ma * __dt)
    in
    Some rluu


  let rl_ux =
    let f ~k:_k ~x:_x ~u:_u = AD.F 0. in
    Some f


  let rl_xx =
    let rlxx ~k ~x:_x ~u:_u =
      let t = AD.Maths.(F (float_of_int k) * __dt) in
      let mu = AD.Maths.(q_start * AD.Mat.eye 2 * sigmoid ((F P.t_prep - t) / F 2E-3)) in
      let mx =
        AD.Maths.(
          sigmoid ((F P.t_prep - t) / F 2E-3)

          * F 0.
          * (transpose P.__c *@ P.__c)
          * __dt)
      in
      let mf1 =
        AD.Maths.concatenate ~axis:1 [| mu; AD.Mat.zeros (n_theta - 2) (P.m + 2) |]
      in
      let mf2 =
        AD.Maths.concatenate
          ~axis:1
          [| AD.Mat.zeros (n_theta - 2) 4; AD.Mat.zeros (n_theta - 2) P.m |]
      in
      let mf3 = AD.Maths.concatenate ~axis:1 [| AD.Mat.zeros P.m n_theta; mx |] in
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
      AD.(Maths.(sum' (dx_p *@ q * dx_p) + (F 0.1 * sum' (dx_vel *@ q * dx_vel))))
    in
    fl


  let fl_x = None
  let fl_xx = None
end

module C_Successive (P : Prms) = struct
  let n_theta = 4
  let pause = 0.7
  let q = Mat.(eye 2 *$ Defaults.q_coeff *$ 0.5) |> AD.pack_arr
  let r = Mat.(eye P.m *$ (P.r_coeff *. 0.5)) |> AD.pack_arr
  let tau_1 = AD.F (target_duration +. P.t_prep)
  let tau_2 = AD.Maths.(tau_1 + F (pause +. target_duration))
  let __dt = AD.F sampling_dt
  let target_1 = P.target_theta.(0) |> AD.pack_arr
  let target_2 = P.target_theta.(1) |> AD.pack_arr
  let tgt_pos_1 = AD.Maths.get_slice [ []; [ 0; 1 ] ] target_1
  let tgt_pos_2 = AD.Maths.get_slice [ []; [ 0; 1 ] ] target_2

  let cost ~u ~x ~k =
    let t = AD.Maths.(F (float_of_int k) * __dt) in
    let thetas = unpack_pos x in
    let dx_p1 = AD.Maths.(tgt_pos_1 - thetas) in
    let dx_p2 = AD.Maths.(tgt_pos_2 - thetas) in
    let c1 =
      AD.Maths.(
        sum' (dx_p1 *@ q * dx_p1)
        * sigmoid ((t - tau_1) / F 50E-3)
        * sigmoid ((tau_1 + F pause - t) / F 50E-3))
    and c2 = AD.Maths.(sum' (dx_p2 *@ q * dx_p2) * sigmoid ((t - tau_2) / F 50E-3)) in
    AD.Maths.((sum' (u *@ r * u) + c1 + c2) * __dt)


  let rl_u =
    let rlu ~k:_k ~x:_x ~u =
      let r = Mat.(eye P.m *$ P.r_coeff) |> AD.pack_arr in
      AD.Maths.(u *@ r * __dt)
    in
    Some rlu


  let rl_uu =
    let rluu ~k:_k ~x:_x ~u:_u =
      let ma = Mat.(eye P.m *$ P.r_coeff) |> AD.pack_arr in
      AD.Maths.(ma * __dt)
    in
    Some rluu


  let rl_ux =
    let f ~k:_k ~x:_x ~u:_u = AD.F 0. in
    Some f


  let rl_xx =
    let rlxx ~k ~x:_x ~u:_u =
      let t = AD.Maths.(F (float_of_int k) * __dt) in
      let mu =
        AD.Maths.(
          AD.Mat.(eye 2)
          * ((sigmoid ((t - tau_1) / F 50E-3) * sigmoid ((tau_1 + F pause - t) / F 50E-3))
            + sigmoid ((t - tau_2) / F 50E-3)))
      in
      let c = AD.Maths.(F Defaults.q_coeff * __dt) in
      let mf =
        AD.Maths.concatenate
          ~axis:1
          [| AD.Maths.(mu * c); AD.Mat.zeros (n_theta - 2) (P.m + 2) |]
      in
      AD.Maths.concatenate ~axis:0 [| mf; AD.Mat.zeros (P.m + 2) P.n |]
    in
    Some rlxx


  let final_cost ~x:_x ~k:_k = AD.F 0.

  let fl_x =
    let flx ~k:_k ~x:_x = AD.F 0. in
    Some flx


  let fl_xx =
    let flx ~k:_k ~x:_x = AD.F 0. in
    Some flx
end
*)
module C_Uncertain (P : Prms) = struct
  (*default setting of qcoeff = 1E-7*)
  let n_theta = 4
  let q = Mat.(eye 2 *$ (Defaults.q_coeff *. 0.5)) |> AD.pack_arr
  let r = Mat.(eye P.m *$ (P.r_coeff *. 0.5)) |> AD.pack_arr
  let tau = AD.F P.t_prep
  let __dt = AD.F sampling_dt
  let target = P.target_theta.(0) |> AD.pack_arr
  let tgt_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] target
  let in_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] initial_theta

  let xstar =
    Mat.(
      (get_slice [ [ 700 ]; [ n_theta; -1 ] ])
        (load_txt "results/inputs_beg/reach_1/traj_700"))
    |> AD.pack_arr

  let __w = AD.pack_arr P.w
  let __a = AD.Maths.((__w - AD.Mat.eye size_net) / tau)
  let __at = AD.Maths.transpose __a
  let __atinv = AD.Maths.inv __at
  let __ainv = AD.Maths.inv __a
  let idt = AD.Mat.eye n
  let a = AD.unpack_arr __a

  let cost ~u ~x ~k =
    let _, x_ac = unpack_full_state x n_theta in
    let t = AD.Maths.(__dt * F (float_of_int k)) in
    let dx = AD.Maths.(x_ac - xstar) in
    let a_of_t = Mat.(a *$ ((sampling_dt *. float k) -. P.t_prep) /$ l2norm' a) in
    let at_of_t =
      Mat.(transpose a *$ ((sampling_dt *. float k) -. P.t_prep) /$ l2norm' a)
    in
    let _expa =
      if t < AD.F P.t_prep then Linalg.D.expm a_of_t |> AD.pack_arr else AD.Mat.zeros n n
    in
    let _expat =
      if t < AD.F P.t_prep then Linalg.D.expm at_of_t |> AD.pack_arr else AD.Mat.zeros n n
    in
    let _ut = AD.Maths.transpose u in
    let _dxt = AD.Maths.transpose dx in
    let _term1 = AD.Maths.(u *@ (idt - _expa) *@ __ainv *@ _dxt) |> AD.Maths.sum' in
    let _term2 = AD.Maths.(__atinv *@ (idt - _expat) *@ _dxt) |> AD.Maths.sum' in
    let _term3 = AD.Maths.(sum' (u *@ r * u)) in
    AD.Maths.((sqr _term1 + sqr _term2 + _term3) * __dt)

  let rl_u = None

  let rl_uu =
    let rluu ~k:_k ~x:_x ~u:_u =
      let ma = Mat.(eye P.m *$ P.r_coeff) |> AD.pack_arr in
      AD.Maths.(ma * __dt)
    in
    Some rluu

  let rl_ux = None

  let rl_xx =
    let rlxx ~k:_k ~x:_x ~u:_u = AD.F 0. in
    Some rlxx

  let final_cost ~x ~k:_k =
    let q = Owl.Mat.(eye 2 *$ Defaults.q_coeff) |> AD.pack_arr in
    let fl =
      let thetas = unpack_pos x in
      let thetas_dot = unpack_vel x in
      let dx_p = AD.Maths.(tgt_pos - thetas)
      and dx_vel = thetas_dot in
      AD.(Maths.(sum' (dx_p *@ q * dx_p) + sum' (dx_vel *@ q * dx_vel)))
    in
    AD.Maths.(F 1. * fl)

  let fl_x = None
  let fl_xx = None
end
(*
module C_Cross (P : Prms) = struct
  (*default setting of qcoeff = 1E-7*)
  let n_theta = 4
  let q = Mat.(eye 2 *$ (Defaults.q_coeff *. 0.5)) |> AD.pack_arr
  let r = Mat.(eye P.m *$ (P.r_coeff *. 0.5)) |> AD.pack_arr
  let t_mat = Mat.(eye 2 *$ (P.t_coeff *. 0.5)) |> AD.pack_arr
  let a_mat = Mat.(eye P.m *$ (Defaults.a_coeff *. 0.5)) |> AD.pack_arr
  let tau = AD.F (P.t_prep +. target_duration)
  let __dt = AD.F sampling_dt
  let target = P.target_theta.(0) |> AD.pack_arr
  let tgt_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] target
  let in_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] initial_theta
  let q_start = AD.F q_start
  let gamma = AD.F P.gamma_exponent
  let cost_function t = AD.Maths.(F 1. + (F 0. * t))

  let cost ~u ~x ~k =
    let thetas = unpack_pos x in
    let vel = unpack_vel x in
    let _, x_state = unpack_full_state x 4 in
    let dx_p = AD.Maths.(tgt_pos - thetas) in
    let dx_vel = vel in
    let dx_start = AD.Maths.(in_pos - thetas) in
    let t = AD.Maths.(__dt * F (float_of_int k)) in
    let torques = AD.Maths.(P.__c *@ transpose x_state) in
    AD.Maths.(
      ((sum' (abs u *@ AD.Mat.ones 200 200 * abs u) * F (0.5 *. r_coeff))
      (*((F r_coeff * AD.Maths.(sum' (pow (abs u) gamma)))*)
      + (sum' (dx_p *@ q * dx_p) * sigmoid ((t - tau) / F 20E-3))
      + (F 0.1 * (sum' (dx_vel *@ q * dx_vel) * sigmoid ((t - tau) / F 20E-3)))
      + (q_start * (sum' (dx_start *@ q * dx_start) * sigmoid ((F P.t_prep - t) / F 2E-4)))
      + (sum' (transpose torques *@ t_mat * transpose torques)
        * sigmoid ((F P.t_prep - t) / F 2E-3)))
      * __dt)


  let rl_u =
    let _rlu ~k:_k ~x:_x ~u =
      AD.Maths.(F r_coeff * __dt * AD.Maths.(u *@ AD.Mat.ones 200 200))
    in
    None


  let rl_uu =
    let _rluu ~k:_k ~x:_x ~u:_u = AD.Maths.(F r_coeff * AD.Mat.ones 200 200 * __dt) in
    None


  let rl_ux =
    let f ~k:_k ~x:_x ~u:_u = AD.F 0. in
    Some f


  let rl_xx =
    let rlxx ~k ~x:_x ~u:_u =
      let t = AD.Maths.(__dt * F (float_of_int k)) in
      let mu =
        AD.Maths.(
          (AD.Mat.eye 2 * AD.F Defaults.q_coeff * __dt * sigmoid ((t - tau) / F 20E-3))
          + (q_start * AD.Mat.eye 2 * sigmoid ((F P.t_prep - t) / F 2E-4) * __dt))
      in
      let mv =
        AD.Maths.(
          __dt * AD.Mat.eye 2 * F Defaults.q_coeff * (F 0.1 * sigmoid ((t - tau) / F 20E-3)))
      in
      let mx =
        AD.Maths.(
          ((sigmoid ((F P.t_prep - t) / F 2E-3)
           * F P.t_coeff
           * (transpose P.__c *@ P.__c))
          + (F Defaults.a_coeff * AD.Mat.eye P.m))
          * __dt)
      in
      let mf1 =
        AD.Maths.concatenate ~axis:1 [| mu; AD.Mat.zeros (n_theta - 2) (P.m + 2) |]
      in
      let mf2 =
        AD.Maths.concatenate
          ~axis:1
          [| AD.Mat.zeros (n_theta - 2) 2; mv; AD.Mat.zeros (n_theta - 2) P.m |]
      in
      let mf3 = AD.Maths.concatenate ~axis:1 [| AD.Mat.zeros P.m n_theta; mx |] in
      AD.Maths.concatenate ~axis:0 [| mf1; mf2; mf3 |]
    in
    Some rlxx


  let final_cost ~x:_x ~k:_k = AD.F 0.

  let fl_x =
    let f ~k:_k ~x:_x = AD.F 0. in
    Some f


  (*let f ~k:_k ~x:_x = AD.F 0. in
    Some f*)

  let fl_xx =
    let f ~k:_k ~x:_x = AD.F 0. in
    Some f
end

module C_Gamma (P : Prms) = struct
  (*default setting of qcoeff = 1E-7*)
  let n_theta = 4
  let q = Mat.(eye 2 *$ (Defaults.q_coeff *. 0.5)) |> AD.pack_arr
  let r = Mat.(eye P.m *$ (P.r_coeff *. 0.5)) |> AD.pack_arr
  let t_mat = Mat.(eye 2 *$ (P.t_coeff *. 0.5)) |> AD.pack_arr
  let a_mat = Mat.(eye P.m *$ (Defaults.a_coeff *. 0.5)) |> AD.pack_arr
  let tau = AD.F (P.t_prep +. target_duration)
  let __dt = AD.F sampling_dt
  let target = P.target_theta.(0) |> AD.pack_arr
  let tgt_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] target
  let in_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] initial_theta
  let q_start = AD.F q_start
  let gamma = AD.F P.gamma_exponent
  let cost_function t = AD.Maths.(F 0. * t)

  let cost ~u ~x ~k =
    let thetas = unpack_pos x in
    let vel = unpack_vel x in
    let _, x_state = unpack_full_state x 4 in
    let dx_p = AD.Maths.(tgt_pos - thetas) in
    let dx_vel = vel in
    let dx_start = AD.Maths.(in_pos - thetas) in
    let t = AD.Maths.(__dt * F (float_of_int k)) in
    let torques = AD.Maths.(P.__c *@ transpose x_state) in
    AD.Maths.(
      ((sum' (pow (abs u) gamma) * F r_coeff)
      (*((F r_coeff * AD.Maths.(sum' (pow (abs u) gamma)))*)
      + (sum' (dx_p *@ q * dx_p) * sigmoid ((t - tau) / F 20E-3))
      + (F 0.1 * (sum' (dx_vel *@ q * dx_vel) * sigmoid ((t - tau) / F 20E-3)))
      + (q_start * (sum' (dx_start *@ q * dx_start) * sigmoid ((F P.t_prep - t) / F 2E-4)))
      + (sum' (transpose torques *@ t_mat * transpose torques)
        * sigmoid ((F P.t_prep - t) / F 2E-3)))
      * __dt)


  let rl_u =
    let _rlu ~k:_k ~x:_x ~u =
      if gamma = AD.F 1.
      then AD.Maths.(F r_coeff * gamma * sum' (signum u))
      else (
        let u_u = AD.unpack_arr u in
        let mat =
          Mat.map
            (fun x ->
              if x < 0.
              then Maths.(neg (pow (neg x) (P.gamma_exponent -. 1.)))
              else Maths.(pow x (P.gamma_exponent -. 1.)))
            u_u
          |> AD.pack_arr
        in
        AD.Maths.(F r_coeff * gamma * AD.Maths.(sum' mat)))
    in
    None


  (*Need to check this




    Some rlu*)

  let rl_uu =
    let _rluu ~k:_k ~x:_x ~u =
      if gamma = AD.F 1.
      then AD.F 0.
      else
        AD.Maths.(
          F r_coeff
          * gamma
          * (gamma - F 1.)
          * AD.Maths.(sum' (pow (abs u) (gamma - F 2.))))
    in
    None


  (*
  let rluu ~k:_k ~x:_x ~u:_u =
   
      AD.Maths.(F r_coeff * (gamma - F 1.) * AD.Maths.(sum' (pow u (gamma - F 1.))))
    in
    Some rluu*)

  let rl_ux =
    let f ~k:_k ~x:_x ~u:_u = AD.F 0. in
    Some f


  let rl_xx =
    let rlxx ~k ~x:_x ~u:_u =
      let t = AD.Maths.(__dt * F (float_of_int k)) in
      let mu =
        AD.Maths.(
          (AD.Mat.eye 2 * AD.F Defaults.q_coeff * __dt * sigmoid ((t - tau) / F 20E-3))
          + (q_start * AD.Mat.eye 2 * sigmoid ((F P.t_prep - t) / F 2E-4) * __dt))
      in
      let mv =
        AD.Maths.(
          __dt * AD.Mat.eye 2 * F Defaults.q_coeff * (F 0.1 * sigmoid ((t - tau) / F 20E-3)))
      in
      let mx =
        AD.Maths.(
          ((sigmoid ((F P.t_prep - t) / F 2E-3)
           * F P.t_coeff
           * (transpose P.__c *@ P.__c))
          + (F Defaults.a_coeff * AD.Mat.eye P.m))
          * __dt)
      in
      let mf1 =
        AD.Maths.concatenate ~axis:1 [| mu; AD.Mat.zeros (n_theta - 2) (P.m + 2) |]
      in
      let mf2 =
        AD.Maths.concatenate
          ~axis:1
          [| AD.Mat.zeros (n_theta - 2) 2; mv; AD.Mat.zeros (n_theta - 2) P.m |]
      in
      let mf3 = AD.Maths.concatenate ~axis:1 [| AD.Mat.zeros P.m n_theta; mx |] in
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
      AD.(Maths.((F 0. * sum' (dx_p *@ q * dx_p)) + (F 0. * sum' (dx_vel *@ q * dx_vel))))
    in
    fl


  let fl_x = None

  (*let f ~k:_k ~x:_x = AD.F 0. in
    Some f*)

  let fl_xx = None
end

module C_Var (P : Prms) = struct
  (*default setting of qcoeff = 1E-7*)
  let n_theta = 4
  let q = Mat.(eye 2 *$ (Defaults.q_coeff *. 0.5)) |> AD.pack_arr
  let r = Mat.(eye P.m *$ (P.r_coeff *. 0.5)) |> AD.pack_arr
  let t_mat = Mat.(eye 2 *$ (P.t_coeff *. 0.5)) |> AD.pack_arr
  let a_mat = Mat.(eye P.m *$ (Defaults.a_coeff *. 0.5)) |> AD.pack_arr
  let tau = AD.F (P.t_prep +. target_duration)
  let __dt = AD.F sampling_dt
  let target = P.target_theta.(0) |> AD.pack_arr
  let tgt_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] target
  let in_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] initial_theta
  let cost_function t = AD.Maths.(F 1. + (F 0. * t))
  let q_start = AD.F q_start

  let power u g =
    AD.pack_flt
      (Mat.sum' (Mat.map (fun x -> Maths.pow x g) (AD.unpack_arr (AD.primal' u))))


  let cost ~u ~x ~k =
    let thetas = unpack_pos x in
    let vel = unpack_vel x in
    let _, x_state = unpack_full_state x 4 in
    let dx_p = AD.Maths.(tgt_pos - thetas) in
    let dx_vel = vel in
    let dx_start = AD.Maths.(in_pos - thetas) in
    let t = AD.Maths.(__dt * F (float_of_int k)) in
    let torques = AD.Maths.(P.__c *@ transpose x_state) in
    let var = AD.Maths.(sum' (sqr (x_state - mean x_state))) in
    AD.Maths.(
      ((sum' (u *@ r * u) / (F 1. + var))
      + sum' (x_state *@ a_mat * x_state)
      + (sum' (dx_p *@ q * dx_p) * sigmoid ((t - tau) / F 20E-3))
      + (F 0.1 * (sum' (dx_vel *@ q * dx_vel) * sigmoid ((t - tau) / F 20E-3)))
      + (q_start * (sum' (dx_start *@ q * dx_start) * sigmoid ((F P.t_prep - t) / F 2E-4)))
      + (sum' (transpose torques *@ t_mat * transpose torques)
        * sigmoid ((F P.t_prep - t) / F 2E-3)))
      * __dt)


  let rl_u = None
  let rl_uu = None
  let rl_ux = None
  let rl_xx = None

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

module C_Dif (P : Prms) = struct
  (*default setting of qcoeff = 1E-7*)
  let n_theta = 4
  let q = Mat.(eye 2 *$ (Defaults.q_coeff *. 0.5)) |> AD.pack_arr
  let r = Mat.(eye P.m *$ (P.r_coeff *. 0.5)) |> AD.pack_arr
  let t_mat = Mat.(eye 2 *$ (P.t_coeff *. 0.5)) |> AD.pack_arr
  let a_mat = Mat.(eye P.m *$ (Defaults.a_coeff *. 0.5)) |> AD.pack_arr
  let tau = AD.F (P.t_prep +. target_duration)
  let __dt = AD.F sampling_dt
  let target = P.target_theta.(0) |> AD.pack_arr
  let tgt_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] target
  let in_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] initial_theta
  let q_start = AD.F q_start
  let gamma = AD.F P.gamma_exponent
  let cost_function t = AD.Maths.(F 0. * t)

  let cost ~u ~x ~k =
    let thetas = unpack_pos x in
    let vel = unpack_vel x in
    let _, x_state = unpack_full_state x 4 in
    let x_s = AD.Maths.get_slice [ []; [ 0; n - 1 ] ] x_state in
    let dx_p = AD.Maths.(tgt_pos - thetas) in
    let dx_vel = vel in
    let dx_start = AD.Maths.(in_pos - thetas) in
    let t = AD.Maths.(__dt * F (float_of_int k)) in
    let torques = AD.Maths.(P.__c *@ transpose x_s) in
    AD.Maths.(
      ((F r_coeff * AD.Maths.(sum' (sqr u)))
      + (sum' (dx_p *@ q * dx_p) * sigmoid ((t - tau) / F 20E-3))
      + (F 0.1 * (sum' (dx_vel *@ q * dx_vel) * sigmoid ((t - tau) / F 20E-3)))
      + (q_start * (sum' (dx_start *@ q * dx_start) * sigmoid ((F P.t_prep - t) / F 2E-4)))
      + (sum' (transpose torques *@ t_mat * transpose torques)
        * sigmoid ((F P.t_prep - t) / F 2E-3)))
      * __dt)


  let rl_u = None

  (*Need to check this
      let rlu ~k:_k ~x:_x ~u =
        AD.Maths.(F r_coeff * AD.Maths.(sum' (pow u (gamma - F 1.))))
      in



    Some rlu*)

  let rl_uu = None

  (*
  let rluu ~k:_k ~x:_x ~u:_u =
   
      AD.Maths.(F r_coeff * (gamma - F 1.) * AD.Maths.(sum' (pow u (gamma - F 1.))))
    in
    Some rluu*)

  let rl_ux =
    let f ~k:_k ~x:_x ~u:_u = AD.F 0. in
    Some f


  let rl_xx = None

  let final_cost ~x ~k:_k =
    let q = Owl.Mat.(eye 2 *$ Defaults.q_coeff) |> AD.pack_arr in
    let fl =
      let thetas = unpack_pos x in
      let thetas_dot = unpack_vel x in
      let dx_p = AD.Maths.(tgt_pos - thetas)
      and dx_vel = thetas_dot in
      AD.(Maths.((F 0. * sum' (dx_p *@ q * dx_p)) + (F 0. * sum' (dx_vel *@ q * dx_vel))))
    in
    fl


  let fl_x = None

  (*let f ~k:_k ~x:_x = AD.F 0. in
    Some f*)

  let fl_xx = None
end

module C_Ddot (P : Prms) = struct
  (*default setting of qcoeff = 1E-7*)
  let n_theta = 4
  let q = Mat.(eye 2 *$ (Defaults.q_coeff *. 0.5)) |> AD.pack_arr
  let r = Mat.(eye P.m *$ (P.r_coeff *. 0.5)) |> AD.pack_arr
  let t_mat = Mat.(eye 2 *$ (P.t_coeff *. 0.5)) |> AD.pack_arr
  let a_mat = Mat.(eye P.m *$ (Defaults.a_coeff *. 0.5)) |> AD.pack_arr
  let tau = AD.F (P.t_prep +. target_duration)
  let __dt = AD.F sampling_dt
  let target = P.target_theta.(0) |> AD.pack_arr
  let tgt_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] target
  let in_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] initial_theta
  let q_start = AD.F q_start
  let gamma = AD.F P.gamma_exponent
  let cost_function t = AD.Maths.(F 0. * t)

  let cost ~u ~x ~k =
    let thetas = unpack_pos x in
    let vel = unpack_vel x in
    let _, x_state = unpack_full_state x 4 in
    let x_s, u_1, u_2 =
      ( AD.Maths.get_slice [ []; [ 0; n - 1 ] ] x_state
      , AD.Maths.get_slice [ []; [ n; n + P.m - 1 ] ] x_state
      , AD.Maths.get_slice [ []; [ n + P.m; -1 ] ] x_state )
    in
    let dx_p = AD.Maths.(tgt_pos - thetas) in
    let dx_vel = vel in
    let dx_start = AD.Maths.(in_pos - thetas) in
    let t = AD.Maths.(__dt * F (float_of_int k)) in
    let torques = AD.Maths.(P.__c *@ transpose x_s) in
    let du_snd = AD.Maths.((u + u_2 - (F 2. * u_1)) / (__dt * __dt)) in
    AD.Maths.(
      (sum' (du_snd *@ r * du_snd)
      + (sum' (dx_p *@ q * dx_p) * sigmoid ((t - tau) / F 20E-3))
      + (F 0.1 * (sum' (dx_vel *@ q * dx_vel) * sigmoid ((t - tau) / F 20E-3)))
      + (q_start * (sum' (dx_start *@ q * dx_start) * sigmoid ((F P.t_prep - t) / F 2E-4)))
      + (sum' (transpose torques *@ t_mat * transpose torques)
        * sigmoid ((F P.t_prep - t) / F 2E-3)))
      * __dt)


  let rl_u = None
  let rl_uu = None

  let rl_ux =
    let f ~k:_k ~x:_x ~u:_u = AD.F 0. in
    Some f


  let rl_xx = None

  let final_cost ~x ~k:_k =
    let q = Owl.Mat.(eye 2 *$ Defaults.q_coeff) |> AD.pack_arr in
    let fl =
      let thetas = unpack_pos x in
      let thetas_dot = unpack_vel x in
      let dx_p = AD.Maths.(tgt_pos - thetas)
      and dx_vel = thetas_dot in
      AD.(Maths.((F 0. * sum' (dx_p *@ q * dx_p)) + (F 0. * sum' (dx_vel *@ q * dx_vel))))
    in
    fl


  let fl_x = None

  (*let f ~k:_k ~x:_x = AD.F 0. in
    Some f*)

  let fl_xx = None
end
(* 
module C_Reduced (P : Prms) = struct
  (*default setting of qcoeff = 1E-7*)
  let n_theta = 4
  let q = Mat.(eye 2 *$ (Defaults.q_coeff *. 0.5)) |> AD.pack_arr
  let r = Mat.(eye P.m *$ (P.r_coeff *. 0.5)) |> AD.pack_arr
  let t_mat = Mat.(eye 2 *$ (Defaults.t_coeff *. 0.5)) |> AD.pack_arr
  let a_mat = Mat.(eye P.m *$ (Defaults.a_coeff *. 0.5)) |> AD.pack_arr
  let tau = AD.F (P.t_prep +. target_duration)
  let __dt = AD.F sampling_dt
  let target = P.target_theta.(0) |> AD.pack_arr
  let tgt_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] target
  let in_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] initial_theta
  let q_start = AD.F q_start
  let gamma = AD.F P.gamma_exponent
  let cost_function t = AD.Maths.(F 0. * t)

  let cost ~u ~x ~k =
    let thetas = unpack_pos x in
    let vel = unpack_vel x in
    let _, x_state = unpack_full_state x 4 in
    let dx_p = AD.Maths.(tgt_pos - thetas) in
    let dx_vel = vel in
    let dx_start = AD.Maths.(in_pos - thetas) in
    let t = AD.Maths.(__dt * F (float_of_int k)) in
    let torques = AD.Maths.(P.__c *@ transpose x_state) in
    let u_tot = AD.Maths.(transpose (top_pc *@ transpose u)) in
    AD.Maths.(
      ((sum' (sqr u_tot) * F (r_coeff /. 2.))
      + sum' (x_state *@ a_mat * x_state)
      + (sum' (dx_p *@ q * dx_p) * sigmoid ((t - tau) / F 20E-3))
      + (F 0.1 * (sum' (dx_vel *@ q * dx_vel) * sigmoid ((t - tau) / F 20E-3)))
      + (q_start * (sum' (dx_start *@ q * dx_start) * sigmoid ((F P.t_prep - t) / F 2E-4)))
      + (sum' (transpose torques *@ t_mat * transpose torques)
        * sigmoid ((F P.t_prep - t) / F 2E-3)))
      * __dt)


  let rl_u =
    (*Need to check this*)
    let _rlu ~k:_k ~x:_x ~u =
      AD.Maths.(
        transpose (transpose top_pc *@ top_pc *@ transpose u) * F P.r_coeff * __dt)
    in
    Some _rlu


  let rl_uu =
    let _rluu ~k:_k ~x:_x ~u:_u =
      AD.Maths.(transpose top_pc *@ top_pc * F P.r_coeff * __dt)
    in
    Some _rluu


  let rl_ux =
    let f ~k:_k ~x:_x ~u:_u = AD.F 0. in
    Some f


  let rl_xx =
    let rlxx ~k ~x:_x ~u:_u =
      let t = AD.Maths.(__dt * F (float_of_int k)) in
      let mu =
        AD.Maths.(
          (AD.Mat.eye 2 * AD.F Defaults.q_coeff * __dt * sigmoid ((t - tau) / F 20E-3))
          + (q_start * AD.Mat.eye 2 * sigmoid ((F P.t_prep - t) / F 2E-4) * __dt))
      in
      let mv =
        AD.Maths.(
          __dt * AD.Mat.eye 2 * F Defaults.q_coeff * (F 0.1 * sigmoid ((t - tau) / F 20E-3)))
      in
      let mx =
        AD.Maths.(
          ((sigmoid ((F P.t_prep - t) / F 2E-3)
           * F Defaults.t_coeff
           * (transpose __c *@ __c))
          + (F Defaults.a_coeff * AD.Mat.eye P.m))
          * __dt)
      in
      let mf1 =
        AD.Maths.concatenate ~axis:1 [| mu; AD.Mat.zeros (n_theta - 2) (P.m + 2) |]
      in
      let mf2 =
        AD.Maths.concatenate
          ~axis:1
          [| AD.Mat.zeros (n_theta - 2) 2; mv; AD.Mat.zeros (n_theta - 2) P.m |]
      in
      let mf3 = AD.Maths.concatenate ~axis:1 [| AD.Mat.zeros P.m n_theta; mx |] in
      AD.Maths.concatenate ~axis:0 [| mf1; mf2; mf3 |]
    in
    Some rlxx


  let final_cost ~x:_k ~k:_k = AD.F 0.
  let fl_x = None

  (*let f ~k:_k ~x:_x = AD.F 0. in
    Some f*)

  let fl_xx = None
end*)

module C_Sparse (P : Prms) = struct
  (*default setting of qcoeff = 1E-7*)
  let n_theta = 4
  let q = Mat.(eye 2 *$ (Defaults.q_coeff *. 0.5)) |> AD.pack_arr
  let r = Mat.(eye P.m *$ (P.r_coeff *. 0.5)) |> AD.pack_arr
  let t_mat = Mat.(eye 2 *$ (Defaults.t_coeff *. 0.5)) |> AD.pack_arr
  let a_mat = Mat.(eye P.m *$ (Defaults.a_coeff *. 0.5)) |> AD.pack_arr
  let tau = AD.F (P.t_prep +. P.t_mov)
  let __dt = AD.F sampling_dt
  let target = P.target_theta.(0) |> AD.pack_arr
  let tgt_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] target
  let in_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] initial_theta
  let cost_function t = AD.Maths.(F 1. + (F 0. * t))
  let q_start = AD.F q_start
  let lambda = AD.F 0.2

  (*AD.Maths.(
    F P.alpha
    * (F 1. + ((F P.beta - F 1.) * sigmoid ((t - F P.t_prep) / F 2E-3))))*)

  let power u g =
    AD.pack_flt
      (Mat.sum' (Mat.map (fun x -> Maths.pow x g) (AD.unpack_arr (AD.primal' u))))


  let cost ~u ~x ~k =
    let thetas = unpack_pos x in
    let vel = unpack_vel x in
    let _, x_state = unpack_full_state x 4 in
    let dx_p = AD.Maths.(tgt_pos - thetas) in
    let dx_vel = vel in
    let dx_start = AD.Maths.(in_pos - thetas) in
    let t = AD.Maths.(__dt * F (float_of_int k)) in
    let torques = AD.Maths.(__c *@ transpose x_state) in
    AD.Maths.(
      ((cost_function t * F r_coeff * sum' (log (F 1. + (sqr u / sqr lambda))))
      + sum' (x_state *@ a_mat * x_state)
      + (sum' (dx_p *@ q * dx_p) * sigmoid ((t - tau) / F 20E-3))
      + (F 0.1 * (sum' (dx_vel *@ q * dx_vel) * sigmoid ((t - tau) / F 20E-3)))
      + ( (sum' (dx_start *@ q * dx_start) * sigmoid ((F P.t_prep - t) / F 2E-4)))
      + (sum' (transpose torques *@ t_mat * transpose torques)
        * sigmoid ((F P.t_prep - t) / F 2E-3)))
      * __dt)


  let rl_u = None
  let rl_uu = None

  let rl_ux =
    let f ~k:_k ~x:_x ~u:_u = AD.F 0. in
    Some f


  let rl_xx =
    let rlxx ~k ~x:_x ~u:_u =
      let t = AD.Maths.(__dt * F (float_of_int k)) in
      let mu =
        AD.Maths.(
          (AD.Mat.eye 2 * AD.F Defaults.q_coeff * __dt * sigmoid ((t - tau) / F 20E-3))

          + ( AD.F Defaults.q_coeff  * AD.Mat.eye 2 * sigmoid ((F P.t_prep - t) / F 2E-4) * __dt))
      in
      let mv =
        AD.Maths.(
          __dt * AD.Mat.eye 2 * F Defaults.q_coeff * (F 0.1 * sigmoid ((t - tau) / F 20E-3)))
      in
      let mx =
        AD.Maths.(
          ((sigmoid ((F P.t_prep - t) / F 2E-3)
           * F Defaults.t_coeff
           * (transpose __c *@ __c))
          + (F Defaults.a_coeff * AD.Mat.eye P.m))
          * __dt)
      in
      let mf1 =
        AD.Maths.concatenate ~axis:1 [| mu; AD.Mat.zeros (n_theta - 2) (P.m + 2) |]
      in
      let mf2 =
        AD.Maths.concatenate
          ~axis:1
          [| AD.Mat.zeros (n_theta - 2) 2; mv; AD.Mat.zeros (n_theta - 2) P.m |]
      in
      let mf3 = AD.Maths.concatenate ~axis:1 [| AD.Mat.zeros P.m n_theta; mx |] in
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
end *)

module C_Sequential (P : Prms) = struct
  (*default setting of qcoeff = 1E-7*)
  let n_theta = 4
  let q = Mat.(eye 2 *$ (q_coeff *. 0.5)) |> AD.pack_arr
  let r = Mat.(eye P.m *$ (P.r_coeff *. 0.5)) |> AD.pack_arr
  let t_mat = Mat.(eye 2 *$ (P.qs_coeff *. 0.5)) |> AD.pack_arr
  let a_mat = Mat.(eye P.m *$ (Defaults.a_coeff *. 0.5)) |> AD.pack_arr
  let tau1 = AD.F (P.t_prep +. P.t_mov)
  let pause = Defaults.pause
  let tau2 = AD.F (P.t_prep +. P.t_mov +. pause)
  let tau3 = AD.F (P.t_prep +. (2. *. P.t_mov) +. pause)
  let __dt = AD.F sampling_dt
  let target1 = P.target_theta.(0) |> AD.pack_arr
  let target2 = P.target_theta.(1) |> AD.pack_arr
  let tgt_pos1 = AD.Maths.get_slice [ []; [ 0; 1 ] ] target1
  let tgt_pos = tgt_pos1
  let tau = tau1
  let target = target1
  let tgt_pos2 = AD.Maths.get_slice [ []; [ 0; 1 ] ] target2
  let in_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] initial_theta
  let cost_function t = AD.Maths.(F 1. + (F 0. * t))
  let q_start = Mat.(eye 2 *$ (P.qs_coeff *. 0.5)) |> AD.pack_arr
  let q_start_coeff = P.qs_coeff

  let power u g =
    AD.pack_flt
      (Mat.sum' (Mat.map (fun x -> Maths.pow x g) (AD.unpack_arr (AD.primal' u))))

  let cost ~u ~x ~k =
    let thetas = unpack_pos x in
    let vel = unpack_vel x in
    let _, x_state = unpack_full_state x 4 in
    let dx_p1 = AD.Maths.(tgt_pos1 - thetas) in
    let dx_p2 = AD.Maths.(tgt_pos2 - thetas) in
    let dx_vel = AD.Maths.(F 1. * vel) in
    let _dx_start = AD.Maths.(in_pos - thetas) in
    let t = AD.Maths.(__dt * F (float_of_int k)) in
    let torques = AD.Maths.(P.__c *@ transpose x_state) in
    let start = AD.Maths.(F 1. * sigmoid ((F P.t_prep - t) / F 2E-4)) in
    let pause_s =
      AD.Maths.(
        F 0.1 / F pause * sigmoid ((t - tau1) / F 2E-4) * sigmoid ((tau2 - t) / F 2E-4))
    in
    let pause_f = AD.Maths.(F 1. * sigmoid ((t - tau3) / F 2E-4)) in
    AD.Maths.(
      ((sum' (u *@ r * u) * cost_function t)
      + ((sum' (dx_p1 *@ q * dx_p1) * F 1. * pause_s)
        + (sum' (dx_vel *@ q * dx_vel) * pause_s * F 1.))
      + (sum' (dx_p2 *@ q * dx_p2) * pause_f)
      + (sum' (dx_vel *@ q * dx_vel) * pause_f)
      + ((F 1. * (sum' (dx_vel *@ q_start * dx_vel) * start))
        + (sum' (_dx_start *@ q_start * _dx_start) * start)
        + (sum' (transpose torques *@ t_mat * transpose torques)
          * sigmoid ((F P.t_prep - t) / F 2E-3))))
      * __dt)

  let rl_u =
    (*Need to check this*)
    let rlu ~k ~x:_x ~u =
      let t = AD.Maths.(__dt * F (float_of_int k)) in
      let c = cost_function t in
      let r = Mat.(eye P.m *$ P.r_coeff) |> AD.pack_arr in
      AD.Maths.(u *@ r * c * __dt)
    in
    Some rlu

  let rl_x =
    let rlx ~k ~x ~u:_u =
      let _, x_state = unpack_full_state x 4 in
      let thetas = unpack_pos x in
      let t = AD.Maths.(__dt * F (float_of_int k)) in
      let start = AD.Maths.(F 1. * sigmoid ((F P.t_prep - t) / F 2E-4)) in
      let pause_s =
        AD.Maths.(
          F 0.1 / F pause * sigmoid ((t - tau1) / F 2E-4) * sigmoid ((tau2 - t) / F 2E-4))
      in
      let pause_f = AD.Maths.(F 1. * sigmoid ((t - tau3) / F 2E-4)) in
      let dx_p1 = AD.Maths.(tgt_pos1 - thetas) in
      let dx_p2 = AD.Maths.(tgt_pos2 - thetas) in
      let __c = AD.pack_arr P.c in
      let dx_vel = AD.Maths.(F 1. * unpack_vel x) in
      let dx_start = AD.Maths.(in_pos - thetas) in
      let r_xstate =
        AD.Maths.(
          x_state
          *@ transpose __c
          *@ t_mat
          *@ __c
          * F 2.
          * sigmoid ((F P.t_prep - t) / F 2E-3))
      in
      let r_xp1 = AD.Maths.(AD.F 2. * neg dx_p1 *@ q * F 1. * pause_s) in
      let r_xp2 = AD.Maths.(AD.F 2. * neg dx_start *@ q_start * start) in
      let r_xp3 = AD.Maths.(AD.F 2. * neg dx_p2 *@ q * pause_f) in
      let r_vx1 = AD.Maths.(dx_vel *@ q * F 2. * ((pause_s * F 1.) + pause_f)) in
      let r_xv2 = AD.Maths.(dx_vel *@ q_start * start * F 2.) in
      AD.Maths.(
        concatenate ~axis:1 [| r_xp1 + r_xp2 + r_xp3; r_vx1 + r_xv2; r_xstate |] * __dt)
    in
    Some rlx

  let rl_uu =
    let rluu ~k ~x:_x ~u:_u =
      let t = AD.Maths.(__dt * F (float_of_int k)) in
      let c = cost_function t in
      let ma = Mat.(eye P.m *$ P.r_coeff) |> AD.pack_arr in
      AD.Maths.(ma * c * __dt)
    in
    Some rluu

  let rl_ux =
    let f ~k:_k ~x:_x ~u:_u = AD.F 0. in
    Some f

  let rl_xx =
    let rlxx ~k ~x:_x ~u:_u =
      let t = AD.Maths.(__dt * F (float_of_int k)) in
      let start = AD.Maths.(sigmoid ((F P.t_prep - t) / F 2E-4)) in
      let pause_s =
        AD.Maths.(
          F 0.1 / F pause * sigmoid ((t - tau1) / F 2E-4) * sigmoid ((tau2 - t) / F 2E-4))
      in
      let pause_f = AD.Maths.(F 1. * sigmoid ((t - tau3) / F 2E-4)) in
      let mu =
        AD.Maths.(
          (AD.Mat.eye 2 * AD.F Defaults.q_coeff * __dt * F 1. * pause_s)
          + (AD.Mat.eye 2 * AD.F Defaults.q_coeff * __dt * pause_f)
          + (AD.F P.qs_coeff * AD.Mat.eye 2 * start * __dt))
      in
      let mv =
        AD.Maths.(
          (__dt * AD.Mat.eye 2 * F q_coeff * ((pause_s * F 1.) + pause_f))
          + (__dt * AD.Mat.eye 2 * F P.qs_coeff * start))
      in
      let mx =
        AD.Maths.(
          ((sigmoid ((F P.t_prep - t) / F 2E-3)
           * (transpose P.__c *@ t_mat * F 2. *@ P.__c))
          + (F Defaults.a_coeff * AD.Mat.eye P.m))
          * __dt)
      in
      let mf1 =
        AD.Maths.concatenate ~axis:1 [| mu; AD.Mat.zeros (n_theta - 2) (P.m + 2) |]
      in
      let mf2 =
        AD.Maths.concatenate
          ~axis:1
          [| AD.Mat.zeros (n_theta - 2) 2
           ; AD.Maths.(F 1. * mv)
           ; AD.Mat.zeros (n_theta - 2) P.m
          |]
      in
      let mf3 = AD.Maths.concatenate ~axis:1 [| AD.Mat.zeros P.m n_theta; mx |] in
      AD.Maths.concatenate ~axis:0 [| mf1; mf2; mf3 |]
    in
    Some rlxx

  let final_cost ~x:_ ~k:_k =
    let fl = AD.F 0. in
    fl

  let fl_x = None

  (*let f ~k:_k ~x:_x = AD.F 0. in
    Some f*)

  let fl_xx = None
end
