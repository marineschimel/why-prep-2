open Owl
module AD = Algodiff.D
module M = Arm.Make (Arm.Defaults)
module T = Core.Time
open Typ
open Lib
open Defaults

let _ = Printexc.record_backtrace true


module Make_4D (P : Prms) (D : C_Prms) = struct
  include P
  module C = D

  let n = P.n
  let m = P.m
  let n_theta = 4
  let rl_x = C.rl_x
    let fl_xx = C.fl_xx
  let fl_x = C.fl_x
  let dt = sampling_dt
  let __dt = AD.F dt

  let __c = P.__c

  let n_steps = int_of_float (P.duration /. dt)


  (*let switch_function t = AD.Maths.(sigmoid ((t - F P.t_prep) / F 2E-3))*)

  let rl_u = 
C.rl_u
  let rl_uu = C.rl_uu
  let rl_ux = C.rl_ux
  let rl_xx = None

  let minv ~x =
    let thetas, _ = unpack_full_state x n_theta in
    let st = Arm.pack_state thetas in
    AD.Linalg.inv (M.inertia st)





 
  let ms ~x =
    let open AD.Maths in
    let _, xs = unpack_full_state x n_theta in
    let thetas = AD.Maths.get_slice [[];[0;3]] x in 
    let st = Arm.pack_state thetas in
    let s = Arm.pack_state thetas in
    let tau = AD.Maths.(cmc*@(phi (__c *@ g (transpose xs)))) |> AD.Maths.transpose in
    let dotdot = M.theta_dot_dot s tau in
    let dotdot1, dotdot2 = AD.Mat.get dotdot 0 0, AD.Mat.get dotdot 1 0 in
    let m12 =
      let cst1 = AD.F M._A2 * sin st.x2 * AD.F (-1.)
      and cst2 = AD.F M._A2 * cos st.x2
      and mat2 =
        [| [| AD.F 0.; neg st.x2_dot * ((F 2. * st.x1_dot) + st.x2_dot) |]
         ; [| AD.F 0.; sqr st.x1_dot |]
        |]
        |> of_arrays
      and mat1 =
        [| [| AD.F 0.; (F 2. * dotdot1) + dotdot2 |]; [| F 0.; dotdot1 |] |] |> of_arrays
      in
      (cst2 * mat2) + (cst1 * mat1)
    in
    let m22 =
      let cst1 = F M._A2 * sin st.x2 in
      let mat1 =
        [| [| AD.F 2. * neg st.x2_dot; AD.F 2. * neg (st.x2_dot + st.x1_dot) |]
         ; [| AD.F 2. * st.x1_dot; AD.F 0. |]
        |]
        |> of_arrays
      in
      (cst1 * mat1) + M._B
    in
    m12, m22


  let _diff_g x =
    let _, xst = unpack_full_state x n_theta in
    AD.Mat.init_2d (n - n_theta) (n - n_theta) (fun i j ->
        if i = j then g_prm ~x:(AD.Mat.get xst 0 i) else F 0.)


  let dyn ~k ~x ~u =
    let _, xs = unpack_full_state x n_theta in
    let thetas = AD.Maths.get_slice [[];[0;3]] x in 
    let xs = AD.Maths.get_slice [ []; [ 0; size_net - 1 ] ] xs in
    let xst = AD.Maths.transpose xs in
    let s = Arm.pack_state thetas in
     let __a = AD.Maths.((((AD.pack_arr P.w)) - (AD.Mat.eye size_net))/F tau) in 
    let t = AD.Maths.(F (float_of_int k)*__dt) in 
    let _switch = if t < AD.F P.t_prep then AD.F 1. else AD.F 1. in 
    let _input_switched =
      let t = AD.Maths.(__dt * F (float_of_int k)) in
      if t < F t_prep then AD.F 1. else AD.F 1.
    in
  let dx =
      AD.Maths.(((__a *@ xst) + (Defaults.__b *@ transpose u * _input_switched)) * __dt)
      |> AD.Maths.transpose
    in

    let tau = AD.Maths.(cmc*@phi(__c *@ g xst)) |> AD.Maths.transpose in
    let dotdot = M.theta_dot_dot s tau in
    let new_thetas =
      AD.Maths.(
        of_arrays
          [| [| s.x1 + (__dt * s.x1_dot * _switch)
              ; s.x2 + (__dt * s.x2_dot * _switch)
              ; s.x1_dot + (__dt * AD.Maths.get_item dotdot 0 0 * _switch)
              ; s.x2_dot + (__dt * AD.Maths.get_item dotdot 1 0 * _switch);
             |]
          |])
    and new_x = AD.Maths.(xs + dx) (*and new_u = u*) in
    AD.Maths.(concatenate ~axis:1 [| new_thetas; new_x (*; new_u *) |])


  let _lin_dyn_x ~k ~x ~u:_u =
    let t = AD.Maths.(F (float_of_int k)*__dt) in 
    let _,x_state = unpack_full_state x n_theta in 
    let _switch = if t < AD.F P.t_prep then AD.F 1. else AD.F 1. 
in
      (*let t = AD.Maths.(__dt * F (float_of_int k)) in
        AD.Maths.(sigmoid ((t - F P.t_prep) / F 2E-3))*)
    
    let nminv = AD.Maths.(neg (minv ~x)) in
    let m21, m22 = ms ~x in
    let __ap = AD.Maths.((((AD.pack_arr P.w)) - (AD.Mat.eye size_net))/F tau) in  
    (*((Defaults.__w *@ diff_g x) - AD.Mat.eye P.m)*)
    let b1 =
      AD.Maths.concatenate ~axis:1 [| AD.Mat.zeros 2 2; AD.Mat.eye 2; AD.Mat.zeros 2 size_net |]
    and b2 =
      AD.Maths.(
        concatenate
          ~axis:1
          [| nminv *@ m21; nminv *@ m22; neg AD.Maths.((nminv *@ (cmc *@ (dphi (__c*@(transpose x_state)))) *@ __c)) |])
    and b3 =
      AD.Maths.concatenate ~axis:1 [| AD.Mat.zeros size_net 2; AD.Mat.zeros size_net 2; __ap |]
    in
    let mat = AD.Maths.(concatenate ~axis:0 [| b1 * _switch; _switch * b2; b3 |]) in
    AD.Maths.((mat * __dt) + AD.Mat.eye n) |> AD.Maths.transpose


  let _lin_dyn_u ~k ~x:_x ~u:_u =
    let _input_switched =
      let t = AD.Maths.(__dt * F (float_of_int k)) in
      if t < F t_prep then AD.F 1. else AD.F 1.
    in
    (*let switch =
        let t = AD.Maths.(__dt * F (float_of_int k)) in
        AD.Maths.(sigmoid ((t - F P.t_prep) / F 2E-3))
      in*)
    let m =
      AD.Maths.concatenate
        ~axis:0
        [| AD.Mat.zeros 2 m
         ; AD.Mat.zeros 2 m
         ; AD.Maths.(Defaults.__b * _input_switched)
        |]
    in
    AD.Maths.(transpose (m * __dt))
  let dyn_x =Some _lin_dyn_x
  let dyn_u = Some _lin_dyn_u
  let running_loss ~k:_k ~x ~u = C.cost ~u ~x ~k:_k
  let final_loss ~k:_k ~x = C.final_cost ~x ~k:_k
end

module Make (P : Prms) (D : C_Prms) = struct
  include P
  module C = D

  let n = P.n
  let m = P.m
  let n_theta = 4
  let rl_x =None
    let fl_xx = C.fl_xx
  let fl_x = C.fl_x
  let dt = sampling_dt
  let __dt = AD.F dt

  let __c = P.__c
  let n_steps = int_of_float (P.duration /. dt)

  (*let switch_function t = AD.Maths.(sigmoid ((t - F P.t_prep) / F 2E-3))*)

  let rl_u = 
C.rl_u
  let rl_uu = C.rl_uu
  let rl_ux = C.rl_ux
  let rl_xx = C.rl_xx

  let minv ~x =
    let thetas, _ = unpack_full_state x n_theta in
    let st = Arm.pack_state thetas in
    AD.Linalg.inv (M.inertia st)



  let ms ~x =
    let open AD.Maths in
    let _, xs = unpack_full_state x n_theta in
    let thetas = AD.Maths.get_slice [[];[0;3]] x in 
    let st = Arm.pack_state thetas in
    let s = Arm.pack_state thetas in
    let tau = AD.Maths.(__c *@ g (transpose xs)) |> AD.Maths.transpose in
    let dotdot = M.theta_dot_dot s tau in
    let dotdot1, dotdot2 = AD.Mat.get dotdot 0 0, AD.Mat.get dotdot 1 0 in
    let m12 =
      let cst1 = AD.F M._A2 * sin st.x2 * AD.F (-1.)
      and cst2 = AD.F M._A2 * cos st.x2
      and mat2 =
        [| [| AD.F 0.; neg st.x2_dot * ((F 2. * st.x1_dot) + st.x2_dot) |]
         ; [| AD.F 0.; sqr st.x1_dot |]
        |]
        |> of_arrays
      and mat1 =
        [| [| AD.F 0.; (F 2. * dotdot1) + dotdot2 |]; [| F 0.; dotdot1 |] |] |> of_arrays
      in
      (cst2 * mat2) + (cst1 * mat1)
    in
    let m22 =
      let cst1 = F M._A2 * sin st.x2 in
      let mat1 =
        [| [| AD.F 2. * neg st.x2_dot; AD.F 2. * neg (st.x2_dot + st.x1_dot) |]
         ; [| AD.F 2. * st.x1_dot; AD.F 0. |]
        |]
        |> of_arrays
      in
      (cst1 * mat1) + M._B
    in
    m12, m22


  let diff_g x =
    let _, xst = unpack_full_state x n_theta in
    AD.Mat.init_2d (n - n_theta) (n - n_theta) (fun i j ->
        if i = j then g_prm ~x:(AD.Mat.get xst 0 i) else F 0.)


  let dyn ~k ~x ~u =
    let _, xs = unpack_full_state x n_theta in
    let thetas = AD.Maths.get_slice [[];[0;3]] x in 
    let xs = AD.Maths.get_slice [ []; [ 0; size_net - 1 ] ] xs in
    let xst = AD.Maths.transpose xs in
    let s = Arm.pack_state thetas in
     let __a x = AD.Maths.((((AD.pack_arr P.w))*@(g x) - (AD.Mat.eye size_net)*@x)/F tau) in 
    let t = AD.Maths.(F (float_of_int k)*__dt) in 
    let _switch = if t < AD.F P.t_prep then AD.F 1. else AD.F 1. in 
    let _input_switched =
      let t = AD.Maths.(__dt * F (float_of_int k)) in
      if t < F t_prep then AD.F 1. else AD.F 1.
    in
  let dx =
      AD.Maths.(((__a xst) + (Defaults.__b *@ transpose u * _input_switched)) * __dt)
      |> AD.Maths.transpose
    in
    let tau = AD.Maths.(__c *@ g xst) |> AD.Maths.transpose in
    let dotdot = M.theta_dot_dot s tau in
    let new_thetas =
      AD.Maths.(
        of_arrays
          [| [| s.x1 + (__dt * s.x1_dot * _switch)
              ; s.x2 + (__dt * s.x2_dot * _switch)
              ; s.x1_dot + (__dt * AD.Maths.get_item dotdot 0 0 * _switch)
              ; s.x2_dot + (__dt * AD.Maths.get_item dotdot 1 0 * _switch);
             |]
          |])
    and new_x = AD.Maths.(xs + dx) (*and new_u = u*) in
    AD.Maths.(concatenate ~axis:1 [| new_thetas; new_x (*; new_u *) |])


  let _lin_dyn_x ~k ~x ~u:_u =
    let t = AD.Maths.(F (float_of_int k)*__dt) in 
    let _switch = if t < AD.F P.t_prep then AD.F 1. else AD.F 1. 
in
      (*let t = AD.Maths.(__dt * F (float_of_int k)) in
        AD.Maths.(sigmoid ((t - F P.t_prep) / F 2E-3))*)
    
    let nminv = AD.Maths.(neg (minv ~x)) in
    let m21, m22 = ms ~x in
    let __ap x = AD.Maths.((((AD.pack_arr P.w))*@(diff_g x) - (AD.Mat.eye size_net))/F tau) in  
    (*((Defaults.__w *@ diff_g x) - AD.Mat.eye P.m)*)
    let b1 =
      AD.Maths.concatenate ~axis:1 [| AD.Mat.zeros 2 2; AD.Mat.eye 2; AD.Mat.zeros 2 size_net |]
    and b2 =
      AD.Maths.(
        concatenate
          ~axis:1
          [| nminv *@ m21; nminv *@ m22; neg AD.Maths.(nminv *@ __c *@ diff_g x) |])
    and b3 =
      AD.Maths.concatenate ~axis:1 [| AD.Mat.zeros size_net 2; AD.Mat.zeros size_net 2; __ap x |]
    in
    let mat = AD.Maths.(concatenate ~axis:0 [| b1 * _switch; _switch * b2; b3 |]) in
    AD.Maths.((mat * __dt) + AD.Mat.eye n) |> AD.Maths.transpose


  let _lin_dyn_u ~k ~x:_x ~u:_u =
    let _input_switched =
      let t = AD.Maths.(__dt * F (float_of_int k)) in
      if t < F t_prep then AD.F 1. else AD.F 1.
    in
    (*let switch =
        let t = AD.Maths.(__dt * F (float_of_int k)) in
        AD.Maths.(sigmoid ((t - F P.t_prep) / F 2E-3))
      in*)
    let m =
      AD.Maths.concatenate
        ~axis:0
        [| AD.Mat.zeros 2 m
         ; AD.Mat.zeros 2 m
         ; AD.Maths.(Defaults.__b * _input_switched)
        |]
    in
    AD.Maths.(transpose (m * __dt))
  let dyn_x = Some _lin_dyn_x
  let dyn_u = Some _lin_dyn_u
  let running_loss ~k:_k ~x ~u = C.cost ~u ~x ~k:_k
  let final_loss ~k:_k ~x = C.final_cost ~x ~k:_k
end
