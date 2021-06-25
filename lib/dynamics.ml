open Base
open Owl
include Dynamics_typ
module AD = Algodiff.D

module Integrate (D : Dynamics_T) = struct
  let integrate ~prms ~task =
    let dyn_k = D.dyn ~task ~theta:prms in
    fun ~n ~u ->
      (* u is a TxN matrix, giving T+1 timesteps,
         but actually the first one is going to give x0 *)
      let n_steps = AD.Mat.row_num u in
      let x0 = AD.Mat.zeros 1 n in
      let us =
        AD.Maths.split ~axis:0 (Array.init n_steps ~f:(fun _ -> 1)) u |> Array.to_list
      in
      let rec dyn k x xs us =
        match us with
        | [] -> List.rev xs
        | u :: unexts ->
          let new_x = dyn_k ~k ~x ~u in
          dyn (k + 1) new_x (new_x :: xs) unexts
      in
      dyn 0 x0 [] us |> Array.of_list |> AD.Maths.concatenate ~axis:0
end

module Arm_Linear = struct
  module P = Owl_parameters.Make (Arm_Linear_P)
  open Arm_Linear_P
  module M = Arm.Make (Arm.Defaults)

  let requires_linesearch = true

  let minv ~x =
    let thetas = AD.Maths.get_slice [ []; [ 0; 3 ] ] x in
    let st = Arm.pack_state thetas in
    AD.Linalg.inv (M.inertia st)


  let ms ~prms ~x =
    let open AD.Maths in
    let xs = AD.Maths.get_slice [ []; [ 4; -1 ] ] x in
    let thetas = AD.Maths.get_slice [ []; [ 0; 3 ] ] x in
    let st = Arm.pack_state thetas in
    let s = Arm.pack_state thetas in
    let c = Owl_parameters.extract prms.c in
    let tau = AD.Maths.(c *@ transpose xs) |> AD.Maths.transpose in
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


  let dyn ~theta ~task =
    let b = Owl_parameters.extract theta.b in
    let c = Owl_parameters.extract theta.c in
    let tau = task.tau in
    let a = Owl_parameters.extract theta.a in
    let a = AD.Maths.(a / AD.F tau) in
    let t_prep = task.t_prep in
    let dt = task.dt in
    let _dt = AD.F dt in
    let n_prep = Float.to_int (t_prep /. dt) in
    fun ~k ~x ~u ->
      let xs = AD.Maths.get_slice [ []; [ 4; -1 ] ] x in
      let thetas = AD.Maths.get_slice [ []; [ 0; 3 ] ] x in
      let xst = AD.Maths.transpose xs in
      let s = Arm.pack_state thetas in
      let dx = AD.Maths.(((xs *@ a) + (u *@ b)) * _dt) in
      let tau = AD.Maths.(c *@ xst) |> AD.Maths.transpose in
      let dotdot = M.theta_dot_dot s tau in
      let new_thetas =
        if k < n_prep
        then thetas
        else
          AD.Maths.(
            of_arrays
              [| [| s.x1 + (_dt * s.x1_dot)
                  ; s.x2 + (_dt * s.x2_dot)
                  ; s.x1_dot + (_dt * AD.Maths.get_item dotdot 0 0)
                  ; s.x2_dot + (_dt * AD.Maths.get_item dotdot 1 0)
                 |]
              |])
      and new_x = AD.Maths.(xs + dx) in
      AD.Maths.(concatenate ~axis:1 [| new_thetas; new_x |])


  let dyn_x =
    (* Marine to check this *)
    let _dyn_x ~theta ~task =
      let tau = task.tau in
      let a = Owl_parameters.extract theta.a in
      let a = AD.Maths.(a / AD.F tau) in
      let at = AD.Maths.transpose a in
      let ms = ms ~prms:theta in
      let m = AD.Mat.row_num a in
      let n = m + 4 in
      let t_prep = task.t_prep in
      let dt = task.dt in
      let _dt = AD.F dt in
      let n_prep = Float.to_int (t_prep /. dt) in
      let c = Owl_parameters.extract theta.c in
      fun ~k ~x ~u:_ ->
        let nminv = AD.Maths.(neg (minv ~x)) in
        let m21, m22 = ms ~x in
        let switch = if k < n_prep then AD.F 0. else AD.F 1. in
        let b1 =
          AD.Maths.concatenate
            ~axis:1
            [| AD.Mat.zeros 2 2; AD.Mat.eye 2; AD.Mat.zeros 2 m |]
        and b2 =
          AD.Maths.(
            concatenate
              ~axis:1
              [| nminv *@ m21; nminv *@ m22; neg AD.Maths.(nminv *@ c) |])
        and b3 =
          AD.Maths.concatenate ~axis:1 [| AD.Mat.zeros m 2; AD.Mat.zeros m 2; at |]
        in
        let mat = AD.Maths.(concatenate ~axis:0 [| b1 * switch; switch * b2; b3 |]) in
        AD.Maths.((mat * _dt) + AD.Mat.eye n) |> AD.Maths.transpose
    in
    Some _dyn_x


  let dyn_u =
    let dyn_u ~theta ~task =
      let b = Owl_parameters.extract theta.b in
      let m = AD.Mat.col_num b in
      let dt = AD.F task.dt in
      let mat =
        AD.Maths.concatenate ~axis:0 [| AD.Mat.zeros 2 m; AD.Mat.zeros 2 m; b |]
      in
      fun ~k:_ ~x:_ ~u:_ -> AD.Maths.(transpose (mat * dt))
    in
    Some dyn_u
end
