open Base
open Owl
include Likelihoods_typ
module AD = Algodiff.D

module End (X : sig
  val label : string
end) =
struct
  module P = Owl_parameters.Make (End_P)
  open End_P

  let requires_linesearch = false
  let label = X.label

  let neg_logp_t ~readout ~prms ~task =
    let t_prep = task.t_prep in
    let target_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] task.target in
    let theta0 = task.theta0 in
    let qs_coeff = Owl_parameters.extract prms.End_P.qs_coeff in
    let g_coeff = Owl_parameters.extract prms.End_P.g_coeff in
    let t_coeff = Owl_parameters.extract prms.End_P.t_coeff in
    let dt = task.dt in
    let n_prep = Float.to_int (t_prep /. dt) in
    let n_mov =
      let d_mov = Float.to_int (task.t_movs.(0) /. dt) in
      n_prep + d_mov
    in
    let c = readout in
    let c_t = AD.Maths.transpose c in
    fun ~k ~z_t ->
      let thetas = AD.Maths.get_slice [ []; [ 0; 3 ] ] z_t in
      let theta_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] thetas in
      let theta_vel = AD.Maths.get_slice [ []; [ 2; 3 ] ] thetas in
      let x_t = AD.Maths.get_slice [ []; [ 4; -1 ] ] z_t in
      if k < n_prep
      then (
        let mu_t = AD.Maths.(x_t *@ c_t) in
        AD.Maths.(
          (F 0.5 * sum' (t_coeff * sqr mu_t))
          + (F 0.5 * qs_coeff * l2norm_sqr' (thetas - theta0))))
      else if k > n_mov
      then
        AD.Maths.(
          F 0.5
          * g_coeff
          * (l2norm_sqr' (theta_pos - target_pos) + (F 0.1 * l2norm_sqr' theta_vel)))
      else AD.F 0.


  let neg_jac_t =
    let _neg_jac_t ~readout ~prms ~task =
      let t_prep = task.t_prep in
      let target_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] task.target in
      let theta0 = AD.Maths.get_slice [ []; [ 0; 1 ] ] task.theta0 in
      let qs_coeff = Owl_parameters.extract prms.qs_coeff in
      let g_coeff = Owl_parameters.extract prms.g_coeff in
      let t_coeff = Owl_parameters.extract prms.t_coeff in
      let dt = task.dt in
      let n_prep = Float.to_int (t_prep /. dt) in
      let n_mov =
        let d_mov = Float.to_int (task.t_movs.(0) /. dt) in
        n_prep + d_mov
      in
      let c = readout in
      let m = AD.Mat.col_num c in
      let c_t = AD.Maths.transpose c in
      fun ~k ~z_t ->
        let thetas = AD.Maths.get_slice [ []; [ 0; 3 ] ] z_t in
        let theta_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] thetas in
        let theta_vel = AD.Maths.get_slice [ []; [ 2; 3 ] ] thetas in
        let x = AD.Maths.get_slice [ []; [ 4; -1 ] ] z_t in
        let r_state =
          if k < n_prep then AD.Maths.(t_coeff * x *@ c_t *@ c) else AD.Mat.zeros 1 m
        in
        let r_xp =
          if k < n_prep
          then AD.Maths.(qs_coeff * (theta_pos - theta0))
          else if k > n_mov
          then AD.Maths.(g_coeff * (theta_pos - target_pos))
          else AD.Mat.zeros 1 2
        in
        let r_xv =
          if k < n_prep
          then AD.Maths.(qs_coeff * theta_vel)
          else if k > n_mov
          then AD.Maths.(F 0.1 * g_coeff * theta_vel)
          else AD.Mat.zeros 1 2
        in
        AD.Maths.(concatenate ~axis:1 [| r_xp; r_xv; r_state |])
    in
    Some _neg_jac_t


  let neg_hess_t =
    let _neg_hess_t ~readout ~prms ~task =
      let c = readout in
      let m = AD.Mat.col_num c in
      let t_prep = task.t_prep in
      let dt = task.dt in
      let t_coeff = Owl_parameters.extract prms.t_coeff in
      let qs_coeff = Owl_parameters.extract prms.qs_coeff in
      let g_coeff = Owl_parameters.extract prms.g_coeff in
      let n_prep = Float.to_int (t_prep /. dt) in
      let n_mov =
        let d_mov = Float.to_int (task.t_movs.(0) /. dt) in
        n_prep + d_mov
      in
      let c_term = AD.Maths.(t_coeff * transpose c *@ c) in
      fun ~k ~z_t:_ ->
        let mu =
          if k < n_prep
          then AD.Maths.(qs_coeff * AD.Mat.eye 2)
          else if k > n_mov
          then AD.Maths.(g_coeff * AD.Mat.eye 2)
          else AD.Mat.zeros 2 2
        in
        let mv =
          if k < n_prep
          then AD.Maths.(qs_coeff * AD.Mat.eye 2)
          else if k > n_mov
          then AD.Maths.(AD.F 0.1 * g_coeff * AD.Mat.eye 2)
          else AD.Mat.zeros 2 2
        in
        let mx = if k < n_prep then c_term else AD.Mat.zeros m m in
        let mf1 = AD.Maths.concatenate ~axis:1 [| mu; AD.Mat.zeros 2 (m + 2) |] in
        let mf2 =
          AD.Maths.concatenate ~axis:1 [| AD.Mat.zeros 2 2; mv; AD.Mat.zeros 2 m |]
        in
        let mf3 = AD.Maths.concatenate ~axis:1 [| AD.Mat.zeros m 4; mx |] in
        AD.Maths.concatenate ~axis:0 [| mf1; mf2; mf3 |]
    in
    Some _neg_hess_t
end

module Soft_End (X : sig
  val label : string
end) =
struct
  module P = Owl_parameters.Make (Soft_End_P)
  open Soft_End_P

  let requires_linesearch = false
  let label = X.label

  let start ~task t =
    let t_prep = AD.F task.t_prep in
    AD.Maths.(sigmoid ((t_prep - t) / F 2E-4))


  let finish ~task t =
    let t_mov = AD.F (task.t_prep +. task.t_movs.(0)) in
    AD.Maths.(sigmoid ((t - t_mov) / F 20E-3))


  let neg_logp_t ~readout ~prms ~task =
    let target_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] task.target in
    let theta0 = task.theta0 in
    let qs_coeff = Owl_parameters.extract prms.qs_coeff in
    let g_coeff = Owl_parameters.extract prms.g_coeff in
    let t_coeff = Owl_parameters.extract prms.t_coeff in
    let dt = task.dt in
    let c = readout in
    let c_t = AD.Maths.transpose c in
    fun ~k ~z_t ->
      let t = AD.F Float.(of_int k *. dt) in
      let thetas = AD.Maths.get_slice [ []; [ 0; 3 ] ] z_t in
      let theta_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] thetas in
      let theta_vel = AD.Maths.get_slice [ []; [ 2; 3 ] ] thetas in
      let mu_t = AD.Maths.(z_t *@ c_t) in
      AD.Maths.(
        (start ~task t
        * ((F 0.5 * sum' (t_coeff * sqr mu_t))
          + (F 0.5 * qs_coeff * l2norm_sqr' (thetas - theta0))))
        + (finish ~task t
          * (F 0.5
            * g_coeff
            * (l2norm_sqr' (theta_pos - target_pos) + (F 0.1 * l2norm_sqr' theta_vel)))))


  let neg_jac_t =
    let _neg_jac_t ~readout ~prms ~task =
      let target_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] task.target in
      let theta0 = AD.Maths.get_slice [ []; [ 0; 1 ] ] task.theta0 in
      let qs_coeff = Owl_parameters.extract prms.qs_coeff in
      let g_coeff = Owl_parameters.extract prms.g_coeff in
      let t_coeff = Owl_parameters.extract prms.t_coeff in
      let dt = task.dt in
      let c = readout in
      let c_t = AD.Maths.transpose c in
      fun ~k ~z_t ->
        let t = AD.F Float.(of_int k *. dt) in
        let thetas = AD.Maths.get_slice [ []; [ 0; 3 ] ] z_t in
        let theta_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] thetas in
        let theta_vel = AD.Maths.get_slice [ []; [ 2; 3 ] ] thetas in
        let x = AD.Maths.get_slice [ []; [ 4; -1 ] ] z_t in
        let r_state = AD.Maths.(start ~task t * (t_coeff * x *@ c_t *@ c)) in
        let r_xp =
          AD.Maths.(
            (start ~task t * (qs_coeff * (theta_pos - theta0)))
            + (finish ~task t * (g_coeff * (theta_pos - target_pos))))
        in
        let r_xv =
          AD.Maths.(
            (start ~task t * (qs_coeff * theta_vel))
            + (finish ~task t * (F 0.1 * g_coeff * theta_vel)))
        in
        AD.Maths.(concatenate ~axis:1 [| r_xp; r_xv; r_state |])
    in
    Some _neg_jac_t


  let neg_hess_t =
    let _neg_hess_t ~readout ~prms ~task =
      let c = readout in
      let m = AD.Mat.col_num c in
      let dt = task.dt in
      let t_coeff = Owl_parameters.extract prms.t_coeff in
      let qs_coeff = Owl_parameters.extract prms.qs_coeff in
      let g_coeff = Owl_parameters.extract prms.g_coeff in
      let c_term = AD.Maths.(t_coeff * transpose c *@ c) in
      fun ~k ~z_t:_ ->
        let t = AD.F Float.(of_int k *. dt) in
        let mu =
          AD.Maths.(
            (start ~task t * (qs_coeff * AD.Mat.eye 2))
            + (finish ~task t * (g_coeff * AD.Mat.eye 2)))
        in
        let mv =
          AD.Maths.(
            (start ~task t * (qs_coeff * AD.Mat.eye 2))
            + (finish ~task t * (AD.F 0.1 * g_coeff * AD.Mat.eye 2)))
        in
        let mx = AD.Maths.(start ~task t * c_term) in
        let mf1 = AD.Maths.concatenate ~axis:1 [| mu; AD.Mat.zeros 2 (m + 2) |] in
        let mf2 =
          AD.Maths.concatenate ~axis:1 [| AD.Mat.zeros 2 2; mv; AD.Mat.zeros 2 m |]
        in
        let mf3 = AD.Maths.concatenate ~axis:1 [| AD.Mat.zeros m 4; mx |] in
        AD.Maths.concatenate ~axis:0 [| mf1; mf2; mf3 |]
    in
    Some _neg_hess_t
end

module Match (X : sig
  val label : string
end) =
struct
  module P = Owl_parameters.Make (Match_P)
  open Match_P

  let requires_linesearch = false
  let label = X.label

  let neg_logp_t ~readout:_ ~prms ~task =
    let q_coeff = Owl_parameters.extract prms.q_coeff in
    let tgt = task.target in
    fun ~k ~z_t ->
      let target_k = AD.Maths.get_slice [ [ k ]; [] ] tgt in
      AD.Maths.(F 0.5 * q_coeff * l2norm_sqr' (z_t - target_k))


  let neg_jac_t =
    let _neg_jac_t ~readout:_ ~prms ~task =
      let q_coeff = Owl_parameters.extract prms.q_coeff in
      let tgt = task.target in
      fun ~k ~z_t ->
        let target_k = AD.Maths.get_slice [ [ k ]; [] ] tgt in
        AD.Maths.(q_coeff * (z_t - target_k))
    in
    Some _neg_jac_t


  let neg_hess_t =
    let _neg_hess_t ~readout:_ ~prms ~task =
      let q_coeff = Owl_parameters.extract prms.q_coeff in
      let tgt = task.target in
      let n = AD.Mat.col_num tgt in
      let mat = AD.Maths.(q_coeff * AD.Mat.eye n) in
      fun ~k:_ ~z_t:_ -> mat
    in
    Some _neg_hess_t
end

module Match_Torques (X : sig
  val label : string
end) =
struct
  module P = Owl_parameters.Make (Match_Torques_P)
  open Match_Torques_P

  let requires_linesearch = false
  let label = X.label

  let neg_logp_t ~readout ~prms ~task =
    let q_coeff = Owl_parameters.extract prms.q_coeff in
    let g_coeff = Owl_parameters.extract prms.g_coeff in
    let t_prep = task.t_prep in
    let dt = task.dt in
    let n_prep = Float.to_int (t_prep /. dt) in
    let c = readout in
    let tgt = task.target in
    let n_mov =
      let d_mov = Float.to_int (task.t_movs.(0) /. dt) in
      n_prep + d_mov
    in
    let ct = AD.Maths.transpose c in
    fun ~k ~z_t ->
      if k < n_prep || k > n_mov
      then (
        let mu_t = AD.Maths.(z_t *@ ct) in
        AD.Maths.(F 0.5 * q_coeff * l2norm_sqr' mu_t))
      else (
        let mu_t = AD.Maths.(z_t *@ ct) in
        let target_k = AD.Maths.get_slice [ [ k - n_prep ]; [] ] tgt in
        AD.Maths.(F 0.5 * g_coeff * l2norm_sqr' (mu_t - target_k)))


  let neg_jac_t =
    None


  let neg_hess_t =
    None
end

module Successive (X : sig
  val label : string
  val phi_x : AD.t -> AD.t
  val d_phi_x : AD.t -> AD.t
  val d2_phi_x : AD.t -> AD.t
end) =
struct
  module P = Owl_parameters.Make (Successive_P)
  open Successive_P

  let requires_linesearch = false
  let label = X.label

  let neg_logp_t ~readout ~prms ~task =
    let n_prep = Float.to_int (task.t_prep /. task.dt) in
    let n_1 = Float.to_int ((task.t_movs.(0) +. task.t_prep) /. task.dt) in
    let pause_0 =
      match task.t_pauses with
      | Some t -> t.(0)
      | None -> 0.
    in
    let n_2 = n_1 + Float.to_int (pause_0 /. task.dt) in
    let n_3 = n_2 + Float.to_int (task.t_movs.(1) /. task.dt) in
    let qs_coeff = Owl_parameters.extract prms.qs_coeff in
    let g_coeff = Owl_parameters.extract prms.g_coeff in
    let t_coeff = Owl_parameters.extract prms.t_coeff in
    let c = readout in
    let c_t = AD.Maths.transpose c in
    let tgt1 =
      AD.Mat.row task.target 0 |> fun z -> AD.Maths.get_slice [ []; [ 0; 1 ] ] z
    in
    let tgt2 =
      AD.Mat.row task.target 1 |> fun z -> AD.Maths.get_slice [ []; [ 0; 1 ] ] z
    in
    let dt = task.dt in
    let _dt = AD.F dt in
    let theta0 = task.theta0 in
    let g_coeff_1 = AD.Maths.(g_coeff * F 0.2 / F pause_0) in
    fun ~k ~z_t ->
      let thetas = AD.Maths.get_slice [ []; [ 0; 3 ] ] z_t in
      let theta_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] thetas in
      let theta_vel = AD.Maths.get_slice [ []; [ 2; 3 ] ] thetas in
      let x_t = AD.Maths.get_slice [ []; [ 4; -1 ] ] z_t in
      let x_t = AD.Maths.(X.phi_x x_t) in
      if k < n_prep
      then (
        let mu_t = AD.Maths.(x_t *@ c_t) in
        AD.Maths.(
          (F 0.5 * sum' (t_coeff * sqr mu_t))
          + (F 0.5 * qs_coeff * l2norm_sqr' (thetas - theta0))))
      else if k > n_1 && k < n_2
      then
        AD.Maths.(
          F 0.5
          * g_coeff_1
          * (l2norm_sqr' (theta_pos - tgt1) + (F 0.1 * l2norm_sqr' theta_vel)))
      else if k > n_3
      then
        AD.Maths.(
          F 0.5
          * g_coeff
          * (l2norm_sqr' (theta_pos - tgt2) + (F 0.1 * l2norm_sqr' theta_vel)))
      else AD.F 0.


  let neg_jac_t =
    let _neg_jac_t ~readout ~prms ~task =
      let n_prep = Float.to_int (task.t_prep /. task.dt) in
      let n_1 = Float.to_int ((task.t_movs.(0) +. task.t_prep) /. task.dt) in
      let pause_0 =
        match task.t_pauses with
        | Some t -> t.(0)
        | None -> 0.
      in
      let n_2 = n_1 + Float.to_int (pause_0 /. task.dt) in
      let n_3 = n_2 + Float.to_int (task.t_movs.(1) /. task.dt) in
      let qs_coeff = Owl_parameters.extract prms.qs_coeff in
      let g_coeff = Owl_parameters.extract prms.g_coeff in
      let t_coeff = Owl_parameters.extract prms.t_coeff in
      let c = readout in
      let m = AD.Mat.col_num c in
      let c_t = AD.Maths.transpose c in
      let tgt1 =
        AD.Mat.row task.target 0 |> fun z -> AD.Maths.get_slice [ []; [ 0; 1 ] ] z
      in
      let tgt2 =
        AD.Mat.row task.target 1 |> fun z -> AD.Maths.get_slice [ []; [ 0; 1 ] ] z
      in
      let dt = task.dt in
      let _dt = AD.F dt in
      let theta0 = task.theta0 |> AD.Maths.get_slice [ []; [ 0; 1 ] ] in
      fun ~k ~z_t ->
        let thetas = AD.Maths.get_slice [ []; [ 0; 3 ] ] z_t in
        let theta_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] thetas in
        let theta_vel = AD.Maths.get_slice [ []; [ 2; 3 ] ] thetas in
        let x = AD.Maths.get_slice [ []; [ 4; -1 ] ] z_t in
        let x_t = X.phi_x x in
        let dphi = AD.Maths.(diagm (X.d_phi_x x)) in
        let r_state =
          if k < n_prep
          then AD.Maths.(t_coeff * x_t *@ c_t *@ c *@ dphi)
          else AD.Mat.zeros 1 m
        in
        let r_xp =
          if k < n_prep
          then AD.Maths.(qs_coeff * (theta_pos - theta0))
          else if k > n_1 && k < n_2
          then AD.Maths.(g_coeff * (theta_pos - tgt1))
          else if k > n_3
          then AD.Maths.(g_coeff * (theta_pos - tgt2))
          else AD.Mat.zeros 1 2
        in
        let r_xv =
          if k < n_prep
          then AD.Maths.(qs_coeff * theta_vel)
          else if k > n_1 && k < n_2
          then AD.Maths.(F 0.1 * g_coeff * theta_vel)
          else if k > n_3
          then AD.Maths.(F 0.1 * g_coeff * theta_vel)
          else AD.Mat.zeros 1 2
        in
        AD.Maths.(concatenate ~axis:1 [| r_xp; r_xv; r_state |])
    in
    Some _neg_jac_t


  let neg_hess_t =
    let _neg_hess_t ~readout ~prms ~task =
      let n_prep = Float.to_int (task.t_prep /. task.dt) in
      let n_1 = Float.to_int ((task.t_movs.(0) +. task.t_prep) /. task.dt) in
      let pause_0 =
        match task.t_pauses with
        | Some t -> t.(0)
        | None -> 0.
      in
      let n_2 = n_1 + Float.to_int (pause_0 /. task.dt) in
      let n_3 = n_2 + Float.to_int (task.t_movs.(1) /. task.dt) in
      let qs_coeff = Owl_parameters.extract prms.qs_coeff in
      let g_coeff = Owl_parameters.extract prms.g_coeff in
      let t_coeff = Owl_parameters.extract prms.t_coeff in
      let c = readout in
      let dt = task.dt in
      let _dt = AD.F dt in
      let m = AD.Mat.col_num c in
      fun ~k ~z_t ->
        let x_t = AD.Maths.get_slice [ []; [ 4; -1 ] ] z_t in
        let dx = AD.Maths.(diagm (X.d_phi_x x_t)) in
        let ddx = X.d2_phi_x x_t in
        let t1 = AD.Maths.(c *@ dx) in
        let t2 = AD.Maths.(diagm (X.phi_x x_t) *@ transpose c *@ c *@ transpose ddx) in
        let c_term = AD.Maths.((t_coeff * (transpose t1 *@ t1)) + t2) in
        let mu =
          if k < n_prep
          then AD.Maths.(qs_coeff * AD.Mat.eye 2)
          else if k > n_1 && k < n_2
          then AD.Maths.(g_coeff * AD.Mat.eye 2)
          else if k > n_3
          then AD.Maths.(g_coeff * AD.Mat.eye 2)
          else AD.Mat.zeros 2 2
        in
        let mv =
          if k < n_prep
          then AD.Maths.(qs_coeff * AD.Mat.eye 2)
          else if k > n_1 && k < n_2
          then AD.Maths.(AD.F 0.1 * g_coeff * AD.Mat.eye 2)
          else if k > n_3
          then AD.Maths.(AD.F 0.1 * g_coeff * AD.Mat.eye 2)
          else AD.Mat.zeros 2 2
        in
        let mx = if k < n_prep then c_term else AD.Mat.zeros m m in
        let mf1 = AD.Maths.concatenate ~axis:1 [| mu; AD.Mat.zeros 2 (m + 2) |] in
        let mf2 =
          AD.Maths.concatenate ~axis:1 [| AD.Mat.zeros 2 2; mv; AD.Mat.zeros 2 m |]
        in
        let mf3 = AD.Maths.concatenate ~axis:1 [| AD.Mat.zeros m 4; mx |] in
        AD.Maths.concatenate ~axis:0 [| mf1; mf2; mf3 |]
    in
    Some _neg_hess_t
end

module End_Phi (X : sig
  val label : string
  val phi_x : AD.t -> AD.t
  val d_phi_x : AD.t -> AD.t
  val d2_phi_x : AD.t -> AD.t
  val speed_end_penalty : float
end) =
struct
  module P = Owl_parameters.Make (End_Phi_P)
  open End_Phi_P
  open X

  let requires_linesearch = false
  let label = X.label

  let neg_logp_t ~readout ~prms ~task =
    let t_prep = task.t_prep in
    let target_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] task.target in
    let theta0 = task.theta0 in
    let qs_coeff = Owl_parameters.extract prms.End_Phi_P.qs_coeff in
    let g_coeff = Owl_parameters.extract prms.End_Phi_P.g_coeff in
    let t_coeff = Owl_parameters.extract prms.End_Phi_P.t_coeff in
    let dt = task.dt in
    let n_prep = Float.to_int (t_prep /. dt) in
    let n_mov =
      let d_mov = Float.to_int (task.t_movs.(0) /. dt) in
      n_prep + d_mov
    in
    let c = readout in
    let c_t = AD.Maths.transpose c in
    let x0 = task.x0 in
    let x0 = AD.Maths.get_slice [ []; [ 4; -1 ] ] x0 in
    let phi_x0 = phi_x x0 in
    fun ~k ~z_t ->
      let thetas = AD.Maths.get_slice [ []; [ 0; 3 ] ] z_t in
      let theta_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] thetas in
      let theta_vel = AD.Maths.get_slice [ []; [ 2; 3 ] ] thetas in
      let x_t = AD.Maths.get_slice [ []; [ 4; -1 ] ] z_t in
      let x_t = phi_x x_t in
      if k < n_prep
      then (
        let mu_t = AD.Maths.(x_t *@ c_t) in
        AD.Maths.(
          (F 0.5 * sum' (t_coeff * sqr mu_t))
          + (F 0.5 * qs_coeff * l2norm_sqr' (thetas - theta0))))
      else if k > n_mov
      then
        AD.Maths.(
          F 0.5
          * g_coeff
          * (l2norm_sqr' (theta_pos - target_pos)
            + (F speed_end_penalty * l2norm_sqr' theta_vel)))
      else AD.F 0.


  let neg_jac_t =
    let _neg_jac_t ~readout ~prms ~task =
      let t_prep = task.t_prep in
      let target_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] task.target in
      let theta0 = AD.Maths.get_slice [ []; [ 0; 1 ] ] task.theta0 in
      let qs_coeff = Owl_parameters.extract prms.qs_coeff in
      let g_coeff = Owl_parameters.extract prms.g_coeff in
      let t_coeff = Owl_parameters.extract prms.t_coeff in
      let dt = task.dt in
      let n_prep = Float.to_int (t_prep /. dt) in
      let n_mov =
        let d_mov = Float.to_int (task.t_movs.(0) /. dt) in
        n_prep + d_mov
      in
      let c = readout in
      let m = AD.Mat.col_num c in
      let c_t = AD.Maths.transpose c in
      let x0 = task.x0 in
      let x0 = AD.Maths.get_slice [ []; [ 4; -1 ] ] x0 in
      let phi_x0 = phi_x x0 in
      fun ~k ~z_t ->
        let thetas = AD.Maths.get_slice [ []; [ 0; 3 ] ] z_t in
        let theta_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] thetas in
        let theta_vel = AD.Maths.get_slice [ []; [ 2; 3 ] ] thetas in
        let x = AD.Maths.get_slice [ []; [ 4; -1 ] ] z_t in
        let x_t = phi_x x in
        let dphi = AD.Maths.(diagm (d_phi_x x)) in
        let r_state =
          if k < n_prep
          then AD.Maths.(t_coeff * x_t *@ c_t *@ c *@ dphi)
          else AD.Mat.zeros 1 m
        in
        let r_xp =
          if k < n_prep
          then AD.Maths.(qs_coeff * (theta_pos - theta0))
          else if k > n_mov
          then AD.Maths.(g_coeff * (theta_pos - target_pos))
          else AD.Mat.zeros 1 2
        in
        let r_xv =
          if k < n_prep
          then AD.Maths.(qs_coeff * theta_vel)
          else if k > n_mov
          then AD.Maths.(F speed_end_penalty * g_coeff * theta_vel)
          else AD.Mat.zeros 1 2
        in
        AD.Maths.(concatenate ~axis:1 [| r_xp; r_xv; r_state |])
    in
    None


  let neg_hess_t =
    let _neg_hess_t ~readout ~prms ~task =
      let c = readout in
      let m = AD.Mat.col_num c in
      let t_prep = task.t_prep in
      let dt = task.dt in
      let t_coeff = Owl_parameters.extract prms.t_coeff in
      let qs_coeff = Owl_parameters.extract prms.qs_coeff in
      let g_coeff = Owl_parameters.extract prms.g_coeff in
      let n_prep = Float.to_int (t_prep /. dt) in
      let n_mov =
        let d_mov = Float.to_int (task.t_movs.(0) /. dt) in
        n_prep + d_mov
      in
      fun ~k ~z_t ->
        let x_t = AD.Maths.get_slice [ []; [ 4; -1 ] ] z_t in
        let dx = AD.Maths.(diagm (d_phi_x x_t)) in
        let ddx = d2_phi_x x_t in
        let t1 = AD.Maths.(c *@ dx) in
        let t2 =
          AD.Maths.(diagm (diagm (phi_x x_t) *@ transpose c *@ c *@ transpose ddx))
        in
        let c_term = AD.Maths.((transpose t1 *@ t1) + t2) in
        let mu =
          if k < n_prep
          then AD.Maths.(qs_coeff * AD.Mat.eye 2)
          else if k > n_mov
          then AD.Maths.(g_coeff * AD.Mat.eye 2)
          else AD.Mat.zeros 2 2
        in
        let mv =
          if k < n_prep
          then AD.Maths.(qs_coeff * AD.Mat.eye 2)
          else if k > n_mov
          then AD.Maths.(AD.F speed_end_penalty * g_coeff * AD.Mat.eye 2)
          else AD.Mat.zeros 2 2
        in
        let mx =
          if k < n_prep
          then AD.Maths.(t_coeff * c_term)
          else if k > n_mov
          then AD.Maths.(F 0. * g_coeff * c_term)
          else AD.Mat.zeros m m
        in
        let mf1 = AD.Maths.concatenate ~axis:1 [| mu; AD.Mat.zeros 2 (m + 2) |] in
        let mf2 =
          AD.Maths.concatenate ~axis:1 [| AD.Mat.zeros 2 2; mv; AD.Mat.zeros 2 m |]
        in
        let mf3 = AD.Maths.concatenate ~axis:1 [| AD.Mat.zeros m 4; mx |] in
        AD.Maths.concatenate ~axis:0 [| mf1; mf2; mf3 |]
    in
    None
end

module Max_Occupancy (X : sig
  val label : string
  val wp : AD.t
  val phi_x : AD.t -> AD.t
  val d_phi_x : AD.t -> AD.t
  val d2_phi_x : AD.t -> AD.t
end) =
struct
  module P = Owl_parameters.Make (Max_Occupancy_P)
  open Max_Occupancy_P

  let requires_linesearch = false
  let label = X.label

  let neg_logp_t ~readout ~prms ~task =
    let n_prep = Float.to_int (task.t_prep /. task.dt) in
    let n_1 = Float.to_int ((task.t_movs.(0) +. task.t_prep) /. task.dt) in
    let pause_0 =
      match task.t_pauses with
      | Some t -> t.(0)
      | None -> 0.
    in
    let n_2 = n_1 + Float.to_int (pause_0 /. task.dt) in
    let n_3 = n_2 + Float.to_int (task.t_movs.(1) /. task.dt) in
    let qs_coeff = Owl_parameters.extract prms.qs_coeff in
    let g_coeff = Owl_parameters.extract prms.g_coeff in
    let t_coeff = Owl_parameters.extract prms.t_coeff in
    let prep_coeff = Owl_parameters.extract prms.prep_coeff in
    let c = readout in
    let c_t = AD.Maths.transpose c in
    let tgt1 =
      AD.Mat.row task.target 0 |> fun z -> AD.Maths.get_slice [ []; [ 0; 1 ] ] z
    in
    let tgt2 =
      AD.Mat.row task.target 1 |> fun z -> AD.Maths.get_slice [ []; [ 0; 1 ] ] z
    in
    let dt = task.dt in
    let _dt = AD.F dt in
    let theta0 = task.theta0 in
    let g_coeff_1 = AD.Maths.(g_coeff * F 0.2 / F pause_0) in
    fun ~k ~z_t ->
      let thetas = AD.Maths.get_slice [ []; [ 0; 3 ] ] z_t in
      let theta_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] thetas in
      let theta_vel = AD.Maths.get_slice [ []; [ 2; 3 ] ] thetas in
      let x_t = AD.Maths.get_slice [ []; [ 4; -1 ] ] z_t in
      let x_t = AD.Maths.(X.phi_x x_t) in
      if k < n_prep
      then (
        let mu_t = AD.Maths.(x_t *@ c_t) in
        AD.Maths.(
          (F 0.5 * sum' (t_coeff * sqr mu_t))
          + (F 0.5 * qs_coeff * l2norm_sqr' (thetas - theta0))))
      else if k > n_1 && k < n_2
      then (
        let proj_x_t = AD.Maths.(x_t *@ X.wp) in
        AD.Maths.(
          (F 0.5
          * g_coeff_1
          * (l2norm_sqr' (theta_pos - tgt1) + (F 0.1 * l2norm_sqr' theta_vel)))
          + (F 0.5 * prep_coeff * l2norm_sqr' proj_x_t / (l2norm_sqr' x_t + F 1E-6))))
      else if k > n_3
      then
        AD.Maths.(
          F 0.5
          * g_coeff
          * (l2norm_sqr' (theta_pos - tgt2) + (F 0.1 * l2norm_sqr' theta_vel)))
      else AD.F 0.


  let neg_jac_t = None
  let neg_hess_t = None
end


module Acquired_Phi (X : sig
  val label : string
  val phi_x : AD.t -> AD.t
  val d_phi_x : AD.t -> AD.t
  val d2_phi_x : AD.t -> AD.t
  val speed_end_penalty : float
end) =
struct
  module P = Owl_parameters.Make (Acquired_Phi_P)
  open Acquired_Phi_P
  open X

  let requires_linesearch = false
  let label = X.label

  let neg_logp_t ~readout ~prms ~task =
    let t_prep = task.t_prep in
    let target_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] task.target in
    let theta0 = task.theta0 in
    let qs_coeff = Owl_parameters.extract prms.Acquired_Phi_P.qs_coeff in
    let g_coeff = Owl_parameters.extract prms.Acquired_Phi_P.g_coeff in
    let t_coeff = Owl_parameters.extract prms.Acquired_Phi_P.t_coeff in
    let rad_thres =  Owl_parameters.extract prms.Acquired_Phi_P.rad_thres in 
    let dt = task.dt in
    let n_prep = Float.to_int (t_prep /. dt) in
    let n_mov =
      let d_mov = Float.to_int (task.t_movs.(0) /. dt) in
      n_prep + d_mov
    in
    let c = readout in
    let c_t = AD.Maths.transpose c in
    let x0 = task.x0 in
    let x0 = AD.Maths.get_slice [ []; [ 4; -1 ] ] x0 in
    let phi_x0 = phi_x x0 in
    fun ~k ~z_t ->
      let thetas = AD.Maths.get_slice [ []; [ 0; 3 ] ] z_t in
      let theta_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] thetas in
      let theta_vel = AD.Maths.get_slice [ []; [ 2; 3 ] ] thetas in
      let x_t = AD.Maths.get_slice [ []; [ 4; -1 ] ] z_t in
      let x_t = phi_x x_t in
      if k < n_prep
      then (
        let mu_t = AD.Maths.(x_t *@ c_t) in
        AD.Maths.(
          (F 0.5 * sum' (t_coeff * sqr mu_t))
          + (F 0.5 * qs_coeff * l2norm_sqr' (thetas - theta0))))
      else if k > n_mov then 
        let rad = (AD.unpack_flt rad_thres) in let curr_rad = 
          (AD.unpack_flt (AD.primal' AD.Maths.(l2norm' (theta_pos - target_pos)))) in 
       if Float.(curr_rad < rad) then 
        AD.Maths.(
          F 0.5
          * g_coeff * (F speed_end_penalty * l2norm_sqr' theta_vel)) else 
            AD.Maths.(
              F 0.5
              * g_coeff * (AD.Maths.(l2norm_sqr' (theta_pos - target_pos))) + F speed_end_penalty * l2norm_sqr' theta_vel)
      else AD.F 0.


  let neg_jac_t =
    None


  let neg_hess_t =
    None
end




module Ramping (X : sig
  val label : string
  val phi_x : AD.t -> AD.t
  val phi_t : AD.t -> AD.t
end) =
struct
  module P = Owl_parameters.Make (Ramping_P)
  open Ramping_P

  let requires_linesearch = false
  let label = X.label

  let neg_logp_t ~readout ~prms ~task =
    let t_prep = task.t_prep in
    let target_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] task.target in
    let theta0 = task.theta0 in
    let qs_coeff = Owl_parameters.extract prms.Ramping_P.qs_coeff in
    let g_coeff = Owl_parameters.extract prms.Ramping_P.g_coeff in
    let t_coeff = Owl_parameters.extract prms.Ramping_P.t_coeff in
    let tau_mov = Owl_parameters.extract prms.Ramping_P.tau_mov in
    let dt = task.dt in
    let n_prep = Float.to_int (t_prep /. dt) in
    let n_mov =
      let d_mov = Float.to_int (task.t_movs.(0) /. dt) in
      n_prep + d_mov
    in
    let c = readout in
    let c_t = AD.Maths.transpose c in
    fun ~k ~z_t ->
      let thetas = AD.Maths.get_slice [ []; [ 0; 3 ] ] z_t in
      let theta_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] thetas in
      let theta_vel = AD.Maths.get_slice [ []; [ 2; 3 ] ] thetas in
      let x_t = AD.Maths.get_slice [ []; [ 4; -1 ] ] z_t in
      let x_t = X.phi_x x_t in
      if k < n_prep
      then (
        let mu_t = AD.Maths.(x_t *@ c_t) in
        AD.Maths.(
          (F 0.5 * sum' (t_coeff * sqr mu_t))
          + (F 0.5 * qs_coeff * l2norm_sqr' (thetas - theta0))))
      else    let t_diff = AD.Maths.(AD.F (Float.(of_int Int.((k - n_prep))*dt))) in  AD.Maths.(g_coeff * X.phi_t ((t_diff / tau_mov)) * l2norm_sqr' (theta_pos - target_pos) )
        (* AD.Maths.(AD.F Float.of_int Int.(k - n_prep)
          * (l2norm_sqr' (theta_pos - target_pos) )) *)


  let neg_jac_t =
    None


  let neg_hess_t =
    None
end



module Successive_Ramping (X : sig
  val label : string
  val phi_x : AD.t -> AD.t
  val phi_t : AD.t -> AD.t
end) =
struct
  module P = Owl_parameters.Make (Successive_Ramping_P)
  open Successive_Ramping_P

  let requires_linesearch = false
  let label = X.label

  let neg_logp_t ~readout ~prms ~task =
    let n_prep = Float.to_int (task.t_prep /. task.dt) in
    let n_1 = Float.to_int ((task.t_movs.(0) +. task.t_prep) /. task.dt) in
    let pause_0 =
      match task.t_pauses with
      | Some t -> t.(0)
      | None -> 0.
    in
    let n_2 = n_1 + Float.to_int (pause_0 /. task.dt) in
    let n_3 = n_2 + Float.to_int (task.t_movs.(1) /. task.dt) in
    let qs_coeff = Owl_parameters.extract prms.qs_coeff in
    let tau_mov_1 =  Owl_parameters.extract prms.tau_mov_1 in
    let tau_mov_2 =  Owl_parameters.extract prms.tau_mov_2 in
    let g_coeff =  Owl_parameters.extract prms.g_coeff in
    let t_coeff = Owl_parameters.extract prms.t_coeff in
    let c = readout in
    let c_t = AD.Maths.transpose c in
    let tgt1 =
      AD.Mat.row task.target 0 |> fun z -> AD.Maths.get_slice [ []; [ 0; 1 ] ] z
    in
    let tgt2 =
      AD.Mat.row task.target 1 |> fun z -> AD.Maths.get_slice [ []; [ 0; 1 ] ] z
    in
    let dt = task.dt in
    let _dt = AD.F dt in
    let theta0 = task.theta0 in
    fun ~k ~z_t -> 
   let thetas = AD.Maths.get_slice [ []; [ 0; 3 ] ] z_t in
      let theta_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] thetas in
      let theta_vel = AD.Maths.get_slice [ []; [ 2; 3 ] ] thetas in
      let x_t = AD.Maths.get_slice [ []; [ 4; -1 ] ] z_t in
      let x_t = AD.Maths.(X.phi_x x_t) in
      if k < n_prep
      then 
        let mu_t = AD.Maths.(x_t *@ c_t) in
        AD.Maths.(
          (F 0.5 * sum' (t_coeff * sqr mu_t) * dt)
          + (F 0.5 * qs_coeff * dt * l2norm_sqr' (thetas - theta0)))
      else
       if k < n_2
      then 
        let t_diff_1 = AD.Maths.(AD.F (Float.(of_int Int.((k - n_1))*dt))) in 
        (AD.Maths.(
          F 0.5
          * g_coeff * dt * X.phi_t ((t_diff_1 / tau_mov_1))
          * (l2norm_sqr' (theta_pos - tgt1) )))
else 
  let t_diff_2 = AD.Maths.(AD.F (Float.(of_int Int.((k - n_2))*dt))) in 
     AD.Maths.(
          F 0.5
          * g_coeff * dt * X.phi_t ((t_diff_2 / tau_mov_2))
          * (l2norm_sqr' (theta_pos - tgt2) ))


  let neg_jac_t =
    None


  let neg_hess_t =
      None
end