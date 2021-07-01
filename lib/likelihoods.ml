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

  let unpack_c ~prms =
    let c = Owl_parameters.extract prms.c in
    match prms.c_mask with
    | None -> c
    | Some cm -> AD.Maths.(c * cm)


  let neg_logp_t ~prms ~task =
    let t_prep = task.t_prep in
    let target_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] task.target in
    let theta0 = task.theta0 in
    let qs_coeff = Owl_parameters.extract prms.qs_coeff in
    let g_coeff = Owl_parameters.extract prms.g_coeff in
      
    let t_coeff = Owl_parameters.extract prms.t_coeff in
    let dt = task.dt in
    let n_prep = Float.to_int (t_prep /. dt) in
    let n_mov = let d_mov = Float.to_int (task.t_mov /. dt) in
    n_prep + d_mov in 
    let c = unpack_c ~prms in
    let c_t = AD.Maths.transpose c in
    fun ~k ~z_t ->
      let thetas = AD.Maths.get_slice [ []; [ 0; 3 ] ] z_t in
      let theta_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] thetas in
      let theta_vel = AD.Maths.get_slice [ []; [ 2; 3 ] ] thetas in
      if k < n_prep
      then (
        let mu_t = AD.Maths.(z_t *@ c_t) in
        AD.Maths.(
          (F 0.5 * sum' (t_coeff * sqr mu_t))
          + (F 0.5 * qs_coeff * l2norm_sqr' (thetas - theta0))))
      else if k > n_mov
      then
        AD.Maths.(
          F 0.5 * g_coeff * (l2norm_sqr' (theta_pos - target_pos) + (F 0.1 * l2norm_sqr' theta_vel)))
      else AD.F 0.


  let neg_jac_t =
    let _neg_jac_t ~prms ~task =
      let t_prep = task.t_prep in
      let target_pos = AD.Maths.get_slice [ []; [ 0; 1 ] ] task.target in
      let theta0 = AD.Maths.get_slice [ []; [ 0; 1 ] ] task.theta0 in
      let qs_coeff = Owl_parameters.extract prms.qs_coeff in
      let g_coeff = Owl_parameters.extract prms.g_coeff in
      let t_coeff = Owl_parameters.extract prms.t_coeff in
      let dt = task.dt in
      let n_prep = Float.to_int (t_prep /. dt) in
      let n_mov = let d_mov = Float.to_int (task.t_mov /. dt) in
      n_prep + d_mov in 
      let c = unpack_c ~prms in
      let c = AD.Maths.get_slice [ []; [ 4; -1 ] ] c in
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
    let _neg_hess_t ~prms ~task =
      let c = unpack_c ~prms in
      let n = AD.Mat.col_num c in
      let c = AD.Maths.get_slice [ []; [ 4; -1 ] ] c in
      let m = n - 4 in
      let t_prep = task.t_prep in
      let dt = task.dt in
      let t_coeff = Owl_parameters.extract prms.t_coeff in
      let qs_coeff = Owl_parameters.extract prms.qs_coeff in
      let g_coeff = Owl_parameters.extract prms.g_coeff in
      let n_prep = Float.to_int (t_prep /. dt) in
      let n_mov = let d_mov = Float.to_int (task.t_mov /. dt) in
      n_prep + d_mov in 
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
