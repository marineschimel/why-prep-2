open Base
include Priors_typ
include Model_typ

module Gaussian = struct
  module P = Owl_parameters.Make (Gaussian_P)
  open Gaussian_P

  let requires_linesearch = false

  let init ?(am = 1.) ~lambda (set : Owl_parameters.setter) =
    { lambda_prep = set (AD.F lambda); lambda_mov = set (AD.F (lambda *. am)) }


  (* returns a column vector *)

  let neg_logp_t ~prms ~task =
    let t_prep = task.t_prep in
    let dt = task.dt in
    let n_prep = Float.to_int (t_prep /. dt) in
    let lp = Owl_parameters.extract prms.lambda_prep in
    let lm = Owl_parameters.extract prms.lambda_mov in
    let sl =
      match task.scale_lambda with
      | None -> AD.F 1.
      | Some sl -> AD.F sl
    in
    fun ~k ~x:_ ~u ->
      let lam = if k < n_prep then AD.Maths.(sl * lp) else AD.Maths.(sl * lm) in
      AD.Maths.(F 0.5 * (lam * l2norm_sqr' u))


  let neg_jac_t =
    let _jac_t ~prms ~task =
      let t_prep = task.t_prep in
      let dt = task.dt in
      let n_prep = Float.to_int (t_prep /. dt) in
      let lp = Owl_parameters.extract prms.lambda_prep in
      let lm = Owl_parameters.extract prms.lambda_mov in
      let sl =
        match task.scale_lambda with
        | None -> AD.F 1.
        | Some sl -> AD.F sl
      in
      fun ~k ~x:_ ~u ->
        let lam = if k < n_prep then AD.Maths.(sl * lp) else AD.Maths.(sl * lm) in
        AD.Maths.(u * lam)
    in
    Some _jac_t


  let neg_hess_t =
    let _hess_t ~prms ~task =
      let t_prep = task.t_prep in
      let dt = task.dt in
      let n_prep = Float.to_int (t_prep /. dt) in
      let lp = Owl_parameters.extract prms.lambda_prep in
      let lm = Owl_parameters.extract prms.lambda_mov in
      let sl =
        match task.scale_lambda with
        | None -> AD.F 1.
        | Some sl -> AD.F sl
      in
      fun ~k ~x:_ ~u ->
        let m = AD.Mat.col_num u in
        let lam = if k < n_prep then AD.Maths.(sl * lp) else AD.Maths.(sl * lm) in
        AD.Maths.(lam * AD.Mat.eye m)
    in
    Some _hess_t
end

module Log_barrier = struct
  module P = Owl_parameters.Make (Log_barrier_P)
  open Log_barrier_P

  let requires_linesearch = true

  let init ?(am = 1.) ~lambda (set : Owl_parameters.setter) =
    { lambda_prep = set (AD.F lambda)
    ; lambda_mov = set (AD.F (lambda *. am))
    ; a = set (AD.F 1.)
    }


  (* returns a column vector *)

  let neg_logp_t ~prms =
    let lambda = Owl_parameters.extract prms.lambda_prep in
    fun ~task:_ ~k:_ ~x:_ ~u ->
      AD.Maths.(lambda * sum' (sqr (neg (log u) + (F 1000. * u))))


  let neg_jac_t = None
  let neg_hess_t = None
end

module Sparse = struct
  module P = Owl_parameters.Make (Sparse_P)
  open Sparse_P

  let requires_linesearch = true

  let init ?(am = 1.) ~lambda (set : Owl_parameters.setter) =
    { lambda_prep = set (AD.F lambda); lambda_mov = set (AD.F (lambda *. am)) }


  let alpha = AD.F 3000.

  (* returns a column vector *)

  let neg_logp_t ~prms ~task =
    let t_prep = task.t_prep in
    let dt = task.dt in
    let n_prep = Float.to_int (t_prep /. dt) in
    let lp = Owl_parameters.extract prms.lambda_prep in
    let lm = Owl_parameters.extract prms.lambda_mov in
    let sl =
      match task.scale_lambda with
      | None -> AD.F 1.
      | Some sl -> AD.F sl
    in
    fun ~k ~x:_ ~u ->
      let lam = if k < n_prep then AD.Maths.(sl * lp) else AD.Maths.(sl * lm) in
      AD.Maths.(
        lam
        * (F 1. / alpha)
        * (log (F 1. + exp (neg (alpha * u))) + log (F 1. + exp (alpha * u))))
      |> AD.Maths.sum'


  let neg_jac_t =
    let _jac_t ~prms ~task =
      let t_prep = task.t_prep in
      let dt = task.dt in
      let n_prep = Float.to_int (t_prep /. dt) in
      let lp = Owl_parameters.extract prms.lambda_prep in
      let lm = Owl_parameters.extract prms.lambda_mov in
      let sl =
        match task.scale_lambda with
        | None -> AD.F 1.
        | Some sl -> AD.F sl
      in
      fun ~k ~x:_ ~u ->
        let lam = if k < n_prep then AD.Maths.(sl * lp) else AD.Maths.(sl * lm) in
        AD.Maths.(u * lam)
    in
    None


  let neg_hess_t =
    let _hess_t ~prms ~task =
      let t_prep = task.t_prep in
      let dt = task.dt in
      let n_prep = Float.to_int (t_prep /. dt) in
      let lp = Owl_parameters.extract prms.lambda_prep in
      let lm = Owl_parameters.extract prms.lambda_mov in
      let sl =
        match task.scale_lambda with
        | None -> AD.F 1.
        | Some sl -> AD.F sl
      in
      fun ~k ~x:_ ~u ->
        let lam = if k < n_prep then AD.Maths.(sl * lp) else AD.Maths.(sl * lm) in
        let m = AD.Mat.col_num u in
        AD.Maths.(lam * AD.Mat.eye m)
    in
    None
end

module Gaussian_Phi (X : sig
  val phi_u : AD.t -> AD.t
  val d_phi_u : AD.t -> AD.t
  val d2_phi_u : AD.t -> AD.t
end) =
struct
  module P = Owl_parameters.Make (Gaussian_Phi_P)
  open Gaussian_Phi_P
  open X

  let requires_linesearch = false

  let init ?(am = 1.) ~lambda (set : Owl_parameters.setter) =
    { lambda_prep = set (AD.F lambda); lambda_mov = set (AD.F (lambda *. am)) }


  (* returns a column vector *)

  let neg_logp_t ~prms ~task =
    let t_prep = task.t_prep in
    let dt = task.dt in
    let n_prep = Float.to_int (t_prep /. dt) in
    let lp = Owl_parameters.extract prms.lambda_prep in
    let lm = Owl_parameters.extract prms.lambda_mov in
    let sl =
      match task.scale_lambda with
      | None -> AD.F 1.
      | Some sl -> AD.F sl
    in
    fun ~k ~x:_ ~u ->
      let lam = if k < n_prep then AD.Maths.(sl * lp) else AD.Maths.(sl * lm) in
      AD.Maths.(F 0.5 * (lam * l2norm_sqr' (phi_u u)))


  let neg_jac_t =
    let _jac_t ~prms ~task =
      let t_prep = task.t_prep in
      let dt = task.dt in
      let n_prep = Float.to_int (t_prep /. dt) in
      let lp = Owl_parameters.extract prms.lambda_prep in
      let lm = Owl_parameters.extract prms.lambda_mov in
      let sl =
        match task.scale_lambda with
        | None -> AD.F 1.
        | Some sl -> AD.F sl
      in
      fun ~k ~x:_ ~u ->
        let lam = if k < n_prep then AD.Maths.(sl * lp) else AD.Maths.(sl * lm) in
        let du = AD.Maths.diagm (d_phi_u u) in
        let u = phi_u u in
        AD.Maths.(u *@ du * lam)
    in
    Some _jac_t


  let neg_hess_t =
    let _hess_t ~prms ~task =
      let t_prep = task.t_prep in
      let dt = task.dt in
      let n_prep = Float.to_int (t_prep /. dt) in
      let lp = Owl_parameters.extract prms.lambda_prep in
      let lm = Owl_parameters.extract prms.lambda_mov in
      let sl =
        match task.scale_lambda with
        | None -> AD.F 1.
        | Some sl -> AD.F sl
      in
      fun ~k ~x:_ ~u ->
        let lam = if k < n_prep then AD.Maths.(sl * lp) else AD.Maths.(sl * lm) in
        let du = AD.Maths.(diagm (d_phi_u u)) in
        let ddu = AD.Maths.(transpose (d2_phi_u u)) in
        let t1 = AD.Maths.(du *@ du) in
        let u = AD.Maths.(diagm (phi_u u)) in
        let t2 = AD.Maths.(diagm (u *@ ddu)) in
        AD.Maths.(lam * (t1 + t2))
    in
    Some _hess_t
end
