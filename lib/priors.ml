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
    let sl = match task.scale_lambda with
    |None -> AD.F 1.
    |Some sl -> AD.F sl in 
    fun ~k ~x:_ ~u ->
      let lam = if k < n_prep then AD.Maths.(sl*lp) else AD.Maths.(sl*lm) in
      AD.Maths.(F 0.5 * (lam * l2norm_sqr' u))


  let neg_jac_t =
    let _jac_t ~prms ~task =
      let t_prep = task.t_prep in
      let dt = task.dt in
      let n_prep = Float.to_int (t_prep /. dt) in
      let lp = Owl_parameters.extract prms.lambda_prep in
      let lm = Owl_parameters.extract prms.lambda_mov in
      let sl = match task.scale_lambda with
      |None -> AD.F 1.
      |Some sl -> AD.F sl in 
      fun ~k ~x:_ ~u ->
        let lam = if k < n_prep then AD.Maths.(sl*lp) else AD.Maths.(sl*lm) in
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
      let sl = match task.scale_lambda with
      |None -> AD.F 1.
      |Some sl -> AD.F sl in 
      fun ~k ~x:_ ~u ->
        let m = AD.Mat.col_num u in
        let lam = if k < n_prep then AD.Maths.(sl*lp) else AD.Maths.(sl*lm) in
        AD.Maths.(lam * AD.Mat.eye m)
    in
    Some _hess_t
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
    let sl = match task.scale_lambda with
    |None -> AD.F 1.
    |Some sl -> AD.F sl in 
    fun ~k ~x:_ ~u ->
      let lam = if k < n_prep then AD.Maths.(sl*lp) else AD.Maths.(sl*lm) in
      AD.Maths.(lam * (F 1. / alpha)* ((log (F 1. + exp (neg (alpha * u)))) + (log (F 1. + exp ((alpha * u)))))) |> AD.Maths.sum' 


  let neg_jac_t =
    let _jac_t ~prms ~task =
      let t_prep = task.t_prep in
      let dt = task.dt in
      let n_prep = Float.to_int (t_prep /. dt) in
      let lp = Owl_parameters.extract prms.lambda_prep in
      let lm = Owl_parameters.extract prms.lambda_mov in
      let sl = match task.scale_lambda with
      |None -> AD.F 1.
      |Some sl -> AD.F sl in 
      fun ~k ~x:_ ~u ->
        let lam = if k < n_prep then AD.Maths.(sl*lp) else AD.Maths.(sl*lm) in
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
      let sl = match task.scale_lambda with
      |None -> AD.F 1.
      |Some sl -> AD.F sl in 
      fun ~k ~x:_ ~u ->
        let lam = if k < n_prep then AD.Maths.(sl*lp) else AD.Maths.(sl*lm) in
        let m = AD.Mat.col_num u in
        AD.Maths.(lam * AD.Mat.eye m)
    in
    None
end 
(* 
module Sparse = struct
  module P = Owl_parameters.Make (Sparse_P)
  open Sparse_P

  let requires_linesearch = true

  let init ?(am = 1.) ~lambda (set : Owl_parameters.setter) =
    { lambda_prep = set (AD.F lambda); lambda_mov = set (AD.F (lambda *. am)) }


  (* returns a column vector *)

  let neg_logp_t ~prms ~task =
    let t_prep = task.t_prep in
    let dt = task.dt in
    let n_prep = Float.to_int (t_prep /. dt) in
    let lp = Owl_parameters.extract prms.lambda_prep in
    let lm = Owl_parameters.extract prms.lambda_mov in
    fun ~k ~x:_ ~u ->
      let lam = if k < n_prep then lp else lm in
      AD.Maths.((lam * (F 1E-1*l2norm_sqr' u + sum' (abs u))))


  let neg_jac_t =
    let _jac_t ~prms ~task =
      let t_prep = task.t_prep in
      let dt = task.dt in
      let n_prep = Float.to_int (t_prep /. dt) in
      let lp = Owl_parameters.extract prms.lambda_prep in
      let lm = Owl_parameters.extract prms.lambda_mov in
      fun ~k ~x:_ ~u ->
        let lam = if k < n_prep then lp else lm in
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
      fun ~k ~x:_ ~u ->
        let m = AD.Mat.col_num u in
        let lam = if k < n_prep then lp else lm in
        AD.Maths.(lam * AD.Mat.eye m)
    in
    None
end *)
