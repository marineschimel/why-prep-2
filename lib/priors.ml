open Base
open Owl
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


  let alpha = AD.F 200.

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

module Student = struct
  module P = Owl_parameters.Make (Student_P)
  open Student_P

  let requires_linesearch = true

  let init
      ?(pin_std = false)
      ?(spatial_std = 1.)
      ?(nu = 10.)
      ~m
      (set : Owl_parameters.setter)
    =
    let spatial_stds = Owl.Mat.create 1 m spatial_std in
    { spatial_stds =
        (if pin_std then Owl_parameters.pinned else set)
          ~above:1E-3
          (AD.pack_arr spatial_stds)
    ; nu = set ~above:2.0 (AD.F nu)
    }


  let spatial_stds ~prms = Owl_parameters.extract prms.spatial_stds
  let kl_to_gaussian = `sampling_based

  let get_eff_prms ~prms =
    let s = Owl_parameters.extract prms.spatial_stds in
    let nu = Owl_parameters.extract prms.nu in
    let sigma = AD.Maths.(sqrt ((nu - F 2.) / nu) * s) in
    nu, sigma


  let neg_logp_t ~prms ~task:_ =
    let nu, sigma = get_eff_prms ~prms in
    let m = AD.Mat.numel sigma in
    let m_half = AD.F Float.(of_int m / 2.) in
    let nu_half = AD.Maths.(F 0.5 * nu) in
    let nu_plus_m_half = AD.Maths.(F 0.5 * (nu + F Float.(of_int m))) in
    let cst0 = Float.(of_int m * log Owl.Const.pi2) in
    let cst =
      let cst1 =
        AD.F
          (Owl.Maths.loggamma (AD.unpack_flt nu_half)
          -. Owl.Maths.loggamma (AD.unpack_flt nu_plus_m_half))
      in
      let cst2 = AD.Maths.(m_half * log (F Const.pi * nu)) in
      let cst3 = AD.Maths.(sum' (log sigma)) in
      AD.Maths.(cst1 + cst2 + cst3)
    in
    fun ~k:_ ~x:_ ~u ->
      let stu =
        let utilde = AD.Maths.(u / sigma) in
        AD.Maths.(
          (F 0. * cst) + (nu_plus_m_half * log (F 1. + (l2norm_sqr' utilde / nu))))
      in
      stu


  let neg_jac_t =
    let jac_t ~prms ~task:_ =
      let nu, sigma = get_eff_prms ~prms in
      let m = AD.Mat.numel sigma in
      let nu_plus_m_half = AD.Maths.(F 0.5 * (nu + F Float.(of_int m))) in
      let sigma2 = AD.Maths.sqr sigma in
      fun ~k:_ ~x:_ ~u ->
        let stu =
          let tmp =
            let utilde = AD.Maths.(u / sigma) in
            AD.Maths.(F 1. + (l2norm_sqr' utilde / nu))
          in
          let tmp' = AD.Maths.(F 2. * u / sigma2 / nu) in
          AD.Maths.(nu_plus_m_half * tmp' / tmp)
        in
        stu
    in
    Some jac_t


  let neg_hess_t =
    let hess_t ~prms ~task:_ =
      let nu, sigma = get_eff_prms ~prms in
      let m = AD.Mat.numel sigma in
      let nu_plus_m_half = AD.Maths.(F 0.5 * (nu + F Float.(of_int m))) in
      let sigma2 = AD.Maths.sqr sigma in
      fun ~k:_ ~x:_ ~u ->
        let stu =
          let u_over_s = AD.Maths.(u / sigma) in
          let tau = AD.Maths.(F 1. + (l2norm_sqr' u_over_s / nu)) in
          let cst = AD.Maths.(F 2. * nu_plus_m_half / nu / sqr tau) in
          let term1 = AD.Maths.(diagm (tau / sigma2)) in
          let term2 = AD.Maths.(F 2. * (transpose u_over_s *@ u_over_s) / nu) in
          AD.Maths.(cst * (term1 - term2))
        in
        stu
    in
    Some hess_t
end
