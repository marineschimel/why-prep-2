open Base
open Owl_parameters
include Model_typ

module type Prior_T = sig
  module P : Owl_parameters.T

  open P

  val requires_linesearch : bool

  val neg_logp_t : prms:p -> task:task -> k:int -> x:AD.t -> u:AD.t -> AD.t

  val neg_jac_t :
    (prms:p -> task:task -> k:int -> x:AD.t -> u:AD.t -> AD.t) option

  val neg_hess_t :
    (prms:p -> task:task -> k:int -> x:AD.t -> u:AD.t -> AD.t) option
end

module Gaussian_P = struct
  type 'a prm = { lambda_prep : 'a; lambda_mov : 'a }

  (* first bin has the interpretation of a rescaling of the std, not the variance *)
  let map ~f x = { lambda_prep = f x.lambda_prep; lambda_mov = f x.lambda_mov }

  let fold ?prefix ~init ~f x =
    let init = f init (x.lambda_prep, with_prefix ?prefix "lambda_prep") in
    f init (x.lambda_mov, with_prefix ?prefix "lambda_mov")
end
