open Base
open Owl_parameters
include Model_typ

module type Likelihood_T = sig
  module P : Owl_parameters.T
  open P

  val requires_linesearch : bool
  val label : string
  val neg_logp_t : prms:p -> task:task -> k:int -> z_t:AD.t -> AD.t
  val neg_jac_t : (prms:p -> task:task -> k:int -> z_t:AD.t -> AD.t) option
  val neg_hess_t : (prms:p -> task:task -> k:int -> z_t:AD.t -> AD.t) option
end

module End_P = struct
  type 'a prm =
    { c : 'a
    ; c_mask : AD.t option
    ; qs_coeff : 'a
    ; t_coeff : 'a
    ; g_coeff : 'a
    }

  let map ~f x =
    { c = f x.c; c_mask = x.c_mask; qs_coeff = f x.qs_coeff; t_coeff = f x.t_coeff ;
    g_coeff = f x.g_coeff }


  let fold ?prefix ~init ~f x =
    let init = f init (x.c, with_prefix ?prefix "c") in
    let init = f init (x.t_coeff, with_prefix ?prefix "t_coeff") in
    f init (x.qs_coeff, with_prefix ?prefix "qs_coeff")
end
