open Base
open Owl_parameters
include Model_typ

module type Prior_T = sig
  module P : Owl_parameters.T
  open P

  val requires_linesearch : bool
  val neg_logp_t : prms:p -> task:task -> k:int -> x:AD.t -> u:AD.t -> AD.t
  val neg_jac_t : (prms:p -> task:task -> k:int -> x:AD.t -> u:AD.t -> AD.t) option
  val neg_hess_t : (prms:p -> task:task -> k:int -> x:AD.t -> u:AD.t -> AD.t) option
end

module Gaussian_P = struct
  type 'a prm =
    { lambda_prep : 'a
    ; lambda_mov : 'a
    }
  [@@deriving accessors ~submodule:A]

  (* first bin has the interpretation of a rescaling of the std, not the variance *)
  let map ~f x = { lambda_prep = f x.lambda_prep; lambda_mov = f x.lambda_mov }

  let fold ?prefix ~init ~f x =
    let init = f init (x.lambda_prep, with_prefix ?prefix "lambda_prep") in
    f init (x.lambda_mov, with_prefix ?prefix "lambda_mov")
end

module Sparse_P = struct
  type 'a prm =
    { lambda_prep : 'a
    ; lambda_mov : 'a
    }
  [@@deriving accessors ~submodule:A]

  (* first bin has the interpretation of a rescaling of the std, not the variance *)
  let map ~f x = { lambda_prep = f x.lambda_prep; lambda_mov = f x.lambda_mov }

  let fold ?prefix ~init ~f x =
    let init = f init (x.lambda_prep, with_prefix ?prefix "lambda_prep") in
    f init (x.lambda_mov, with_prefix ?prefix "lambda_mov")
end

module Log_barrier_P = struct
  type 'a prm =
    { lambda_prep : 'a
    ; lambda_mov : 'a
    ; a : 'a
    }
  [@@deriving accessors ~submodule:A]

  (* first bin has the interpretation of a rescaling of the std, not the variance *)
  let map ~f x = { lambda_prep = f x.lambda_prep; lambda_mov = f x.lambda_mov; a = f x.a }

  let fold ?prefix ~init ~f x =
    let init = f init (x.lambda_prep, with_prefix ?prefix "lambda_prep") in
    let init = f init (x.lambda_mov, with_prefix ?prefix "lambda_mov") in
    f init (x.a, with_prefix ?prefix "a")
end

module Gaussian_Phi_P = struct
  type 'a prm =
    { lambda_prep : 'a
    ; lambda_mov : 'a
    }
  [@@deriving accessors ~submodule:A]

  (* first bin has the interpretation of a rescaling of the std, not the variance *)
  let map ~f x = { lambda_prep = f x.lambda_prep; lambda_mov = f x.lambda_mov }

  let fold ?prefix ~init ~f x =
    let init = f init (x.lambda_prep, with_prefix ?prefix "lambda_prep") in
    f init (x.lambda_mov, with_prefix ?prefix "lambda_mov")
end

module Student_P = struct
  type 'a prm =
    { spatial_stds : 'a
    ; nu : 'a
    }
  [@@deriving accessors ~submodule:A]

  let map ~f x = { spatial_stds = f x.spatial_stds; nu = f x.nu }

  let fold ?prefix ~init ~f x =
    let init = f init (x.spatial_stds, with_prefix ?prefix "spatial_stds") in
    f init (x.nu, with_prefix ?prefix "nu")
end

module Gaussian_B_P = struct
  type 'a prm =
    { b : 'a
    ; lambda_prep : 'a
    ; lambda_mov : 'a
    }
  [@@deriving accessors ~submodule:A]

  (* first bin has the interpretation of a rescaling of the std, not the variance *)
  let map ~f x = { b = f x.b; lambda_prep = f x.lambda_prep; lambda_mov = f x.lambda_mov }

  let fold ?prefix ~init ~f x =
    let init = f init (x.lambda_prep, with_prefix ?prefix "lambda_prep") in
    let init = f init (x.b, with_prefix ?prefix "b") in
    f init (x.lambda_mov, with_prefix ?prefix "lambda_mov")
end
