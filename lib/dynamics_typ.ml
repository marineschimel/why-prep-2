open Base
open Owl_parameters
include Model_typ

module type Dynamics_T = sig
  module P : Owl_parameters.T

  val requires_linesearch : bool
  val dyn : readout:AD.t -> theta:P.p -> task:task -> k:int -> x:AD.t -> u:AD.t -> AD.t

  val dyn_x
    : (readout:AD.t -> theta:P.p -> task:task -> k:int -> x:AD.t -> u:AD.t -> AD.t) option

  val dyn_u
    : (readout:AD.t -> theta:P.p -> task:task -> k:int -> x:AD.t -> u:AD.t -> AD.t) option
end

module Arm_Linear_P = struct
  type 'a prm =
    { a : 'a
    ; b : 'a
    ; c : 'a
    }
  [@@deriving accessors ~submodule:A]

  let map ~f x = { a = f x.a; b = f x.b; c = f x.c }

  let fold ?prefix ~init ~f x =
    let init = f init (x.a, with_prefix ?prefix "a") in
    f init (x.b, with_prefix ?prefix "b")
end

module Linear_P = struct
  type 'a prm =
    { a : 'a
    ; b : 'a
    }
  [@@deriving accessors ~submodule:A]

  let map ~f x = { a = f x.a; b = f x.b }

  let fold ?prefix ~init ~f x =
    let init = f init (x.a, with_prefix ?prefix "a") in
    f init (x.b, with_prefix ?prefix "b")
end

module Arm_Plus_P = struct
  type 'a prm =
    { a : 'a
    ; b : 'a
    }
  [@@deriving accessors ~submodule:A]

  let map ~f x = { a = f x.a; b = f x.b }

  let fold ?prefix ~init ~f x =
    let init = f init (x.a, with_prefix ?prefix "a") in
    f init (x.b, with_prefix ?prefix "b")
end
