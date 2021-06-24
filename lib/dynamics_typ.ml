open Base
open Owl_parameters
include Model_typ

module type Dynamics_T = sig
  module P : Owl_parameters.T
  open P

  val requires_linesearch : bool
  val dyn : theta:p -> task:task -> k:int -> x:AD.t -> u:AD.t -> AD.t
  val dyn_x : (theta:p -> task:task -> k:int -> x:AD.t -> u:AD.t -> AD.t) option
  val dyn_u : (theta:p -> task:task -> k:int -> x:AD.t -> u:AD.t -> AD.t) option
end

module Arm_Linear_P = struct
  type 'a prm =
    { a : 'a
    ; b : 'a
    ; c : 'a
    }

  let map ~f x = { a = f x.a; b = f x.b; c = f x.c }

  let fold ?prefix ~init ~f x =
    let init = f init (x.a, with_prefix ?prefix "a") in
    f init (x.b, with_prefix ?prefix "b")
end
