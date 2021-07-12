open Owl

(* no need for the nested parameter type, just basic prms *)
include module type of Model_typ
open Priors
open Dynamics
open Likelihoods

module ILQR (U : Prior_T) (D : Dynamics_T) (L : Likelihood_T) : sig
  module G : module type of Owl_parameters.Make (Generative_P.Make (U.P) (D.P) (L.P))

  val solve
    :  ?arm:bool
    -> ?u_init:Mat.mat
    -> ?single_run:bool
    -> n:int
    -> m:int
    -> prms:G.p
    -> task
    -> AD.t * AD.t

  val run : ustars:Mat.mat -> n:int -> m:int -> prms:G.p -> task -> AD.t * AD.t
end
