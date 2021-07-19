open Owl

(* no need for the nested parameter type, just basic prms *)
include module type of Model_typ
open Priors
open Dynamics
open Likelihoods

module ILQR (U : Prior_T) (D : Dynamics_T) (L : Likelihood_T) : sig
  module G : module type of Owl_parameters.Make (Generative_P.Make (U.P) (D.P) (L.P))
  module F : module type of Owl_parameters.Make (Full_P.Make (U.P) (D.P) (L.P))

  val solve
    :  ?arm:bool
    -> ?u_init:Mat.mat
    -> ?single_run:bool
    -> n:int
    -> m:int
    -> prms:F.p
    -> task
    -> AD.t * AD.t * AD.t

  val run : ustars:Mat.mat -> n:int -> m:int -> prms:F.p -> task -> AD.t * AD.t

  val train
    :  ?max_iter:int
    -> ?recycle_u:bool
    -> ?save_progress_to:int * int * string
    -> ?eta:[ `constant of float | `of_iter of int -> float ]
    -> loss:(u_init:'b option -> prms:F.p -> 'c -> AD.t * 'b)
    -> init_prms:F.p
    -> 'c array
    -> F.p
end
