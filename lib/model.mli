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
    :  ?save:(AD.t list -> unit)
    -> ?u_init:Mat.mat
    -> ?single_run:bool
    -> ?opt:bool
    -> n:int
    -> m:int
    -> x0:AD.t
    -> prms:F.p
    -> task
    -> AD.t * AD.t * AD.t

  val run : ustars:Mat.mat -> n:int -> m:int -> x0:AD.t -> prms:F.p -> task -> AD.t * AD.t

  val train
    :  ?max_iter:int
    -> ?recycle_u:bool
    -> ?save_progress_to:int * int * string
    -> ?eta:[ `constant of float | `of_iter of int -> float ]
    -> ?in_each_iteration:(prms:F.p -> int -> unit)
    -> loss:(u_init:'a option -> prms:F.p -> 'b -> AD.t * 'a)
    -> init_prms:F.p
    -> 'b array
    -> F.p
end
