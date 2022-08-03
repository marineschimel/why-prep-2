(* open Owl 
include module type of Model_typ
open Base
module AD = Algodiff.D
open Priors
open Dynamics
open Likelihoods
open Model *)

(* module D = Owl_parameters.Make
val q_cov_chol: ?epsilon:float -> AD.t list -> AD.t array
val eval_ml :
           ?epsilon:float ->
           mu_u:AD.t ->
           q_cov_chol:AD.t array ->
           n_samples:int ->
           prms:(Priors.Gaussian_Prior.P.p, D0.P.p Owl_parameters.tag,
                 AD.t Owl_parameters.tag L0.P.prm, AD.t Owl_parameters.tag)
                Full_P.prm_ ->
           task:task -> AD.t *)
