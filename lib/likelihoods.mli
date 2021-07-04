include module type of Likelihoods_typ

module End (X : sig
  val label : string
end) : Likelihood_T with type 'a P.prm = 'a End_P.prm

module Soft_End (X : sig
  val label : string
end) : Likelihood_T with type 'a P.prm = 'a Soft_End_P.prm
