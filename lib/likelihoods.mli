include module type of Likelihoods_typ

module End (X : sig
  val label : string
end) : Likelihood_T with type 'a P.prm = 'a End_P.prm