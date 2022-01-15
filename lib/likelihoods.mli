include module type of Likelihoods_typ

module End (X : sig
  val label : string
end) : Likelihood_T with type 'a P.prm = 'a End_P.prm

module Soft_End (X : sig
  val label : string
end) : Likelihood_T with type 'a P.prm = 'a Soft_End_P.prm

module Match (X : sig
  val label : string
end) : Likelihood_T with type 'a P.prm = 'a Match_P.prm

module Match_Torques (X : sig
  val label : string
end) : Likelihood_T with type 'a P.prm = 'a Match_Torques_P.prm

module Successive (X : sig
  val label : string
  val phi_x : AD.t -> AD.t
  val d_phi_x : AD.t -> AD.t
  val d2_phi_x : AD.t -> AD.t
end) : Likelihood_T with type 'a P.prm = 'a Successive_P.prm

module End_Phi (X : sig
  val label : string
  val phi_x : AD.t -> AD.t
  val d_phi_x : AD.t -> AD.t
  val d2_phi_x : AD.t -> AD.t
  val speed_end_penalty : float
end) : Likelihood_T with type 'a P.prm = 'a End_Phi_P.prm
