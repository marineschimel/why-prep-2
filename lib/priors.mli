include module type of Priors_typ

module Gaussian : sig
  include Prior_T with type 'a P.prm = 'a Gaussian_P.prm

  val init : ?am:float -> lambda:float -> Owl_parameters.setter -> P.p
end

module Log_barrier : sig
  include Prior_T with type 'a P.prm = 'a Log_barrier_P.prm

  val init : ?am:float -> lambda:float -> Owl_parameters.setter -> P.p
end

module Sparse : sig
  include Prior_T with type 'a P.prm = 'a Sparse_P.prm

  val init : ?am:float -> lambda:float -> Owl_parameters.setter -> P.p
end
