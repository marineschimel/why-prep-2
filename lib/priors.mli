include module type of Priors_typ

module Gaussian : sig
  include Prior_T with type 'a P.prm = 'a Gaussian_P.prm

  val init : ?am:float -> lambda:float -> Owl_parameters.setter -> P.p
end

module Gaussian_Double : sig
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

module Gaussian_Phi (X : sig
  val phi_u : AD.t -> AD.t
  val d_phi_u : AD.t -> AD.t
  val d2_phi_u : AD.t -> AD.t
end) : sig
  include Prior_T with type 'a P.prm = 'a Gaussian_Phi_P.prm

  val init : ?am:float -> lambda:float -> Owl_parameters.setter -> P.p
end

module Student : sig
  include Prior_T with type 'a P.prm = 'a Student_P.prm

  val init
    :  ?pin_std:bool
    -> ?spatial_std:float
    -> ?nu:float
    -> m:int
    -> Owl_parameters.setter
    -> P.p
end

module Gaussian_B : sig
  include Prior_T with type 'a P.prm = 'a Gaussian_B_P.prm

  val init : ?am:float -> lambda:float -> Owl_parameters.setter -> P.p
end

module Gaussian_Prior : sig
  include Prior_T with type 'a P.prm = 'a Gaussian_P.prm

  val init : ?am:float -> lambda:float -> Owl_parameters.setter -> P.p
end