include module type of Dynamics_typ
module AD = Owl.Algodiff.D

module Integrate (D : Dynamics_T) : sig
  val integrate : readout:AD.t -> prms:D.P.p -> task:task -> n:int -> u:AD.t -> AD.t
end

module Arm_Linear : sig
  include Dynamics_T with type 'a P.prm = 'a Arm_Linear_P.prm
end

module Linear : sig
  include Dynamics_T with type 'a P.prm = 'a Linear_P.prm
end

module Arm_Plus (X : sig
  val phi_x : AD.t -> AD.t
  val d_phi_x : AD.t -> AD.t
  val phi_u : AD.t -> AD.t
  val d_phi_u : AD.t -> AD.t
end) : Dynamics_T with type 'a P.prm = 'a Arm_Plus_P.prm

module Arm_Discrete (X : sig
  val phi_x : AD.t -> AD.t
  val d_phi_x : AD.t -> AD.t
  val phi_u : AD.t -> AD.t
  val d_phi_u : AD.t -> AD.t
end) : Dynamics_T with type 'a P.prm = 'a Arm_Plus_P.prm


module Arm_Gated_Plus (X : sig
  val phi_x : AD.t -> AD.t
  val d_phi_x : AD.t -> AD.t
  val phi_u : AD.t -> AD.t
  val d_phi_u : AD.t -> AD.t
end) : Dynamics_T with type 'a P.prm = 'a Arm_Plus_P.prm


module Nonlinear (X : sig
  val phi_x : AD.t -> AD.t
  val d_phi_x : AD.t -> AD.t
  val phi_u : AD.t -> AD.t
  val d_phi_u : AD.t -> AD.t
end) : Dynamics_T with type 'a P.prm = 'a Arm_Plus_P.prm
