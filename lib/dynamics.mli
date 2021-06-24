include module type of Dynamics_typ

module AD = Owl.Algodiff.D

module Integrate (D : Dynamics_T) : sig
  val integrate : prms:D.P.p -> task:task -> n:int -> u:AD.t -> AD.t
end

module Arm_Linear : sig
  include Dynamics_T with type 'a P.prm = 'a Arm_Linear_P.prm
end
