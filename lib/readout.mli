include module type of Readout_typ

module Readout : sig
  module P : Owl_parameters.T with type 'a prm = 'a Readout_P.prm

  val init : Owl_parameters.setter -> int -> P.p
end
