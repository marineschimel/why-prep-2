open Owl
include Readout_typ
module AD = Algodiff.D

module Readout = struct
  module P = Owl_parameters.Make (Readout_P)
  open Readout_P

  let init (set : Owl_parameters.setter) n = { c = set (AD.Mat.gaussian ~sigma:0.1 2 n) }
end
