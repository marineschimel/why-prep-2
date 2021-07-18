open Owl_parameters

module Readout_P = struct
  type 'a prm = { c : 'a } [@@deriving accessors ~submodule:A]

  let map ~f x = { c = f x.c }
  let fold ?prefix ~init ~f x = f init (x.c, with_prefix ?prefix "readout")
end
