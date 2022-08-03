module AD = Owl.Algodiff.D
open Readout

(** Data type *)
open Owl_parameters

type task =
  { target : AD.t
  ; theta0 : AD.t
  ; t_prep : float
  ; x0 : AD.t
  ; t_movs : float array
  ; t_hold : float option
  ; t_pauses : float array option (*give an array of pause times during all the reaches*)
  ; scale_lambda : float option
  ; dt : float
  ; tau : float
  }

module Generative_P = struct
  type ('a, 'b, 'c) prm =
    { prior : 'a
    ; dynamics : 'b
    ; likelihood : 'c
    }
  [@@deriving accessors ~submodule:A]

  module Make (U : Owl_parameters.T) (D : Owl_parameters.T) (L : Owl_parameters.T) =
  struct
    type ('a, 'b, 'c) prm_ = ('a, 'b, 'c) prm
    type 'a prm = ('a U.prm, 'a D.prm, 'a L.prm) prm_

    let map ~f prms =
      { prior = U.map ~f prms.prior
      ; dynamics = D.map ~f prms.dynamics
      ; likelihood = L.map ~f prms.likelihood
      }

    let fold ?prefix ~init ~f prms =
      let w = with_prefix ?prefix in
      let init = U.fold ~prefix:(w "prior") ~init ~f prms.prior in
      let init = D.fold ~prefix:(w "dynamics") ~init ~f prms.dynamics in
      L.fold ~prefix:(w "likelihood") ~init ~f prms.likelihood
  end
end

module Full_P = struct
  type ('a, 'b, 'c, 'd) prm_ =
    { generative : ('a, 'b, 'c) Generative_P.prm
    ; readout : 'd Readout_P.prm
    }
  [@@deriving accessors ~submodule:A]

  module Make (U : Owl_parameters.T) (D : Owl_parameters.T) (L : Owl_parameters.T) =
  struct
    module G = Generative_P.Make (U) (D) (L)
    module R = Readout_P

    type 'a prm = ('a U.prm, 'a D.prm, 'a L.prm, 'a) prm_

    let map ~f prms =
      { generative = G.map ~f prms.generative; readout = R.map ~f prms.readout }

    let fold ?prefix ~init ~f prms =
      let w = with_prefix ?prefix in
      R.fold ~prefix:(w "readout") ~init ~f prms.readout
      |> fun init -> G.fold ~prefix:(w "generative") ~init ~f prms.generative
  end
end
