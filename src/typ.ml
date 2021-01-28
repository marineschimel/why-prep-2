open Owl
module AD = Algodiff.D

type t = k:int -> x:AD.t -> u:AD.t -> AD.t
type s = k:int -> x:AD.t -> AD.t

module type Prms = sig
  val qs_coeff : float
  val r_coeff : float
  val t_prep : float
  val t_mov : float
  val cost : string
  val gamma_exponent : float
  val target_theta : Owl.Mat.mat array
  val aug : int
  val n : int
  val m : int
  val saving_dir : string -> string
  val duration : float
  val w : Mat.mat
  val c : Mat.mat
  val __c : AD.t
end

module type JPrms = sig
  val x : Mat.mat
end

module type C_Prms = sig
  val n_theta : int
  val q : AD.t
  val t_mat : AD.t
  val a_mat : AD.t
  val tau : AD.t
  val __dt : AD.t
  val target : AD.t
  val tgt_pos : Owl_algodiff.D.t
  val in_pos : Owl_algodiff.D.t
  val q_start : AD.t
  val cost_function : AD.t -> AD.t
  val cost : u:Owl_algodiff.D.t -> x:AD.t -> k:int -> Owl_algodiff.D.t
  val rl_u : t option
  val rl_x : t option
  val rl_uu : t option
  val rl_ux : (k:'a -> x:'b -> u:'c -> AD.t) option
  val rl_xx : (k:int -> x:AD.t -> u:AD.t -> Owl_algodiff.D.t) option
  val final_cost : x:Owl_algodiff.D.t -> k:'a -> Owl_algodiff.D.t
  val fl_x : s option
  val fl_xx : s option
end
