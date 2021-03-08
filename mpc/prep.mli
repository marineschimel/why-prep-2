open Typ

module Make (P : Prms) (D : C_Prms) : sig
  val n : int

  val m : int

  val dyn : k:int -> x:AD.t -> u:AD.t -> AD.t

  val final_loss : Ilqr.Default.final_loss

  val running_loss : Ilqr.Default.running_loss

  val n_steps : int

  val dyn_x : Ilqr.Default.t option

  val dyn_u : Ilqr.Default.t option

  val rl_uu : t option

  val rl_xx : Ilqr.Default.t option

  val rl_ux : Ilqr.Default.t option

  val rl_u : t option

  val rl_x : Ilqr.Default.t option

  val fl_xx : Ilqr.Default.s option

  val fl_x : Ilqr.Default.s option
end
