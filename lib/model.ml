open Owl
open Priors
open Dynamics
open Base
open Likelihoods
include Model_typ

(*iLQR model build : takes in the likelihood and the input penalty *)

module ILQR (U : Prior_T) (D : Dynamics_T) (L : Likelihood_T) = struct
  module G = Owl_parameters.Make (Generative_P.Make (U.P) (D.P) (L.P))

  let primal' = G.map ~f:(Owl_parameters.map AD.primal')
  let linesearch = U.requires_linesearch || D.requires_linesearch || L.requires_linesearch

  (* n : dimensionality of state space; m : input dimension *)
  let solve ?u_init ?(single_run = false) ~n ~m ~prms task =
    let open Generative_P in
    let module M = struct
      type theta = G.p

      let primal' = primal'

      let n_steps =
        match task.t_hold with
        | Some th -> Float.to_int ((task.t_prep +. task.t_mov +. th) /. task.dt)
        | None -> Float.to_int ((task.t_prep +. task.t_mov) /. task.dt)


      let cost ~theta =
        let cost_lik = L.neg_logp_t ~prms:theta.likelihood ~task in
        let cost_u = U.neg_logp_t ~prms:theta.prior ~task in
        fun ~k ~x ~u ->
          let cost_lik = cost_lik ~k ~z_t:x in
          let cost_u = cost_u ~k ~x ~u in
          AD.Maths.(cost_u + cost_lik)


      let m = m
      let n = n

      let rl_u =
        Option.map U.neg_jac_t ~f:(fun neg_jac_t ~theta ->
            neg_jac_t ~prms:theta.prior ~task)


      let rl_x =
        Option.map L.neg_jac_t ~f:(fun neg_jac_t ~theta ->
            let neg_jac_t = neg_jac_t ~prms:theta.likelihood ~task in
            fun ~k ~x ~u:_ -> neg_jac_t ~k ~z_t:x)


      let rl_xx =
        Option.map L.neg_hess_t ~f:(fun neg_hess_t ~theta ->
            let neg_hess_t = neg_hess_t ~prms:theta.likelihood ~task in
            fun ~k ~x ~u:_ -> neg_hess_t ~k ~z_t:x)


      let rl_uu =
        Option.map U.neg_hess_t ~f:(fun neg_hess_t ~theta ->
            neg_hess_t ~prms:theta.prior ~task)


      let rl_ux = Some (fun ~theta:_ ~k:_ ~x:_ ~u:_ -> AD.Mat.zeros m n)
      let final_cost ~theta:_ ~k:_ ~x:_ = AD.F 0.

      let fl_x =
        let z = AD.Mat.zeros 1 n in
        Some (fun ~theta:_ ~k:_ ~x:_ -> z)


      let fl_xx =
        let z = AD.Mat.zeros n n in
        Some (fun ~theta:_ ~k:_ ~x:_ -> z)


      let dyn ~theta =
        let dyna = D.dyn ~theta:theta.dynamics ~task in
        fun ~k ~x ~u -> dyna ~k ~x ~u


      let dyn_x = Option.map D.dyn_x ~f:(fun d ~theta -> d ~theta:theta.dynamics ~task)
      let dyn_u = Option.map D.dyn_u ~f:(fun d ~theta -> d ~theta:theta.dynamics ~task)
      let running_loss = cost
      let final_loss = final_cost
    end
    in
    let module IP = Dilqr.Default.Make (M) in
    let stop_ilqr loss ~prms =
      let x0, theta =
        AD.Maths.concatenate ~axis:1 [| task.theta0; AD.Mat.zeros 1 m |], prms
      in
      let cprev = ref 1E9 in
      fun k us ->
        let c = loss ~theta x0 us in
        let pct_change = Float.(abs (c -. !cprev) /. !cprev) in
        cprev := c;
        Stdio.printf "\n loss %f || Iter %i \n%!" c k;
        if single_run then k >= 0 else Float.(pct_change < 1E-4)
    in
    let us =
      match u_init with
      | None -> List.init M.n_steps ~f:(fun _ -> AD.Mat.zeros 1 m)
      | Some us ->
        List.init M.n_steps ~f:(fun k -> AD.pack_arr (Mat.get_slice [ [ k ] ] us))
    in
    (*
          u0        u1  u2 ......   uT
          x0 = 0    x1  x2 ......   xT xT+1
      *)
    let tau =
      IP.ilqr
        ~linesearch
        ~stop:(stop_ilqr IP.loss ~prms)
        ~us
        ~x0:(AD.Maths.concatenate ~axis:1 [| task.theta0; AD.Mat.zeros 1 m |])
        ~theta:prms
        ()
    in
    let tau = AD.Maths.reshape tau [| M.n_steps + 1; -1 |] in
    ( AD.Maths.get_slice [ [ 0; -1 ]; [ 0; n - 1 ] ] tau
    , AD.Maths.get_slice [ [ 0; -1 ]; [ n; -1 ] ] tau )


  let run ~ustars ~n ~m ~prms task =
    solve ~u_init:ustars ~single_run:true ~n ~m ~prms task
end
