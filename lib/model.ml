open Owl
open Priors
open Dynamics
open Base
open Likelihoods
include Model_typ
include Readout_typ

(*iLQR model build : takes in the likelihood and the input penalty *)

module ILQR (U : Prior_T) (D : Dynamics_T) (L : Likelihood_T) = struct
  module G = Owl_parameters.Make (Generative_P.Make (U.P) (D.P) (L.P))
  module F = Owl_parameters.Make (Full_P.Make (U.P) (D.P) (L.P))

  let primal' = F.map ~f:(Owl_parameters.map AD.primal')
  let linesearch = U.requires_linesearch || D.requires_linesearch || L.requires_linesearch

  (* n : dimensionality of state space; m : input dimension *)
  let solve ?(arm = true) ?u_init ?(single_run = false) ~n ~m ~prms task =
    let open Full_P in
    let module M = struct
      type theta = F.p

      let primal' = primal'

      let n_steps =
        match task.t_pauses with
        | Some tps ->
          let n_tgts = Array.length tps in
          let sum_pauses = Mat.sum' (Mat.of_arrays [| tps |]) in
          (match task.t_hold with
          | Some th ->
            Float.to_int
              ((task.t_prep +. (Float.of_int n_tgts *. task.t_mov) +. sum_pauses +. th)
              /. task.dt)
          | None -> Float.to_int ((task.t_prep +. task.t_mov) /. task.dt))
        | None ->
          (match task.t_hold with
          | Some th -> Float.to_int ((task.t_prep +. task.t_mov +. th) /. task.dt)
          | None -> Float.to_int ((task.t_prep +. task.t_mov) /. task.dt))


      let cost ~theta =
        let cost_lik =
          L.neg_logp_t
            ~readout:(Owl_parameters.extract theta.readout.c)
            ~prms:theta.generative.likelihood
            ~task
        in
        let cost_u = U.neg_logp_t ~prms:theta.generative.prior ~task in
        fun ~k ~x ~u ->
          let cost_lik = cost_lik ~k ~z_t:x in
          let cost_u = cost_u ~k ~x ~u in
          AD.Maths.(cost_u + cost_lik)


      let m = m
      let n = n

      let rl_u =
        Option.map U.neg_jac_t ~f:(fun neg_jac_t ~theta ->
            neg_jac_t ~prms:theta.generative.prior ~task)


      let rl_x =
        Option.map L.neg_jac_t ~f:(fun neg_jac_t ~theta ->
            let neg_jac_t =
              neg_jac_t
                ~readout:(Owl_parameters.extract theta.readout.c)
                ~prms:theta.generative.likelihood
                ~task
            in
            fun ~k ~x ~u:_ -> neg_jac_t ~k ~z_t:x)


      let rl_xx =
        Option.map L.neg_hess_t ~f:(fun neg_hess_t ~theta ->
            let neg_hess_t =
              neg_hess_t
                ~readout:(Owl_parameters.extract theta.readout.c)
                ~prms:theta.generative.likelihood
                ~task
            in
            fun ~k ~x ~u:_ -> neg_hess_t ~k ~z_t:x)


      let rl_uu =
        Option.map U.neg_hess_t ~f:(fun neg_hess_t ~theta ->
            neg_hess_t ~prms:theta.generative.prior ~task)


      let rl_ux = Some (fun ~theta:_ ~k:_ ~x:_ ~u:_ -> AD.Mat.zeros m n)
      let final_cost ~theta:_ ~k:_ ~x:_ = AD.F 0.

      let fl_x =
        let z = AD.Mat.zeros 1 n in
        Some (fun ~theta:_ ~k:_ ~x:_ -> z)


      let fl_xx =
        let z = AD.Mat.zeros n n in
        Some (fun ~theta:_ ~k:_ ~x:_ -> z)


      let dyn ~theta =
        let dyna =
          D.dyn
            ~readout:(Owl_parameters.extract theta.readout.c)
            ~theta:theta.generative.dynamics
            ~task
        in
        fun ~k ~x ~u -> dyna ~k ~x ~u


      let dyn_x =
        Option.map D.dyn_x ~f:(fun d ~theta ->
            d
              ~readout:(Owl_parameters.extract theta.readout.c)
              ~theta:theta.generative.dynamics
              ~task)


      let dyn_u =
        Option.map D.dyn_u ~f:(fun d ~theta ->
            d
              ~readout:(Owl_parameters.extract theta.readout.c)
              ~theta:theta.generative.dynamics
              ~task)


      let running_loss = cost
      let final_loss = final_cost
    end
    in
    let x0 =
      if arm
      then AD.Maths.concatenate ~axis:1 [| task.theta0; AD.Mat.zeros 1 m |]
      else AD.Mat.zeros 1 m
    in
    let module IP = Dilqr.Default.Make (M) in
    let stop_ilqr loss ~prms =
      let x0, theta = x0, prms in
      let cprev = ref 1E9 in
      fun k us ->
        let c = loss ~theta x0 us in
        let pct_change = Float.(abs (c -. !cprev) /. !cprev) in
        cprev := c;
        Stdio.printf "\n loss %f || pct_change %f || Iter %i \n%!" c pct_change k;
        if single_run then k >= 0 else Float.(pct_change < 1E-2)
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
      IP.ilqr ~linesearch ~stop:(stop_ilqr IP.loss ~prms) ~us ~x0 ~theta:prms ()
    in
    let tau = AD.Maths.reshape tau [| M.n_steps + 1; -1 |] in
    ( AD.Maths.get_slice [ [ 0; -1 ]; [ 0; n - 1 ] ] tau
    , AD.Maths.get_slice [ [ 0; -1 ]; [ n; -1 ] ] tau
    , IP.differentiable_loss ~theta:prms tau )


  let run ~ustars ~n ~m ~prms task =
    let a, b, _ = solve ~u_init:ustars ~single_run:true ~n ~m ~prms task in
    a, b


  let train
      ?max_iter
      ?(recycle_u = true)
      ?in_each_iteration
      ?save_progress_to
      ?eta
      ~loss
      ~init_prms
      data
    =
    let n_trials = Array.length data in
    assert (Int.(n_trials % C.n_nodes = 0));
    (* make sure all workers have the same data *)
    let data = C.broadcast data in
    (* make sure all workers have different random seeds *)
    C.self_init_rng ();
    let module Packer = Owl_parameters.Packer () in
    let handle = F.pack (module Packer) init_prms in
    let theta, lbound, ubound = Packer.finalize () in
    let theta = AD.unpack_arr theta in
    let us_init = Array.create ~len:n_trials None in
    let adam_loss theta gradient =
      Stdlib.Gc.full_major ();
      let theta = C.broadcast theta in
      let loss, g =
        Array.foldi
          data
          ~init:(0., Arr.(zeros (shape theta)))
          ~f:(fun i (accu_loss, accu_g) datai ->
            if Int.(i % C.n_nodes = C.rank)
            then (
              let open AD in
              let theta = make_reverse (Arr (Owl.Mat.copy theta)) (AD.tag ()) in
              let prms = F.unpack handle theta in
              (* we sample parameters at the outer level of the ELBOBO *)
              let u_init = us_init.(i) in
              let loss, mu_u = loss ~u_init ~prms datai in
              if recycle_u then us_init.(i) <- Some mu_u;
              reverse_prop (F 1.) loss;
              accu_loss +. unpack_flt loss, Owl.Mat.(accu_g + unpack_arr (adjval theta)))
            else accu_loss, accu_g)
      in
      let loss = Mpi.reduce_float loss Mpi.Float_sum 0 Mpi.comm_world in
      Mpi.reduce_bigarray g gradient Mpi.Sum 0 Mpi.comm_world;
      loss
    in
    let stop iter current_loss =
      Option.iter in_each_iteration ~f:(fun do_this ->
          let prms = F.unpack handle (AD.pack_arr theta) in
          let u_init = Array.map us_init ~f:(fun _ -> None) in
          do_this ~u_init ~prms iter);
      C.root_perform (fun () ->
          Stdio.printf "\r[%05i]%!" iter;
          Option.iter save_progress_to ~f:(fun (loss_every, prms_every, prefix) ->
              let kk = Int.((iter - 1) / loss_every) in
              if Int.((iter - 1) % prms_every) = 0
              then (
                let prefix = Printf.sprintf "%s_%i" prefix kk in
                let prms = F.unpack handle (AD.pack_arr theta) in
                Misc.save_bin (prefix ^ ".params.bin") prms;
                F.save_to_files ~prefix ~prms);
              if Int.((iter - 1) % loss_every) = 0
              then (
                Stdio.printf "\r[%05i] %.4f%!" iter current_loss;
                let z = [| [| Float.of_int kk; current_loss |] |] in
                Mat.(save_txt ~append:true (of_arrays z) ~out:(prefix ^ ".loss")))));
      match max_iter with
      | Some mi -> iter > mi
      | None -> false
    in
    let _ = Adam.min ?eta ?lb:lbound ?ub:ubound ~stop adam_loss theta in
    theta |> AD.pack_arr |> F.unpack handle
end
