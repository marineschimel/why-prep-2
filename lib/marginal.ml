(* open Owl 
include Model_typ
open Base
module AD = Algodiff.D
open Priors
open Dynamics
open Likelihoods
open Model
open C
(* 
1) define the single sample elbo function 
2) get the covariance and mean of the approx posterior
3) draw samples from it 
4) compute logsumexp of elbo for each sample and average over samples

*)


let dt = 2E-3
let m = 200
(* module L0 = Likelihoods.Ramping (struct
  let label = "output"
  let phi_x x = x
  let phi_t t = AD.Maths.(t ** F 2.)
end)


module D0 = Dynamics.Arm_Plus (struct
  let phi_x x = x
  let d_phi_x x = AD.Maths.(F 0.5 * (F 1. + F 0.* x))
  let phi_u x = x
  let d_phi_u _ = AD.Mat.eye m
end) *)


let expand0 x = AD.Maths.reshape x (Array.append [| 1 |] (AD.shape x))

let sample_from_q ~mu ~q_cov_chol ~n_samples =
  (*cov = List of blocks of the covariance*)
  let m =  (AD.Arr.shape mu).(1) in 
  let mu = AD.Maths.reshape mu [| 1; (AD.Arr.shape mu).(0);m |] in
  let z =
    Array.map q_cov_chol ~f:(fun ell_t ->
        (* let z_t =
          AD.Linalg.linsolve ~trans:true ~typ:`l ell_t (AD.Mat.gaussian m n_samples)
        in *)
        let eta = AD.Mat.gaussian m n_samples in 
        let z_t = AD.Maths.(ell_t*@eta) in 
        let z_t = AD.Maths.transpose z_t in
        AD.Maths.reshape z_t [| 1; n_samples; m |])
    |> AD.Maths.concatenate ~axis:0
  in
  let z = AD.Maths.transpose ~axis:[| 1; 0; 2 |] z in
  let eta = AD.Arr.gaussian (AD.shape mu) in 
  AD.Maths.(mu + z)


let q_cov_chol ?(epsilon = 0.) inv_q_covs = 
    (*inv_q_covs is a list with each element being the hessian of that block of inputs*)
    inv_q_covs
    |> Array.of_list
    |> Array.map ~f:(fun x ->
            assert (Owl.Linalg.D.is_posdef (AD.unpack_arr x));
            let n = AD.Mat.row_num x in 
            let y = AD.Maths.(AD.Linalg.inv x + F epsilon * AD.Mat.eye n) in AD.Linalg.chol y |> AD.Maths.transpose)
            (* inv_q_cov = U^T U = LL^T *)
            (* AD.Linalg.chol x |> AD.Maths.transpose) *)


let log_lik_term ~prms ~dyn_prms ~readout ~task u f_int (module L0 : Likelihood_T) = 
  let u = AD.Maths.reshape u [|(AD.shape u).(1); (AD.shape u).(2)|] in 
  let xs = f_int ~prms:dyn_prms ~task ~n:0 ~u in 
  let xs = AD.unpack_arr xs in 
  let n_prep = Int.(of_float (task.t_prep /.dt))
  in 
  let le, lm = Analysis_funs.cost_x
  ~f:(fun k x ->
    let c =
      L0.neg_logp_t
        ~readout 
        ~prms
        ~task
        ~k
        ~z_t:(AD.pack_arr xs)
    in
    AD.unpack_flt c)
  ~n_prep
  xs in AD.pack_flt Float.(le +. lm)

let log_p ~prms ~task us = 
  let n_prep = Int.(of_float (task.t_prep /.dt)) in 
  let _,_,l = let f k u =
  Priors.Gaussian_Prior.neg_logp_t
    ~prms:prms.Full_P.generative.prior
    ~task
    ~k
    ~x:(AD.pack_arr u)
    ~u:(AD.pack_arr u)
in
Analysis_funs.cost_u ~f:(fun k x -> 
AD.unpack_flt (f k x)) ~n_prep us
in l |> AD.pack_flt

let logq us mu_u q_cov_chol =
  let us = AD.pack_arr us in 
    let q_cov0 = q_cov_chol.(0) in
    let m_ = AD.Mat.col_num q_cov0 in
    let t_ = Array.length q_cov_chol in
    let m = Float.of_int m_ in
    let t = Float.of_int t_ in
    let cst = Float.(m * t * log Const.pi2) in
    let log_det_term =
      Array.fold q_cov_chol ~init:(AD.F 0.) ~f:(fun accu x ->
          AD.Maths.(accu + (F 2. * sum' (log (diag x)))))
    in
    let u_s = AD.shape us in
    assert (Array.length u_s = 3);
    (*us is of size 1 x M x T*)
    let n_samples = u_s.(0) in
    let du = AD.Maths.(us - expand0 mu_u) in
    (* quadratic term: 
            assuming vec is stacking columns, du = vec(dU) and dU = is T x N
              du^t ((S^t S)⊗(T^t T))^{-1} du
          =  du^t ((S^{-1} S^{-t})⊗(T^{-1} T^{-t})) du 
          =  du^t (S^{-1}⊗T^{-1}) (S^{-t}⊗T^{-t}) du
          =  || (S^{-t}⊗T^{-t}) du ||^2 
          =  || vec(T^{-t} dU S^{-1}) ||^2 *)
    let quadratic_term =
      Array.foldi q_cov_chol ~init:(AD.F 0.) ~f:(fun t accu ell_t ->
          let du =
            AD.Maths.get_slice [ []; [ t ] ] du
            |> fun z -> AD.Maths.reshape z [| -1; m_ |]
            (*now it's n_samples x 1 x m*)
          in
          let _ = assert AD.Maths.(ell_t == AD.Maths.tril ell_t) in 
          let du_ell =  AD.Linalg.linsolve ~trans:true ~typ:`l ell_t du in
          AD.Maths.(accu + l2norm_sqr' du_ell))
    in
    AD.Maths.(
      F (-0.5) * ((F Float.(of_int n_samples) * (F cst + log_det_term)) + quadratic_term))  


let eval_ml ?(epsilon=0.) ~mu_u ~q_cov_chol ~n_samples ~prms ~task (module D0 : Dynamics_T) (module L0 : Likelihoods_T)= 
  let module I = Dynamics.Integrate (D0) in 
  let integrate = I.integrate ~readout:(Owl_parameters.extract prms.Full_P.readout.c) in 
  let u_samples = sample_from_q ~mu:mu_u ~q_cov_chol ~n_samples |> AD.unpack_arr 
 in let u_samples = Arr.split ~axis:0 (Array.init n_samples ~f:(fun _ ->1)) u_samples 
in let v = Array.foldi u_samples ~init:[] ~f:(fun i acc u -> 
    if Int.(i % C.n_nodes = C.rank) then 
    let llik = 
      log_lik_term ~prms:(prms.Full_P.generative.likelihood)  ~dyn_prms:(Owl_parameters.extract prms.Full_P.generative.dynamics) ~readout:(Owl_parameters.extract prms.Full_P.readout.c) ~task (AD.pack_arr u) integrate (module L0) in
    let lp = log_p ~prms ~task u 
in let lq = logq u mu_u q_cov_chol 
in AD.Maths.(llik + lp - lq)::acc
else acc)
|> C.gather
in 
C.broadcast' (fun () ->
      let v = v |> Array.to_list |> List.concat |> Array.of_list |> fun z -> AD.Maths.of_arrays [|z|]
in AD.Maths.(log_sum_exp' v - log (AD.F (Float.of_int n_samples))))
     *)
