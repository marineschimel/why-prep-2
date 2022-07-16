open Owl
module AD = Algodiff.D
open Lib
open Owl_parameters
let _ = Printexc.record_backtrace true
let n_lr = Cmdargs.(get_int "-n_lr" |> force ~usage:"-n_lr")
let n = 200
let m = n_lr
let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")
let data_dir = "/home/mmcs3/rds/rds-t2-cs156-T7o4pEA8QoU/mmcs3/random_monkeys_random_lambda_2E-6/ramping/seed_5_200_2"
let n_reaches = 8
let phi_x x = AD.Maths.relu x
let d_phi_x x = AD.Maths.(F 0.5 * (F 1. + signum x))
let d2_phi_x x = AD.Maths.(diagm (F 0. * x))
let c = AD.Mat.gaussian ~sigma:(1./.Maths.sqrt ((2.*.(Float.of_int n)))) n 2
let target = Array.init n_reaches (fun i -> Mat.load_txt (Printf.sprintf "%s/torques_%i_0" data_dir i) |> fun z -> Arr.reshape z [|1; (Arr.shape z).(0); (Arr.shape z).(1)|]) |> Arr.concatenate ~axis:0 |> AD.pack_arr

let placeholder_task x0 = Model.
  { t_prep = 0.
  ; x0
  ; t_movs = [| 0.4 |]
  ; dt = 2E-3
  ; t_hold = Some 0.2
  ; t_pauses = None
  ; scale_lambda = None
  ; target = target
  ; theta0 = x0
  ; tau = 150E-3
  }
(*Load the torques to be matched + sample a random readout*)



let n_steps = (AD.Arr.shape target).(1)
(*fit w and initial conditio,*)
module D = Dynamics.Nonlinear (struct
  let phi_x x = phi_x x
  let d_phi_x x = d_phi_x x
  let phi_u x = x
  let d_phi_u _ = AD.Mat.eye m
end)
module R = Readout

let readout = R.Readout_P.{ c = (pinned : setter) c } 
module Prms = struct
  type 'a t = { x0 : 'a; u : 'a ; v : 'a} [@@deriving prms]
end

open Prms
module P = Owl_opt_lbfgs.Make (Prms)

let integrate ~readout ~prms ~x0 ~task ~n_steps =
    let dyn_k = D.dyn ~readout ~task ~theta:prms in
    fun ~n ->
      let us =
        List.init n_steps (fun _ -> AD.Mat.zeros 1 n)
      in
      let rec dyn k x xs us=
        match us with
        | [] -> List.rev xs
        | u :: unexts ->
          let new_x = dyn_k ~k ~x ~u in
          dyn (k + 1) new_x (new_x :: xs) unexts
      in
      dyn 0 x0 [] us |> Array.of_list |> fun z -> Array.map (fun x -> let n_out = AD.Mat.col_num readout in let x = AD.Maths.(x*@readout) in AD.Arr.reshape x [|8; 1; n_out|]) z |> AD.Maths.concatenate ~axis:1


let cost prms = 
  let x0 = prms.x0 in 
  let w = AD.Maths.(prms.u*@prms.v) in
  let dyn_prms =  Dynamics.Arm_Plus_P.{ a = (learned : setter) (AD.Maths.(transpose w))
  ; b = (pinned : setter) ((AD.Mat.eye n))
  ; bias = (pinned : setter) ((AD.Mat.zeros 1 n))
  }
in 
let traj = integrate ~readout:c ~prms:dyn_prms ~x0 ~task:(placeholder_task x0) ~n ~n_steps in
AD.Maths.l2norm_sqr' AD.Maths.(traj - target) 

let losses = [] 

let final_prms =
  let prms0 = { u = AD.Mat.gaussian ~sigma:0.1 n n_lr; v = AD.Mat.gaussian ~sigma:0.1 n_lr n; x0 =  AD.Mat.gaussian ~sigma:0.1 n_reaches n} in
let stop =
    let fv_prev = ref 1E9 in
    fun s ->
      let k = P.iter s in
      let fv = P.(fv s) in
      let pct_change = abs_float ((fv -. !fv_prev) /. !fv_prev) in
      fv_prev := fv;
      [fv] |> Array.of_list |> fun z -> Mat.of_array z (-1) 1 |> fun x -> Mat.save_txt ~append:true ~out:(Printf.sprintf "%s/losses" dir) x;
      Printf.printf "\r iter %i | fv %f | pct change %f %!" k fv pct_change;
      pct_change < 1E-5 && Int.(k > 200)
  in
  let f prms =
    let f = cost prms in
    f
  in
  let s0 = P.init ~f ~prms0 () in
  let sf = P.min ~stop s0 in
  let prms = P.prms sf in
  prms

let _ = losses |> Array.of_list |> fun z -> Mat.of_array z (-1) 1 |> fun x -> Mat.save_txt ~out:(Printf.sprintf "%s/losses" dir) x


let eigenvalues m =
  let v = Linalg.D.eigvals m in
  let re = Dense.Matrix.Z.re v in
  let im = Dense.Matrix.Z.im v in
  Mat.(concat_horizontal (transpose re) (transpose im))


let _ = 
  let w = AD.Maths.(final_prms.u*@final_prms.v) in 
  let eigs = eigenvalues (AD.unpack_arr w) in 
  Mat.save_txt ~out:(Printf.sprintf "%s/u" dir) (AD.unpack_arr (final_prms.u)); 
  Mat.save_txt ~out:(Printf.sprintf "%s/v" dir) (AD.unpack_arr (final_prms.v)); 
  Mat.save_txt ~out:(Printf.sprintf "%s/eigs" dir) (eigs); 
  Mat.save_txt ~out:(Printf.sprintf "%s/x0" dir) (AD.unpack_arr (final_prms.x0));
  Mat.save_txt ~out:(Printf.sprintf "%s/c" dir) (AD.unpack_arr c) 


let w = AD.Maths.(final_prms.u*@final_prms.v)

let dyn_prms =  Dynamics.Arm_Plus_P.{ a = (learned : setter) (AD.Maths.(transpose w))
  ; b = (pinned : setter) ((AD.Mat.eye n))
  ; bias = (pinned : setter) ((AD.Mat.zeros 1 n))
  }

let torques = integrate ~readout:c ~prms:dyn_prms ~x0:final_prms.x0 ~task:(placeholder_task final_prms.x0) ~n ~n_steps

let _ = 
  let torques i = AD.Maths.get_slice [[i]] torques |> AD.unpack_arr |> fun x -> Mat.reshape x [|-1;2|]
in Array.init n_reaches (fun i -> let t = torques i in Mat.save_txt ~out:(Printf.sprintf "%s/torques_%i" dir i) t)

let traj = integrate ~readout:(AD.Mat.eye n) ~prms:dyn_prms ~x0:final_prms.x0 ~task:(placeholder_task final_prms.x0) ~n ~n_steps

let _ = 
  let traj i = AD.Maths.get_slice [[i]] traj |> AD.unpack_arr |> fun x -> Mat.reshape x [|-1;n|]
in Array.init n_reaches (fun i -> let t = traj i in Mat.save_txt ~out:(Printf.sprintf "%s/rates_%i" dir i) t) 