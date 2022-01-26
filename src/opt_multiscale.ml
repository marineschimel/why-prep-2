open Owl
module AD = Algodiff.D
open Lib
open Base
open Accessor.O

let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")
let reuse = Cmdargs.(get_string "-reuse")
let in_dir s = Printf.sprintf "%s/%s" dir s
let dt = 5E-3
let t = Float.(to_int (0.4 /. dt))
let n_samples = 32
let two_n_samples = 2 * n_samples
let rational_quadratic ~tau dt = AD.Maths.(F 1. / (F 1. + sqr (dt / F tau)))
let unpack_fun f x = AD.unpack_flt (f (AD.pack_flt x))

let g t1 t2 kernel =
  let dt = Mat.(transpose t1 - t2) in
  let gg = Mat.map (unpack_fun kernel) dt in
  let gg' = Mat.map (unpack_fun (AD.diff kernel)) dt in
  let g'g' = Mat.map (unpack_fun (AD.diff (AD.diff kernel))) dt in
  Mat.(concat_vh [| [| gg; gg' |]; [| neg gg'; neg g'g' |] |])


let norm_c0 = AD.F 1.

let samples kernel =
  let t_max = 0.5 in
  let ts = Mat.linspace 0. t_max t in
  let ts_red = Mat.of_array [| 0.; t_max |] 1 2 in
  let gyy = g ts ts kernel in
  let gxx = g ts_red ts_red kernel in
  let gxy = g ts_red ts kernel in
  let gxx_inv =
    let u, s, _ = Linalg.D.svd gxx in
    Mat.(u / (s +$ 1E-6) *@ transpose u)
  in
  let cov = Mat.(gyy - (transpose gxy *@ gxx_inv *@ gxy)) in
  let u, s, _ = Linalg.D.svd cov in
  Mat.(u * sqrt s *@ gaussian (numel s) two_n_samples)


let t1, t2 =
  C.broadcast' (fun () ->
      let torques_1 =
        samples (rational_quadratic ~tau:0.1)
        |> Mat.get_slice [ [ 0; t - 1 ]; [ 0; n_samples - 1 ] ]
      in
      let _ = Mat.save_txt ~out:(in_dir "torque_x") torques_1 in
      let torques_2 =
        samples (rational_quadratic ~tau:0.1)
        |> Mat.get_slice [ [ 0; t - 1 ]; [ n_samples; -1 ] ]
      in
      let _ = Mat.save_txt ~out:(in_dir "torque_y") torques_2 in
      ( Arr.split ~axis:1 (Array.init n_samples ~f:(fun _ -> 1)) torques_1
      , Arr.split ~axis:1 (Array.init n_samples ~f:(fun _ -> 1)) torques_2 ))


let t_prep = 0.3
let lambda_prep = 5E-3
let lambda_mov = 5E-3
let n_out = 2
let m = 50
let n = 50
let tau = 150E-3
let n_output = 2
let theta0 = Mat.of_arrays [| [| 0.174533; 2.50532; 0.; 0. |] |] |> AD.pack_arr
let c = AD.Mat.gaussian ~sigma:0.001 n_out m
let x0 = C.broadcast' (fun () -> AD.Mat.zeros 1 m)

let tasks, test_tasks =
  match reuse with
  | Some dir ->
    ( Misc.read_bin (Printf.sprintf "%s/tasks" dir)
    , Misc.read_bin (Printf.sprintf "%s/test_tasks" dir) )
  | None ->
    let tasks =
      C.broadcast' (fun () ->
          Array.init n_samples ~f:(fun i ->
              let target = Mat.(t1.(i) @|| t2.(i)) in
              let target =
                let zsz = Mat.zeros 20 2 in
                Mat.(target @= zsz)
              in
              (*I do this to make sure that with rounding etc we don't get some indices issues*)
              Model.
                { t_prep = Stats.uniform_rvs ~a:0.2 ~b:0.3
                ; t_mov = 0.4
                ; dt
                ; x0
                ; t_hold = Some 0.1
                ; scale_lambda = None
                ; target = AD.pack_arr target
                ; t_pauses = None
                ; theta0
                ; tau = 150E-3
                }))
    in
    let test_tasks =
      C.broadcast' (fun () ->
          Array.init n_samples ~f:(fun i ->
              let target = Mat.(t1.(i) @|| t2.(i)) in
              let target =
                let zsz = Mat.zeros 20 2 in
                Mat.(target @= zsz)
              in
              (*I do this to make sure that with rounding etc we don't get some indices issues*)
              Model.
                { t_prep = 0.4
                ; t_mov = 0.4
                ; dt
                ; x0
                ; t_hold = Some 0.1
                ; scale_lambda = None
                ; target = AD.pack_arr target
                ; t_pauses = None
                ; theta0
                ; tau = 150E-3
                }))
    in
    tasks, test_tasks


let _ = C.root_perform (fun () -> Misc.save_bin (in_dir "tasks") tasks)
let _ = C.root_perform (fun () -> Misc.save_bin (in_dir "test_tasks") test_tasks)
let _ = C.print (Printf.sprintf "array len : %i %!" (Array.length tasks))

module U = Priors.Gaussian
module D = Dynamics.Linear

module L = Likelihoods.Match_Torques (struct
  let label = "output"
end)

module R = Readout

let init_prms =
  let open Owl_parameters in
  C.broadcast' (fun () ->
      let likelihood =
        Likelihoods.Match_Torques_P.{ q_coeff = (pinned : setter) (AD.F 1.) }
      in
      let dynamics =
        Dynamics.Linear_P.
          { a = (learned : setter) (AD.Mat.zeros m m)
          ; b = (pinned : setter) (AD.Mat.eye m)
          }
      in
      let prior =
        Priors.Gaussian_P.
          { lambda_prep = (pinned : setter) (AD.F lambda_prep)
          ; lambda_mov = (pinned : setter) (AD.F lambda_mov)
          }
      in
      let readout = R.Readout_P.{ c = (learned : setter) c } in
      let generative = Model.Generative_P.{ prior; dynamics; likelihood } in
      Model.Full_P.{ generative; readout })


(* let init_prms = Misc.read_bin (in_dir "progress_1501.params.bin") *)

module I = Model.ILQR (U) (D) (L)

let renorm_c c =
  let open Owl_parameters in
  let c = extract c in
  let norm_c = AD.Maths.l2norm' c in
  (learned : setter) AD.Maths.(c / norm_c * norm_c0)


let loss ~u_init ~prms t =
  let prms =
    C.broadcast
      (Accessor.map (Model.Full_P.A.readout @> Readout.Readout_P.A.c) prms ~f:renorm_c)
  in
  let prms = C.broadcast prms in
  let _, us, l =
    match u_init with
    | Some _u ->
      let _ = C.print "Some u_init  : " in
      I.solve ~opt:true ~n ~m ~x0 ~prms t
    | None -> I.solve ~opt:true ~n ~m ~x0 ~prms t
  in
  let a = Owl_parameters.extract prms.generative.dynamics.a in
  AD.Maths.((l + (F 0. * l2norm_sqr' a)) / F (Float.of_int n_samples)), AD.unpack_arr us


let save_results suffix prms tasks =
  Array.iteri tasks ~f:(fun i t ->
      if Int.(i % C.n_nodes = C.rank)
      then (
        let file s = Printf.sprintf "%s.%s_%i" suffix s i in
        let xs, us, _ = I.solve ~opt:true ~x0 ~n ~m ~prms t in
        let c = Owl_parameters.extract prms.readout.c |> AD.unpack_arr in
        let xs, us = AD.unpack_arr xs, AD.unpack_arr us in
        let torques = Mat.(xs *@ transpose c) in
        Owl.Mat.save_txt ~out:(file "torques") torques;
        Owl.Mat.save_txt ~out:(file "xs") xs;
        Owl.Mat.save_txt ~out:(file "us") us))


let final_prms =
  let in_each_iteration ~prms k =
    if Int.(k % 30 = 0) then save_results (in_dir "train") prms tasks;
    if Int.(k % 30 = 0)
    then (
      let _ = C.print "testing..." in
      save_results (in_dir "test") prms test_tasks;
      if Int.(k % 10 = 0)
      then
        C.root_perform (fun () ->
            let a = Owl_parameters.extract prms.generative.dynamics.a in
            let eigs = Linalg.D.eigvals (AD.unpack_arr a) in
            let er, ei = Dense.Matrix.Z.(re eigs, im eigs) in
            Mat.(
              save_txt ~out:(in_dir (Printf.sprintf "eigs_%i" k)) (transpose (er @= ei)))))
  in
  I.train
    ~max_iter:3000
    ~loss
    ~recycle_u:true
    ~eta:(`of_iter (fun k -> Float.(0.01 / sqrt (of_int k /. 1.5))))
    ~init_prms
    ~in_each_iteration
    ~save_progress_to:(1, 500, in_dir "progress")
    tasks