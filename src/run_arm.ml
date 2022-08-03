open Owl
module AD = Algodiff.D
module M = Arm.Make (Arm.Defaults)
open Lib
open Base

let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")
let data_dir = Cmdargs.(get_string "-data" |> force ~usage:"-data [dir in which data is]")
let in_dir s = Printf.sprintf "%s/%s" dir s
let in_data_dir s = Printf.sprintf "%s/%s" data_dir s
let t_prep = 0.3
let dt = 1E-3
let lambda_prep = 1E-6
let lambda_mov = 1.
let n_out = 2
let n = 44
let m = 40
let tau = 150E-3
let n_output = 2
let n_targets = 144
let theta0 = Mat.of_arrays [| [| 0.174533; 2.50532; 0.; 0. |] |] |> AD.pack_arr
let t_preps = [| 0.3 |]
let c = Mat.gaussian ~sigma:0.5 n_out m
let _ = C.root_perform (fun () -> Mat.save_txt ~out:(in_data_dir "c") c)

let targets =
  C.broadcast' (fun () ->
      let m1 =
        Array.init 48 ~f:(fun i ->
            let radius = 0.08 in
            let theta = Float.(of_int i * Const.pi / 48. * 2.) in
            let x = Maths.(radius *. cos theta) in
            let y = Maths.(radius *. sin theta) in
            (* let x = Stats.uniform_rvs ~a:(Maths.neg radius) ~b:radius in
          let y = Stats.uniform_rvs ~a:(Maths.neg radius) ~b:radius in *)
            let t1, t2 = Misc.pos_to_angles x (y +. 0.199112) in
            [| t1; t2; 0.; 0. |])
        |> Mat.of_arrays
      in
      let m2 =
        Array.init 48 ~f:(fun i ->
            let radius = 0.12 in
            let theta = Float.(of_int i * Const.pi / 48. * 2.) in
            let x = Maths.(radius *. cos theta) in
            let y = Maths.(radius *. sin theta) in
            (* let x = Stats.uniform_rvs ~a:(Maths.neg radius) ~b:radius in
          let y = Stats.uniform_rvs ~a:(Maths.neg radius) ~b:radius in *)
            let t1, t2 = Misc.pos_to_angles x (y +. 0.199112) in
            [| t1; t2; 0.; 0. |])
        |> Mat.of_arrays
      in
      let m3 =
        Array.init 48 ~f:(fun i ->
            let radius = 0.15 in
            let theta = Float.(of_int i * Const.pi / 48. * 2.) in
            let x = Maths.(radius *. cos theta) in
            let y = Maths.(radius *. sin theta) in
            (* let x = Stats.uniform_rvs ~a:(Maths.neg radius) ~b:radius in
      let y = Stats.uniform_rvs ~a:(Maths.neg radius) ~b:radius in *)
            let t1, t2 = Misc.pos_to_angles x (y +. 0.199112) in
            [| t1; t2; 0.; 0. |])
        |> Mat.of_arrays
      in
      Mat.(m1 @= m2 @= m3))

let _ = C.root_perform (fun () -> Mat.save_txt ~out:(in_dir "target_thetas") targets)
let target i = Mat.row targets i

let tasks =
  Array.init
    (Array.length t_preps * n_targets)
    ~f:(fun i ->
      let n_time = i / n_targets in
      let n_target = Int.rem i n_targets in
      Model.
        { t_prep = t_preps.(n_time)
        ; t_mov = 0.3
        ; dt
        ; t_hold = Some 0.2
        ; t_pauses = None
        ; scale_lambda = None
        ; target = AD.pack_arr (target n_target)
        ; theta0
        ; tau = 150E-3
        })

let save_results suffix xs us =
  let file s = Printf.sprintf "%s/%s_%s" dir s suffix in
  let thetas, xs, us =
    ( Mat.get_slice [ []; [ 0; 3 ] ] (AD.unpack_arr xs)
    , Mat.get_slice [ []; [ 4; -1 ] ] (AD.unpack_arr xs)
    , AD.unpack_arr us )
  in
  let hands =
    let h =
      AD.Mat.map_by_row
        (fun x -> Arm.unpack_state_diff (M.hand_of (Arm.pack_state x)))
        (AD.pack_arr thetas)
    in
    AD.unpack_arr h
  in
  let subsamp x =
    let x = Mat.get_slice [ [ 0; -200 ]; [] ] x in
    let n_steps = Mat.row_num x in
    let bin = 10 in
    let y =
      Array.init (n_steps / 10) ~f:(fun i ->
          let j = i * bin in
          Mat.row x j)
    in
    Mat.concatenate ~axis:0 y
  in
  Owl.Mat.save_txt ~out:(file "thetas") thetas;
  Owl.Mat.save_txt ~out:(file "hands") hands;
  Owl.Mat.save_txt ~out:(file "xs") xs;
  Owl.Mat.save_txt ~out:(file "us") us;
  Owl.Mat.save_txt ~out:(file "subsamp_us") (subsamp us);
  Owl.Mat.save_txt ~out:(file "subsamp_xs") (subsamp xs);
  Owl.Mat.save_txt ~out:(file "subsamp_hands") (subsamp hands)

let save_prms suffix prms = Misc.save_bin (Printf.sprintf "%s/prms_%s" dir suffix) prms
let save_task suffix task = Misc.save_bin (Printf.sprintf "%s/prms_%s" dir suffix) task

module U = Priors.Gaussian
module D = Dynamics.Arm_Linear
module R = Readout

module L = Likelihoods.End (struct
  let label = "output"
end)

let prms =
  let open Owl_parameters in
  C.broadcast' (fun () ->
      let likelihood =
        Likelihoods.End_P.
          { c =
              (pinned : setter)
                (AD.pack_arr (Mat.load_txt (Printf.sprintf "%s/c" data_dir)))
          ; c_mask = None
          ; qs_coeff = (pinned : setter) (AD.F 1.)
          ; t_coeff = (pinned : setter) (AD.F 0.5)
          ; g_coeff = (pinned : setter) (AD.F 1.)
          }
      in
      let dynamics =
        Dynamics.Arm_Linear_P.
          { a =
              (pinned : setter)
                (AD.pack_arr
                   (Mat.transpose
                      Mat.(load_txt (Printf.sprintf "%s/w_rec_40" data_dir) - eye m)))
          ; b = (pinned : setter) (AD.Mat.eye m)
          }
      in
      let prior =
        Priors.Gaussian_P.
          { lambda_prep = (pinned : setter) (AD.F lambda_prep)
          ; lambda_mov = (pinned : setter) (AD.F lambda_mov)
          }
      in
      let readout = R.Readout_P.{ c = (pinned : setter) (AD.pack_arr c) } in
      let generative = Model.Generative_P.{ prior; dynamics; likelihood } in
      Model.Full_P.{ generative; readout })

module I = Model.ILQR (U) (D) (L)

let _ =
  let _ = save_prms "" prms in
  let x0 = AD.Maths.concatenate ~axis:1 [| theta0; AD.Maths.(F 0. * AD.Mat.ones 1 m) |] in
  Array.mapi tasks ~f:(fun i t ->
      if Int.(i % C.n_nodes = C.rank)
      then (
        let n_target = Int.rem i n_targets in
        let t_prep = Float.to_int (1000. *. t.t_prep) in
        let xs, us, _ = I.solve ~x0 ~n ~m ~prms t in
        save_results (Printf.sprintf "%i_%i" n_target t_prep) xs us;
        save_task (Printf.sprintf "%i_%i" n_target t_prep) t))
