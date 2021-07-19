open Owl
module AD = Algodiff.D
module M = Arm.Make (Arm.Defaults)
open Lib
open Base

let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")
let data_dir = Cmdargs.(get_string "-data" |> force ~usage:"-data [dir in which data is]")
let in_dir s = Printf.sprintf "%s/%s" dir s
let in_data_dir s = Printf.sprintf "%s/%s" data_dir s
let n_targets = 21
(* let targets = C.broadcast' (fun () -> Mat.load_txt (Printf.sprintf "%s/target_thetas" data_dir))  *)
let targets = (Array.init n_targets ~f:(fun _ -> 
  [|Stats.uniform_rvs ~a:(-0.287147) ~b:0.764927; 
  Stats.uniform_rvs ~a:2.09031 ~b:2.88352;
  0.;0.|])) |> Mat.of_arrays
let _ = C.root_perform (fun ()-> Mat.save_txt ~out:(in_dir "target_thetas") targets)
let target i = Mat.row targets i
let t_prep = 0.3
let dt = 5E-3
let lambda_prep = 1E-6
let lambda_mov = 1E-6
let n_out = 2
let n = 204
let m = 200
let tau = 150E-3
let n_output = 2
let w = C.broadcast' (fun () -> Mat.(load_txt (Printf.sprintf "%s/w_rec" data_dir)))


let theta0 = Mat.of_arrays [| [| 0.174533; 2.50532; 0.; 0. |] |] |> AD.pack_arr
let t_preps = [| 0.; 0.05; 0.3|]

let c = Mat.gaussian ~sigma:0.1 n_out m
let _ = C.root_perform (fun () -> Mat.save_txt ~out:(in_data_dir "c") c)
let c = C.broadcast' (fun () -> AD.pack_arr (Mat.load_txt (Printf.sprintf "%s/c" data_dir)))


let tasks =
  Array.init
    (Array.length t_preps * n_targets)
    ~f:(fun i ->
      let n_time = i / n_targets in
      let n_target = Int.rem i n_targets in
      Model.
        { t_prep = t_preps.(n_time)
        ; t_mov = 0.4
        ; dt
        ; t_hold = Some 0.1
        ; scale_lambda = None
        ; target = AD.pack_arr (target n_target)
        ; t_pauses = None
        ; theta0
        ; tau = 150E-3
        })

let _ = C.print (Printf.sprintf "array len : %i %!" (Array.length tasks))
let save_results suffix xs us =
  let file s = Printf.sprintf "%s/%s_%s" dir s suffix in
  let thetas, xs, us =
    ( Mat.get_slice [ []; [ 0; 3 ] ] (AD.unpack_arr xs)
    , Mat.get_slice [ []; [ 4; -1 ] ] (AD.unpack_arr xs)
    , AD.unpack_arr us )
  in
  Owl.Mat.save_txt ~out:(file "thetas") thetas;
  Owl.Mat.save_txt ~out:(file "xs") xs;
  Owl.Mat.save_txt ~out:(file "us") us


let save_prms suffix prms = Misc.save_bin (Printf.sprintf "%s/prms_%s" dir suffix) prms
let save_task suffix task = Misc.save_bin (Printf.sprintf "%s/prms_%s" dir suffix) task

module U = Priors.Gaussian

module D = Dynamics.Arm_Plus (struct
  let phi_x x = x
  let d_phi_x _ = AD.Mat.eye 200
  let phi_u x = x
  let d_phi_u u = AD.Maths.(diag (F 0.5 *(signum u + F 1.)))
end)

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
              (learned : setter)
                c
          ; c_mask = None
          ; qs_coeff = (pinned : setter) (AD.F 1.)
          ; t_coeff = (pinned : setter) (AD.F 0.5)
          ; g_coeff = (pinned : setter) (AD.F 1.)
          }
      in
      let dynamics =
        Dynamics.Arm_Plus_P.
          { a =
          (pinned : setter)
          (AD.pack_arr 
           Mat.(transpose w - eye m))
          ; b = (pinned : setter) (AD.Mat.eye m)
          }
      in
      let prior =
        Priors.Gaussian_P.
          { lambda_prep = (pinned : setter) (AD.F lambda_prep)
          ; lambda_mov = (pinned : setter) (AD.F lambda_mov)
          }
      in
      let readout =
        R.Readout_P.
          { c =
              (learned : setter)
                c
          }
      in
      let generative = Model.Generative_P.{ prior; dynamics; likelihood } in
      Model.Full_P.{ generative; readout })


module I = Model.ILQR (U) (D) (L)


let loss ~u_init ~prms t =
  let _, us, l =
    match u_init with
    | Some u -> let _ = C.print "Some u_init %!" in I.solve ~u_init:u ~n ~m ~prms t
    | None -> I.solve ~n ~m ~prms t
  in
  l, AD.unpack_arr us


let final_prms =
  I.train
    ~max_iter:2000
    ~loss
    ~init_prms:prms
    ~save_progress_to:(1, 100, in_dir "progress")
    tasks



let _ =
  Mat.save_txt ~out:"c_test" (AD.unpack_arr (Owl_parameters.extract final_prms.readout.c))
