open Owl
module AD = Algodiff.D
module M = Arm.Make (Arm.Defaults)
open Lib
open Base

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


(* Setting up the parameters/directories
*)
let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")


let lambda = 0.00001

let scale_mov = Cmdargs.(get_float "-scale_mov" |> default 1.)
  
let scale_prep =  Cmdargs.(get_float "-scale_prep" |> default 1.)
(* 
let rad_w = Cmdargs.(get_float "-rad" |> force ~usage:"-rad") *)

let nonnormal =  Cmdargs.(check "-nonnormal")

(*let pc =  Cmdargs.(get_int "-pc" |> force ~usage:"-pc")*)

let skew =  Cmdargs.(check "-skew")

let results_dir = dir

let nc_angle =  Cmdargs.(get_int "-nc" |> force ~usage:"-nc")

let c_angle = Maths.(Const.pi/.8.*.Float.of_int nc_angle)
let rads = [|4.0|]
let save_all = Cmdargs.(check "-save_all")
let in_dir s = Printf.sprintf "%s/%s" dir s
let n_targets = 1

let target = Mat.of_array [| 20.; 0.|] 1 (-1)(*distance and velocity along the readout axis (1D)*)

let phi_t t = AD.Maths.(t ** F 2.)
let phi_x x = x
let d_phi_x x = AD.Maths.(F 0.5 * (F 1. + F 0. * x))
let d2_phi_x x = AD.Maths.(diagm (F 0. * x))
let link_f x = phi_x x
let dt = 2E-3
let lambda_prep = lambda *. scale_prep
let lambda_mov = lambda *. scale_mov
let n_out = 1
let _n = 4
let m = 2
let tau = 150E-3
let n_output = 1

let theta0 = Mat.of_arrays [| [| 0.; 0. |] |] |> AD.pack_arr

let x0 = AD.Mat.zeros 2 1
let t_preps = [| 0.|]

(* let w = C.broadcast' (fun () -> if nonnormal then Mat.of_arrays [|[|0.; 0.|]; [|rad_w; 0.|]|] else if skew then Mat.of_arrays [|[|0.; Float.(neg rad_w)|]; [|rad_w; 0.|]|] else 
  Mat.of_arrays [|[|rad_w; 0.|]; [|0.; rad_w|]|]) *)
(* let w = C.broadcast' (fun () -> Mat.of_arrays [|[|10.; -10.|]; [|10.; -10.|]|]) *)

(* let c = C.broadcast' (fun () ->  if Int.(pc = 1) then AD.Mat.of_arrays [|[|1.;0.|]|] else AD.Mat.of_arrays [|[|0.;1.|]|]) *)

let eigenvalues m =
  let v = Linalg.D.eigvals m in
  let re = Dense.Matrix.Z.re v in
  let im = Dense.Matrix.Z.im v in
  Mat.(concat_horizontal (transpose re) (transpose im))

(* 
let _ = C.root_perform (fun () -> Mat.save_txt ~out:(in_dir "eigs") (eigenvalues w); 
Mat.save_txt ~out:(in_dir "w") w) *)
let u0 = phi_x x0
let norm_u0 = AD.Maths.(l2norm_sqr' u0)

(* let _ =
  C.root_perform (fun () ->
      Mat.save_txt
        ~out:(Printf.sprintf "%s/nullspace" dir)
        (AD.unpack_arr AD.Maths.(c *@ u0))) *)


(* let _ =
  C.root_perform (fun () ->
      Mat.save_txt ~out:(Printf.sprintf "%s/c" dir) (AD.unpack_arr c)) *)


let b = Mat.eye m


let _ = Stdio.printf "N targets : %i %!" n_targets

(* let baseline_input =
  C.broadcast' (fun () -> 
      AD.Maths.(neg ((AD.pack_arr w *@ link_f x0) - x0)) |> AD.Maths.transpose)
 *)
let x0 =
        AD.Maths.concatenate [| x0 ; AD.Maths.transpose theta0|] ~axis:0 |> AD.Maths.transpose

        
let tasks =
  Array.init
    (n_targets * 40)
    ~f:(fun i ->
      let n_time = 0 in 
      let x0 = AD.Mat.gaussian ~sigma:Float.(of_int i/.40.) 2 1 in 
      let x0 =
              AD.Maths.concatenate [| x0 ; AD.Maths.transpose theta0|] ~axis:0 |> AD.Maths.transpose in 
      let n_target = Int.rem i n_targets in
      i, Model.
        { t_prep = t_preps.(0)
        ; x0
        ; t_movs = [| 0.4 |]
        ; dt
        ; t_hold = Some 0.2
        ; t_pauses = None
        ; scale_lambda = None
        ; target = AD.pack_arr (target)
        ; theta0
        ; tau = 150E-3
        })

let _ = Stdio.printf "Size of tasks : %i %!" (Array.length tasks)
let save_prms suffix prms = Misc.save_bin (Printf.sprintf "%s/prms_%s" dir suffix) prms
let save_task suffix task = Misc.save_bin (Printf.sprintf "%s/prms_%s" dir suffix) task
let epsilon = 1E-1

module U = Priors.Gaussian

module D0 = Dynamics.Integrator (struct
  let phi_x x = phi_x x
  let d_phi_x x = d_phi_x x
  let phi_u x = x
  let d_phi_u _ = AD.Mat.eye m
end)


module L0 = Likelihoods.Ramping_Integrator (struct
  let label = "output"
  let phi_x x = phi_x x
  let phi_t t = phi_t t
end)

let t_tot = 0.6
module R = Readout

let c = C.broadcast' (fun () -> let rot_c = Mat.of_arrays [|[|Maths.cos c_angle; Maths.(neg (sin c_angle))|]; [|Maths.(sin c_angle); Maths.(cos c_angle)|]|]
in let c0 = Mat.of_array [|1.;0.|] (-1) 1 in 
 Mat.(rot_c*@c0) |> Mat.transpose) 

let dt_scaling = Float.(dt/.1E-3)
let r = 4.0
let prms =
  let open Owl_parameters in
  C.broadcast' (fun () ->
        let w = if nonnormal then Mat.of_arrays [|[|0.; 0.|]; [|r; 0.|]|] else if skew then Mat.of_arrays [|[|0.; Float.(neg r)|]; [|r; 0.|]|] else 
          Mat.of_arrays [|[|r; 0.|]; [|0.; r|]|]
        in let c = c |> AD.pack_arr in 
        let baseline_input =
          AD.Mat.zeros 1 2 in 
      let likelihood =
        Likelihoods.Ramping_P.
          { c = (pinned : setter) c
          ; c_mask = None
          ; qs_coeff = (pinned : setter) (AD.F Float.( dt_scaling))
          ; t_coeff = (pinned : setter) (AD.F Float.(dt_scaling))
          ; g_coeff = (pinned : setter) (AD.F Float.(dt_scaling))
          ; tau_mov = (pinned : setter) (AD.F Float.(t_tot))
          }
      in
      let dynamics =
        Dynamics.Arm_Plus_P.
          { a = (pinned : setter) (AD.pack_arr Mat.(transpose w))
          ; b = (pinned : setter) (AD.pack_arr b)
          ; bias = (pinned : setter) baseline_input
          }
      in
      let prior =
        Priors.Gaussian_P.
          { lambda_prep = (pinned : setter) (AD.F Float.(lambda_prep*.dt_scaling))
          ; lambda_mov = (pinned : setter) (AD.F Float.(lambda_mov*.dt_scaling))
          }
      in
      let readout = R.Readout_P.{ c = (pinned : setter) c } in
      let generative = Model.Generative_P.{ prior; dynamics; likelihood } in
      Model.Full_P.{ generative; readout })


module I = Model.ILQR (U) (D0) (L0)

let summary_tasks =
  Array.init (Array.length t_preps) ~f:(fun _ ->
      Array.init n_targets ~f:(fun _ -> Mat.ones 1 1, false))


let get_idx t =
  let _ = Stdio.printf "ts are %f %f %!" t t_preps.(0) in
  let idx, _ = Array.findi_exn t_preps ~f:(fun _ tp -> Float.(t = tp)) in
  idx

(* let sum_mode, diff_mode = Mat.of_array [|Float.(1./.sqrt(2.)); Float.(1./.sqrt(2.))|] (-1) 1,
Mat.of_array [|Float.(1./.sqrt(2.)); Float.(-1./.sqrt(2.))|] (-1) 1 *)

let sum_mode, diff_mode = Mat.of_array [|1.;0.|] (-1) 1,
Mat.of_array [|0.;1.|] (-1) 1


let projs w xs suffix = 
  let a = Mat.(w - eye m) in 
  let obs_gramian =
    Linalg.D.lyapunov Mat.(transpose a) Mat.(neg (transpose c) *@ c)
  in let ctr_gramian = Linalg.D.lyapunov a Mat.(eye m)
in let obs_modes, _,_ = Linalg.D.svd obs_gramian
in let ctr_modes, _,_ = Linalg.D.svd ctr_gramian
in let proj_ctr =  Mat.(xs*@ctr_modes)
in let proj_obs =  Mat.(xs*@obs_modes)
in let proj_sum_diff =  Mat.((xs*@sum_mode)@||(xs*@diff_mode)) |> fun z -> Mat.save_txt (in_dir (Printf.sprintf "proj_modes_%s" suffix)) z
in Mat.((transpose obs_modes*@sum_mode)@||(transpose obs_modes*@diff_mode)) |> fun z -> Mat.save_txt (in_dir (Printf.sprintf "obs_modes_%s" suffix)) z;
Mat.((transpose ctr_modes*@sum_mode)@||(transpose ctr_modes*@diff_mode)) |> fun z -> Mat.save_txt (in_dir (Printf.sprintf "ctr_modes_%s" suffix)) z;
Mat.save_txt (in_dir (Printf.sprintf "proj_obs_%s" suffix)) proj_obs; 
Mat.save_txt (in_dir (Printf.sprintf "proj_ctr_%s" suffix)) proj_ctr


let save_results suffix xs us quus n_target n_prep task =
  let file s = Printf.sprintf "%s/%s_%s" dir s suffix in
  let _ = Stdio.printf "nprep is %i %!" n_prep in
  let xs = AD.unpack_arr xs in
  let us = AD.unpack_arr us in
  let _ = Stdio.printf "%i %i %!" (Mat.row_num xs) (Mat.col_num xs) in 
  let xs, thetas, us =
  (Mat.get_slice [ []; [ 0; -3 ] ]xs), Mat.get_slice [ []; [ -2; -1 ] ] xs, us
  in
  let x0 = Mat.get_slice [ []; [ 0; -3 ] ] (AD.unpack_arr x0) in
  let rates = AD.unpack_arr (link_f (AD.pack_arr xs)) in
  let ue_prep, ue_mov, ue_tot =
    Analysis_funs.cost_u ~f:(fun _ x -> Mat.l2norm_sqr' x) ~n_prep us
  in
  let t_prep = Float.of_int n_prep *. dt in
    Owl.Mat.save_txt ~out:(file "thetas") thetas;
    Owl.Mat.save_txt ~out:(file "xs") xs;
    Owl.Mat.save_txt ~out:(file "us") us;
    Owl.Mat.save_txt ~out:(file "rates") rates;
    Owl.Mat.save_txt ~out:(file "eff_us") (AD.unpack_arr (link_f (AD.pack_arr us)));
    Owl.Mat.save_txt
      ~out:(file "torques")
      Mat.(
        (rates *@ transpose (c)))




let () =
  let x0 = x0 in
  let _ = save_prms "" prms in
    Array.iteri tasks ~f:(fun i (_, t) ->
      let n_target = 0 in
        if Int.(i % C.n_nodes = C.rank)
        then (
          let n_prep = Float.to_int (t.t_prep /. dt) in
          let t_prep_int = Float.to_int (1000. *. t.t_prep) in
          let xs, us, l, quus,_ =
            I.solve ~single_run:true ~u_init:Mat.(gaussian ~sigma:0. 2001 m) ~n:(m + 2) ~m ~x0:t.x0 ~prms t
          in
          let prep_idx = let us_prep = AD.Maths.get_slice [[0;n_prep - 1 ]] us in 
          let us_mov = AD.Maths.get_slice [[n_prep;-1]] us
        in AD.Maths.(l2norm_sqr' us_prep / l2norm_sqr' us_mov)|> AD.unpack_flt 
      in 
          let () =
          save_results
            (Printf.sprintf "%i" i)
            xs
            us quus
            n_target
            n_prep
            t; 
            let xs = AD.Maths.get_slice [[]; [0;-3]] xs in 
            let open Model_typ in 
            let open Dynamics in 
            let w = (Owl_parameters.extract prms.Full_P.generative.dynamics.Arm_Plus_P.a) |> AD.unpack_arr in 
            projs w (AD.unpack_arr xs) (Printf.sprintf "x_%i" i); projs w (AD.unpack_arr us)  (Printf.sprintf "u_%i" i)
        in
        if save_all
        then (
          let _ = Stdio.printf "success %i %i" n_target t_prep_int in 
          Mat.save_txt
            ~out:(in_dir (Printf.sprintf "loss_%i" i))
            (Mat.of_array [| AD.unpack_flt l |] 1 (-1));
          Mat.save_txt
            ~out:(in_dir (Printf.sprintf "prep_idx_%i" i))
            (Mat.of_array [|  prep_idx |] 1 (-1));
          save_task (Printf.sprintf "%i" i) t)))
