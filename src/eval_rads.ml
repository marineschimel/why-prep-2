open Owl
module AD = Algodiff.D
module M = Arm.Make (Arm.Defaults)
open Lib
open Base

(* Setting up the parameters/directories
*)
let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")
let data_dir = Cmdargs.(get_string "-data" |> force ~usage:"-data [dir in which data is]")

let lambda =
  Cmdargs.(get_float "-lambda" |> force ~usage:"-lambda [dir in which data is]")

let scale_mov = Cmdargs.(get_float "-scale_mov" |> default 1.)
let scale_prep = Cmdargs.(get_float "-scale_prep" |> default 1.)

let results_dir =
  Cmdargs.(get_string "-rdir" |> force ~usage:"-rdir [where to save the key results]")

let save_all = Cmdargs.(check "-save_all")
let soc = Cmdargs.(check "-soc")
let sim_soc = Cmdargs.(check "-sim_soc")
let skew = Cmdargs.(check "-skew")
let rdn = Cmdargs.(check "-rdn")
let discrete = Cmdargs.(check "-discrete")
let rad_c = Cmdargs.(get_float "-rad_c" |> default 0.5)
let rad_w = Cmdargs.(get_float "-rad_w" |> default 0.5)
let tau_mov = Cmdargs.(get_float "-tau_mov" |> default 600.)
let t_coeff = Cmdargs.(get_float "-t_coeff" |> default 1.)
let sa = Cmdargs.(get_float "-sa" |> default 0.8)
let exponent = Cmdargs.(get_int "-exponent" |> force ~usage:"exponent")
let n_lr = Cmdargs.(get_int "-n_lr" |> default 1)
let seed = Cmdargs.(get_int "-seed" |> default 1)
let in_dir s = Printf.sprintf "%s/%s" dir s
let in_data_dir s = Printf.sprintf "%s/%s" data_dir s
let n_targets = 1
let n_target = Cmdargs.(get_int "-n_target" |> default 1)
let net_name = if soc then "soc" else if skew then "skew" else "rdn"

let targets =
  C.broadcast' (fun () ->
      Array.init n_targets ~f:(fun _ ->
          let i = n_target in
          let radius = 0.12 in
          let theta = Float.(of_int i *. Const.pi *. 2. /. of_int n_targets) in
          let x = Maths.(radius *. cos theta) in
          let y = Maths.(radius *. sin theta) in
          (* let x = Stats.uniform_rvs ~a:(Maths.neg radius) ~b:radius in
      let y = Stats.uniform_rvs ~a:(Maths.neg radius) ~b:radius in *)
          let t1, t2 = Misc.pos_to_angles x (y +. 0.199112) in
          Mat.of_array [| t1; t2; 0.; 0. |] 1 (-1)))

let hand_targets =
  C.broadcast' (fun () ->
      Array.init n_targets ~f:(fun _ ->
          let i = n_target in
          let radius = 0.12 in
          let theta = Float.(of_int i *. Const.pi *. 2. /. of_int n_targets) in
          let x = Maths.(radius *. cos theta) in
          let y = Maths.(radius *. sin theta) in
          Mat.of_array [| x; y |] 1 (-1)))

let _ =
  C.root_perform (fun () ->
      Mat.save_txt ~out:(in_dir "targets") (Mat.concatenate ~axis:0 targets);
      Mat.save_txt ~out:(in_dir "hand_targets") (Mat.concatenate ~axis:0 hand_targets))

let beta = AD.F 1E-2
let exponent = AD.F (Float.of_int exponent)
let phi_t t = AD.Maths.(t ** exponent)
let phi_x x = AD.Maths.relu x
let d_phi_x x = AD.Maths.(F 0.5 * (F 1. + signum x))
let d2_phi_x x = AD.Maths.(diagm (F 0. * x))

(* let phi_x x = AD.Maths.(F 5. * tanh x)
let d_phi_x x = AD.Maths.(F 5. * (F 1. - F 2. * (sqr (tanh x))))
let d2_phi_x x = AD.Maths.(F (-10.) * tanh x * (d_phi_x x)) *)
(* let phi_x x = AD.Maths.(beta * AD.requad (x / beta))
let d_phi_x x = AD.Maths.(AD.d_requad (x / beta))
let d2_phi_x x = AD.Maths.(F 1. / beta * AD.d2_requad (x / beta)) *)
let link_f x = phi_x x
let target i = targets.(i)
let dt = 2E-3
let lambda_prep = lambda *. scale_prep
let lambda_mov = lambda *. scale_mov
let n_out = 2
let _n = 204
let m = 200
let tau = 150E-3
let n_output = 2

let _ =
  Mat.save_txt
    ~out:(in_dir "prms")
    (Mat.of_array
       [| tau; lambda_prep; lambda_mov; dt; AD.unpack_flt beta; rad_c |]
       1
       (-1))

let theta0 = Mat.of_arrays [| [| 0.174533; 2.50532; 0.; 0. |] |] |> AD.pack_arr
let t_preps = [| 0.3 |]

(* let radii = [|0.;0.05;0.1;0.15;0.2;0.25;0.3;0.35;0.4;0.45;0.5;0.55;0.6;0.65;0.7;0.75;0.8;0.85;0.9;0.95;1.0;1.05|] *)
let radii =
  if soc
  then
    [| 0.3
     ; 0.5
     ; 0.7
     ; 0.8
     ; 1.0
     ; 1.1
     ; 1.3
     ; 1.4
     ; 1.5
     ; 1.6
     ; 1.7
     ; 1.8
     ; 1.9
     ; 2.1
     ; 2.2
     ; 2.7
     ; 3.2
     ; 3.7
     ; 1.0
     ; 1.5
     ; 2.0
     ; 2.5
     ; 3.0
     ; 3.5
     ; 4.0
     ; 4.5
     ; 5.0
     ; 5.5
     ; 6.0
     ; 6.5
     ; 10.0
     ; 1.3
     ; 3.1
     ; 3.4
     ; 3.6
     ; 3.8
     ; 4.1
     ; 4.2
     ; 4.3
     ; 4.6
     ; 4.9
     ; 5.2
     ; 5.3
     ; 5.4
     ; 5.7
     ; 5.9
     ; 6.1
     ; 6.3
     ; 6.7
     ; 6.9
    |]
    (* [|0.3; 0.5; 0.7; 0.8; 1.0; 1.1; 1.3; 1.4; 1.6; 1.7; 1.8; 1.9; 2.1|] *)
    (*   
[|0.2; 0.4; 0.6; 0.9; 1.3; 2.2; 2.7; 3.2; 3.7; 1.0; 1.5; 2.0; 2.5; 3.0; 3.5; 4.0; 4.5; 5.0; 5.5; 6.0; 6.5; 10.0; 1.3; 3.1; 3.4; 3.6; 3.8; 4.1; 4.2;4.3; 4.6; 4.9; 5.2; 5.3; 5.4; 5.7; 5.9; 6.1; 6.3; 6.7; 6.9|] *)
    (*4.0; 4.5; 5.0; 5.5; 6.0; 6.5; 7.0; 10.0|]*)
  else if skew
  then
    [| 0.3
     ; 0.5
     ; 0.7
     ; 0.8
     ; 1.0
     ; 1.1
     ; 1.5
     ; 1.7
     ; 1.8
     ; 2.0
     ; 2.3
     ; 2.5
     ; 2.9
     ; 3.1
     ; 3.3
     ; 3.5
     ; 3.8
     ; 4.1
     ; 4.3
     ; 4.7
     ; 4.9
     ; 5.1
     ; 5.3
     ; 5.4
     ; 5.5
     ; 6.0
     ; 6.2
     ; 6.5
     ; 6.7
     ; 7.0
     ; 7.3
     ; 7.5
     ; 7.7
     ; 8.0
     ; 8.2
     ; 8.5
     ; 8.7
     ; 8.9
     ; 9.1
     ; 9.3
     ; 9.5
     ; 9.7
     ; 10.0
     ; 10.5
     ; 10.7
     ; 11.0
     ; 11.5
     ; 12.0
     ; 12.5
    |]
    (* ; 13.0
    ; 13.5
    ; 14.0
    ; 14.5
    ; 15.0 *)
    (* [|0.1; 0.2; 0.4; 0.6; 0.9; 1.3; 2.2; 2.7; 3.2; 3.7; 0.5; 0.8; 1.0; 1.5; 2.0; 2.5; 3.0; 3.5; 4.0; 4.5; 5.0; 5.5; 6.0; 6.2;  6.5;6.7; 7.0; 7.3; 7.5; 7.7; 8.0; 8.2; 8.5; 8.7; 8.9; 9.1; 9.3; 9.5; 9.7; 10.0; 10.5; 10.7; 11.0; 11.5; 12.0; 12.5; 13.0; 13.5;  14.0; 14.5; 15.0|]  *)
  else [| 0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9; 1.0; 1.1 |]

let c =
  C.broadcast' (fun () -> AD.Mat.gaussian ~sigma:Float.(rad_c / sqrt (of_int m)) 2 m)

let x0 = C.broadcast' (fun () -> AD.Maths.(F 0.5 * AD.Mat.uniform ~a:5. ~b:15. m 1))

let eigenvalues m =
  let v = Linalg.D.eigvals m in
  let re = Dense.Matrix.Z.re v in
  let im = Dense.Matrix.Z.im v in
  Mat.(concat_horizontal (transpose re) (transpose im)), re, im

let u0 = phi_x x0
let norm_u0 = AD.Maths.(l2norm_sqr' u0)
let c = AD.Maths.(c - (c *@ u0 *@ transpose u0 / norm_u0))
let b = Mat.eye m

(*structure : get inputs from task, then from the same SOC just generating inputs, 
then from yet another SOC
Maybe do the same thing using 200 - 100 - 50 
what about one population receiving modulated inputs about the target? 
*)
(* x(t+1)- x(t) = Wx(t) + baseline => 0 when x = -W^(-1)*baseline  *)
(* let baseline_input = AD.Mat.ones 1 m *)

(* let x0 =
  let winv = AD.Linalg.linsolve (AD.pack_arr w) (AD.Mat.eye m) in
  AD.Maths.(neg (winv *@ transpose baseline_input)) |> AD.Maths.transpose *)

(* let x0 =
  C.broadcast' (fun () ->
      AD.Maths.transpose (AD.pack_arr (Mat.load_txt (Printf.sprintf "%s/x0" dir)))) *)

let _ = Stdio.printf "N targets : %i %!" n_targets

module U = Priors.Gaussian

module D0 =
(* if discrete then Dynamics.Arm_Discrete (struct
  let phi_x x = phi_x x
  let d_phi_x x = d_phi_x x
  let phi_u x = x
  let d_phi_u _ = AD.Mat.eye m
end) else  *)
Dynamics.Arm_Plus (struct
  let phi_x x = phi_x x
  let d_phi_x x = d_phi_x x
  let phi_u x = x
  let d_phi_u _ = AD.Mat.eye m
end)

module L0 = Likelihoods.Ramping (struct
  let label = "output"
  let phi_x x = phi_x x
  let phi_t t = phi_t t
end)

let t_tot = 0.6

module R = Readout

let dt_scaling = Float.(dt /. 1E-3 *. 1000.)

let make_prms radius =
  let w =
    if soc
    then
      Mat.(
        load_txt (Printf.sprintf "%s/socs/w_rad_%.1f_sa_%.2f_%i" data_dir radius sa seed))
    else (
      let m = Mat.gaussian ~sigma:Float.(radius /. sqrt (of_int m)) m m in
      if skew then Mat.(((m - transpose m) /$ 2.) + (sa $* Mat.(eye 200))) else m)
  in
  (* let baseline_input =
    let a_discrete =
      Linalg.D.(expm Mat.(transpose w *$ Float.(dt /. tau))) |> AD.pack_arr
    in
    let x0 = AD.Maths.transpose x0 in
    AD.Maths.(neg (link_f x0 *@ a_discrete) + x0) *)
  let baseline_input =
    AD.Maths.(neg ((AD.pack_arr w *@ link_f x0) - x0)) |> AD.Maths.transpose
  in
  let open Owl_parameters in
  let likelihood =
    Likelihoods.Ramping_P.
      { c = (pinned : setter) c
      ; c_mask = None
      ; qs_coeff = (pinned : setter) (AD.F Float.(dt_scaling))
      ; t_coeff = (pinned : setter) (AD.F Float.(t_coeff *. dt_scaling))
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
      { lambda_prep = (pinned : setter) (AD.F Float.(lambda_prep *. dt_scaling))
      ; lambda_mov = (pinned : setter) (AD.F Float.(lambda_mov *. dt_scaling))
      }
  in
  let readout = R.Readout_P.{ c = (pinned : setter) c } in
  let generative = Model.Generative_P.{ prior; dynamics; likelihood } in
  Model.Full_P.{ generative; readout }, w

let x0 =
  AD.Maths.concatenate [| AD.Maths.transpose theta0; x0 |] ~axis:0 |> AD.Maths.transpose

let tasks =
  C.broadcast' (fun () ->
      let tasks = [] in
      let rec add_task i t =
        if Int.(i < Array.length radii)
        then (
          let n_rad = i in
          let n_target = n_target in
          let new_t =
            try
              let prms, w = make_prms radii.(n_rad) in
              ( n_target
              , Model.
                  { t_prep = t_preps.(0)
                  ; x0
                  ; t_movs = [| 0.4 |]
                  ; dt
                  ; t_hold = Some 0.2
                  ; t_pauses = None
                  ; scale_lambda = None
                  ; target = AD.pack_arr (target 0)
                  ; theta0
                  ; tau = 150E-3
                  }
              , radii.(n_rad)
              , w
              , prms )
              :: t
            with
            | e ->
              Stdio.printf "%s failed rad %f %!" (Exn.to_string e) radii.(n_rad);
              t
          in
          add_task Int.(i + 1) new_t)
        else t
      in
      let all_tasks = add_task 0 tasks in
      Array.of_list all_tasks)

let _ = Stdio.printf "Size of tasks : %i %!" (Array.length tasks)
let epsilon = 1E-1

module I = Model.ILQR (U) (D0) (L0)

let summary_tasks =
  Array.init (Array.length t_preps) ~f:(fun _ ->
      Array.init n_targets ~f:(fun _ -> Mat.zeros 1 1, false))

let get_idx t =
  let _ = Stdio.printf "ts are %f %f %!" t t_preps.(0) in
  let idx, _ = Array.findi_exn t_preps ~f:(fun _ tp -> Float.(t = tp)) in
  idx

let () =
  let x0 = x0 in
  Array.iteri tasks ~f:(fun i (_, t, rad, w, prms) ->
      if Int.(i % C.n_nodes = C.rank)
      then (
        try
          let norm_w = Mat.l2norm' w in
          let _ = Stdio.printf "%i %f %!" n_target rad in
          let _ =
            if Int.(n_target = 0)
            then (
              Mat.save_txt
                ~out:
                  (in_dir
                     (Printf.sprintf "w_%s_%.1f_%.1f_%i_%i" net_name rad sa seed n_target))
                w;
              Mat.save_txt
                ~out:
                  (in_dir
                     (Printf.sprintf "c_%s_%.1f_%.1f_%i_%i" net_name rad sa seed n_target))
                (AD.unpack_arr c))
            else Stdio.printf "tree %!"
          in
          let n_prep = Float.to_int (t.t_prep /. dt) in
          let t_prep_int = Float.to_int (1000. *. t.t_prep) in
          let _ = Stdio.printf "abc %!" in
          let xs, us, l, quus, _ =
            I.solve ~u_init:Mat.(gaussian ~sigma:0. 5001 m) ~n:(m + 4) ~m ~x0 ~prms t
          in
          let _ = Stdio.printf "abc %!" in
          let eigs, re, im = eigenvalues w in
          let sa = Mat.max' re in
          let max_eig = Mat.max' (Mat.l2norm ~axis:1 eigs) in
          let us = AD.unpack_arr us in
          let us_prep = Mat.get_slice [ [ 0; n_prep - 2 ] ] us in
          let us_mov = Mat.get_slice [ [ n_prep - 1; -1 ] ] us in
          let pi = Mat.(l2norm' us_prep) /. Mat.(l2norm' us_mov) in
          let res =
            Mat.of_arrays [| [| rad; norm_w; pi; AD.unpack_flt l; max_eig; sa |] |]
          in
          Mat.save_txt
            ~out:
              (in_dir (Printf.sprintf "%s_%i_%.1f_%.1f_%i" net_name n_target rad sa seed))
            res
        with
        | e -> Stdio.printf "%s" (Exn.to_string e)))
