open Owl
module AD = Algodiff.D
open Arm
module M = Arm.Make (Arm.Defaults)

(* open Misc *)

let sampling_dt = 1E-3

let tau = 150E-3

let size_net =
  Cmdargs.(get_int "-size" |> force ~usage:"-size [size of network]")

let targets = Mat.load_txt "data/target_thetas"

let size_inputs =
  Cmdargs.(get_int "-b_size" |> force ~usage:"-b_size [size of inputs]")

let __b = AD.Maths.(AD.Mat.eye size_inputs / F tau)

(* let __c = Mat.(0.01$*(Mat.load_txt "reach_1/c")) *)

let __c = Mat.gaussian ~sigma:0.1 2 200

let n = size_net

let r0 = AD.F 0.1

let lambda = AD.F 0.2

let g x = x

let g_prm ~x:_ = AD.F 1.

(* 
let g x = AD.Maths.relu x 
let g_prm ~x = AD.Maths.(F 0.5 *(signum x + F 1.)) *)

let initial_theta =
  Mat.of_arrays [| [| 0.174533; 2.50532; 0.; 0. |] |] |> AD.pack_arr

let target_duration = 0.4

let pause = 0.05

let r_coeff = 1E-1 /. float size_net

let q_coeff = 1.

let q_start = 0.

let a_coeff = 0.

let t_coeff = 0.

let t_mat =
  AD.Mat.of_arrays
    [|
      [| r_coeff *. 1000. /. float size_net; 0. |];
      [| 0.; 1000. *. r_coeff /. float size_net /. 3. |];
    |]

let new_targets n duration =
  let module Task = struct
    let angle_i = 0.

    let angle_f = 360.

    let reach_angles ?(angle_i = angle_i) ?(angle_f = angle_f) n_angles =
      Mat.linspace angle_i angle_f n_angles |> Mat.to_array
  end in
  let reach_angles = Task.reach_angles n in
  Array.iter
    (fun reach_angle ->
      let hand_trajectory, _ =
        M.straight_reach ~tau:0.140 ~dt:sampling_dt ~duration
          ~angle:(Arm.deg reach_angle) ~radius:0.12
      in
      hand_trajectory |> Arm.unpack_sequence
      |> Mat.get_fancy [ R [ -1 ]; R [] ]
      |> Mat.save_txt ~append:true ~out:"data/target_hands")
    reach_angles

let unpack_full_state x n_theta =
  let thetas = AD.Maths.get_slice [ []; [ 0; n_theta - 1 ] ] x
  and xs =
    AD.Maths.get_slice [ []; [ n_theta; -1 ] ] x
    (*and x_switch = AD.Maths.get_slice [ []; [ -1 ] ] x |> AD.Maths.sum' *)
  in
  (thetas, xs)

let unpack_pos x =
  let thetas = AD.Maths.get_slice [ []; [ 0; 1 ] ] x in
  thetas

let unpack_vel x =
  let thetas = AD.Maths.get_slice [ []; [ 2; 3 ] ] x in
  thetas

let stack ?(axis = 0) xs =
  let shp = Owl.Arr.shape xs.(0) in
  let ndim = Array.length shp + 1 in
  let axis = Owl_utils.adjust_index axis ndim in
  let new_shp =
    Array.init ndim (fun i ->
        if i < axis then shp.(i) else if i = axis then 1 else shp.(i - 1))
  in
  let y =
    Array.map
      (fun x ->
        let shp' = Owl.Arr.shape x in
        if shp' <> shp then
          failwith "stack: ndarrays in [xs] must all have the same shape";
        Owl.Arr.reshape x new_shp)
      xs
  in
  Owl.Arr.concatenate ~axis y

let nonnormality a =
  let a, _, _ = Linalg.D.schur a in
  let module Z = Dense.Matrix.Z in
  let _, evals = Linalg.D.eig a in
  let re, im = (Z.re evals, Z.im evals) in
  1. -. ((Mat.l2norm_sqr' re +. Mat.l2norm_sqr' im) /. Mat.l2norm_sqr' a)

let a_sym_psi radius =
  let sigma = Maths.(sqrt (radius /. float n)) in
  let x = Mat.gaussian ~sigma n n in
  Mat.((((x + transpose x) /$ 2.) - eye n) /$ tau)

let a_random radius =
  let sigma = Maths.(sqrt (radius /. float n)) in
  let x = Mat.gaussian ~sigma n n in
  Mat.((x - eye n) /$ tau)

let a_antisym_psi radius =
  let sigma = Maths.(sqrt (2. *. radius /. float n)) in
  let x = Mat.gaussian ~sigma n n in
  Mat.((((x - transpose x) /$ 2.) - eye n) /$ tau)

let obs_gramian a c =
  Linalg.D.lyapunov Mat.(transpose a) Mat.(neg (transpose c) *@ c)

(* let cmc = AD.Mat.of_arrays [|[|1.;-2.;0.;0.|]; [|0.;0.;2.;-1.|]|] *)

let cmc =
  AD.Mat.of_arrays
    [|
      [|
        0.2;
        -0.3;
        0.8;
        -1.;
        0.1;
        -0.01;
        1.;
        -0.5;
        0.2;
        0.5;
        -0.5;
        0.1;
        0.;
        0.;
        0.;
        0.;
        0.;
        0.;
        0.;
        0.;
        0.;
        0.;
        0.;
        0.;
      |];
      [|
        0.;
        0.;
        0.;
        0.;
        0.;
        0.;
        0.;
        0.;
        0.;
        0.;
        0.;
        0.;
        2.;
        -1.;
        0.1;
        0.8;
        -2.;
        0.3;
        -0.12;
        0.8;
        0.1;
        -0.3;
        0.5;
        -0.8;
      |];
    |]

let phi x = AD.Maths.(log (exp x + F 1.))

let dphi y =
  let k = AD.Mat.row_num y in
  AD.Mat.init_2d k k (fun i j ->
      if i = j then
        let x = AD.Mat.get y i 0 in
        AD.Maths.(exp x / (exp x + F 1.))
      else F 0.)

(* let phi x = AD.Maths.(log (exp x + F 1.))

   let dphi y = let k = AD.Mat.row_num y in
   AD.Mat.init_2d k k (fun i j ->
       if i = j then (let x = AD.Mat.get y i 0 in  AD.Maths.(exp x/(exp x + F 1.))) else F 0.) *)
(* let dphi y = let k = AD.Mat.row_num y in
   AD.Mat.init_2d k k (fun i j ->
       if i = j then (if ((AD.Mat.get y i 0) < AD.F 0.) then AD.F 0. else F 1.) else F 0.) *)
