(*Similarity transforms don't change the transfer function of the system : therefore, every SOC should be able to be transformed into a system wuth a normal A that prepares just as much?? In this file I'll compute all the matrices obtained from the similarity transform and see how the internal dynamics compare... we expecct that the optimal inputs and the output should be the same as for the original system (I think?) or at least that for the initial set of inputs we'll recover the right output*)
open Owl
module AD = Algodiff.D
module M = Arm.Make (Arm.Defaults)
open Lib
module T = Core.Time
open Defaults
open Printf
open Misc
module A = Analysis_funs


let transform w b c = 
  let w_tilde, t = transform w in let b_tilde = Mat.((inv t)*@b) in let c_tilde = Mat.(c*@t)
in w_tilde, b_tilde, c_tilde

let _ = Printexc.record_backtrace true
let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")

let in_dir s = Printf.sprintf "%s/%s" dir s

(* let num_init = Cmdargs.(get_int "-init" |> force ~usage:"-init [init number]") *)
let net = Cmdargs.(get_string "-net" |> force ~usage:"-net [dir to save in]")
 
let evaluate
    ?target:_target
    ?qs_coeff:_t_coeff
    ?r_coeff:_r_coeff
    ?t_prep:_t_prep
    ?gamma:_gamma
    ?cost:_cost
    ?x0:_x0
    ?weighing_pm:_wpm
    ~t_mov
    ~c
    ~w
    ~b
    ~annealing
    subdir 
  =
  let module PT = struct
    let qs_coeff =
      match _t_coeff with
      | Some a -> a
      | None   -> 1000.


    let cost =
      match _cost with
      | Some a -> a
      | None   -> "running"


    let r_coeff =
      match _r_coeff with
      | Some a -> Defaults.r_coeff*.a
      | None   -> Defaults.r_coeff


    let t_prep =
      match _t_prep with
      | Some a -> a
      | None   -> 0.


    let gamma_exponent =
      match _gamma with
      | Some a -> a
      | None   -> 2.


    let x0 =
      match _x0 with
      | Some x0 -> x0
      | None    ->
        AD.Maths.(
          concatenate
            ~axis:1
            [| initial_theta
             ; AD.Mat.zeros 1 n (*; Mat.of_array [| -0.5 |] 1 (-1) |> AD.pack_arr*)
            |])


    let target_theta =
      match _target with
      | Some a -> a
      | None   -> [| Mat.of_arrays [| [| -0.287147; 2.39155; 0.; 0. |] |] |]

    let weighing_pm =
        match _wpm with
        | Some (wp,wm) -> (wp,wm)
        | None   -> (1.,1.)
    let b = AD.pack_arr b
    let w = w
    let n = size_net+4
    let m = size_inputs
    let aug = 0
    let t_mov = t_mov
    let duration = t_prep +. t_mov +. 0.2
    let saving_dir = Printf.sprintf "%s/%s/%s" dir subdir
    let c = c 
    let __c = AD.pack_arr c
  end
  in
  let module P = Prep.Make (PT) (Costs.C_Running (PT)) in
  let module I = Ilqr.Default.Make (P) in
  let __w = AD.pack_arr PT.w in 
  let angles_to_x x0 us =
    let traj = I.trajectory x0 us in
    AD.Mat.map_by_row (fun x -> Arm.unpack_state_diff (M.hand_of (Arm.pack_state x))) traj
  in
  let x0 = PT.x0 in
  let us = match annealing with 
  |(false,_)  ->  List.init P.n_steps (fun _ -> AD.Mat.zeros 1 size_inputs)
  |(true,t_pred) ->  let _ = Printf.printf "%f %f %!" t_pred PT.t_prep in let n_pred_prep = int_of_float (t_pred/.sampling_dt) in 
  let _n_prep = int_of_float (PT.t_prep/.sampling_dt) in 
  let _ = Printf.printf "%i = i %!" n_pred_prep in 
  let inpts = Mat.load_txt (PT.saving_dir (Printf.sprintf "results_us_%i"
  n_pred_prep)) |> AD.pack_arr in 
  let _ = Printf.printf "n_steps_pred : %i %i new : %i %!" (AD.Mat.row_num inpts) (AD.Mat.col_num inpts) P.n_steps in 
  let delta_n = _n_prep - n_pred_prep +1 in 
 List.init P.n_steps (fun i -> if i<= delta_n then AD.Mat.zeros 1 size_inputs
  else AD.Maths.get_slice [[i - delta_n - 1];[]] inpts)
  in
  I.trajectory x0 us |> AD.unpack_arr |> Mat.save_txt ~out:(PT.saving_dir "traj0");
  angles_to_x x0 us |> AD.unpack_arr |> Mat.save_txt ~out:(PT.saving_dir "x0");
  let t0 = T.now () in
  let stop =
    let cprev = ref 1E9 in
    fun k us ->
      let c = I.loss x0 us in
      let pct_change = abs_float (c -. !cprev) /. !cprev in
      let traj_ad = I.trajectory x0 us in
      let traj = traj_ad |> AD.unpack_arr in
      let inputs = us |> Array.of_list |> AD.Maths.concatenate ~axis:0 |> AD.unpack_arr in
      let _ = Printf.printf "%i %i %i %i %!" (AD.Mat.row_num __w) (AD.Mat.col_num __w)
      (AD.Mat.row_num traj_ad) (AD.Mat.col_num traj_ad) in 
      let recurrent =
          AD.Maths.(
            transpose (__w  *@ g (transpose (get_slice [ [ 0; -2 ]; [ 4; pred PT.n ] ] traj_ad))))
          |> AD.unpack_arr
        in
      if k mod 1 = 0
      then (
        Printf.printf "iter %i | cost %f | pct change %f\n%!" k c pct_change;
        cprev := c;
        Mat.save_txt
          ~out:(PT.saving_dir (sprintf "traj_%i" (int_of_float (PT.t_prep *. 1000.))))
          traj;
        Mat.save_txt
          ~out:
            (PT.saving_dir
               (sprintf "floored_traj_%i" (int_of_float (PT.t_prep *. 1000.))))
          (Mat.map (fun y -> if y > 0. then y else 0.) traj);
        Mat.save_txt
          ~out:(PT.saving_dir (sprintf "torques_%i" (int_of_float (PT.t_prep *. 1000.))))
          (Mat.concatenate
             ~axis:0
             (Mat.map_rows
                (fun x -> Mat.(transpose (PT.c *@ (AD.unpack_arr (g (AD.pack_arr (transpose x)))))))
                (Mat.get_slice [ []; [ 4; PT.n - 1 ] ] traj)));
        AD.unpack_arr (angles_to_x x0 us)
        |> Mat.save_txt
             ~out:(PT.saving_dir (sprintf "hands_%i" (int_of_float (PT.t_prep *. 1000.))));
        inputs
        |> Mat.save_txt
             ~out:
               (PT.saving_dir
                  (sprintf "results_us_%i" (int_of_float (PT.t_prep *. 1000.))));
        Mat.(inputs*@(transpose (AD.unpack_arr PT.b)))
                  |> Mat.save_txt
                       ~out:
                         (PT.saving_dir
                            (sprintf "transformed_us_%i" (int_of_float (PT.t_prep *. 1000.))));
        (* Mat.(inputs *@ transpose v)
        |> Mat.save_txt
             ~out:
               (PT.saving_dir
                  (sprintf "results_us_scaled_%i" (int_of_float (PT.t_prep *. 1000.)))); *)
        (* Mat.save_txt
          ~out:
            (PT.saving_dir (sprintf "traj_scaled_%i" (int_of_float (PT.t_prep *. 1000.))))
          Mat.(get_slice [ []; [ 4; -1 ] ] traj *@ transpose v); *)
        let t = T.now () in
        let dt = T.diff t t0 in
        Printf.printf "Time iter %i | %s | %!" k (T.Span.to_string dt));
      (* Mat.(inputs + recurrent)
        |> Mat.save_txt
             ~out:
               (PT.saving_dir
                  (sprintf "effective_us_%i" (int_of_float (PT.t_prep *. 1000.)))); *)
      recurrent
        |> Mat.save_txt
             ~out:
               (PT.saving_dir
                  (sprintf "recurrent_us_%i" (int_of_float (PT.t_prep *. 1000.))));
      if pct_change < 1E-3
      then (
        Mat.(
          save_txt
            ~append:true
            ~out:(PT.saving_dir "loss_time")
            (Mat.of_array
               [| PT.t_prep; PT.r_coeff; I.loss x0 us; nonnormality w |]
               1
               (-1)));
        let ip, im = A.inputs ~t_prep:PT.t_prep inputs in
        let cp, cm = A.cost_inputs ~t_prep:PT.t_prep inputs in
        Mat.save_txt
          ~append:true
          ~out:(PT.saving_dir "input_distrib")
          (Mat.of_arrays [| [| PT.t_prep; ip; im |] |]);
        Mat.save_txt
          ~append:true
          ~out:(PT.saving_dir "sqr_input_distrib")
          (Mat.of_arrays [| [| PT.t_prep; cp; cm |] |]);
        Mat.(
          save_txt
            ~append:true
            ~out:(PT.saving_dir "sum_us")
            (Mat.of_array
               [| PT.t_prep
                ; PT.r_coeff
                ; inputs |> Mat.map (fun x -> Maths.abs x *. r_coeff) |> Mat.sum'
               |]
               1
               (-1)));
        Mat.(
          save_txt
            ~append:true
            ~out:(PT.saving_dir "accuracy")
            (Mat.of_array [| PT.r_coeff; PT.t_prep; A.accuracy_theta traj |] 1 (-1))));
      pct_change < 1E-3
  in
  let final_us = I.learn ~stop x0 us in
  let energy =
    List.fold_left
      (fun acc x ->
        let x = AD.unpack_arr x in
        acc +. (Mat.(sqr x |> sum') *. sampling_dt))
      0.
      final_us
  in
  energy


let targets = Mat.load_txt "data/target_thetas"

let run_run =

    
let _dir_rad =  "big_soc" in 

  (* let dir_rad = (Printf.sprintf "random/rad_%i" (int_of_float (radius*.1000.))) in 

  let w = Mat.load_txt (Printf.sprintf "results_c/%s/size_%i" dir_rad) in   *)
  
  (* let c = Mat.load_txt "results_smal" in  *)
 (* let c = Mat.load_txt (Printf.sprintf "results_small/%s/new_c" _dir_rad) in   *)
 (*let _tpreds =  [| 0.1;0.2;0.3;0.4;0.5;0.6;0.7;0.8 |]
 in let null_c = Linalg.D.null c in 
 let _ = Printf.printf "%i %i %!" (Mat.row_num null_c) (Mat.col_num null_c) in 
let _x_init = Mat.(null_c *@ (gaussian ~sigma:0.2 198 1)) |> Mat.transpose |> AD.pack_arr in*) 
 Array.map (fun _ ->
  (* let w = Mat.load_txt (Printf.sprintf "results_c/%s/alpha_%i/w" net _alpha) in
  let c = Mat.load_txt (Printf.sprintf "results_c/%s/alpha_%i/c" net _alpha) in *)
  (* let w = Mat.load_txt (Printf.sprintf "results_bc_r2/%s/size_%i/w" net size_inputs) in
  let c = Mat.load_txt (Printf.sprintf "results_bc_r2/%s/size_%i/c" net size_inputs) in *)
  let _w = Mat.load_txt ("results/reach_1/w") in 
  let w = make_random 0.8 in 
  let c = Mat.load_txt ("results/reach_1/c") in 
  let b = ( Defaults.__b) in 
   (*let w,_b,_c = transform w (AD.unpack_arr b) c in *)
   let b = AD.unpack_arr ( Defaults.__b) in 
  let _ = Mat.save_txt ~out:(in_dir "similarity/w") w;  Mat.save_txt ~out:(in_dir "similarity/c") c; Mat.save_txt ~out:(in_dir "similarity/b") b in 
  (* let c = Mat.load_txt (Printf.sprintf "results_bc_r2/%s/size_%i/c" net size_inputs) in *)
    let x0 =
    AD.Maths.(concatenate ~axis:1 [| initial_theta; AD.Mat.zeros 1 n |])
  in Array.mapi
    (fun _ t ->
      evaluate
        ~x0
        ~t_mov:0.4
        ~t_prep:t
        ~gamma:2.
        ~target:[| Mat.row targets 0|]
        ~c
        ~w
        ~b
        ~r_coeff:0.1
        ~qs_coeff:1.
        ~annealing:(false,0.)
        ~weighing_pm:(1.,1.)
        (Printf.sprintf "similarity/random_net"))
        (* (Printf.sprintf "%s/compound_%i%i" net i j)) *)
       
    [|0.4|] )
    [|0|]