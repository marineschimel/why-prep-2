open Owl
open Lib
open Defaults

let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")

let in_dir s = Printf.sprintf "%s/%s" dir s
(*Q(t) = max_t' ||(\dot{x}(t)-\dot{x}(t')||^2/(||(\dot{x}(t)-\dot{x}(t')||^2+\epsilon)*)


let epsilon = 0.1*. (Mat.(get_slice [[];[4;-1]] (Mat.load_txt (in_dir "traj_400"))) |> fun z -> Mat.(sum' (var ~axis:0 z)))
let dt = sampling_dt

let sampling_window = 1


let tangling t traj_1 traj_2 dt =
  let epsilon =  0.1*. (traj_1 |>  fun z -> Mat.(sum' (var ~axis:0 z))) in 
  let get_row i y = Mat.get_slice [ [ i ]; [] ] y in
  let state_t = Mat.get_slice [ [ t ]; [] ] traj_1
  and dx_t =
    Mat.((get_slice [ [ t]; [] ] traj_1 - get_slice [ [ pred t ]; [] ] traj_1) /$ dt)
  and dx_2 =
    Mat.((get_slice [ [ 1; -1 ]; [] ] traj_2 - get_slice [ [ 0; -2 ]; [] ] traj_2) /$ dt)
  in
  let entangled =
    Mat.mapi_rows
      (fun i r ->
        let dx' = get_row i dx_2 in
        let top = Mat.l2norm_sqr' Mat.(dx' - dx_t)
        and bottom = Mat.l2norm_sqr' Mat.(r - state_t) +. epsilon in
        top /. bottom)
      (Mat.get_slice [ [ 1; -1 ]; [] ] traj_2)
  in
  let max_idx = Utils.Array.max_i entangled in
  Mat.max' (Mat.of_array entangled 1 (-1)), max_idx


(* || Tests || *)

let activity ?reach:_reach time =
  let _reach =
    match _reach with
    | Some a -> a
    | None   -> 1
  in
  let idces = 
  List.init (400/sampling_window) (fun i -> time + i*sampling_window) in 
  let x = Mat.get_fancy
    [ L idces; R [ 4; -1 ] ]
    (Mat.load_txt (Printf.sprintf "%s/traj_%i" dir time))
  in Mat.(x/$(max' x -. min' x +. 0.05)) |> fun z -> Mat.(z - mean ~axis:0 z)
  


let muscles ?reach:_reach time =
  let _reach =
    match _reach with
    | Some a -> a
    | None   -> 1
  in
  let idces = 
    List.init (400/sampling_window) (fun i -> i*sampling_window) in 
  let x = Mat.get_fancy
      [ L idces; R [0;1] ]
  (Mat.load_txt (Printf.sprintf "%s/torques_%i" dir time))
  in Mat.(x/$(Mat.max' x -. Mat.min' x)) |> fun z -> Mat.(z - mean ~axis:0 z)


let _ =
  let tang, i = tangling 1 (activity 400) (activity 400) dt in
  Printf.printf "%f %i \n %!" tang i


  let pca_proj x dim = let cov = Mat.((transpose x)*@x) in let modes,_,_ = Linalg.D.svd cov
  in let top_modes = Mat.get_slice [[];[0;dim]] modes in Mat.(x*@top_modes)


let save_tangling t_prep =
  let nr = Mat.row_num (activity t_prep) - 1 in
  let traj_tang =  Array.init 1 (fun idx -> let idx = idx in let j = idx/4 in let i = idx mod 4 in Mat.(
    init nr 1 (fun k ->
        let tang, _ = tangling (succ k) (pca_proj (activity ~reach:(succ i) t_prep) 8) (pca_proj (activity ~reach:(succ j) t_prep) 8) dt in
        tang))) |> Mat.concatenate ~axis:1 in 
            let mot_tang = 
               Array.init 1 (fun idx -> let j = idx/4 in let i = idx mod 4 in let _ = Printf.printf "%i %i %!" i j in Mat.(
                init nr 1 (fun k ->
                    let tang, _ = tangling (succ k) ((muscles ~reach:(succ i) t_prep)) ((muscles ~reach:(succ j) t_prep)) dt in
                    tang))) |> Mat.concatenate ~axis:1 in        
  Mat.save_txt
    ~out:((in_dir (Printf.sprintf "entanglement/traj_%i" t_prep))) traj_tang
   ;
  Mat.save_txt
    ~out:(in_dir (Printf.sprintf "entanglement/torques_%i" t_prep)) mot_tang;
  Mat.save_txt ~out:(in_dir "both_entanglements") Mat.((traj_tang@||mot_tang))


let _ = save_tangling 600


let _ = let m_mov = Mat.load_txt "results/weighing_pm/w_1000_1/entanglement/traj_600" |> Mat.max'
                in let m_mixed = Mat.load_txt "results/weighing_pm/w_1_1/entanglement/traj_600" |> Mat.max'
              in let m_prep = Mat.load_txt "results/weighing_pm/w_1_1000/entanglement/traj_600" |> Mat.max'
            in Mat.save_txt ~out:"results/weighing_pm/comp_entanglement" (Mat.transpose (Mat.of_arrays [|[|m_mov;m_mixed;m_prep|]|]))
(*do across conditions if I can!
Compute entanglement for projection of activity/compare w full state
get results for multiple netwtorks with inputs at different times and make a plot? 
For visualization make PCA plot across conditions
 *)
