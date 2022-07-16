open Owl
open Base

let _ = Printexc.record_backtrace true
let typ_net = "random_0.9"
let seed = 1

let dir_skew =
  "/home/mmcs3/rds/rds-t2-cs133-hh9aMiOkJqI/mmcs3/random_monkeys_skew_lambda_1E-6/rad_0.9/seed_1"


let dir_rdn =
  "/home/mmcs3/rds/rds-t2-cs133-hh9aMiOkJqI/mmcs3/random_monkeys_random_lambda_1E-6/rad_0.9/seed_1"


let dir_soc =
  "/home/mmcs3/rds/rds-t2-cs133-hh9aMiOkJqI/mmcs3/random_monkeys_lambda_1E-7/seed_6"


let lambda = 1E-6
let t_preps = [| 0.; 50.; 100. |]

let loss_over_t i d =
  Array.map t_preps ~f:(fun t ->
      let t_int = Float.to_int t in
      let loss = Mat.load_txt (Printf.sprintf "%s/loss_%i_%i" d i t_int) |> Mat.sum' in
      let u_cost_prep, u_cost_mov, u_cost_tot =
        Mat.load_txt (Printf.sprintf "%s/u_cost_%i_%i" d i t_int)
        |> fun x -> Mat.get x 0 0, Mat.get x 0 1, Mat.get x 0 2
      in
      let tgt_cost, torques_cost =
        Mat.load_txt (Printf.sprintf "%s/task_cost_%i_%i" d i t_int)
        |> fun x -> Mat.get x 0 1, Mat.get x 0 0
      in
      [| t; loss; u_cost_tot; tgt_cost; torques_cost; u_cost_prep; u_cost_mov |])
  |> Mat.of_arrays
  |> fun z -> Arr.reshape z [| 1; Array.length t_preps; -1 |]


let rdn_summary =
  Array.fold
    (Array.init 8 ~f:(fun i -> i))
    ~init:[]
    ~f:(fun accu i ->
      try
        let l = loss_over_t i dir_rdn in
        l :: accu
      with
      | _ -> accu)
  |> Array.of_list
  |> Arr.concatenate ~axis:0
  |> fun v ->
  let m, va =
    Arr.mean ~axis:0 v
    |> fun z ->
    ( Arr.reshape z [| Array.length t_preps; -1 |]
    , Arr.std ~axis:0 v |> fun z -> Arr.reshape z [| Array.length t_preps; -1 |] )
  in
  Mat.(m @|| va)


let _ =
  Mat.save_txt
    ~out:
      (Printf.sprintf
         "/home/mmcs3/rds/rds-t2-cs133-hh9aMiOkJqI/mmcs3/why_prep_summaries/%s_%f_%i"
         typ_net
         lambda
         seed)
    rdn_summary
