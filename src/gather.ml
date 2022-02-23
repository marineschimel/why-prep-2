open Owl
open Base
open Lib
open Misc

(* let dir = "/home/mmcs3/rds/hpc-work/_results/why_prep/baselines_larger_prep_discrete"
let in_dir s = Printf.sprintf "%s/%s" dir s  *)

(* let _ = Misc.read_bin (in_dir "all_trials") *)
let _ =
  Array.init 1520 ~f:(fun i ->
      try
        let x = Mat.load_txt (in_dir (Printf.sprintf "rates_%i" i)) in
        let y = Mat.load_txt (in_dir (Printf.sprintf "us_%i" i)) in 
        assert (Int.(Mat.row_num x = 551));
        let n_steps = Mat.row_num y in
        Some (Arr.reshape y [| 1; n_steps; -1 |])
      with
      | e ->
        Stdio.printf "%s" (Exn.to_string e);
        None)
  |> fun a ->
  Array.fold a ~init:[] ~f:(fun acc x ->
      match x with
      | Some m -> m :: acc
      | None -> acc)
  |> Array.of_list
  |> Arr.concatenate ~axis:0
  |> Misc.save_bin (in_dir "all_inputs_nonlinear")
