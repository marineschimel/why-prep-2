open Owl

let mat =
  Arr.load_npy
    "/home/mmcs3/rds/hpc-work/_results/why_prep/baselines_larger_prep_2/all_trials"

let _ =
  Stdio.printf "%i %i %i %!" (Arr.shape mat).(0) (Arr.shape mat).(1) (Arr.shape mat).(2)
