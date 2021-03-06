open Owl
module A = Analysis

let _ = Printexc.record_backtrace true

let rinv =
  Cmdargs.(get_int "-rinv" |> force ~usage:"-rinv [inverse of r_coeff]")

let noises = [| 0.000; 0.001; 0.010; 0.100; 1.000; 2.000; 10.000 |]

let tprep = 300

let thetas_in, traj_in =
  let x =
    Mat.load_txt (Printf.sprintf "results/noise/soc_r_%i/traj_%i" rinv tprep)
  in
  (Mat.get_slice [ []; [ 0; 1 ] ] x, Mat.get_slice [ []; [ 4; -1 ] ] x)

let gather_loss =
  Array.map
    (fun n ->
      let noise = Mat.of_arrays [| [| n |] |] in
      Mat.(
        noise
        @|| Mat.(
              load_txt
                (Printf.sprintf "results/noise/soc_r_%i/level_%.3f/loss_time"
                   rinv n))))
    noises

let comp_neural_traj =
  Array.map
    (fun n ->
      let traj =
        Mat.get_slice [ []; [ 4; -1 ] ]
          (Mat.load_txt
             (Printf.sprintf "results/noise/soc_r_%i/level_%.3f/traj_%i" rinv n
                tprep))
      in
      Mat.of_arrays [| [| n; Mat.(l2norm' (traj - traj_in)) |] |])
    noises

let comp_theta_traj =
  Array.map
    (fun n ->
      let thetas =
        Mat.get_slice [ []; [ 0; 1 ] ]
          (Mat.load_txt
             (Printf.sprintf "results/noise/soc_r_%i/level_%.3f/traj_%i" rinv n
                tprep))
      in
      Mat.of_arrays [| [| n; Mat.(l2norm' (thetas - thetas_in)) |] |])
    noises

let _ =
  Mat.save_txt
    ~out:(Printf.sprintf "results/noise/soc_r_%i/losses" rinv)
    (Mat.concatenate ~axis:0 gather_loss)
