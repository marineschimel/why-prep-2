open Owl

let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")

let _ =
  Array.init 56 (fun i ->
      let _ = Printf.printf "%i \n %!" i in
      let x = Mat.load_txt (Printf.sprintf "%s/qtuu_%i" dir i) in
      let _, s, _ = Linalg.D.svd x in
      Printf.printf "%f \n!" (Mat.min' s))
