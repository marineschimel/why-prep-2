open Owl
module Z = Dense.Matrix.Z

let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir to save in]")
let file = Cmdargs.(get_string "-f" |> force ~usage:"-f [file we're fetching]")
let mat = Mat.load_txt (Printf.sprintf "%s/%s" dir file)
let _, eigvals = Linalg.D.eig mat

let _ =
  Mat.save_txt
    ~out:(Printf.sprintf "%s/%s_eigs" dir file)
    Mat.(transpose (Z.re eigvals @= Z.im eigvals))
