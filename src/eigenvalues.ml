open Owl
module Z = Dense.Matrix.Z

let mat =
  let m = Mat.load_txt "data/w_end" in
  Mat.(m - eye 200)


let _, eigvals = Linalg.D.eig mat

let _ =
  Mat.save_txt
    ~out:"eigenvalues/eig_end_to_end"
    Mat.(transpose (Z.re eigvals @= Z.im eigvals))


let _ =
  let u, _ = Mat.load_txt "data/u_1", Mat.load_txt "data/v_1" in
  Printf.printf "%f %!" (Mat.l2norm_sqr' u /. 200.)
