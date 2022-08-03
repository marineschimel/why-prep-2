open Owl
module AD = Algodiff.D

let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir]")
let in_dir = Printf.sprintf "%s/%s" dir

let us i t_prep =
  Mat.load_txt (in_dir (Printf.sprintf "us_%i_%i" i Float.(to_int (1000. *. t_prep))))

let xs i t_prep =
  Mat.load_txt (in_dir (Printf.sprintf "xs_%i_%i" i Float.(to_int (1000. *. t_prep))))

(*general function to perform the regression from inputs to activity, we might want do this either using the same regression matrix across movements, or using a single one every time
We have the option to do it either using least squares, or using optimization (might be useful if the matrices are too big)
We return A and the residuals *)

(*here us = TxM and xs = TxN, a = NxM*)
let regress ?(analytical = true) us xs =
  let a =
    if analytical
    then (
      let xinv = Linalg.D.linsolve Mat.(transpose xs *@ xs) Mat.(transpose xs) in
      Mat.(xinv *@ us))
    else
      let module Prms = struct
        type 'a t = { m : 'a } [@@deriving prms]
      end
      in
      let us, xs = AD.pack_arr us, AD.pack_arr xs in
      let open Prms in
      let prms0 =
        { m = AD.Mat.gaussian ~sigma:0.01 (AD.Mat.col_num xs) (AD.Mat.col_num us) }
      in
      let module O = Owl_opt_lbfgs.D.Make (Prms) in
      let f _ prms = AD.Maths.(l2norm' (us - (xs *@ prms.m))) in
      let s = O.init ~prms0 () in
      let stop fv _ = fv < 1E-4 in
      let _ = O.min ~f ~stop s in
      let prms = O.prms s in
      AD.unpack_arr prms.m
  in
  a, Mat.(us - (xs *@ a))

let n = 7
let t_prep = 0.3
let us, xs = us n t_prep, xs n t_prep
let a, res = regress us xs

let _ =
  Mat.save_txt
    ~out:(in_dir (Printf.sprintf "residuals_%i_%i" n (Float.to_int (1000. *. t_prep))))
    res;
  Mat.save_txt
    ~out:(in_dir (Printf.sprintf "a_%i_%i" n (Float.to_int (1000. *. t_prep))))
    a
