open Owl
module Z = Dense.Matrix.Z

let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir]")
let in_dir s = Printf.sprintf "%s/%s" dir s

let save_bin filename m =
  let output = open_out filename in
  Marshal.to_channel output m [ Marshal.No_sharing ];
  close_out output


(* reads whatever was saved using [save_bin] *)
let read_bin filename =
  if not Sys.(file_exists filename)
  then failwith Printf.(sprintf "%s does not exist" filename)
  else (
    let input = open_in filename in
    let m = Marshal.from_channel input in
    close_in input;
    m)


let print_dim_2d x =
  let shp = Owl.Arr.shape x in
  Printf.printf "%i, %i\n%!" shp.(0) shp.(1)


let print_dim_3d x =
  let shp = Owl.Arr.shape x in
  Printf.printf "%i, %i, %i\n%!" shp.(0) shp.(1) shp.(2)


let stack ?(axis = 0) xs =
  let shp = Owl.Arr.shape xs.(0) in
  let ndim = Array.length shp + 1 in
  let axis = Owl_utils.adjust_index axis ndim in
  let new_shp =
    Array.init ndim (fun i ->
        if i < axis then shp.(i) else if i = axis then 1 else shp.(i - 1))
  in
  let y =
    Array.map
      (fun x ->
        let shp' = Owl.Arr.shape x in
        if shp' <> shp
        then failwith "stack: ndarrays in [xs] must all have the same shape";
        Owl.Arr.reshape x new_shp)
      xs
  in
  Owl.Arr.concatenate ~axis y


let print_algodiff_dim x =
  let shp = Algodiff.D.shape x in
  Printf.printf "%i, %i\n%!" shp.(0) shp.(1)


(* participation_ratio [x] does PCA and calculate (sum_i lambda_i)^2 /. (sum_i lambda_i^2) 
x has dimensions n x n_samples  *)
let participation_ratio x =
  let x = Mat.(x - mean ~axis:1 x) in
  let _, s, _ = Linalg.D.svd x in
  let v = Mat.sqr s in
  Maths.sqr Mat.(sum' v) /. Mat.(sum' (sqr v))


let transform a =
  let n = Mat.row_num a in
  let v, lam = Linalg.D.eig a in
  let data = List.init n (fun i -> Z.get lam 0 i, Z.col v i) in
  let rec reorder data accu =
    match data with
    | []                         -> List.rev accu
    | [ (_, v) ]                 -> List.rev (Z.re v :: accu)
    | (l1, v1) :: (l2, v2) :: tl ->
      if Complex.(norm (sub (conj l1) l2)) < 1E-8
      then reorder tl (Z.re v1 :: Z.im v1 :: accu)
      else reorder ((l2, v2) :: tl) (Z.re v1 :: accu)
  in
  let v = reorder data [] |> Array.of_list |> Mat.concatenate ~axis:1 in
  Mat.(inv v *@ a *@ v), v
