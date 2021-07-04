open Owl

let xs = Mat.load_txt "test_c/xs_5_200"
let c = Mat.load_txt "data/c"
let torques = Mat.(xs *@ transpose c)
let _ = Mat.save_txt ~out:"test_c/torques_5_200" torques
