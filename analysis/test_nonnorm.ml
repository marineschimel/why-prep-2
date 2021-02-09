open Owl
open Lib
open Defaults 
let _ = Printexc.record_backtrace true

let file = Cmdargs.(get_string "-file" |> force ~usage:"-file [file we're looking at]")

let _ = Printf.printf "%f %!" (nonnormality (Mat.load_txt (Printf.sprintf "%s" file)))

let _ = (Array.map (fun (i,j)-> let x = Mat.load_txt (Printf.sprintf "results/weighing_pm/w_%i_%i/loss_time" (int_of_float i) (int_of_float j))
in let ratio = Mat.of_array [|j/.i;j/.i;j/.i|] (-1) 1 in Mat.((ratio@||x))) [|(1000.,1.);(500.,1.);(100.,1.);(10.,1.);(5.,1.);(1.,1.);(1.,5.);(1.,10.);(1.,100.);(1.,500.);(1.,1000.);
|])
|> Mat.concatenate ~axis:0 |> fun z -> Mat.save_txt ~out:"results/weighing_pm/summary" z 

let filter_t t = let x = Mat.load_txt "results/weighing_pm/summary" in (Mat.filter_rows (fun r ->  Mat.get r 0 1 = t ) x) |> fun y -> Mat.get_fancy [L (Array.to_list y);R [0;-1]] x |> fun l -> Mat.save_txt ~out:(Printf.sprintf "results/weighing_pm/summary_%i" (int_of_float (1000.*.t))) l

let _ = filter_t 0.4; filter_t 0.6 


(*fun z -> Mat.get_slice [(Array.to_list z);[]] x |> Mat.save_txt ~out:(Printf.sprintf "results/weighing_pm/summary_%i" (int_of_float (1000.*.t)))