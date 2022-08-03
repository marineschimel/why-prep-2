open Owl
open Lib
module AD = Algodiff.D
open Defaults

let mean_t = 0.
let expected_t = mean_t
let go_t = 0.
let max_t = 0.1
let t_mov = 0.4
let duration = go_t +. t_mov +. 0.2
let get_step t = int_of_float (t /. sampling_dt)
let step_exp = get_step expected_t
let step_go = get_step go_t
let step_max = get_step max_t
let end_mov = get_step (t_mov +. go_t)
let last_step = end_mov + 200
let n_horizon = 800
let n_planned = 40
