open Bigarray

external sort
  :  (float, float32_elt, c_layout) Array2.t
  -> (int32, int32_elt, c_layout) Array1.t
  -> unit
  = "ml_ndsort"

let sort input =
  let output = Array1.create Int32 c_layout (Array2.dim1 input) in
  sort input output;
  output
;;
