open Bigarray
open Matplotlib

let () =
  let size = 100_000 in
  let input = Array2.init Float32 c_layout size 3 (fun _ _ -> Random.float 2.0 -. 1.0) in
  let ranks = Ndsort.sort input in
  let get_point i = input.{i, 0}, input.{i, 1}, input.{i, 2} in
  let get_front n =
    let res = ref [] in
    let n = Int32.of_int n in
    for i = 0 to size - 1 do
      let v = Array1.get ranks i in
      if Int32.equal v n then res := get_point i :: !res
    done;
    List.rev !res |> Array.of_list
  in
  let fronts = [ get_front 0; get_front 1; get_front 2; get_front 3 ] in
  let fig = Fig.create () in
  let ax = Fig.add_subplot_3d fig ~nrows:1 ~ncols:1 ~index:1 in
  Ax3d.set_xlabel ax "max drawback";
  Ax3d.set_ylabel ax "loss";
  Ax3d.set_zlabel ax "cost";
  Ax3d.grid ax true;
  List.iteri
    (fun idx front ->
      Printf.printf "Front %d = %d\n" idx (Array.length front);
      Ax3d.scatter
        ax
        ~c:
          (match idx with
          | 0 -> Green
          | 1 -> Blue
          | 2 -> Orange
          | 3 -> Red
          | _ -> failwith "too many fronts")
        front)
    fronts;
  Mpl.show ()
;;
