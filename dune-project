(lang dune 2.7)

(name ndsort)

(formatting
 (enabled_for ocaml))

(authors "\197\129ukasz Kurowski")

(maintainers "crackcomm@gmail.com")

(generate_opam_files true)

(package
 (name ndsort)
 (synopsis "OCaml bindings to parallel non-dominated sorting")
 (description "Parallel implementations of the Best Order Sort algorithm for Non-Dominated Sorting")
 (depends
  (matplotlib :with-test)
  (ocaml
   (>= 4.07.0))))
