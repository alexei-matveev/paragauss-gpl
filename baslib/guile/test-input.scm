((tasks (task "gradients") (dipole #t))
 (main_options (spin_restricted #t) (relativistic "false"))
 (symmetry_group (point_group "C2V"))
 (unique_atom_number (n_unique_atoms 2))
 (unique_atom
  (name "O")
  (z 8.0)
  (n_equal_atoms 1))
 (0.0 0.0 0.0) ; literal entries
 (unique_atom (name "H") (z 1.0) (n_equal_atoms 2))
 (1.429960 0.000000 -1.107190) ; literal entries
 (xc_control (xc "hf"))
 (grid (sym_reduce #t) (weight_grads #t))
 (gridatom (nrad 80) (nang 291))
 (gridatom (nrad 80) (nang 291))
 (basis "nwchem" "o" "6-31g")
 (basis "nwchem" "h" "6-31g"))

