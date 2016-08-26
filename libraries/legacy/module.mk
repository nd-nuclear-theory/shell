$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h := pair_indexing
module_units_cpp-h := shell_2body shell_indexing_nlj shell_radial_nl shell_xform shell_separable

# omit preliminary lj-space work: shell_indexing_lj shell_radial_nlj
# omit nl orbitals: shell_indexing_nl

# module_units_f := 
module_programs_cpp := pair_indexing_test shell_indexing_test

# module_programs_f :=
# module_generated :=

################################################################
# library creation flag
################################################################

$(eval $(library))

$(eval $(end-module))
