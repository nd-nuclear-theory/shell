$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h :=
module_units_cpp-h :=
# module_units_f :=
module_programs_cpp := h2stat h2conv h2mixer me2jfilter me2j2h2
# module_programs_cpp_test += make-radial-dummy make-radial-permutation

# module_programs_f :=
# module_generated :=

################################################################
# library creation flag
################################################################

## $(eval $(library))

$(eval $(end-module))
