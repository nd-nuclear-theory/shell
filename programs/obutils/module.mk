$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h :=
module_units_cpp-h :=
# module_units_f :=
module_programs_cpp := obmixer obscalc-ob ew-gen
module_programs_cpp += orbital-gen natorb-gen
module_programs_cpp += radial-gen radial-xform radial-compose
module_programs_cpp += obsolve obdme-conv orbital-extract
# module_programs_cpp_test := obdme-compare_test

# module_programs_f :=
# module_generated :=

################################################################
# library creation flag
################################################################

## $(eval $(library))

$(eval $(end-module))
