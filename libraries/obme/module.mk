$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h := obme_operator
module_units_cpp-h := obme_io radial obme
# module_units_f :=
module_programs_cpp_test := obme_io_test obme_test

# module_programs_f :=
# module_generated :=

################################################################
# library creation flag
################################################################

$(eval $(library))

$(eval $(end-module))
