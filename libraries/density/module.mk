$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h :=
module_units_cpp-h := obdme_io
# module_units_f :=
module_programs_cpp_test := obdme_io_test

# module_programs_f :=
# module_generated :=

################################################################
# library creation flag
################################################################

$(eval $(library))

$(eval $(end-module))
