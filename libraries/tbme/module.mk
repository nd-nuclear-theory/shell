$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h :=
module_units_cpp-h := oscillator_parameters h2_io tbme_separable tbme_xform tbme_mapping
# module_units_f :=
module_programs_cpp_test := tbme_separable_test tbme_mapping_test
## module_programs_cpp += moshinsky_bracket_table

# module_programs_f :=
# module_generated :=

module_extras := h2_io_v0.cpp h2_io_v15099.cpp

################################################################
# library creation flag
################################################################

$(eval $(library))

################################################################
# custom dependencies
################################################################

h2_io: h2_io_v0.cpp h2_io_v15099.cpp


$(eval $(end-module))
