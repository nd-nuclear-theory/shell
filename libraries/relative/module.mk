$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h := lenpic_constants
module_units_cpp-h := relative_me relative_xform jpv_io pn_io
module_units_cpp-h += lenpic_relcm_me lenpic_relative_me
# module_units_f :=
module_programs_cpp_test := relative_me_test

# module_programs_f :=
# module_generated :=

################################################################
# library creation flag
################################################################

$(eval $(library))

$(eval $(end-module))
