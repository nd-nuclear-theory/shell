$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h := 
module_units_cpp-h := oscillator_parameters h2_io tbme_separable tbme_xform two_body_mapping
# module_units_f := 
module_programs_cpp := tbme_separable_test two_body_mapping_test
## module_programs_cpp += moshinsky_bracket_table

# module_programs_f :=
# module_generated :=

################################################################
# library creation flag
################################################################

$(eval $(library))

$(eval $(end-module))
