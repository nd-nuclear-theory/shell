$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h := 
module_units_cpp-h := h2_io oscillator_parameters
# module_units_f := 
module_programs_cpp :=
## module_programs_cpp += moshinsky_bracket_table

# module_programs_f :=
# module_generated :=

################################################################
# library creation flag
################################################################

$(eval $(library))

$(eval $(end-module))
