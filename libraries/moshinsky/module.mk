$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h := 
module_units_cpp-h := moshinsky_bracket
# module_units_f := 
module_programs_cpp := moshinsky_test
# module_programs_cpp += moshinsky_bracket_table moshinsky

# module_programs_f :=
# module_generated :=

################################################################
# library creation flag
################################################################

$(eval $(library))

$(eval $(end-module))
