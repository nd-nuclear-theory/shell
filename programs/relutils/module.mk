$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h := 
module_units_cpp-h :=
# module_units_f := 
module_programs_cpp := relative-gen jpv2rel jpvstat
module_programs_cpp += moshinsky 

# module_programs_f :=
# module_generated :=

################################################################
# library creation flag
################################################################

## $(eval $(library))

$(eval $(end-module))
