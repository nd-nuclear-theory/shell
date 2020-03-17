$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h := 
module_units_cpp-h :=
# module_units_f := 
module_programs_cpp := jpv2rel jpvstat relative-gen relative-filter
module_programs_cpp += moshinsky 
module_programs_cpp += pn2rel

# module_programs_f :=
# module_generated :=

################################################################
# library creation flag
################################################################

## $(eval $(library))

$(eval $(end-module))
