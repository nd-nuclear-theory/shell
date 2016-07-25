$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h := 
module_units_cpp-h := construct_relative
# module_units_f := 
module_programs_cpp := writerel jpv2rel

# module_programs_f :=
# module_generated :=

################################################################
# library creation flag
################################################################

$(eval $(library))

$(eval $(end-module))
