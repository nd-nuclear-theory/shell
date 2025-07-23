$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h :=
module_units_cpp-h :=
# module_units_f :=
module_programs_cpp := smwf-convert smwf-truncate
module_programs_cpp_test := # group_read introduces parallel_hashmap dependency

# module_programs_f :=
# module_generated :=

$(eval $(end-module))
