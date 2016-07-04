$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h := arithmetic vector_tuple memoizer profiling
module_units_cpp-h :=
# module_units_f := 
module_programs_cpp := test_arithmetic test_vector_tuple test_memoizer
# module_programs_f :=
# module_generated :=

################################################################
# library creation flag
################################################################

$(eval $(library))

$(eval $(end-module))
