$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h := readWavefunction
module_units_cpp-h :=
# module_units_f :=
module_programs_cpp :=  nomixer
module_programs_cpp += nomixer_V0
module_programs_cpp += readWavefunction
# module_programs_cpp_test += 

# module_programs_f :=
# module_generated :=

################################################################
# library creation flag
################################################################

## $(eval $(library))

$(eval $(end-module))
