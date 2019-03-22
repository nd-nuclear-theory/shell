$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h :=
module_units_cpp-h :=
# module_units_f :=
module_programs_cpp := orbital-gen radial-gen
module_programs_cpp += radial-xform natorb-gen radial-compose
module_programs_cpp += radial-extract
# module_programs_cpp += radial-scale radial-stat

# module_programs_f :=
# module_generated :=

################################################################
# library creation flag
################################################################

## $(eval $(library))

$(eval $(end-module))
