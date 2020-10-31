$(eval $(begin-module))

################################################################
# variable definitions
################################################################

CPPFLAGS += -DUSE_DAEJEON16

################################################################
# unit definitions
################################################################

# module_units_h :=
module_units_cpp-h := Daejeon16_wrapper
module_units_f := Daejeon16_public_v2
# module_programs_cpp_test :=

# module_programs_f :=
# module_generated :=

################################################################
# library creation flag
################################################################

$(eval $(library))

$(eval $(end-module))
