$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h := no_filename_retrieval
module_units_cpp-h := mfdn_file_processing

# module_units_f := 
module_programs_cpp := nogen noradial noxform
## module_programs_cpp += h2gen_no h2xform_no  # TODO

# omit for now: noradial_iter noxform_iter

# module_programs_f :=
# module_generated :=

################################################################
# library creation flag
################################################################

$(eval $(library))

$(eval $(end-module))
