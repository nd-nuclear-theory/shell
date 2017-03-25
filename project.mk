################################################################
# project name
################################################################

# This name is used for the distribution tar file.

project_name := shell

################################################################
# modules -- list of directories in which to search
# for module.mk include files
################################################################

# Caution: The order in which modules are declares is important, since
# it is also used in linking.  The object/library file for the
# *caller* must precede the object/library file for the *callee* for
# successful linking.

################
# programs
################

modules += programs/h2utils programs/radialutils

# legacy programs -- DEPRECATED
##modules += programs/h2utils_legacy

################
# libraries
################

modules += libraries/relative libraries/moshinsky
modules += libraries/tbme libraries/radial libraries/density

# legacy libraries -- DEPRECATED
##modules += libraries/no libraries/legacy

# additional libraries -- cloned as submodule

modules += libraries/basis  # ordering note: basis depends on am and mcutils
modules += libraries/am libraries/mcutils libraries/spline libraries/cppformat

################################################################
# extras -- list of extra files to be included
# in distribution tar file
################################################################

extras :=

################################################################
# additional project-specific make settings and rules
################################################################

# Gnu Scientific Library
LDLIBS += -lgsl
LDLIBS += -lgslcblas
CPPFLAGS += -DHAVE_INLINE

# spline submodule
CPPFLAGS += -DSPLINE_NO_FANCY_INTEGRATION

# program algorithm choices
#   map vs. hash for space lookup in basis library
CPPFLAGS += -DINDEXING_HASH

# optimization mode
CPPFLAGS += -O3

# debugging mode
CXXFLAGS += -g
