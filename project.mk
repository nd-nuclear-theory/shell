################################################################
# project name
################################################################

project_name := shellutils

################################################################
# modules -- list of directories in which to search 
# for module.mk include files
################################################################

# programs

modules += programs/h2utils

# libraries

# Caution: Order is important since used also in linking.  Caller must
# precede callee.

modules += libraries/moshinsky libraries/relative libraries/tbme libraries/no

modules += libraries/legacy  # DEPRECATED

# additional libraries -- cloned as submodule
#
# Dependency: basis <- am
modules += libraries/basis libraries/am libraries/mcutils libraries/spline libraries/cppformat


#programs
## modules += programs/interactions

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



# debugging mode
CXXFLAGS += -g
