################################################################
# project name
################################################################

project_name := shellutils

################################################################
# modules -- list of directories in which to search 
# for module.mk include files
################################################################

# libraries

# Caution: Order is important since used also in linking.  Caller must
# precede callee.

modules := libraries/moshinsky libraries/relative libraries/legacy

# additional libraries -- static copy (perhaps to be converted to submodule)
modules += libraries/mcpp

# additional libraries -- cloned as submodule
#
# Dependency: basis <- am
modules += libraries/basis libraries/am libraries/cppformat


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

# program algorithm choices
#   map vs. hash for space lookup in basis library
CPPFLAGS += -DINDEXING_HASH

# debugging mode
CXXFLAGS += -g
