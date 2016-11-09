################################################################
# project name
################################################################

project_name := shellutils

################################################################
# modules -- list of directories in which to search
# for module.mk include files
################################################################

# programs

modules += programs/h2utils programs/radialutils

# legacy programs -- DEPRECATED
## modules += programs/h2utils_legacy

# libraries

# Caution: Order is important since used also in linking.  Caller must
# precede callee.

modules += libraries/moshinsky libraries/relative libraries/tbme
modules += libraries/radial

# legacy libraries -- DEPRECATED
## modules += libraries/no libraries/legacy

# additional libraries -- cloned as submodule
#
# Dependency: basis depends on am and mcutils
modules += libraries/basis libraries/am libraries/mcutils libraries/spline libraries/cppformat

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
