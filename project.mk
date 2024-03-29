################################################################
# project name
################################################################

# This name is used for the distribution tar file.

project_name := shell

install_prefix := $(install_prefix)/shell

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

modules += programs/radialutils programs/relutils programs/h2utils
modules += programs/obutils

# legacy programs -- DEPRECATED
##modules += programs/h2utils_legacy

################
# libraries
################

modules += contrib/Daejeon16
modules += libraries/relative libraries/moshinsky
modules += libraries/tbme libraries/obme libraries/density libraries/analytic
modules += libraries/quadrature

# legacy libraries -- DEPRECATED
##modules += libraries/no libraries/legacy

# additional libraries -- cloned as submodule

modules += libraries/basis  # ordering note: basis depends on am and mcutils
modules += libraries/am libraries/mcutils libraries/spline libraries/fmt

################################################################
# extras -- list of extra files to be included
# in distribution tar file
################################################################

extras :=

################################################################
# additional project-specific make settings and rules
################################################################

# contrib
src_library_dir += ./libraries ./contrib

# gsl
LDLIBS += -lgsl
LDLIBS += -lgslcblas
CPPFLAGS += -DHAVE_INLINE

# compile git information into executables
CPPFLAGS += -D'VCS_REVISION="$(vcs-git)"'

# basis submodule
#   map vs. hash for space lookup in basis library
CPPFLAGS += -DBASIS_HASH

# spline submodule
#   disable integration routines requiring later versions of gsl
CPPFLAGS += -DSPLINE_NO_FANCY_INTEGRATION
