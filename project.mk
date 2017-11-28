################################################################
# project name
################################################################

# This name is used for the distribution tar file.

project_name := shell

install_prefix := $(install_prefix)/shell

$(eval $(vcs-git))

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
modules += programs/obscalc

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

# gsl
LDLIBS += -lgsl
LDLIBS += -lgslcblas
CPPFLAGS += -DHAVE_INLINE

# basis submodule
#   map vs. hash for space lookup in basis library
CPPFLAGS += -DBASIS_HASH

# mcutils submodule
#   allow legacy global access to variables now wrapped in mcutils namespace
CPPFLAGS += -DMCUTILS_ALLOW_LEGACY_GLOBAL

# spline submodule
#   disable integration routines requiring later versions of gsl
CPPFLAGS += -DSPLINE_NO_FANCY_INTEGRATION
