# generic intel+ompi base definitions
#
# need to augment with MKL, ACML, or LAPACK
#
# 11/20/16 (mac): Created.

# fortran
FC := mpif90
FFLAGS := -qopenmp 
FFLAGS += -O3

# C
CC := mpicc
CFLAGS := $(FFLAGS)
