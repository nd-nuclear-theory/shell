# generic gnu+ompi base definitions
#
# need to augment with MKL, ACML, or LAPACK
#
# 11/20/16 (mac): Created.

# fortran
FC := mpif90
FFLAGS := -fopenmp 
FFLAGS += -O3 -ffast-math -funroll-loops

# C
CC := mpicc
CFLAGS := $(FFLAGS)
