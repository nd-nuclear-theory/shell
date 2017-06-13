# generic intel+impi base definitions
#
# need to augment with MKL, ACML, or LAPACK
#
# 11/20/16 (mac): Created.
# 06/13/17 (pjf): Switch to Intel MPI (included with ifort).

# fortran
FC := mpiifort
FFLAGS := -qopenmp
FFLAGS += -O3

# C
CC := mpiicc
CFLAGS := $(FFLAGS)
