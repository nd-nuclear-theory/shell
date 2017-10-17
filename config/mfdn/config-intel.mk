# generic intel+impi base definitions
#
# need to augment with MKL, ACML, or LAPACK
#
# 11/20/16 (mac): Created.
# 06/13/17 (pjf): Switch to Intel MPI (included with ifort).
# 10/17/17 (pjf): Add DUCTTAPE_POSTPROCESSOR flag

# fortran
FC := mpiifort
FFLAGS := -qopenmp
FFLAGS += -O3
FFLAGS += -cpp


ifdef DUCTTAPE_POSTPROCESSOR
  FFLAGS += -DDUCTTAPE_POSTPROCESSOR
endif

F90FLAGS := $(FFLAGS)

# C
CC := mpiicc
CFLAGS := $(FFLAGS)

# linking
LDLIBS = -mkl
