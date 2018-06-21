include config/config-intel.mk

# use Intel compilers
FC := mpifort
CC := mpicc

FFLAGS += -ip
F90FLAGS += -ip

# keep binaries separate by target architecture
## install_prefix := install/$(CRAY_CPU_TARGET)
