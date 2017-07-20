include config/config-gnu.mk

# use Cray programming environment compilers
FC := ftn
CC := cc


# keep binaries separate by target architecture
## install_prefix := install/$(CRAY_CPU_TARGET)

# use MKL
LDLIBS = -mkl
