include config/config-intel.mk

# use Cray programming environment compilers
FC := ftn
CC := cc


# keep binaries separate by target architecture
## install_prefix := install/$(CRAY_CPU_TARGET)
