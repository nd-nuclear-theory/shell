include config/config-gcc.mk

# use Cray programming environment compilers
FC := ftn
CXX := CC

# keep binaries separate by target architecture
install_prefix := install/$(CRAY_CPU_TARGET)