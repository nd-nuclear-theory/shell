include config/config-gnu.mk

FFLAGS += -L$(HOME)/local/opt/OpenBLAS/lib  # temporary path

LDLIBS += -lopenblas
