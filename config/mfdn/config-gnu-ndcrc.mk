# 11/20/16 (mac): Based on ACML linking notes at NDCRC 9/22/16
#   for use with module ompi/1.10.2-gcc-4.9.2.

include config/config-gnu.mk

LDLIBS += /afs/crc.nd.edu/x86_64_linux/scilib/acml/current/gfortran/gfortran64/lib/libacml.a
