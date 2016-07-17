################################
# project configuration
################################

libname = mcpp

# modules -- header-only
modules_h = arithmetic vector_tuple memoizer profiling

# modules -- header-plus-object 
modules_ho = parsing

# programs
programs = test_arithmetic test_vector_tuple 
programs += test_memoizer
CC := $(CXX)

################################
# common definitions
################################

COMMON_MAKE_DIR ?= .
include $(COMMON_MAKE_DIR)/common.mk

################################
# special dependencies
################################

# program linking
CC := $(CXX)

# external libraries
LDFLAGS += -L../am
LDLIBS += -lam

CXXFLAGS += -std=c++11
