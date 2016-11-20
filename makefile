################################################################
#
# makefile -- hybrid C++/FORTRAN project makefile
#
# Language: GNU make
#
# Created by Mark A. Caprio, University of Notre Dame, for SU3NCSM project.
# 2/24/11 (mac): Originated.
# 3/10/11 (mac): Initial development completed.
# 4/25/11 (mac): Addition of search_dirs_include and search_dirs_lib
#   configuration variables.  Addition of optional environment variable 
#   HYBRID_MAKE_DIR for config.mk location.
# 5/02/11 (tdyt): Launch prepare_install_directory script during install.
# 11/5/11 (mac): Add installation script hook install_script.
# 7/4/13 (mac): Dereference symlinks when creating distribution tarball.
# 9/16/16 (mac):
#   - Add defaults for config.mk variables.
#   - Add "make install_lib" to optionally also install libraries.
# 11/20/16 (mac): Restore mkdir in install.
################################################################

################################################################
# configuration include files
################################################################
#
# config.mk -- describes the external library and compiler configuration
# 
# This file contains the information which is likely to be different
# for each machine on which the project is built.
#
# --------
#
# Variables to set up in config.mk:
#
#   search_prefix -- list of one or more directory trees into which
#     external libraries have been installed, e.g., $HOME/local,
#     instead of or in addition to the compiler's default paths (e.g.,
#     /usr/local)
#
#   search_dirs_include -- additional directory to search for include
#     files *only* (useful if a source library is installed with a
#     nonconventional tree structure)
#
#   search_dirs_lib -- additional directory to search for library
#     files *only* (useful if a source library is installed with a
#     nonconventional tree structure)
#
#   fortran_libs -- the libraries required for a C++ project to
#     properly link when including FORTRAN objects
#
#   fortran_flags -- linking flags required for a C++ project to 
#     include FORTRAN objects
#
#   MPICXX -- (optional) the MPI C++ command, if your project uses MPI
#     (DEPRECATED)
#
# These are all initialized to a null string, so config.mk can
# append to them (with +=) if preferred.
#
#   install_

#
# --------
#
# This file would also normally define the compile commands to use
# (FC and CXX) and any machine specific flags for these compilers.
#
################################################################
#
# project.mk -- defines the contents of the project 
#
# This file contains the information on the overall structure of the
# project and compiler/linking configuration for the whole project
# (but is machine independent).
#
# --------
#
# Variables to define in project.mk:
#
# project_name -- project name, e.g., for use in the tar file name
#   when make produces a source tar distribution
# modules -- list of "modules", i.e., subdirectories each of which 
#   a module.mk include file
# extras -- list of extra files or directories to be bundled in the 
#   source tar distribution (e.g., README)
# install_prefix (default: ".") -- install prefix, i.e., binaries installed to 
#   install_prefix/bin (e.g., "$(current-dir)/.." for install from
#   project_dir/src to project_dir/bin) 
# install_script -- script to be run during "make install", after binaries 
#   have been installed to binary directory
# list-mpi-programs-cpp -- (optional) *function* which determines which 
#   program files should be compiled or linked with MPI; will be 
#   evaluated *after* all module.mk files have been read 
# list-mpi-objects-cpp -- (optional) *function* which determines which 
#   C++ object files should be compiled with MPI; will be evaluated 
#   *after* all module.mk files have been read 
#
# --------
#
# This file would also normally define various compiler and linker
# flags common to the whole project.
#
################################################################
#
# module.mk (for each project "module", i.e., subdirectory):
#
# --------
#
# File should begin with $(eval $(begin-module))
#
# --------
#
# Variables defining the compilation units in this module (and other
# date files to be generated):
#
# module_units_h -- units consisting of file.h only
# module_units_cpp-h -- units consisting of file.cpp + file.h 
#   for compilation to file.o and inclusion in library
# module_units_f -- units consisting of file.f only
# module_programs_cpp -- programs consisting of file.cpp to be compiled and linked
# module_programs_f -- programs consisting of file.f to be compiled and linked
# module_generated -- generated files created by rules defined in module.mk files 
#
# These variables also determine the *source* files to be included in the source tar.  
#
# --------
#
# $(eval $(library)) -- should be included if the object files are to
#   be assembled into a library archive (named after the module
#   directory)
# 
# --------
#
# This file would also normally define any dependencies and target-specific
# variable assignments.  
#
# Notes:
# (1) All filenames must be relative to the top of the source tree, 
#   i.e., SU3NCSM/src.
# (2) The current directory for this module can be referenced
#   as $(current-dir), in the target and prerequisite names.
# (3) Predefined prerequisites:
#    All object files from units_cpp-h are already dependent on both 
#      their .cpp and .h file.
#    All object files from units_f are already dependent on 
#      their .f file.
#
# EX: $(current-dir)/myunit.o : $(current-dir)/anotherunit.h libraries/somelibrary/somelibrary.h
#    Note $(current-dir)/myunit.h already assumed by makefile.
#
# (4) However, $(current-dir) should not appear in the rule to build the target
#    since evaluation of the rules is delayed, and $(current-dir) will not longer
#    hold the correct directory name at run time.
#
# EX: $(current-dir)/data-file :  $(current-dir)/data-generator
#     [TAB] cd $(dir $@); ./data-generator
#    Note we have used $(dir $@) rather than the (inappropriate) $(current-dir)
#    to recover the module's directory in the build command.
#
# --------
#
# File should end with $(eval $(end-module))
#
################################################################
#
# Environment variables
#
# HYBRID_MAKE_DIR (optional) -- If set, config.mk will be sought
# in this directory instead of the current working directory.
#
################################################################

################################################################
# declare default target
################################################################

.PHONY: splash
splash:

################################################################
################################################################
#
# macro definitions
#
################################################################
################################################################

################################################################
# string and path utilities
################################################################

#($call last-word,list)
#  returns last word of list
define last-word
$(word $(words $1),$1)
endef

#($call strip-trailing-slash,list)
#  strips any trailing slash from each word of list
#  after Mecklenburg p. 68
define strip-trailing-slash
$(patsubst %/,%,$1)
endef

#($call sandwich,prefix,list,suffix)
#  adds both prefix and suffix to each entry in list
define sandwich
$(patsubst %,$(1)%$(3),$(2))
endef

#$(current-dir)
#  returns subdirectory from which project makefile was invoked
#  relative to make directory
#  after Mecklenburg p. 142
define current-dir
$(strip $(call strip-trailing-slash,$(dir $(call last-word,$(MAKEFILE_LIST)))))
endef

#$(get-exe-ext)
#  returns system-specific executable extension
#  assuming COMSPEC is defined only under MS Windows
ifdef COMSPEC
  get-exe-ext := .exe
else
  get-exe-ext :=
endif

################################################################
# accumulation
################################################################

#$(eval $(begin-module))
#  Clear all module "local" variables
#  Note: eval needed for multi-line macro at top level of file
define begin-module
  module_units_h :=
  module_units_cpp-h :=
  module_units_f :=
  module_programs_cpp :=
  module_programs_f :=
  module_generated :=
endef

#$(eval $(library))
#  creates definitions so that local module forms a library
#  and defines dependency on object files
define library
  $(eval module_library := lib$(notdir $(current-dir)))
  $(eval module_library_ar_name := $(call sandwich,$(current-dir)/,$(module_library),.a))
  $(eval module_library_units := $(module_units_cpp-h) $(module_units_f))
  $(eval module_library_objects := $(call sandwich,$(current-dir)/,$(module_library_units),.o))
  libraries += $(addprefix $(current-dir)/,$(module_library))
  library_journal += $(module_library) - $(module_library_units) ...
  $(foreach obj,$(module_library_objects),$(eval $(module_library_ar_name): $(obj)) )
endef

#$(eval $(end-module))
#   accumulate local unit information into global lists
#  Note: eval needed for multi-line macro at top level of file
define end-module
  units_h += $(addprefix $(current-dir)/,$(module_units_h))
  units_cpp-h += $(addprefix $(current-dir)/,$(module_units_cpp-h))
  units_f += $(addprefix $(current-dir)/,$(module_units_f))
  programs_cpp += $(addprefix $(current-dir)/,$(module_programs_cpp))
  programs_f += $(addprefix $(current-dir)/,$(module_programs_f))
  generated += $(addprefix $(current-dir)/,$(module_generated))
endef

# Debugging note: The dependency declaration
#   $(module_library_ar_name): $(module_library_objects)
# yields virtual memory overflow errors (GNU make 3.80 Cygwin) for long 
# lists of object files (more than ~5).  This is due to a known bug in GNU 
# make 3.80, fixed in 3.81:
#   http://stackoverflow.com/questions/2428506/workaround-for-gnu-make-3-80-eval-bug
#   When $(eval) evaluates a line that is over 193 characters, Make
#   crashes with a "Virtual Memory Exhausted" error.
# The workaround applied here is to add the dependencies one by on in a foreach 
# *without* enclosing the foreach in an eval.  


################################################################
################################################################
#
# project configuration
#
################################################################
################################################################

################################################################
# read in system-specific configuration
################################################################

search_prefix :=
search_dirs_include :=
search_dirs_lib :=
fortran_libs := 
fortran_flags :=
install_prefix := .

HYBRID_MAKE_DIR ?= .
include $(HYBRID_MAKE_DIR)/config.mk

################################################################
# read in project details
################################################################

include project.mk

################################################################
# define basic directories
################################################################

# assumes tree:
# <project_base>/src
# <project_base>/src/libraries   
# <project_base>/src/programs
# <project_base>/bin
# with relative paths taken relative to <project_base>/src

project_base := $(dir $(PWD))
bin_dir := $(project_base)/bin
src_library_dir := ./libraries
#src_programs := $(src)/programs

################################################################
# include search path setup
################################################################

# search path for external library include files
#   from config.mk
CPPFLAGS += $(call sandwich,-I,$(search_prefix),/include)
vpath %.h $(call sandwich,,$(search_prefix),/include)
CPPFLAGS += $(call sandwich,-I,$(search_dirs_include),)
vpath %.h $(call sandwich,,$(search_dirs_include),)

# search path for internal library include files
CPPFLAGS += -I$(src_library_dir)
vpath %.h $(src_library_dir)

################################################################
# library search path setup
################################################################

# search path for required external libraries
#   from config.mk
LDFLAGS += $(call sandwich,-L,$(search_prefix),/lib)
LDFLAGS += $(call sandwich,-L,$(search_dirs_lib),)

# search path for internal library archive files
#   which will be specified using libxxxx.a rather than -lxxxx syntax
vpath %.a $(src_library_dir)

################################################################
# FORTRAN linking setup
################################################################

# linking options for linking to FORTRAN
#   from config.mk
LDLIBS += $(fortran_libs)
LDFLAGS += $(fortran_flags)


################################################################
################################################################
#
# module processing
#
################################################################
################################################################

################################################################
# accumulator list declarations
################################################################

# units_h -- units consisting of file.h only
units_h :=

# units_cpp-h -- units consisting of file.cpp + file.h 
#   for compilation to file.o and inclusion in library
units_cpp-h :=

# units_f -- units consisting of file.f only
units_f :=

# libraries -- libraries consisting of file.a to be used in program linking
libraries :=

# programs_cpp -- programs consisting of file.cpp to be compiled and linked
programs_cpp :=

# programs_f -- programs consisting of file.f to be compiled and linked
programs_f :=

# generated -- generated files created by rules defined in module.mk files 
generated := 

# diagnostic accumulators
library_journal :=
debug_output := 

################################################################
# iteration over modules
################################################################

module_files := $(addsuffix /module.mk,$(modules))
include $(module_files)

################################################################
# deduced filenames
################################################################

cpp_ext := .cpp
h_ext := .h
f_ext := .F
o_ext := .o
archive_ext := .a
binary_ext := $(get-exe-ext)

sources_cpp := $(addsuffix $(cpp_ext),$(units_cpp-h) $(programs_cpp)) 
sources_h := $(addsuffix $(h_ext),$(units_cpp-h) $(units_h)) 
sources_f := $(addsuffix $(f_ext),$(units_f) $(programs_f)) 
sources := $(sources_cpp) $(sources_h) $(sources_f)

makefiles := $(MAKEFILE_LIST)

objects := $(addsuffix $(o_ext),$(units_cpp-h) $(units_f) $(programs_cpp) $(programs_f)) 

archives := $(addsuffix $(archive_ext),$(libraries))

programs := $(programs_cpp) $(programs_f)
executables := $(addsuffix $(binary_ext),$(programs))

################################################################
################################################################
#
# targets
#
################################################################
################################################################

################################################################
# rules and dependencies
################################################################

# make all executables dependent upon all project libraries
#   Note: There is then no need to add the libraries to LDLIBS,
#   since they will appear in the dependencies argument to the linker.
ifneq "$(strip $(programs))" ""
$(programs): $(archives)
endif

# make all units_cpp-h object files dependent on the header file
#   using static pattern rule
ifneq "$(strip $(units_cpp-h))" ""
$(addsuffix .o,$(units_cpp-h)): %.o: %.h
endif

# object library rule
ARFLAGS = r
%.a:
	$(RM) $@           
	$(AR) $(ARFLAGS) $@ $^

################################################################
# linker
################################################################

# link C++ programs using C++ compiler
ifneq "$(strip $(programs_cpp))" ""
$(programs_cpp): CC := $(CXX)
endif

# link FORTRAN programs using FORTRAN compiler
ifneq "$(strip $(programs_f))" ""
$(programs_f): CC=$(FC)
endif

# compile and link MPI C++ programs using MPICXX wrapper
#   Functions $(list-mpi-programs-cpp) and $(list-mpi-objects-cpp) for deciding which 
#   files need MPI for this project must be defined in project.mk.
#   EX: 
#     list-mpi-programs-cpp = $(filter %MPI,$(programs_cpp))
#     list-mpi-objects-cpp = $(addsuffix $(o_ext),$(list-mpi-programs-cpp))
#   Note the use of *delayed* assignment (=), since programs_cpp
#   has not been populated yet when project.mk is included.

ifneq "$(strip $(list-mpi-programs-cpp))" ""
$(list-mpi-objects-cpp): CXX := $(MPICXX)
$(list-mpi-programs-cpp): CXX := $(MPICXX)
$(list-mpi-programs-cpp): CC := $(MPICXX)
endif

################################################################
# shorthand library/program/generated targets
################################################################
# allow target to be specified without qualifying path
#   e.g., "make libfoo" instead of "make libraries/foo/libfoo.a"

$(foreach target,$(libraries),$(eval .PHONY : $(notdir $(target))) )
$(foreach target,$(libraries),$(eval $(notdir $(target)) : $(addsuffix .a,$(target)) ))

$(foreach target,$(programs),$(eval .PHONY : $(notdir $(target))) )
$(foreach target,$(programs),$(eval $(notdir $(target)) : $(target)) )

$(foreach target,$(generated),$(eval .PHONY : $(notdir $(target))) )
$(foreach target,$(generated),$(eval $(notdir $(target)) : $(target)) )


################################################################
# diagnostic output target
################################################################

splash:
	@echo "makefile -- hybrid C++/FORTRAN project                     "
	@echo "M. A. Caprio, University of Notre Dame			  "
	@echo
	@echo $(project_name) make information
	@echo
	@echo "Working directory:" $(PWD)
	@echo
	@echo "Modules:" $(modules)
	@echo
	@echo "Libraries:" $(notdir $(libraries))
	@echo
	@echo $(library_journal)
	@echo
	@echo "Programs:" $(notdir $(programs))
	@echo
	@echo "Generated:" $(notdir $(generated))
	@echo
	@echo $(debug_output)
	@echo "To build the project, run \"make all\"."
	@echo "Or for further instructions, run \"make help\"."

################################################################
# help target
################################################################

.PHONY: help
help:
	@echo
	@echo "makefile -- hybrid C++/FORTRAN project                     "
	@echo "M. A. Caprio, University of Notre Dame			  "
	@echo "								  "
	@echo "Syntax:							  "
	@echo "  make <target>						  "
	@echo "								  "
	@echo "Targets:							  "
	@echo "  (none) -- a summary of the project is displayed	  "
	@echo "  all  -- all libraries, programs, and generated files 	  "
	@echo "  libraries -- just the libraries			  "
	@echo "  programs -- just the programs				  "
	@echo "  generated -- just the generated files (e.g., data files) "
	@echo "  distrib [tag=<version>] -- make source tarball		  "
	@echo "    Default tag is YYMMDD date.				  "
	@echo "  clean -- delete binaries and generated files		  "
	@echo "								  "
	@echo "  <library> -- shorthand for full path to library	  "
	@echo "    EX: Use shorthand target libmylib for xxxx/xxxx/libmylib.a.	  "
	@echo "  <program> -- shorthand for full path to program	  "
	@echo "    EX: Use shorthand target myprog for xxxx/xxxx/myprog.            "

# ending message 

.PHONY: finished
finished:
# Displayed since otherwise it is confusing...  Otherwise, if make all
# is invoked with all targets already made, successful termination
# would leave us with the help message `To build the project, run
# "make all"', which looks like a failure, rather than an indication
# that the project was actually already successfully built.
	@echo "								  "
	@echo "Full-project build completed.                              "

################################################################
# general targets
################################################################

.PHONY: all
all: splash libraries programs generated finished

.PHONY: libraries
libraries: $(archives)

.PHONY: programs 
programs: $(programs)

.PHONY: generated
generated: $(generated)

################################################################
# install
################################################################

install_dir_bin := $(install_prefix)/bin
install_dir_include := $(install_prefix)/include
install_dir_lib := $(install_prefix)/lib
MKDIR := mkdir --parents 

# 11/20/16 (mac): The "mkdir --parents" *should* be redundant to
# "install -D", but it seems to be necessary on the ndcrc.

.PHONY: install_bin
install_bin: programs
	@echo Installing binaries to $(install_dir_bin)...
	$(MKDIR) $(install_dir_bin)
	install -D $(executables) --target-directory=$(install_dir_bin)

.PHONY: install_include
install_include: ${sources_h}
	@echo Installing includes to $(install_dir_include)...
	@echo WARNING: not yet supported
##	$(MKDIR) $(install_dir_lib)
##	install -D ${sources_h} --target-directory=$(install_dir_lib)
##	@ $(foreach source,$(sources_h),echo $(source); )


.PHONY: install_lib
install_lib: libraries
	@echo Installing libraries to $(install_dir_lib)...
	$(MKDIR) $(install_dir_lib)
	install -D $(archives) --target-directory=$(install_dir_lib) --mode=u=rw,go=r

.PHONY: install
##install: install_bin install_include install_lib
install: install_bin
	$(install_script) 


################################################################
# source tarball
################################################################

# construct last name of current directory (e.g., "src")
pwd_tail := $(notdir $(CURDIR))

# construct tar filename

tag ?= $(shell date +%y%m%d)
tarball = $(project_name)-$(tag).tgz

# construct list of items to include in distribution
tar_constituents = $(sources) $(makefiles) $(extras)

# Rule for tarball
# To put in source directory instead of parent: $(pwd_tail)/$(tarball)

.PHONY: distrib
distrib: 
	@ echo Making source tarball $(tarball)...
	@ cd ..; tar --dereference -zcvf $(tarball) $(addprefix $(pwd_tail)/,$(tar_constituents))
	@ ls -Fla ../$(tarball)

# Make alias "distribution"
.PHONY: distribution
distribution: distrib 

################################################################
# cleanup
################################################################

.PHONY: clean
clean: 
	$(RM) $(objects) $(archives) $(executables) $(generated)

.PHONY: distclean
distclean: clean
