""" ncsm_config.py

  Originated as auxiliary file to ncsm.py
  5/17/15 (mac): Modify to rely upon environment variables.
  11/18/15 (mac): Add logging of QSUBM environment variables.

  Last modified 5/17/15 (mac).
"""

import os


# Note: Use get method, rather than indexing, to allow None to be
# quietly used if the directory is not defined in the environment.


# directory for radial integral files
# 
# These are lightweight files which can easily be hosted in a home
# filesystem.

data_dir_radial = os.environ.get("QSUBM_NCSM_DATA_DIR_RADIAL")
print("QSUBM_NCSM_DATA_DIR_RADIAL",data_dir_radial)

# directory for partitioning files
# 
# Partitioning files are not strictly *necessary* for running MFDn,
# but they are critical for optimizing large MFDn runs.  These are
# lightweight files which can easily be hosted in a home filesystem.

data_dir_partitioning = os.environ.get("QSUBM_NCSM_DATA_DIR_PARTITIONING")
print("QSUBM_NCSM_DATA_DIR_PARTITIONING",data_dir_partitioning)

# directory for h2 files
# 
# This is where TBME files will be sought by h2utils.  Specifically,
# h2utils will search in named *subdirectories* of this directory.
# These are heavyweight files, but ones which should persist for later
# use.  They might be stored in scratch filesystems (only in clusters where
# scratch systems are not automatically purged) or in a shared project
# filespace (where available).

data_dir_h2 = os.environ.get("QSUBM_NCSM_DATA_DIR_H2")
print("QSUBM_NCSM_DATA_DIR_H2",data_dir_h2)

# directory for saving new archives of h2 files
#
# May typically be the same as data_dir_h2, or might be on scratch space.

data_dir_h2_archive = os.environ.get("QSUBM_NCSM_DATA_DIR_H2_ARCHIVE")
print("QSUBM_NCSM_DATA_DIR_H2_ARCHIVE",data_dir_h2_archive)

# directory for saving new archives of results from MFDn
#
# These range from lightweight to midweight, depending on the size of
# the one-body density output files.  These might be stored in a home
# directory system (to insure persistence) or a scratch system (if
# care is taken to archive them to persistent storage before they are
# automatically purged), or perhaps even shared project space (so
# collaborators can access them).

data_dir_results_archive = os.environ.get("QSUBM_NCSM_DATA_DIR_RESULTS_ARCHIVE")
print("QSUBM_NCSM_DATA_DIR_RESULTS_ARCHIVE",data_dir_results_archive)

# MFDn executable directory
#
# This is generally the base directory.  Specific versions of MFDn
# (e.g., "version14beta06") will be in subdirectories.

exec_dir_mfdn = os.environ.get("QSUBM_NCSM_EXEC_DIR_MFDN")
print("QSUBM_NCSM_EXEC_DIR_MFDN",exec_dir_mfdn)

# h2utils executable directory & shortcuts
#
# This is the directory where the executables themselves will be found.

exec_dir_h2utils = os.environ.get("QSUBM_NCSM_EXEC_DIR_H2UTILS")
if (exec_dir_h2utils is not None):
    h2gen = os.path.join(exec_dir_h2utils,"h2gen")
    h2xform = os.path.join(exec_dir_h2utils,"h2xform")
    h2mixer = os.path.join(exec_dir_h2utils,"h2mixer")

# emcalc executable directory & shortcuts
#
# This is the directory where the executables themselves will be found.

exec_dir_emcalc = os.environ.get("QSUBM_NCSM_EXEC_DIR_emcalc")
if (exec_dir_emcalc is not None):
    emcalc = os.path.join(exec_dir_emcalc,"emcalc")


# old code (for temporary reference)
## if (os.environ["HOST"] == "mac03"):  
##     # mac03 Cygwin
##     data_dir_radial = os.path.join(os.environ["HOME"],"research","data","radial-basis")
##     data_dir_h2 = os.path.join(os.environ["HOME"],"research","data","h2")
## elif ("NERSC_HOST" in os.environ):
##     # NERSC
##     data_dir_radial = os.path.join(os.environ["HOME"],"data","radial-basis")
##     data_dir_partitioning = os.path.join(os.environ["HOME"],"data","partitioning")
##     data_dir_h2 = os.path.join(os.environ["m94"],"personal","mcaprio","data", "h2")
##     data_dir_h2_archive = os.path.join(os.environ["m94"],"personal","mcaprio","data", "h2")
##     if (os.environ["NERSC_HOST"] == "hopper"):
##         mfdn_dir = os.path.join(os.environ["HOME"],"mfdn","hopper")
##     elif (os.environ["NERSC_HOST"] == "edison"):
##         mfdn_dir = os.path.join(os.environ["HOME"],"mfdn","edison")
##     data_dir_results_archive = os.path.join(os.environ["HOME"],"results")
## elif (("SGE_CELL" in os.environ) and (os.environ["SGE_CELL"] == "crc")):
##     # ND CRC
##     data_dir_radial = os.path.join(os.environ["SCRATCH"],"data","radial-basis")
##     data_dir_partitioning = os.path.join(os.environ["HOME"],"data","partitioning")
##     data_dir_h2 = os.path.join(os.environ["SCRATCH"],"data","h2")
##     data_dir_h2_archive = os.path.join(os.environ["SCRATCH"],"data","h2")
##     mfdn_dir = os.path.join(os.environ["HOME"],"mfdn")
##     data_dir_results_archive = os.path.join(os.environ["HOME"],"archive","results")

