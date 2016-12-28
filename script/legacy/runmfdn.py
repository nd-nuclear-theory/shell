#!/usr/bin/env python3

"""runmfdn00.py -- running h2mixer + MFDn

    Example runs of h2mixer (for setup) followed by MFDn, for 4He.

    Instructions are provided below:

        Setup

        Phase 0

        Phase 1

        Common problems

    ----------------------------------------------------------------

    Setup

    Prerequisites: You should already have installed and built the
    h2utils codes before you attempt this example.  You should also
    have installed and built MFDn.  This example currently assumes you
    are running MFDn Version 14 Beta 06, but the scripting can/will be
    updated to accommodate other versions.

    Environment: You will first need to set several environment
    variables, to identify directories where various data and
    executable files are found.  These are documented in
    mfdn_config.py.  The basic setup might look something like:

        setenv QSUBM_NCSM_EXEC_DIR_MFDN ${HOME}/codes/${NERSC_HOST}/mfdn
        setenv QSUBM_NCSM_EXEC_DIR_H2UTILS ${HOME}/codes/${NERSC_HOST}/h2utils/programs/h2utils
        setenv QSUBM_NCSM_DATA_DIR_H2 ${SCRATCH}/data/h2
        setenv QSUBM_NCSM_DATA_DIR_PARTITIONING ${HOME}/data/partitioning
        setenv QSUBM_NCSM_DATA_DIR_RESULTS_ARCHIVE ${HOME}/results


    Example TBME input files: This example requires certain two-body
    matrix element input files.  They are contained in
    mfdn-example-data.tgz.  Untar this archive into the directory you
    identified above as QSUBM_NCSM_DATA_DIR_H2.

        cd ${QSUBM_NCSM_DATA_DIR_H2}
        tar xvf mfdn-example-data.tgz


    Taking a look at the tasks defined by this example run: You can
    try this out on the front end:

       qsubm mfdn00 --toc

    You should see an output like the following:

         Run: runmfdn00
         Thu Jun 25 17:40:46 2015
         Tasks: 10
         Archive phases: 2
         0000 Nmax02 - - Z2-N2-JISP16-1-hw20.000-aL20-Nmax02-Mj0.0-lan500
         0001 Nmax02 - - Z2-N2-JISP16-1-hw25.000-aL20-Nmax02-Mj0.0-lan500
         0002 Nmax04 - - Z2-N2-JISP16-1-hw20.000-aL20-Nmax04-Mj0.0-lan500
         0003 Nmax04 - - Z2-N2-JISP16-1-hw25.000-aL20-Nmax04-Mj0.0-lan500
         0004 Nmax06 - - Z2-N2-JISP16-1-hw20.000-aL20-Nmax06-Mj0.0-lan500
         0005 Nmax06 - - Z2-N2-JISP16-1-hw25.000-aL20-Nmax06-Mj0.0-lan500
         0006 Nmax08 - - Z2-N2-JISP16-1-hw20.000-aL20-Nmax08-Mj0.0-lan500
         0007 Nmax08 - - Z2-N2-JISP16-1-hw25.000-aL20-Nmax08-Mj0.0-lan500
         0008 Nmax10 - - Z2-N2-JISP16-1-hw20.000-aL20-Nmax10-Mj0.0-lan500
         0009 Nmax10 - - Z2-N2-JISP16-1-hw25.000-aL20-Nmax10-Mj0.0-lan500

    Each row indicates:
 
         task pool status0 status1 descriptor

    Notice that each of the ten "tasks" has been assigned to a "pool"
    (this allows us to group similar runs by name, say, "Nmax02", and
    tell the script "I want you to execute the next available Nmax02
    run", for instance, as we shall do below).

    ----------------------------------------------------------------

    Phase 0

    This "phase" of each task is where we construct the Hamiltonian
    (and other operator) two-body matrix element (TBME) input files
    for MFDn.  This is done by the serial code h2mixer, from the
    h2utils package (it is assumed that you have installed and built
    these programs).

    You can run this phase on the front end.  Select the pool you want
    to run (well, to keep the example data files small, only example
    input files for Nmax02 have been provided):

        qsubm mfdn00 --pool=Nmax02 --phase=0

    This should run phase 0 of both tasks in the pool "Nmax02" (these two tasks
    are for hw=20 and hw=25, respectively).  You should see output
    like the following:


         ----------------------------------------------------------------
         task 0 phase 0
         Z2-N2-JISP16-1-hw20.000-aL20-Nmax02-Mj0.0-lan500
         ----------------------------------------------------------------
         Redirecting to /home/mcaprio/research/runs/scratch/runmfdn00/output/task-0000-0.out
         (Task time: 0.20 sec)

         Task timing: elapsed 0.210001, remaining 3599.79, last task 0.200001, required 0.220001

         ----------------------------------------------------------------
         task 1 phase 0
         Z2-N2-JISP16-1-hw25.000-aL20-Nmax02-Mj0.0-lan500
         ----------------------------------------------------------------
         Redirecting to /home/mcaprio/research/runs/scratch/runmfdn00/output/task-0001-0.out
         (Task time: 0.12 sec)

         Task timing: elapsed 0.330001, remaining 3599.67, last task 0.12, required 0.132

         Found no available tasks in pool Nmax02
         ----------------------------------------------------------------
         End script
         Thu Jun 25 17:54:43 2015
         ----------------------------------------------------------------

    If you look at the table of contents now, you should see that
    these phases are flagged as completed (their status will show as
    "X"):

        Run: runmfdn00
        Thu Jun 25 17:56:41 2015
        Tasks: 10
        Archive phases: 2
        0000 Nmax02 X - Z2-N2-JISP16-1-hw20.000-aL20-Nmax02-Mj0.0-lan500
        0001 Nmax02 X - Z2-N2-JISP16-1-hw25.000-aL20-Nmax02-Mj0.0-lan500
        0002 Nmax04 - - Z2-N2-JISP16-1-hw20.000-aL20-Nmax04-Mj0.0-lan500
        0003 Nmax04 - - Z2-N2-JISP16-1-hw25.000-aL20-Nmax04-Mj0.0-lan500
        0004 Nmax06 - - Z2-N2-JISP16-1-hw20.000-aL20-Nmax06-Mj0.0-lan500
        0005 Nmax06 - - Z2-N2-JISP16-1-hw25.000-aL20-Nmax06-Mj0.0-lan500
        0006 Nmax08 - - Z2-N2-JISP16-1-hw20.000-aL20-Nmax08-Mj0.0-lan500
        0007 Nmax08 - - Z2-N2-JISP16-1-hw25.000-aL20-Nmax08-Mj0.0-lan500
        0008 Nmax10 - - Z2-N2-JISP16-1-hw20.000-aL20-Nmax10-Mj0.0-lan500
        0009 Nmax10 - - Z2-N2-JISP16-1-hw25.000-aL20-Nmax10-Mj0.0-lan500

    Recovering from a mistake: Suppose, instead the setup run failed.
    Perhaps you had forgotten to set one of the required environment
    variables.  You may either see "locked" ("L") or "failed" ("F") as
    the status of the failed task:

         Run: runmfdn00
         Thu Jun 25 18:08:18 2015
         Tasks: 10
         Archive phases: 2
         0000 Nmax02 F - Z2-N2-JISP16-1-hw20.000-aL20-Nmax02-Mj0.0-lan500
         0001 Nmax02 - - Z2-N2-JISP16-1-hw25.000-aL20-Nmax02-Mj0.0-lan500
         0002 Nmax04 - - Z2-N2-JISP16-1-hw20.000-aL20-Nmax04-Mj0.0-lan500
         0003 Nmax04 - - Z2-N2-JISP16-1-hw25.000-aL20-Nmax04-Mj0.0-lan500
         0004 Nmax06 - - Z2-N2-JISP16-1-hw20.000-aL20-Nmax06-Mj0.0-lan500
         0005 Nmax06 - - Z2-N2-JISP16-1-hw25.000-aL20-Nmax06-Mj0.0-lan500
         0006 Nmax08 - - Z2-N2-JISP16-1-hw20.000-aL20-Nmax08-Mj0.0-lan500
         0007 Nmax08 - - Z2-N2-JISP16-1-hw25.000-aL20-Nmax08-Mj0.0-lan500
         0008 Nmax10 - - Z2-N2-JISP16-1-hw20.000-aL20-Nmax10-Mj0.0-lan500
         0009 Nmax10 - - Z2-N2-JISP16-1-hw25.000-aL20-Nmax10-Mj0.0-lan500

    As long as this flag is set, the task will not be eligible for
    execution.  To clear all "lock" or "fail" flags:

         qsubm mfdn00 --unlock

    Alternatively, you might want to run the setup phase in a queue
    for serial jobs.  Here is an example for NERSC:
 
        module load serial
        qsubm mfdn00 serial 10 --pool=Nmax02 --phase=0
        module unload serial

    Or on the "short" queue on Notre Dame's cluster

        qsubm mfdn00 short 10 --pool=Nmax02 --phase=0

    ----------------------------------------------------------------

    Phase 1

    This is where we run MFDn.

    Here is an example run on NERSC:

        qsubm mfdn00 debug 10 --pool=Nmax02 --phase=1 --width=15 --opt="-m ae"

    Here are some example runs at higher Nmax on the compute nodes at
    the ND CRC:

        qsubm mfdn00 short 60              --pool=Nmax04 --phase=1 --width=6 --nodesize=8 --start=0 --limit=1 --opt="-m ae"
        qsubm mfdn00 "*@@dqcneh_253GHz" 60 --pool=Nmax06 --phase=1 --width=6 --nodesize=8 --start=0 --limit=1 --opt="-m ae"
        qsubm mfdn00 "*@@ivybridge" 60     --pool=Nmax08 --phase=1 --width=15 --nodesize=16 --start=0 --limit=1 --opt="-m ae"
 

    ----------------------------------------------------------------

    Common issues

    (1) Depending how the queue submission is set up on your system,
    and depending on the qsubm "local" file you have set up, you might
    or might not need a "shebang" line like the one at the start of
    this file:

        #!/usr/bin/env python3

    If it is needed, it should be set so that the node responsible for
    running submission scripts (this may or may not be the compute
    node) can use it to find Python 3.

    (2) Programming environment module: Depending on the compilers you
    used for the MFDn executable, it might be necessary to either: (1)
    make sure you have an appropriate module loaded at the time you
    submit the job or (2) include a module load command
    (mcscript.module) as in the commented-out example line in the code
    below.  This may also be true for h2utils.


    ----------------------------------------------------------------

    Technical documentation on provenance of example data files:
  
        run0229/r2k2-HO-1.000-HO-1.000-3-4-r2.bin
        run0229/r2k2-HO-1.000-HO-1.000-3-4-r1r2.bin
        run0229/r2k2-HO-1.000-HO-1.000-3-4-k2.bin
        run0229/r2k2-HO-1.000-HO-1.000-3-4-k1k2.bin
        run0177/JISP16-HO-1.000-HO-1.000-ob-13-3-4-20.bin
        run0177/JISP16-HO-1.000-HO-1.000-ob-13-3-4-25.bin
        run0177/VC-HO-1.000-HO-1.000-ob-13-3-4-20.bin

    ----------------------------------------------------------------

    Last modified 6/30/15 (mac).

"""

import math
import sys
import os
import shutil
import subprocess

import mcscript
import ncsm
import ncsm_config
import mfdn_v14_b06

################################################################
# module loads
################################################################

## mcscript.module(["load","mpich/3.0.4-pgi"])

################################################################
# general run configuration
################################################################

ncsm.set_run_search_lists(
    radial = [
        "131114"
    ],
    r2k2 = [
        "mfdn-example-data",  # minimal data needed for this example
        "run0229",  # HO/CS including halo   
        "run0244", "run0250"  # hybrid
    ],
    xform = [
        "mfdn-example-data",  # minimal data needed for this example
        "run0177", # standard HO hw 2.5 mesh   <<<< actually needed for this example
        "run0245", "run0251"  # hybrid
    ]
)

##################################################################
# task list builder
##################################################################

# nuclide parameters
nuclide = (2,2)
Nv = 1 # target Hamiltonian valence shell 

# Hamiltonian parameters
interaction = "JISP16" 
coulomb = True
aLawson = 20
hwLawson_list = [ None ]

# basis parameters
hw_int_list = None  # None = use hw for hw_int
hw_int_coul = 20.
radial_basis_pair_list = [ ("HO","HO")  ]
traditional_ho = True
## hw_log_range = mcscript(10,40,8)  
## hw_list = mcscript.log_range(*hw_log_range)
##hw_list = mcscript.value_range(10,40,5)
hw_list = [20.,25.]
beta_ratio_log_range = (1.0,1.0,1)
beta_ratio_list = mcscript.log_range(*beta_ratio_log_range,power=0.5)
Nmax_start = 2
Nmax_limit = 10
Nstep = 2
Nmax_list = mcscript.value_range(Nmax_start,Nmax_limit,2)
xform_cutoff_list = [
    ("ob", N1b_cut)
    for N1b_cut in mcscript.value_range(13,13,2)
    ]  

# diagonalization
Mj = 0
eigenvectors = 10
initial_vector = -2
lanczos = 500
tolerance = 0

# obdme parameters
obdme_reference_state_list = [(0.0,0,1)]
obdme_multipolarity = 2 # WARNING: reducing below 2 (on HO runs) crashes MFDn

# version parameters
mfdn_executable = "version14-beta06-newmake/xmfdn-h2-lan"
mfdn_wrapper = mfdn_v14_b06.call_mfdn_h2

# generate task list
tasks = [
    {
        # nuclide parameters
        "nuclide" : nuclide,

        #Hamiltonian parameters
        "interaction" : interaction,
        "coulomb" : coulomb,
        "aLawson" : aLawson,
        "hwLawson" : mcscript.auto_value(hwLawson,hw),
        
        # basis parameters
        "hw_int" : hw_int,  # hw of source interaction
        "hw_int_coul" : hw_int_coul,  # hw of coulomb interaction
        "hw" : hw,  # logical hw, for scaling basis
        "basis" :  # logical basis, to be scaled by hw
        (
            radial_basis_p,
            1.,
            radial_basis_n,
            beta_ratio
            ),
        "scaled_basis" :  # deduced basis used in task descriptor, etc.
        (
            radial_basis_p,
            ncsm.beta_for_hw(hw_int,hw),
            radial_basis_n,
            ncsm.beta_for_hw(hw_int,hw)*beta_ratio
            ),
        "xform_cutoff" : xform_cutoff,
        "xform_cutoff_coul" : xform_cutoff,
        "Nmax" : Nmax,
        "Nstep" : Nstep,
        "Nv" : Nv,

        # diagonalization parameters
        "Mj" : Mj,
        "eigenvectors" : eigenvectors,
        "initial_vector" : initial_vector,   
        "lanczos" : lanczos,
        "tolerance" : tolerance,

        # obdme parameters
        "obdme_reference_state_list" : obdme_reference_state_list,
        "obdme_multipolarity" : obdme_multipolarity,
        "traditional_ho" : traditional_ho,

        # two-body parameters
        "obs-R20K20" : False,
        "obs-am-sqr" : True,

        # version parameters
        "mfdn_executable" : mfdn_executable,
        "mfdn_wrapper" : mfdn_wrapper

        }    
    # many-body-truncation
    for Nmax in Nmax_list
    # unit basis
    for (radial_basis_p, radial_basis_n) in radial_basis_pair_list
    for beta_ratio in beta_ratio_list
    # scaling
    for hw in hw_list
    # Lawson
    for hwLawson in hwLawson_list 
    # source xform
    for hw_int in mcscript.auto_value(hw_int_list,[hw])
    for xform_cutoff in xform_cutoff_list
]

##################################################################
# task postprocessing functions
##################################################################

def task_pool (current_task):
    """Defines the pool name for the given task, used to identify sets of
    similar tasks.

    """

    pool = "Nmax{Nmax:02d}".format(**current_task)
    return pool

def task_mask (current_task):
    """ Defines whether or not the given task is eligible for
    execution.
    """
    ## allow = mcscript.approx_equal(current_task["hw"],20.,0.1)
    allow = True
    return allow

##################################################################
# master loop
##################################################################

mcscript.task.init(
    tasks,
    task_descriptor=ncsm.task_descriptor_format_6, 
    task_pool=task_pool,
    task_mask=task_mask,
    phase_handler_list=[
        ncsm.task_handler_h2mixer,
        ncsm.task_handler_mfdn_h2
        ],
    archive_handler_list=[
        ncsm.archive_handler_mfdn_res_only,
        ncsm.archive_handler_mfdn_archive
        ]
    )

################################################################
# termination
################################################################

mcscript.termination()
