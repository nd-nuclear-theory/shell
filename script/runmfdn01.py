""" runmfdn01.py -- h2mixer+mfdn mscript test run

rm runmfdn01/flags/*
qsubm --here mfdn01 --pool=test --limit=1 --noredirect 


    Mark A. Caprio
    University of Notre Dame

    - 12/14/16 (mac): Created.
      
"""

import mcscript
import mfdn_h2

# load mcscript
import mcscript
mcscript.init()

##################################################################
# build task list
##################################################################

params = {
    # stored radial
    "Nmax_orb" : 2,  # for target space
    "Nmax_orb_int" : 20,  # for overlaps from interaction tbmes
    "orbitals_filename" : "orbitals.dat",  # for target space
    "orbitals_int_filename" : "orbitals-int.dat",  # for overlaps from interaction tbmes
    "radial_me_filename_pattern" : "radial-me-{}{}.dat",  # "{}{}" will be replaced by {"r1","r2","k1","k2"}
    "radial_olap_filename" : "radial-olap.dat",
    # stored input TBMEs -- TO GENERALIZE
    "VNN_filename" : "run0164-ob-9/JISP16-ob-9-20.bin",
    "VC_filename" : "run0164-ob-9/VC-ob-9-20.bin",
    # scaling/transformation/basis parameters
    "hw" : 20,
    "target_truncation" : ("tb",2),
    "hw_int" : 20,
    "hw_c" : 20,
    "xform_enabled" : False,
    "xform_truncation" : ("ob",9),
    # Hamiltonian parameters
    "A" : 4,   # atomic mass number of nucleus (for kinematic operators)
    "a" : 20,  # coefficient on Ncm in Lawson term
    "hw_cm" : None,  # hw for Ncm in Lawson term; if None, defaults to "hw"
    "use_coulomb" : True,
    # output file format
    "h2_format" : 0  # h2 formats: 0, 15099
}

# generate task list

tasks = [
    params
]

##################################################################
# implementation functions for doing a "hello world" task
#
# For a more complicated application, you would separate these out
# into their own module.
##################################################################

def task_descriptor(current_task):
    """ Return task descriptor for hello task.
    """

    return "test"

##################################################################
# task list entry annotation functions
##################################################################

def task_pool (current_task):
    """ Create task pool identifier.
    """
    
    return "test"

##################################################################
# master loop
##################################################################

mcscript.task.init(
    tasks,
    task_descriptor=task_descriptor,
    task_pool=task_pool,
    phase_handler_list=[mfdn_h2.task_handler_ho]
    )

################################################################
# termination
################################################################

mcscript.termination()
