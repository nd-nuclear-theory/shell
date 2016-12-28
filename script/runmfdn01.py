""" runmfdn01.py -- h2mixer+mfdn mscript test run

rm runmfdn01/flags/*
qsubm --here mfdn01 --pool=test --limit=1 --noredirect 

qsubm --here mfdn01

    Mark A. Caprio
    University of Notre Dame

    - 12/14/16 (mac): Created.
      
"""

import mcscript
import mfdn

# initialize mcscript
mcscript.init()

##################################################################
# build task list
##################################################################

mfdn.configuration.interaction_run_list = ["run0164-ob-9"]

task = {
    # nuclide parameters
    "nuclide" : (2,2),

    # Hamiltonian parameters
    "interaction" : "JISP16",
    "use_coulomb" : True,
    "a_cm" : 20.,
    "hw_cm" : None,

    # input TBME parameters
    "truncation_int" : ("ob",9),
    "hw_int" : 20.,
    "truncation_coul" : ("ob",9),
    "hw_coul" : 20.,

    # basis parameters
    "basis_mode" : mfdn.k_basis_mode_direct,
    "hw" : 20.,

    # transformation parameters
    "xform_truncation_int" : None,
    "xform_truncation_coul" : None,
    "hw_coul_rescaled" : None,
    "target_truncation" : None,

    # traditional oscillator many-body truncation
    "ho_truncation" : True,
    "Nv" : 0,
    "Nmax" : 2,
    "many_body_truncation" : "Nmax",
    "Nstep" : 2,

    # diagonalization parameters
    "Mj" : 0,
    "eigenvectors" : 2,
    "initial_vector" : -2,
    "lanczos" : 200,
    "tolerance" : 1e-6,

    # obdme parameters
    ## "hw_for_trans" : 20,
    "obdme_multipolarity" : 2,
    "obdme_reference_state_list" : [(0,0,1)],
    "save_obdme" : True,

    # two-body observables
    ## "obs_basename_list" : ["tbme-rrel2","tbme-Ncm"],
    "observable_sets" : ["H-components","am-sqr"],

    # version parameters
    "h2_format" : 0,
    "mfdn_executable" : "mfdn-v14-beta06-newmake/xmfdn-h2-lan"

}

## mfdn.configuration.interaction_filename("JISP16-ob-9-20.bin")

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

## mcscript.task.init(
##     tasks,
##     task_descriptor=task_descriptor,
##     task_pool=task_pool,
##     phase_handler_list=[mfdn_h2.task_handler_ho]
##     )

mfdn.set_up_orbitals(task)
mfdn.set_up_radial_analytic(task)
mfdn.generate_tbme(task)
mfdn.run_mfdn_v14(task)

################################################################
# termination
################################################################

mcscript.termination()
