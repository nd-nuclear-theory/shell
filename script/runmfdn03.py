""" runmfdn02.py

    See runmfdn.txt for description.

    Mark A. Caprio
    University of Notre Dame

    - 12/29/16 (mac): Created.
    - 1/29/17 (pjf): Updated for new truncation_mode parameter.

"""

import mcscript
import mfdn

# initialize mcscript
mcscript.init()

##################################################################
# build task list
##################################################################

mfdn.configuration.interaction_run_list = [
    "run0164-JISP16-ob-9",
    "run0164-JISP16-ob-13",
    "run0164-JISP16-tb-10",
    "run0164-JISP16-tb-20",
    "run0306-N2LOopt500",  # up to tb-20
    "runvc0083-Daejeon16-ob-13"
]

task = {
    # nuclide parameters
    "nuclide" : (2,2),

    # Hamiltonian parameters
    "interaction" : "JISP16",
    "use_coulomb" : True,
    "a_cm" : 20.,
    "hw_cm" : None,

    # input TBME parameters
    "truncation_int" : ("ob",13),
    "hw_int" : 40.,
    "truncation_coul" : ("ob",13),
    "hw_coul" : 20.,

    # basis parameters
    "basis_mode" : mfdn.k_basis_mode_dilated,
    "hw" : 20.,

    # transformation parameters
    "xform_truncation_int" : ("ob",13),
    "xform_truncation_coul" : None,
    "hw_coul_rescaled" : None,
    "target_truncation" : None,

    # traditional oscillator many-body truncation
    "truncation_mode" : mfdn.k_truncation_mode_ho,
    "truncation_parameters": {
        "Nv" : 0,
        "Nmax" : 2,
        "many_body_truncation" : "Nmax",
        "Nstep" : 2,
        },

    # diagonalization parameters
    "Mj" : 0,
    "eigenvectors" : 2,
    "initial_vector" : -2,
    "lanczos" : 200,
    "tolerance" : 1e-6,
    "partition_filename" : None,

    # obdme parameters
    ## "hw_for_trans" : 20,
    "obdme_multipolarity" : 2,
    "obdme_reference_state_list" : [(0,0,1)],
    "save_obdme" : True,

    # two-body observables
    ## "observable_sets" : ["H-components","am-sqr"],
    "observable_sets" : ["H-components"],

    # version parameters
    "h2_format" : 0,
    "mfdn_executable" : "mfdn-v14-beta06-newmake/xmfdn-h2-lan"

}

################################################################
# run control
################################################################

# add task descriptor field (needed for filenames)
task["descriptor"] = mfdn.task_descriptor_7(task)

mfdn.set_up_orbitals(task)
mfdn.set_up_radial(task)
mfdn.generate_tbme(task)
mfdn.run_mfdn_v14_b06(task)
mfdn.save_mfdn_v14_output(task)

##################################################################
# task control
##################################################################

## mcscript.task.init(
##     tasks,
##     task_descriptor=mfdn.task_descriptor_7,
##     task_pool=task_pool,
##     phase_handler_list=[mfdn_h2.task_handler_ho]
##     )

################################################################
# termination
################################################################

mcscript.termination()
