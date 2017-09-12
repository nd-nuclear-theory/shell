""" runmfdn07.py

    See runmfdn.txt for description.

    Mark A. Caprio
    University of Notre Dame

    - 06/08/17 (pjf): Created, copied from runmfd01; switch to MFDn v15.
    - 07/31/17 (pjf): Set MFDn driver module in task dictionary.
    - 08/11/17 (pjf): Update for split single-particle and many-body truncation modes.
"""

import mcscript
import mfdn
import mfdn.mfdn_v15

# initialize mcscript
mcscript.init()

##################################################################
# build task list
##################################################################

mfdn.environ.environ.interaction_run_list = [
    "run0164-JISP16-ob-9",
    "run0164-JISP16-ob-13",
    "run0164-JISP16-tb-10",
    "run0164-JISP16-tb-20",
    "run0306-N2LOopt500",  # up to tb-20
    "runvc0083-Daejeon16-ob-13"
]

task = {
    # nuclide parameters
    "nuclide": (2, 2),

    # Hamiltonian parameters
    "interaction": "JISP16",
    "use_coulomb": True,
    "a_cm": 20.,
    "hw_cm": None,

    # input TBME parameters
    "truncation_int": ("tb", 10),
    "hw_int": 20.,
    "truncation_coul": ("tb", 10),
    "hw_coul": 20.,

    # basis parameters
    "basis_mode": mfdn.modes.BasisMode.kDirect,
    "hw": 20.,

    # transformation parameters
    "xform_truncation_int": None,
    "xform_truncation_coul": None,
    "hw_coul_rescaled": None,
    "target_truncation": None,

    # traditional oscillator many-body truncation
    "sp_truncation_mode": mfdn.modes.SingleParticleTruncationMode.kNmax,
    "mb_truncation_mode": mfdn.modes.ManyBodyTruncationMode.kNmax,
    "truncation_parameters": {
        "Nv": 0,
        "Nmax": 2,
        "Nstep": 2,
        },

    # diagonalization parameters
    "Mj": 0,
    "eigenvectors": 2,
    "initial_vector": -2,
    "lanczos": 200,
    "tolerance": 1e-6,
    "partition_filename": None,

    # obdme parameters
    ## "hw_for_trans": 20,
    "obdme_multipolarity": 2,
    "obdme_reference_state_list": [(0, 0, 1)],
    "save_obdme": True,

    # two-body observables
    ## "observable_sets": ["H-components","am-sqr"],
    "observable_sets": ["H-components"],

    # version parameters
    "h2_format": 15099,
    "mfdn_executable": "v15-beta00/xmfdn-h2-lan",
    "mfdn_driver": mfdn.mfdn_v15,

}

################################################################
# run control
################################################################

# add task descriptor metadata field (needed for filenames)
task["metadata"] = {
    "descriptor": mfdn.descriptors.task_descriptor_7(task)
    }

mfdn.radial.set_up_interaction_orbitals(task)
mfdn.radial.set_up_orbitals(task)
mfdn.radial.set_up_radial_analytic(task)
mfdn.tbme.generate_tbme(task)
mfdn.mfdn_v15.run_mfdn(task)
mfdn.mfdn_v15.save_mfdn_output(task)
# mfdn.handlers.task_handler_oscillator(task)

##################################################################
# task control
##################################################################

## mcscript.task.init(
##     tasks,
##     task_descriptor=mfdn.descriptors.task_descriptor_7,
##     task_pool=task_pool,
##     phase_handler_list=[mfdn_h2.task_handler_ho]
##     )

################################################################
# termination
################################################################

mcscript.termination()
