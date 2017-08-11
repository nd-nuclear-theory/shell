""" runmfdn02.py

    See runmfdn.txt for description.

    Mark A. Caprio
    University of Notre Dame

    - 12/29/16 (mac): Created.
    - 01/29/17 (pjf): Updated for new truncation_mode parameter.
    - 06/03/17 (pjf): Updated for new scripting.
    - 07/31/17 (pjf): Set MFDn driver module in task dictionary.
    - 08/11/17 (pjf): Update for split single-particle and many-body truncation modes.
"""

import mcscript
import mfdn
import mfdn.mfdn_v14

# initialize mcscript
mcscript.init()

##################################################################
# build task list
##################################################################

mfdn.config.environ.interaction_run_list = [
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
    "truncation_int": ("ob", 13),
    "hw_int": 40.,
    "truncation_coul": ("ob", 13),
    "hw_coul": 20.,

    # basis parameters
    "basis_mode": mfdn.config.BasisMode.kDilated,
    "hw": 20.,

    # transformation parameters
    "xform_truncation_int": ("ob", 13),
    "xform_truncation_coul": None,
    "hw_coul_rescaled": None,
    "target_truncation": None,

    # traditional oscillator many-body truncation
    "sp_truncation_mode": mfdn.config.SingleParticleTruncationMode.kNmax,
    "mb_truncation_mode": mfdn.config.ManyBodyTruncationMode.kNmax,
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
    "h2_format": 0,
    "mfdn_executable": "mfdn-v14-beta06-newmake/xmfdn-h2-lan",
    "mfdn_driver": mfdn.mfdn_v14,

}

################################################################
# run control
################################################################

# add task descriptor metadata field (needed for filenames)
task["metadata"] = {
    "descriptor": mfdn.descriptors.task_descriptor_7(task)
    }

mfdn.radial.set_up_orbitals_ho(task)
mfdn.radial.set_up_radial_analytic(task)
mfdn.tbme.generate_tbme(task)
mfdn.mfdn_v14.run_mfdn(task)
mfdn.mfdn_v14.save_mfdn_output(task)

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
