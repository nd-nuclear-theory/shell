""" runmfdn10.py

    See runmfdn.txt for description.

    Patrick J. Fasano
    University of Notre Dame

    - 09/27/17 (pjf): Created, copied from runmfd07.
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

tasks = [{
    # nuclide parameters
    "nuclide": (2, 6),

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
        "Nv": 1,
        "Nmax": 8,
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
    # "obdme_reference_state_list": [(0, 0, 1)],
    "save_obdme": True,

    # two-body observables
    ## "observable_sets": ["H-components","am-sqr"],
    "observable_sets": ["H-components"],
    "tb_observables": [],

    # wavefunction storage
    "save_wavefunctions": True,

    # version parameters
    "h2_format": 15099,
    "mfdn_executable": "v15-beta01/xmfdn-h2-lan",
    "mfdn_driver": mfdn.mfdn_v15
}]

##################################################################
# task dictionary postprocessing functions
##################################################################

def task_pool(current_task):
    pool = "Nmax{truncation_parameters[Nmax]:02d}-Mj{Mj:3.1f}".format(**current_task)
    return pool


##################################################################
# task control
##################################################################

mcscript.task.init(
    tasks,
    task_descriptor=mfdn.descriptors.task_descriptor_c1,
    task_pool=task_pool,
    phase_handler_list=[
        mfdn.handlers.task_handler_dimension,
        mfdn.handlers.task_handler_nonzeros
        ]
    )

################################################################
# termination
################################################################

mcscript.termination()
