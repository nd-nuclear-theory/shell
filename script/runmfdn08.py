""" runmfdn08.py

    See runmfdn.txt for description.

    Mark A. Caprio, Patrick J. Fasano
    University of Notre Dame

    - 08/01/17 (pjf): Created, copied from runmfd07; switch to MFDn v15 b01.
    - 08/11/17 (pjf): Update for split single-particle and many-body truncation modes.
"""

import mcscript
import mfdn
import mfdn.mfdn_v15

try:
    mcscript.control.module(["swap", "craype-haswell", "craype-mic-knl"])
    mcscript.control.module(["load", "craype-hugepages2M"])
    mcscript.control.module(["list"])
except:
    print("problem with modules")

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

# hw -- linear mesh
hw_range = (15, 25, 5)
hw_list = mcscript.utils.value_range(*hw_range)

# hw -- log mesh
## hw_log_range = (5, 40, 8)
## hw_list = mcscript.utils.log_range(*hw_log_range)

# Nmax
Nmax_range = (2, 20, 2)
Nmax_list = mcscript.utils.value_range(*Nmax_range)

tasks = [{
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
    "basis_mode": mfdn.modes.BasisMode.kDilated,
    "hw": hw,

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
        "Nmax": Nmax,
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
    "obdme_reference_state_list": None,
    "save_obdme": True,

    # two-body observables
    ## "observable_sets": ["H-components","am-sqr"],
    "observable_sets": ["H-components"],

    # version parameters
    "h2_format": 15099,
    "mfdn_executable": "mfdn-v15-beta01/xmfdn-h2-lan",
    "mfdn_driver": mfdn.mfdn_v15,

    }
    for Nmax in Nmax_list
    for hw in hw_list
]

################################################################
# run control
################################################################

# add task descriptor metadata field (needed for filenames)
## task["metadata"] = {
##     "descriptor": mfdn.descriptors.task_descriptor_7(task)
##     }

## mfdn.radial.set_up_orbitals_ho(task)
## mfdn.radial.set_up_radial_analytic(task)
## mfdn.tbme.generate_tbme(task)
## mfdn.mfdn_v15.run_mfdn(task)
## mfdn.mfdn_v15.save_mfdn_output(task)

## mfdn.handlers.task_handler_oscillator(task)

def task_pool(current_task):
    pool = "Nmax{truncation_parameters[Nmax]:02d}-Mj{Mj:3.1f}".format(**current_task)
    return pool

##################################################################
# task control
##################################################################

mcscript.task.init(
    tasks,
    task_descriptor=mfdn.descriptors.task_descriptor_7,
    task_pool=task_pool,
    phase_handler_list=[mfdn.handlers.task_handler_oscillator]
    )

################################################################
# termination
################################################################

mcscript.termination()
