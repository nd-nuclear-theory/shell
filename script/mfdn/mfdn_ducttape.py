"""mfdn_ducttape.py -- driver module for MFDn v15 duct-tape postprocessor.

Patrick Fasano
University of Notre Dame

- 10/05/17 (pjf): Copied/stripped from mfdn_v15.py.
- 10/18/17 (pjf): Use separate work directory for each postfix.
"""
import os
import glob
import collections

import mcscript

from . import modes, environ
from .mfdn_v15 import truncation_setup_functions


def run_mfdn(task, postfix=""):
    """Generate input file and execute MFDn version 15 beta 00 duct-tape postprocessor.

    Arguments:
        task (dict): as described in module docstring
        run_mode (modes.MFDnRunMode): run mode for MFDn
        postfix (string, optional): identifier to add to generated files

    Raises:
        mcscript.exception.ScriptError: if MFDn output not found
    """
    # create work directory if it doesn't exist yet (-p)
    work_dir = "work{:s}".format(postfix)
    mcscript.call(["mkdir", "-p", work_dir])

    # inputlist namelist dictionary
    inputlist = collections.OrderedDict()
    # tbo: two-body observable namelist
    obslist = collections.OrderedDict()

    # nucleus
    inputlist["Nprotons"], inputlist["Nneutrons"] = task["nuclide"]

    # Mj
    inputlist["TwoMj"] = int(2*task["Mj"])

    # single-particle orbitals
    inputlist["orbitalfile"] = environ.filenames.orbitals_filename(postfix)

    # truncation mode
    truncation_setup_functions[task["mb_truncation_mode"]](task, inputlist)

    if (task["basis_mode"] in {modes.BasisMode.kDirect, modes.BasisMode.kDilated}):
        inputlist["hbomeg"] = float(task["hw"])

    # diagonalization parameters
    inputlist["neivals"] = int(task["eigenvectors"])

    # obdme: parameters
    inputlist["obdme"] = True
    obslist["max2K"] = int(2*task["obdme_multipolarity"])

    # construct transition observable input if reference states given
    if task.get("obdme_reference_state_list") is not None:
        # obdme: validate reference state list
        #
        # guard against pathetically common mistakes
        for (J, g, i) in task["obdme_reference_state_list"]:
            # validate integer/half-integer character of angular momentum
            twice_J = int(2*J)
            if (twice_J % 2) != (sum(task["nuclide"]) % 2):
                raise ValueError("invalid angular momentum for reference state")
            # validate grade
            # if (-1)**g != inputlist["parity"]:
            #     raise ValueError("invalid parity for reference state")

        # obdme: construct input
        inputlist["nrefstates"] = len(task["obdme_reference_state_list"])
        obslist["ref2J"] = []
        obslist["refseq"] = []
        for (J, g, i) in task["obdme_reference_state_list"]:
            obslist["ref2J"].append(int(2*J))
            obslist["refseq"].append(i)

    # generate MFDn input file
    mcscript.utils.write_namelist(
        os.path.join(work_dir, "mfdn.input"),
        input_dict={"inputlist": inputlist, "obslist": obslist}
    )

    # enter work directory
    os.chdir(work_dir)

    # invoke MFDn
    mcscript.call(
        [
            environ.environ.mfdn_filename(task["ducttape_executable"])
        ],
        mode=mcscript.CallMode.kHybrid,
        check_return=True
    )

    # test for basic indications of success
    if not os.path.exists("mfdn.out"):
        raise mcscript.exception.ScriptError("mfdn.out not found")
    if not os.path.exists("mfdn.res"):
        raise mcscript.exception.ScriptError("mfdn.res not found")

    # leave work directory
    os.chdir("..")
