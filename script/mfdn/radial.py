"""handlers.py -- task handlers for MFDn runs.

Patrick Fasano
University of Notre Dame

- 03/22/17 (pjf): Created, split from __init__.py.
- 04/06/17 (pjf): Correctly reference config submodule (mfdn.config).
- 04/07/17 (pjf): Update for mcscript namespace changes.
- 04/11/17 (pjf): Fix broken imports.
- 06/03/17 (pjf): Remove explicit references to natural orbitals from bulk of
    scripting.
- 06/07/17 (pjf): Clean up style.
- 08/11/17 (pjf): Use new TruncationModes.
"""
import math

import mcscript
import mcscript.exception

from . import config

##################################################################
# traditional ho run
##################################################################


def set_up_orbitals_ho(task, postfix=""):
    """Set up source and target orbitals for MFDn run.

    Arguments:
        task (dict): as described in module docstring
        postfix (string, optional): identifier to add to generated files
    """
    # validate truncation mode
    if task["sp_truncation_mode"] is not config.SingleParticleTruncationMode.kNmax:
        raise ValueError("expecting truncation_mode to be {} but found {truncation_mode}".format(config.SingleParticleTruncationMode.kNmax, **task))

    # generate orbitals -- interaction bases
    mcscript.call(
        [
            config.environ.shell_filename("orbital-gen"),
            "--oscillator",
            "{truncation_int[1]:d}".format(**task),
            "{:s}".format(config.filenames.orbitals_int_filename(postfix))
        ]
    )
    if (task["use_coulomb"]):
        mcscript.call(
            [
                config.environ.shell_filename("orbital-gen"),
                "--oscillator",
                "{truncation_coul[1]:d}".format(**task),
                "{:s}".format(config.filenames.orbitals_coul_filename(postfix))
            ]
        )

    # generate orbitals -- target basis
    truncation_parameters = task["truncation_parameters"]
    if task["mb_truncation_mode"] == config.ManyBodyTruncationMode.kNmax:
        Nmax_orb = truncation_parameters["Nmax"] + truncation_parameters["Nv"]
    elif task["mb_truncation_mode"] == config.ManyBodyTruncationMode.kFCI:
        Nmax_orb = truncation_parameters["Nmax"]
    mcscript.call(
        [
            config.environ.shell_filename("orbital-gen"),
            "--oscillator",
            "{Nmax_orb:d}".format(Nmax_orb=Nmax_orb),
            "{:s}".format(config.filenames.orbitals_filename(postfix))
        ]
    )
    mcscript.call(
        [
            config.environ.shell_filename("radial-gen"),
            "--identity",
            config.filenames.orbitals_filename(postfix),
            config.filenames.radial_xform_filename(postfix)
        ]
    )


def set_up_orbitals_natorb(task, source_postfix, target_postfix):
    """Set up natural orbitals for MFDn run.

    Arguments:
        task (dict): as described in module docstring
        source_postfix (string): identifier for source of natural orbital information
        target_postfix (string): identifier to add to generated files

    Limitation: Currently only supports harmonic oscillator style
    truncation.
    """
    # validate truncation mode
    if task["sp_truncation_mode"] is not config.SingleParticleTruncationMode.kNmax:
        raise ValueError("expecting truncation_mode to be {} but found {truncation_mode}".format(config.SingleParticleTruncationMode.kNmax, **task))

    # validate natural orbitals enabled
    if not task.get("natural_orbitals"):
        raise mcscript.exception.ScriptError("natural orbitals are not enabled")

    mcscript.call(
        [
            config.environ.shell_filename("natorb-gen"),
            config.filenames.orbitals_filename(source_postfix),
            config.filenames.natorb_info_filename(source_postfix),
            config.filenames.natorb_obdme_filename(source_postfix),
            config.filenames.natorb_xform_filename(target_postfix),
            config.filenames.orbitals_filename(target_postfix)
        ]
    )

# def set_up_radial(task):
#     """Generate radial integrals and overlaps for MFDn run.
#
#     Operation mode may in general be direct oscillator, dilated
#     oscillator, or generic.
#
#     Arguments:
#         task (dict): as described in module_docstring
#
#     """
#
#     # get natural orbital iteration; None should be treated the same as 0
#     natural_orbital_iteration = task.get("natorb_iteration")
#     if (natural_orbital_iteration in {None,0}):
#         return set_up_radial_analytic(task)
#     elif (natural_orbital_iteration > 0):
#         return set_up_radial_natorb(task)


def set_up_radial_analytic(task, postfix=""):
    """Generate radial integrals and overlaps by integration for MFDn run
    in analytic basis.

    Operation mode may in general be direct oscillator, dilated
    oscillator, or generic (TODO).

    Arguments:
        task (dict): as described in module docstring
        postfix (string, optional): identifier to add to generated files
    """
    # validate basis mode
    if (task["basis_mode"] not in {config.BasisMode.kDirect, config.BasisMode.kDilated}):  # no config.BasisMode.kGeneric yet
        raise ValueError("invalid basis mode {basis_mode}".format(**task))

    # basis radial code -- expected by radial_utils codes
    basis_radial_code = "oscillator"  # TO GENERALIZE: if not oscillator basis

    # generate radial integrals
    for operator_type in ["r", "k"]:
        for power in [1, 2]:
            mcscript.call(
                [
                    config.environ.shell_filename("radial-gen"),
                    "--kinematic",
                    "{:s}".format(operator_type),
                    "{:d}".format(power),
                    basis_radial_code,
                    config.filenames.orbitals_filename(postfix),
                    config.filenames.radial_me_filename(postfix, operator_type, power)
                ],
                mode=mcscript.CallMode.kSerial
            )

    # generate radial overlaps -- generate trivial identities if applicable
    if (task["basis_mode"] in {config.BasisMode.kDirect}):
        mcscript.call(
            [
                config.environ.shell_filename("radial-gen"),
                "--identity",
                config.filenames.orbitals_int_filename(postfix),
                config.filenames.orbitals_filename(postfix),
                config.filenames.radial_olap_int_filename(postfix)
            ],
            mode=mcscript.CallMode.kSerial
        )
    else:
        b_ratio = math.sqrt(task["hw_int"]/task["hw"])
        mcscript.call(
            [
                config.environ.shell_filename("radial-gen"),
                "--overlaps",
                "{:g}".format(b_ratio),
                basis_radial_code,
                config.filenames.orbitals_int_filename(postfix),
                config.filenames.orbitals_filename(postfix),
                config.filenames.radial_olap_int_filename(postfix)
            ],
            mode=mcscript.CallMode.kSerial
        )
    if (task["use_coulomb"]):
        if (task["basis_mode"] in {config.BasisMode.kDirect, config.BasisMode.kDilated}):
            mcscript.call(
                [
                    config.environ.shell_filename("radial-gen"),
                    "--identity",
                    config.filenames.orbitals_coul_filename(postfix),
                    config.filenames.orbitals_filename(postfix),
                    config.filenames.radial_olap_coul_filename(postfix)
                ],
                mode=mcscript.CallMode.kSerial
            )
        else:
            if task.get("hw_coul_rescaled") is None:
                b_ratio = 1
            else:
                b_ratio = math.sqrt(task["hw_coul_rescaled"]/task["hw"])
            mcscript.call(
                [
                    config.environ.shell_filename("radial-gen"),
                    "--overlaps",
                    "{:g}".format(b_ratio),
                    basis_radial_code,
                    config.filenames.orbitals_coul_filename(postfix),
                    config.filenames.orbitals_filename(postfix),
                    config.filenames.radial_olap_coul_filename(postfix)
                ],
                mode=mcscript.CallMode.kSerial
            )


def set_up_radial_natorb(task, source_postfix, target_postfix):
    """Generate radial integrals and overlaps by transformation for MFDn run in natural orbital basis.

    Operation mode must be generic.

    Arguments:
        task (dict): as described in module docstring
        source_postfix (str): postfix for old basis
        target_postfix (str): postfix for new basis
    """
    # validate natural orbitals enabled
    if not task.get("natural_orbitals"):
        raise mcscript.exception.ScriptError("natural orbitals are not enabled")

    # compose radial transform
    mcscript.call(
        [
            config.environ.shell_filename("radial-compose"),
            config.filenames.radial_xform_filename(source_postfix),
            config.filenames.natorb_xform_filename(target_postfix),
            config.filenames.radial_xform_filename(target_postfix)
        ]
    )

    # compose interaction transform
    mcscript.call(
        [
            config.environ.shell_filename("radial-compose"),
            config.filenames.radial_olap_int_filename(source_postfix),
            config.filenames.natorb_xform_filename(target_postfix),
            config.filenames.radial_olap_int_filename(target_postfix)
        ]
    )

    # compose Coulomb transform
    mcscript.call(
        [
            config.environ.shell_filename("radial-compose"),
            config.filenames.radial_olap_coul_filename(source_postfix),
            config.filenames.natorb_xform_filename(target_postfix),
            config.filenames.radial_olap_coul_filename(target_postfix)
        ],
        mode=mcscript.CallMode.kSerial
    )

    # transform radial integrals
    for operator_type in ["r", "k"]:
        for power in [1, 2]:
            mcscript.call(
                [
                    config.environ.shell_filename("radial-xform"),
                    config.filenames.orbitals_filename(target_postfix),
                    config.filenames.natorb_xform_filename(target_postfix),
                    config.filenames.radial_me_filename(source_postfix, operator_type, power),
                    config.filenames.radial_me_filename(target_postfix, operator_type, power)
                ],
                mode=mcscript.CallMode.kSerial
            )
