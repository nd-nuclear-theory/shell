import math

import mcscript

from . import config

##################################################################
# traditional ho run
##################################################################

def set_up_orbitals(task):
    """Set up source and target orbitals for MFDn run.

    Arguments:
        task (dict): as described in module docstring

    Limitation: Currently only supports harmonic oscillator style
    truncation.
    """

    # validate truncation mode
    if (task["truncation_mode"] is not config.TruncationMode.kHO):
        raise ValueError("expecting truncation_mode to be {} but found {truncation_mode}".format(config.TruncationMode.kHO,**task))

    natural_orbital_iteration = task.get("natorb_iteration")

    # generating orbitals for a non-natural orbitals run is the same
    # as generating orbitals for the 0th natural orbital run
    # generate orbitals and basis xform -- target basis
    if (natural_orbital_iteration in {None,0}):
        # generate orbitals -- interaction bases
        mcscript.call(
            [
                config.environ.shell_filename("orbital-gen"),
                "--oscillator",
                "{truncation_int[1]:d}".format(**task),
                "{:s}".format(config.filenames.orbitals_int_filename(natural_orbital_iteration))
            ]
        )
        if (task["use_coulomb"]):
            mcscript.call(
                [
                    config.environ.shell_filename("orbital-gen"),
                    "--oscillator",
                    "{truncation_coul[1]:d}".format(**task),
                    "{:s}".format(config.filenames.orbitals_coul_filename(natural_orbital_iteration))
                ]
            )

        # generate orbitals -- target basis
        if (task["truncation_mode"] is config.TruncationMode.kHO):
            truncation_parameters = task["truncation_parameters"]
            if (truncation_parameters["many_body_truncation"]=="Nmax"):
                Nmax_orb = truncation_parameters["Nmax"] + truncation_parameters["Nv"]
            elif (truncation_parameters["many_body_truncation"]=="FCI"):
                Nmax_orb = truncation_parameters["Nmax"]
            mcscript.call(
                [
                    config.environ.shell_filename("orbital-gen"),
                    "--oscillator",
                    "{Nmax_orb:d}".format(Nmax_orb=Nmax_orb),
                    "{:s}".format(config.filenames.orbitals_filename(natural_orbital_iteration))
                ]
            )
            mcscript.call(
                [
                    config.environ.shell_filename("radial-gen"),
                    "--identity",
                    config.filenames.orbitals_filename(natural_orbital_iteration),
                    config.filenames.radial_xform_filename(natural_orbital_iteration)
                ]
            )
        else:
            # TODO implement non-HO truncation
            raise mcscript.ScriptError("expecting ho_truncation to be True but found {ho_truncation}".format(**task))
    else:
        mcscript.call(
            [
                config.environ.shell_filename("natorb-gen"),
                config.filenames.orbitals_filename(natural_orbital_iteration-1),
                config.filenames.natorb_info_filename(natural_orbital_iteration-1),
                config.filenames.natorb_obdme_filename(natural_orbital_iteration-1),
                config.filenames.natorb_xform_filename(natural_orbital_iteration),
                config.filenames.orbitals_filename(natural_orbital_iteration)
            ]
        )

def set_up_radial(task):
    """Generate radial integrals and overlaps for MFDn run.

    Operation mode may in general be direct oscillator, dilated
    oscillator, or generic.

    Arguments:
        task (dict): as described in module_docstring

    """

    # get natural orbital iteration; None should be treated the same as 0
    natural_orbital_iteration = task.get("natorb_iteration")
    if (natural_orbital_iteration in {None,0}):
        return set_up_radial_analytic(task)
    elif (natural_orbital_iteration > 0):
        return set_up_radial_natorb(task)

def set_up_radial_analytic(task):
    """Generate radial integrals and overlaps by integration for MFDn run
    in analytic basis.

    Operation mode may in general be direct oscillator, dilated
    oscillator, or generic (TODO).

    Arguments:
        task (dict): as described in module docstring

    """

    # validate basis mode
    if (task["basis_mode"] not in {config.BasisMode.kDirect,config.BasisMode.kDilated}):  # no config.BasisMode.kGeneric yet
        raise ValueError("invalid basis mode {basis_mode}".format(**task))

    # get natural orbital iteration
    natural_orbital_iteration = task.get("natorb_iteration")

    # basis radial code -- expected by radial_utils codes
    basis_radial_code = "oscillator"  # TO GENERALIZE: if not oscillator basis

    # generate radial integrals
    for operator_type in ["r","k"]:
        for power in [1,2]:
            mcscript.call(
                [
                    config.environ.shell_filename("radial-gen"),
                    "--kinematic",
                    "{:s}".format(operator_type),
                    "{:d}".format(power),
                    basis_radial_code,
                    config.filenames.orbitals_filename(natural_orbital_iteration),
                    config.filenames.radial_me_filename(natural_orbital_iteration, operator_type, power)
                ],
                mode = mcscript.call.serial
            )

    # generate radial overlaps -- generate trivial identities if applicable
    if (task["basis_mode"] in {config.BasisMode.kDirect}):
        mcscript.call(
            [
                config.environ.shell_filename("radial-gen"),
                "--identity",
                config.filenames.orbitals_int_filename(natural_orbital_iteration),
                config.filenames.orbitals_filename(natural_orbital_iteration),
                config.filenames.radial_olap_int_filename(natural_orbital_iteration)
            ],
            mode = mcscript.call.serial
        )
    else:
        b_ratio = math.sqrt(task["hw_int"]/task["hw"])
        mcscript.call(
            [
                config.environ.shell_filename("radial-gen"),
                "--overlaps",
                "{:g}".format(b_ratio),
                basis_radial_code,
                config.filenames.orbitals_int_filename(natural_orbital_iteration),
                config.filenames.orbitals_filename(natural_orbital_iteration),
                config.filenames.radial_olap_int_filename(natural_orbital_iteration)
            ],
            mode = mcscript.call.serial
        )
    if (task["use_coulomb"]):
        if (task["basis_mode"] in {config.BasisMode.kDirect,config.BasisMode.kDilated}):
            mcscript.call(
                [
                    config.environ.shell_filename("radial-gen"),
                    "--identity",
                    config.filenames.orbitals_coul_filename(natural_orbital_iteration),
                    config.filenames.orbitals_filename(natural_orbital_iteration),
                    config.filenames.radial_olap_coul_filename(natural_orbital_iteration)
                ],
                mode = mcscript.call.serial
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
                    config.filenames.orbitals_coul_filename(natural_orbital_iteration),
                    config.filenames.orbitals_filename(natural_orbital_iteration),
                    config.filenames.radial_olap_coul_filename(natural_orbital_iteration)
                ],
                mode = mcscript.call.serial
            )

def set_up_radial_natorb(task):
    """Generate radial integrals and overlaps by transformation for MFDn run
    in natural orbital basis.

    Operation mode must be generic.

    Arguments:
        task (dict): as described in module docstring

    """

    # get natural orbital iteration
    natural_orbital_iteration = task.get("natorb_iteration")

    if (natural_orbital_iteration in {None,0}):
        raise mcscript.ScriptError(
            "cannot set up natural orbital radial files in natorb_iteration {}".format(natural_orbital_iteration)
        )

    # compose radial transform
    mcscript.call(
        [
            config.environ.shell_filename("radial-compose"),
            config.filenames.radial_xform_filename(natural_orbital_iteration-1),
            config.filenames.natorb_xform_filename(natural_orbital_iteration),
            config.filenames.radial_xform_filename(natural_orbital_iteration)
        ]
    )

    # compose interaction transform
    mcscript.call(
        [
            config.environ.shell_filename("radial-compose"),
            config.filenames.radial_olap_int_filename(natural_orbital_iteration-1),
            config.filenames.natorb_xform_filename(natural_orbital_iteration),
            config.filenames.radial_olap_int_filename(natural_orbital_iteration)
        ]
    )

    # compose Coulomb transform
    mcscript.call(
        [
            config.environ.shell_filename("radial-compose"),
            config.filenames.radial_olap_coul_filename(natural_orbital_iteration-1),
            config.filenames.natorb_xform_filename(natural_orbital_iteration),
            config.filenames.radial_olap_coul_filename(natural_orbital_iteration)
        ]
    )

    # transform radial integrals
    for operator_type in ["r","k"]:
        for power in [1,2]:
            mcscript.call(
                [
                    config.environ.shell_filename("radial-xform"),
                    config.filenames.orbitals_filename(natural_orbital_iteration),
                    config.filenames.radial_xform_filename(natural_orbital_iteration),
                    config.filenames.radial_me_filename(0, operator_type, power),
                    config.filenames.radial_me_filename(natural_orbital_iteration, operator_type, power)
                ],
                mode = mcscript.call.serial
            )
