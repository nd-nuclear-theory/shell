import math

import mcscript.utils

from . import utils, config, operators

def generate_tbme(task):
    """Generate TBMEs for MFDn run.

    Operation mode may in general be direct oscillator, dilated
    oscillator, or generic (TODO).

    Arguments:
        task (dict): as described in module docstring

    """

    # validate basis mode
    if (task["truncation_mode"] is not config.TruncationMode.kHO):
        raise ValueError("expecting truncation_mode to be {} but found {ho_truncation}".format(config.TruncationMode.kHO,**task))

    # extract parameters for convenience
    natural_orbital_iteration = task.get("natorb_iteration")
    A = sum(task["nuclide"])
    a_cm = task["a_cm"]
    hw = task["hw"]
    hw_cm = task["hw_cm"]
    if (hw_cm is None):
        hw_cm = hw
    hw_coul = task["hw_coul"]
    hw_coul_rescaled = task["hw_coul_rescaled"]
    if (hw_coul_rescaled is None):
        hw_coul_rescaled = hw
    xform_truncation_int = task.get("xform_truncation_int")
    if (xform_truncation_int is None):
        xform_truncation_int = task["truncation_int"]
    xform_truncation_coul = task.get("xform_truncation_coul")
    if (xform_truncation_coul is None):
        xform_truncation_coul = task["truncation_coul"]

    # accumulate h2mixer targets
    targets = {}

    # target: Hamiltonian
    if (task.get("hamiltonian")):
        targets["tbme-H"] = task["hamiltonian"]
    else:
        targets["tbme-H"] = operators.Hamiltonian(
            A=A, hw=hw, a_cm=a_cm, bsqr_intr=hw_cm/hw,
            use_coulomb=task["use_coulomb"], bsqr_coul=hw_coul_rescaled/hw_coul
        )

    # accumulate observables
    if (task.get("observables")):
        targets.update(task.get("observables"))

    # target: radius squared
    if ("tbme-rrel2" not in targets.keys()):
        targets["tbme-rrel2"] = operators.rrel2(A, hw)
    # target: Ncm
    if ("tbme-Ncm" not in targets.keys()):
        targets["tbme-Ncm"] = operators.Ncm(A, hw/hw_cm)

    # optional observable sets
    # Hamiltonian components
    if ("H-components" in task["observable_sets"]):
        # target: Trel (diagnostic)
        targets["tbme-Trel"] = operators.Trel(A, hw)
        # target: VNN (diagnostic)
        targets["tbme-VNN"] = operators.VNN()
        # target: VC (diagnostic)
        if (task["use_coulomb"]):
            targets["tbme-VC"] = operators.VC(hw_coul_rescaled/hw_coul)
    # squared angular momenta
    if ("am-sqr" in task["observable_sets"]):
        targets["tbme-L"] = operators.L()
        targets["tbme-Sp"] = operators.Sp()
        targets["tbme-Sn"] = operators.Sn()
        targets["tbme-S"] = operators.S()
        targets["tbme-J"] = operators.J()

    # get set of required sources
    required_sources = set()
    required_sources.update(*[op.keys() for op in targets.values()])

    # accumulate h2mixer input lines
    lines = []

    # initial comment
    lines.append("# task: {}".format(task))
    lines.append("")

    # global mode definitions
    target_truncation = task["target_truncation"]
    if (target_truncation is None):
        # automatic derivation
        if (task["truncation_mode"] is config.TruncationMode.kHO):
            truncation_parameters = task["truncation_parameters"]
            if (truncation_parameters["many_body_truncation"]=="Nmax"):
                # important: truncation of orbitals file, one-body
                # truncation of interaction file, and MFDn
                # single-particle shells (beware 1-based) must agree
                N1_max = truncation_parameters["Nv"]+truncation_parameters["Nmax"]
                N2_max = 2*truncation_parameters["Nv"]+truncation_parameters["Nmax"]
                target_weight_max = utils.weight_max_string((N1_max,N2_max))
            elif (truncation_parameters["many_body_truncation"]=="FCI"):
                N1_max = truncation_parameters["Nmax"]
                target_weight_max = utils.weight_max_string(("ob",N1_max))
        else:
            # calculation of required weight_max will require external program for occupation counting
            pass
    else:
        # given value
        target_weight_max = target_truncation
    lines.append("set-target-indexing {orbitals_filename} {target_weight_max}".format(
        orbitals_filename=config.filenames.orbitals_filename(natural_orbital_iteration),
        target_weight_max=target_weight_max,
        **task
    ))
    lines.append("set-target-multipolarity 0 0 0")
    lines.append("set-output-format {h2_format}".format(**task))
    lines.append("set-mass {A}".format(A=A,**task))
    lines.append("")

    # radial operator inputs
    for operator_type in ["r","k"]:
        for power in [1,2]:
            radial_me_filename = config.filenames.radial_me_filename(natural_orbital_iteration, operator_type, power)
            lines.append("define-radial-operator {} {} {}".format(operator_type,power,radial_me_filename))
    lines.append("")

    # sources: h2mixer built-ins
    builtin_sources = (operators.kinematic_operator_set | operators.angular_momentum_operator_set)
    for source in sorted(builtin_sources & required_sources):
        lines.append("define-source operator " + source)
    lines.append("")

    # sources: VNN
    if ("VNN" in required_sources):
        VNN_filename = config.environ.interaction_filename(
            "{}-{}-{:g}.bin".format(
                task["interaction"],
                mcscript.utils.dashify(task["truncation_int"]),
                task["hw_int"]
            )
        )
        if (task["basis_mode"]==config.BasisMode.kDirect and natural_orbital_iteration in {None,0}):
            lines.append("define-source input VNN {VNN_filename}".format(VNN_filename=VNN_filename,**task))
        else:
            xform_weight_max_int = utils.weight_max_string(xform_truncation_int)
            lines.append("define-source xform VNN {VNN_filename} {xform_weight_max_int} {radial_olap_int_filename}".format(
                VNN_filename=VNN_filename,
                xform_weight_max_int=xform_weight_max_int,
                radial_olap_int_filename=config.filenames.radial_olap_int_filename(natural_orbital_iteration),
                **task
            ))

    # sources: Coulomb
    #
    # Note: This is the "unscaled" Coulomb, still awaiting the scaling
    # factor from dilation.
    if ("VC_unscaled" in required_sources):
        VC_filename = config.environ.interaction_filename(
            "{}-{}-{:g}.bin".format(
                "VC",
                mcscript.utils.dashify(task["truncation_coul"]),
                task["hw_coul"]
            )
        )
        if (task["basis_mode"] in {config.BasisMode.kDirect,config.BasisMode.kDilated} and natural_orbital_iteration in {None,0}):
            lines.append("define-source input VC_unscaled {VC_filename}".format(VC_filename=VC_filename,**task))
        else:
            xform_weight_max_coul = utils.weight_max_string(xform_truncation_coul)
            lines.append("define-source xform VC_unscaled {VC_filename} {xform_weight_max_coul} {radial_olap_coul_filename}".format(
                VC_filename=VC_filename,
                xform_weight_max_coul=xform_weight_max_coul,
                radial_olap_coul_filename=config.filenames.radial_olap_coul_filename(natural_orbital_iteration),
                **task
            ))

    lines.append("")


    # targets: generate h2mixer input
    for (basename, operator) in targets.items():
        lines.append("define-target work/"+basename+".bin")
        for (source, coefficient) in operator.items():
            lines.append("  add-source {:s} {:e}".format(source, coefficient))
        lines.append("")

    # ensure terminal line
    lines.append("")

    # diagnostic: log input lines to file
    #
    # This is purely for easy diagnostic purposes, since lines will be
    # fed directly to h2mixer as stdin below.
    mcscript.utils.write_input(config.filenames.h2mixer_filename(natural_orbital_iteration),input_lines=lines,verbose=False)

    # create work directory if it doesn't exist yet (-p)
    mcscript.call(["mkdir","-p","work"])

    # invoke h2mixer
    mcscript.call(
        [
            config.environ.shell_filename("h2mixer")
        ],
        input_lines=lines,
        mode = mcscript.call.serial
    )
