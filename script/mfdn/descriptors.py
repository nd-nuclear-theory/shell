"""descriptors.py -- task descriptors for MFDn runs.

Patrick Fasano
University of Notre Dame

- 3/22/17 (pjf): Created, split from __init__.py.
- 4/7/17 (pjf): Reference new config submodule.
- 6/3/17 (pjf): Switch natorb_iteration -> natural_orbitals.
"""
import mcscript.exception
import mcscript.utils

from . import config


################################################################
# task descriptor for h2mixer + mfdn run
################################################################

def task_descriptor_7(task):
    """Task descriptor format 7

        Overhaul for new h2utils scripting:
        - Strip back down to basic form for oscillator-like runs only.
        - Adjust some field labels.
        - Add tolerance.
    """
    if(
        task["truncation_mode"] is config.TruncationMode.kHO
        and
        task["basis_mode"] in {config.BasisMode.kDirect, config.BasisMode.kDilated}
    ):
        # traditional oscillator run
        template_string = (
            "Z{nuclide[0]}-N{nuclide[1]}-{interaction}-coul{coulomb_flag:d}"
            "-hw{hw:06.3f}"
            "-a_cm{a_cm:g}"
            "-Nmax{Nmax:02d}{mixed_parity_indicator}{fci_indicator}-Mj{Mj:03.1f}"
            "-lan{lanczos:d}-tol{tolerance:.1e}"
            "{natural_orbital_indicator}"
            )
    else:
        raise mcscript.exception.ScriptError("mode not supported by task descriptor")

    truncation_parameters = task["truncation_parameters"]
    if (truncation_parameters["many_body_truncation"] == "fci"):
        fci_indicator = "-fci"
    else:
        fci_indicator = ""
    mixed_parity_indicator = mcscript.utils.ifelse(truncation_parameters["Nstep"] == 1, "x", "")
    coulomb_flag = int(task["use_coulomb"])
    natural_orbital_indicator = mcscript.utils.ifelse(task.get("natural_orbitals"), "-natorb", "")

    descriptor = template_string.format(
        coulomb_flag=coulomb_flag,
        mixed_parity_indicator=mixed_parity_indicator,
        fci_indicator=fci_indicator,
        natural_orbital_indicator=natural_orbital_indicator,
        **mcscript.utils.dict_union(task, truncation_parameters)
        )

    return descriptor
