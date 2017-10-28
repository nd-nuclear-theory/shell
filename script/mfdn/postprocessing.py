"""postprocessing.py -- perform postprocessing and observable calculation tasks

Patrick Fasano
University of Notre Dame

- 10/25/17 (pjf): Created.
"""
import os
import glob
import re

import mcscript

from . import modes, environ, utils


def generate_em(task, postfix=""):
    """Generate electromagnetic matrix elements.

    Arguments:
        task (dict): as described in module docstring
        postfix (string, optional): identifier to add to generated files
    """

    # accumulate em-gen input lines
    lines = []

    # set up orbitals
    lines += [
        "set-indexing {:s}".format(environ.filenames.orbitals_filename(postfix)),
        "set-basis-scale-factor {:e}".format(utils.oscillator_length(task["hw"])),
        ]

    for (operator_type, order) in task.get("ob_observables", []):
        lines.append("define-radial-source {:s}".format(
            environ.filenames.radial_me_filename(postfix, operator_type, order)
            ))
        for species in ["p", "n"]:
            if operator_type == "E":
                lines.append(
                    "define-target E {order:d} {species:s} {output_filename:s}".format(
                        order=order, species=species,
                        output_filename=environ.filenames.observable_me_filename(postfix, operator_type, order, species)
                        )
                    )
            elif operator_type == "M":
                lines.append(
                    "define-target Dl {order:d} {species:s} {output_filename:s}".format(
                        order=order, species=species,
                        output_filename=environ.filenames.observable_me_filename(postfix, "Dl", order, species)
                        )
                    )
                lines.append(
                    "define-target Ds {order:d} {species:s} {output_filename:s}".format(
                        order=order, species=species,
                        output_filename=environ.filenames.observable_me_filename(postfix, "Ds", order, species)
                        )
                    )

    # ensure trailing line
    lines.append("")

    # write input file
    mcscript.utils.write_input(
        environ.filenames.emgen_filename(postfix),
        input_lines=lines,
        verbose=False
        )

    # invoke em-gen
    mcscript.call(
        [
            environ.environ.shell_filename("em-gen")
        ],
        input_lines=lines,
        mode=mcscript.CallMode.kSerial
    )

def evaluate_ob_observables(task, postfix=""):
    """Evaluate one-body observables with obscalc-ob.

    Arguments:
        task (dict): as described in module docstring
        postfix (string, optional): identifier to add to generated files
    """

    work_dir = "work{:s}".format(postfix)

    # accumulate obscalc-ob input lines
    lines = []

    # initial comment
    lines.append("# task: {}".format(task))
    lines.append("")

    # indexing setup
    lines += [
        "set-indexing {:s}".format(environ.filenames.orbitals_filename(postfix)),
        "set-robdme-info {:s}".format(os.path.join(work_dir, "mfdn.rppobdme.info")),
        "set-output-file {:s}".format(environ.filenames.obscalc_ob_res_filename(postfix)),
        ]

    # set up operators
    for (operator_type, order) in task.get("ob_observables", []):
        lines.append("define-radial-source {:s}".format(
            environ.filenames.radial_me_filename(postfix, operator_type, order)
            ))
        for species in ["p", "n"]:
            if operator_type == "M":
                # convenience definition for M observable
                lines.append(
                    "define-operator Dl({:s}) {:s}".format(
                        species,
                        environ.filenames.observable_me_filename(postfix, "Dl", order, species)
                        )
                    )
                lines.append(
                    "define-operator Ds({:s}) {:s}".format(
                        species,
                        environ.filenames.observable_me_filename(postfix, "Ds", order, species)
                        )
                    )
            else:
                lines.append(
                    "define-operator {:s}({:s}) {:s}".format(
                        operator_type,
                        species,
                        environ.filenames.observable_me_filename(postfix, operator_type, order, species)
                        )
                    )

    # get filenames for static densities and extract quantum numbers
    filenames = glob.glob(os.path.join(work_dir, "mfdn.statrobdme.*"))
    regex = re.compile(
        # directory prefix
        r"{}".format(os.path.join(work_dir, "")) +
        # prolog
        r"mfdn\.statrobdme"
        # sequence number
        r"\.seq(?P<seq>\d{3})"
        # 2J
        r"\.2J(?P<twoJ>\d{2})"
        # parity (v14 only)
        r"(\.p(?P<g>\d))?"
        # n
        r"\.n(?P<n>\d{2})"
        # 2T
        r"\.2T(?P<twoT>\d{2})"
        )
    conversions = {
        "seq": int,
        "twoJ": int,
        "g": lambda x: int(x) if x is not None else 0,
        "n": int,
        "twoT": int
        }
    statrobdme_files = []
    for filename in filenames:
        match = regex.match(filename)
        if match is None:
            print(regex)
            raise ValueError("bad statrobdme filename: {}".format(filename))
        info = match.groupdict()

        # convert fields
        for key in info:
            conversion = conversions[key]
            info[key] = conversion(info[key])

        statrobdme_files.append(mcscript.utils.dict_union(info, {"filename": filename}))

    # sort states by sequence number
    statrobdme_files.sort(key=lambda item: item["seq"])
    for statrobdme_file in statrobdme_files:
        lines.append(
            "define-static-densities {twoJ:d} {g:d} {n:d} {filename:s}".format(**statrobdme_file)
            )

    # define-transition-densities 2Jf gf nf 2Ji gi fi robdme_info_filename robdme_filename
    # get filenames for static densities and extract quantum numbers
    filenames = glob.glob(os.path.join(work_dir, "mfdn.robdme.*"))
    regex = re.compile(
        r"{}".format(os.path.join(work_dir, "")) +
        # prolog
        r"mfdn\.robdme"
        # final sequence number
        r"\.seq(?P<seqf>\d{3})"
        # final 2J
        r"\.2J(?P<twoJf>\d{2})"
        # final parity (v14 only)
        r"(\.p(?P<gf>\d))?"
        # final n
        r"\.n(?P<nf>\d{2})"
        # final 2T
        r"\.2T(?P<twoTf>\d{2})"
        # initial sequence number
        r"\.seq(?P<seqi>\d{3})"
        # initial 2J
        r"\.2J(?P<twoJi>\d{2})"
        # initial parity (v14 only)
        r"(\.p(?P<gi>\d))?"
        # initial n
        r"\.n(?P<ni>\d{2})"
        # inital 2T
        r"\.2T(?P<twoTi>\d{2})"
        )
    conversions = {
        "seqf": int,
        "twoJf": int,
        "gf": lambda x: int(x) if x is not None else 0,
        "nf": int,
        "twoTf": int,
        "seqi": int,
        "twoJi": int,
        "gi": lambda x: int(x) if x is not None else 0,
        "ni": int,
        "twoTi": int
        }
    robdme_files = []
    for filename in filenames:
        match = regex.match(filename)
        if match is None:
            raise ValueError("bad statrobdme filename format")
        info = match.groupdict()

        # convert fields
        for key in info:
            conversion = conversions[key]
            info[key] = conversion(info[key])

        if "gf" not in info:
            info["gf"] = 0
        if "gi" not in info:
            info["gi"] = 0

        robdme_files.append(mcscript.utils.dict_union(info, {"filename": filename}))

    # sort by sequence number of final state, then sequence number of initial state
    robdme_files.sort(key=lambda item: (item["seqf"], item["seqi"]))
    for robdme_file in robdme_files:
        lines.append(
            "define-transition-densities {twoJf:d} {gf:d} {nf:d} {twoJi:d} {gi:d} {ni:d} {filename:s}".format(**robdme_file)
            )

    # ensure trailing line
    lines.append("")

    # write input file
    mcscript.utils.write_input(
        environ.filenames.obscalc_ob_filename(postfix),
        input_lines=lines,
        verbose=False
        )

    # invoke em-gen
    mcscript.call(
        [
            environ.environ.shell_filename("obscalc-ob")
        ],
        input_lines=lines,
        mode=mcscript.CallMode.kSerial
    )
