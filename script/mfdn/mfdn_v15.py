"""mfdn_v15.py -- driver module for MFDn v15.

Patrick Fasano
University of Notre Dame

- 06/08/17 (pjf): Initial version of MFDn v15 scripting, built from mfdn_v14.py.
- 06/22/17 (pjf): Update references to mcscript.exception.ScriptError.
- 07/31/17 (pjf): Create work directory if nonexistent when running mfdn.
- 08/11/17 (pjf):
  + Use new TruncationModes.
  + Fix FCI truncation.
- 08/25/17 (pjf): Implement general truncation:
  + Always use orbital file, even for Nmax runs.
  + Break out different many-body truncations into their own functions.
- 08/26/17 (pjf): Add parity constraints for general truncations.
"""
import os
import glob
import collections

import mcscript

from . import config


def set_up_Nmax_truncation(task, inputlist):
    """Generate Nmax truncation inputs for MFDn v15.

    Arguments:
        task(dict): as described in module docstring
        inputlist (dict): MFDn v15 inputlist
    """
    # sanity check
    if task["sp_truncation_mode"] is not config.SingleParticleTruncationMode.kNmax:
        raise ValueError("expecting sp_truncation_mode to be {} but found {sp_truncation_mode}".format(config.SingleParticleTruncationMode.kNmax, **task))
    if task["mb_truncation_mode"] is not config.ManyBodyTruncationMode.kNmax:
        raise ValueError("expecting mb_truncation_mode to be {} but found {mb_truncation_mode}".format(config.ManyBodyTruncationMode.kNmax, **task))

    truncation_parameters = task["truncation_parameters"]

    # many-body truncation: Nmin, Nmax, deltaN
    if (truncation_parameters["Nstep"] == 2):
        inputlist["Nmin"] = truncation_parameters["Nmax"] % 2
    else:
        inputlist["Nmin"] = 1
    inputlist["Nmax"] = int(truncation_parameters["Nmax"])
    inputlist["deltaN"] = int(truncation_parameters["Nstep"])


def set_up_WeightMax_truncation(task, inputlist):
    """Generate weight max truncation inputs for MFDn v15.

    Arguments:
        task(dict): as described in module docstring
        inputlist (dict): MFDn v15 inputlist
    """
    # sanity check
    if task["mb_truncation_mode"] is not config.ManyBodyTruncationMode.kWeightMax:
        raise ValueError("expecting sp_truncation_mode to be {} but found {sp_truncation_mode}".format(config.ManyBodyTruncationMode.kWeightMax, **task))

    truncation_parameters = task["truncation_parameters"]

    inputlist["WTmax"] = truncation_parameters["mb_weight_max"]
    inputlist["parity"] = truncation_parameters["parity"]


def set_up_FCI_truncation(task, inputlist):
    """Generate FCI truncation inputs for MFDn v15.

    Arguments:
        task(dict): as described in module docstring
        inputlist (dict): MFDn v15 inputlist
    """
    # sanity check
    if task["mb_truncation_mode"] is not config.ManyBodyTruncationMode.kFCI:
        raise ValueError("expecting sp_truncation_mode to be {} but found {sp_truncation_mode}".format(config.ManyBodyTruncationMode.kFCI, **task))

    truncation_parameters = task["truncation_parameters"]

    # maximum weight of an orbital is either Nmax or sp_weight_max
    if task["sp_truncation_mode"] is config.SingleParticleTruncationMode.kNmax:
        max_sp_weight = truncation_parameters["Nmax"]
        parity = (-1)**(truncation_parameters["Nmax"] % truncation_parameters["Nstep"])
    else:
        max_sp_weight = truncation_parameters["sp_weight_max"]
        parity = truncation_parameters.get("parity", 0)

    inputlist["WTmax"] = sum(task["nuclide"])*max_sp_weight
    inputlist["parity"] = int(parity)


truncation_setup_functions = {
    config.ManyBodyTruncationMode.kNmax: set_up_Nmax_truncation,
    config.ManyBodyTruncationMode.kWeightMax: set_up_WeightMax_truncation,
    config.ManyBodyTruncationMode.kFCI: set_up_FCI_truncation
}


def run_mfdn(task, postfix=""):
    """Generate input file and execute MFDn version 15 beta 00.

    Arguments:
        task (dict): as described in module docstring
        postfix (str, optional):
    Raises:
        mcscript.exception.ScriptError: if MFDn output not found
    """
    # create work directory if it doesn't exist yet (-p)
    mcscript.call(["mkdir", "-p", "work"])

    # inputlist namelist dictionary
    inputlist = collections.OrderedDict()

    # nucleus
    inputlist["Nprotons"], inputlist["Nneutrons"] = task["nuclide"]

    # Mj
    inputlist["TwoMj"] = int(2*task["Mj"])

    # single-particle orbitals
    inputlist["orbitalfile"] = config.filenames.orbitals_filename(postfix)
    mcscript.call([
        "cp", "--verbose",
        config.filenames.orbitals_filename(postfix),
        os.path.join("work", config.filenames.orbitals_filename(postfix))
    ])

    # truncation mode
    truncation_setup_functions[task["mb_truncation_mode"]](task, inputlist)

    if (task["basis_mode"] in {config.BasisMode.kDirect, config.BasisMode.kDilated}):
        inputlist["hbomeg"] = float(task["hw"])

    # diagonalization parameters
    inputlist["neivals"] = int(task["eigenvectors"])
    inputlist["maxits"] = int(task["lanczos"])
    inputlist["tol"] = float(task["tolerance"])

    # Hamiltonian input
    inputlist["TBMEfile"] = "tbme-H"

    # tbo: two-body observable namelist
    obslist = {}

    # tbo: collect tbo names
    obs_basename_list = ["tbme-rrel2"]
    if ("H-components" in task["observable_sets"]):
        obs_basename_list += ["tbme-Trel", "tbme-Ncm", "tbme-VNN"]
        if (task["use_coulomb"]):
            obs_basename_list += ["tbme-VC"]
    if ("am-sqr" in task["observable_sets"]):
        obs_basename_list += ["tbme-L", "tbme-Sp", "tbme-Sn", "tbme-S", "tbme-J"]
    if ("observables" in task):
        obs_basename_list += list(task["observables"].keys())

    # tbo: log tbo names in separate file to aid future data analysis
    mcscript.utils.write_input("tbo_names{:s}.dat".format(postfix), input_lines=obs_basename_list)

    # tbo: count number of observables
    num_obs = len(obs_basename_list)
    if (num_obs > 32):
        raise mcscript.exception.ScriptError("Too many observables for MFDn v15")

    inputlist["numTBops"] = num_obs
    obslist["TBMEoperators"] = obs_basename_list

    # obdme: parameters
    inputlist["obdme"] = True
    obslist["max2K"] = int(2*task["obdme_multipolarity"])

    # construct transition observable input if reference states given
    if task.get("obdme_reference_state_list") is not None:
        # obdme: validate reference state list
        #
        # guard against pathetically common mistakes
        for (J, g_rel, i) in task["obdme_reference_state_list"]:
            # validate integer/half-integer character of angular momentum
            twice_J = int(2*J)
            if ((twice_J % 2) != (sum(task["nuclide"]) % 2)):
                raise ValueError("invalid angular momentum for reference state")
            # validate grade (here taken relative to natural grade)
            # if ((g_rel != (truncation_parameters["Nmax"] % 2)) and (truncation_parameters["Nstep"] != 1)):
            #     raise ValueError("invalid parity for reference state")

        # obdme: construct input
        inputlist["nrefstates"] = len(task["obdme_reference_state_list"])
        obslist["ref2J"] = []
        obslist["refseq"] = []
        for (J, g_rel, i) in task["obdme_reference_state_list"]:
            obslist["ref2J"].append(int(2*J))
            obslist["refseq"].append(i)

    # generate MFDn input file
    mcscript.utils.write_namelist(
        os.path.join("work", "mfdn.input"), input_dict={"inputlist": inputlist, "obslist": obslist}
    )

    # mcscript.utils.write_input("work/mfdn.dat",input_lines=lines)

    # import partitioning file
    if (task["partition_filename"] is not None):
        if (not os.path.exists(task["partition_filename"])):
            raise mcscript.exception.ScriptError("partition file not found")
        mcscript.call(["cp", "--verbose", task["partition_filename"], "work/mfdn_partitioning.info"])

    # enter work directory
    os.chdir("work")

    # invoke MFDn
    mcscript.call(
        [
            config.environ.mfdn_filename(task["mfdn_executable"])
        ],
        mode=mcscript.CallMode.kHybrid,
        check_return=True
    )

    # test for basic indications of success
    if (not os.path.exists("mfdn.out")):
        raise mcscript.exception.ScriptError("mfdn.out not found")
    if (not os.path.exists("mfdn.res")):
        raise mcscript.exception.ScriptError("mfdn.res not found")

    # leave work directory
    os.chdir("..")


def save_mfdn_output(task, postfix=""):
    """Generate input file and execute MFDn version 14 beta 06.

    Arguments:
        task (dict): as described in module docstring

    Raises:
        mcscript.exception.ScriptError: if MFDn output not found
    """
    # save quick inspection copies of mfdn.{res,out}
    natural_orbitals = task.get("natural_orbitals")
    descriptor = task["metadata"]["descriptor"]
    print("Saving basic output files...")
    res_filename = "{:s}-mfdn15-{:s}{:s}.res".format(mcscript.parameters.run.name, descriptor, postfix)
    mcscript.call(["cp", "--verbose", "work/mfdn.res", res_filename])
    out_filename = "{:s}-mfdn15-{:s}{:s}.out".format(mcscript.parameters.run.name, descriptor, postfix)
    mcscript.call(["cp", "--verbose", "work/mfdn.out", out_filename])

    # save OBDME files for next natural orbital iteration
    if natural_orbitals:
        print("Saving OBDME files for natural orbital generation...")
        obdme_info_filename = "mfdn.rppobdme.info"
        mcscript.call(
            [
                "cp", "--verbose",
                "work/{}".format(obdme_info_filename),
                config.filenames.natorb_info_filename(postfix)
            ]
        )
        obdme_filename = glob.glob("work/mfdn.statrobdme.seq{:03d}*".format(task["natorb_base_state"]))
        mcscript.call(
            [
                "cp", "--verbose",
                obdme_filename[0],
                config.filenames.natorb_obdme_filename(postfix)
            ]
        )

    # save full archive of input, log, and output files
    print("Saving full output files...")
    # logging
    archive_file_list = [
        config.filenames.h2mixer_filename(postfix),
        "tbo_names{:s}.dat".format(postfix)
        ]
    # orbital information
    archive_file_list += [
        config.filenames.orbitals_int_filename(postfix),
        config.filenames.orbitals_filename(postfix),
        ]
    # transformation information
    archive_file_list += [
        config.filenames.radial_xform_filename(postfix),
        # config.filenames.radial_me_filename(postfix, operator_type, power),
        config.filenames.radial_olap_int_filename(postfix),
        ]
    # Coulomb information:
    if task["use_coulomb"]:
        archive_file_list += [
            config.filenames.orbitals_coul_filename(postfix),
            config.filenames.radial_olap_coul_filename(postfix),
        ]
    # natural orbital information
    if natural_orbitals:
        archive_file_list += [
            config.filenames.natorb_info_filename(postfix),
            config.filenames.natorb_obdme_filename(postfix),
            ]
        # glob for natural orbital xform
        archive_file_list += glob.glob(config.filenames.natorb_xform_filename(postfix))
    # MFDn output
    archive_file_list += [
        "work/mfdn.input", "work/mfdn.out", "work/mfdn.res",
        "work/mfdn_partitioning.generated"
    ]
    if os.path.isfile("work/mfdn_sp_orbitals.info"):
        archive_file_list += ["work/mfdn_sp_orbitals.info"]
    # renamed versions
    archive_file_list += [out_filename, res_filename]
    # MFDN obdme
    if (task["save_obdme"]):
        archive_file_list += glob.glob("work/*obdme*")
    # generate archive (outside work directory)
    archive_filename = "{:s}-mfdn15-{:s}{:s}.tgz".format(mcscript.parameters.run.name, descriptor, postfix)
    mcscript.call(
        [
            "tar", "zcvf", archive_filename,
            "--transform=s,work/,,",
            "--show-transformed"
        ] + archive_file_list
    )

    # copy results out (if in multi-task run)
    if (mcscript.task.results_dir is not None):
        mcscript.call(
            [
                "cp",
                "--verbose",
                res_filename, out_filename, archive_filename,
                "--target-directory={}".format(mcscript.task.results_dir)
            ]
        )

    # cleanup of wave function files
    scratch_file_list = glob.glob("work/*")
    mcscript.call(["rm", "-vf"] + scratch_file_list)
