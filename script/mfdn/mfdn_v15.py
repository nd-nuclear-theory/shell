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
- 09/12/17 (pjf): Update for config -> modes + environ split.
- 09/22/17 (pjf): Take "observables" as list of tuples instead of dict.
- 09/24/17 (pjf):
  + Archive wavefunction files in separate archive.
  + Create tar files with minimal directory structure for easy inflation.
- 09/27/17 (pjf): Allow for counting modes with run_mode.
- 10/05/17 (pjf): Add save_mfdn_output_out_only.
- 10/18/17 (pjf):
  + Use separate work directory for each postfix.
  + Factor out extract_natural_orbitals().
- 10/25/17 (pjf): Rename "observables" to "tb_observables".
"""
import os
import glob
import collections

import mcscript

from . import modes, environ


def set_up_Nmax_truncation(task, inputlist):
    """Generate Nmax truncation inputs for MFDn v15.

    Arguments:
        task(dict): as described in module docstring
        inputlist (dict): MFDn v15 inputlist
    """
    # sanity check
    if task["sp_truncation_mode"] is not modes.SingleParticleTruncationMode.kNmax:
        raise ValueError("expecting sp_truncation_mode to be {} but found {sp_truncation_mode}".format(modes.SingleParticleTruncationMode.kNmax, **task))
    if task["mb_truncation_mode"] is not modes.ManyBodyTruncationMode.kNmax:
        raise ValueError("expecting mb_truncation_mode to be {} but found {mb_truncation_mode}".format(modes.ManyBodyTruncationMode.kNmax, **task))

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
    if task["mb_truncation_mode"] is not modes.ManyBodyTruncationMode.kWeightMax:
        raise ValueError("expecting sp_truncation_mode to be {} but found {sp_truncation_mode}".format(modes.ManyBodyTruncationMode.kWeightMax, **task))

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
    if task["mb_truncation_mode"] is not modes.ManyBodyTruncationMode.kFCI:
        raise ValueError("expecting sp_truncation_mode to be {} but found {sp_truncation_mode}".format(modes.ManyBodyTruncationMode.kFCI, **task))

    truncation_parameters = task["truncation_parameters"]

    # maximum weight of an orbital is either Nmax or sp_weight_max
    if task["sp_truncation_mode"] is modes.SingleParticleTruncationMode.kNmax:
        max_sp_weight = truncation_parameters["Nmax"]
        parity = (-1)**(truncation_parameters["Nmax"] % truncation_parameters["Nstep"])
    else:
        max_sp_weight = truncation_parameters["sp_weight_max"]
        parity = truncation_parameters.get("parity", 0)

    inputlist["WTmax"] = sum(task["nuclide"])*max_sp_weight
    inputlist["parity"] = int(parity)


truncation_setup_functions = {
    modes.ManyBodyTruncationMode.kNmax: set_up_Nmax_truncation,
    modes.ManyBodyTruncationMode.kWeightMax: set_up_WeightMax_truncation,
    modes.ManyBodyTruncationMode.kFCI: set_up_FCI_truncation
}


def run_mfdn(task, run_mode=modes.MFDnRunMode.kNormal, postfix=""):
    """Generate input file and execute MFDn version 15 beta 00.

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

    # run mode
    inputlist["IFLAG_mode"] = int(run_mode)

    # nucleus
    inputlist["Nprotons"], inputlist["Nneutrons"] = task["nuclide"]

    # Mj
    inputlist["TwoMj"] = int(2*task["Mj"])

    # single-particle orbitals
    inputlist["orbitalfile"] = environ.filenames.orbitals_filename(postfix)
    mcscript.call([
        "cp", "--verbose",
        environ.filenames.orbitals_filename(postfix),
        os.path.join(work_dir, environ.filenames.orbitals_filename(postfix))
    ])

    # truncation mode
    truncation_setup_functions[task["mb_truncation_mode"]](task, inputlist)

    if run_mode is modes.MFDnRunMode.kNormal:
        if (task["basis_mode"] in {modes.BasisMode.kDirect, modes.BasisMode.kDilated}):
            inputlist["hbomeg"] = float(task["hw"])

        # diagonalization parameters
        inputlist["neivals"] = int(task["eigenvectors"])
        inputlist["maxits"] = int(task["lanczos"])
        inputlist["tol"] = float(task["tolerance"])
        if task.get("reduce_solver_threads"):
            inputlist["reduce_solver_threads"] = task["reduce_solver_threads"]

        # Hamiltonian input
        inputlist["TBMEfile"] = "tbme-H"

        # tbo: collect tbo names
        obs_basename_list = ["tbme-rrel2"]
        if ("H-components" in task["observable_sets"]):
            obs_basename_list += ["tbme-Trel", "tbme-Ncm", "tbme-VNN"]
            if (task["use_coulomb"]):
                obs_basename_list += ["tbme-VC"]
        if ("am-sqr" in task["observable_sets"]):
            obs_basename_list += ["tbme-L", "tbme-Sp", "tbme-Sn", "tbme-S", "tbme-J"]
        if ("isospin" in task["observable_sets"]):
            obs_basename_list += ["tbme-T"]
        if ("tb_observables" in task):
            obs_basename_list += [basename for (basename, operator) in task["tb_observables"]]

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
        os.path.join(work_dir, "mfdn.input"),
        input_dict={"inputlist": inputlist, "obslist": obslist}
    )

    # import partitioning file
    if (task["partition_filename"] is not None):
        if (not os.path.exists(task["partition_filename"])):
            raise mcscript.exception.ScriptError("partition file not found")
        mcscript.call([
            "cp", "--verbose",
            task["partition_filename"],
            os.path.join(work_dir, "mfdn_partitioning.info")
            ])

    # enter work directory
    os.chdir(work_dir)

    # invoke MFDn
    mcscript.call(
        [
            environ.environ.mfdn_filename(task["mfdn_executable"])
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


def extract_natural_orbitals(task, postfix=""):
    """Extract OBDME files for subsequent natural orbital iterations.

    Arguments:
        task (dict): as described in module docstring
        postfix (string, optional): identifier to add to generated files
    """
    # save OBDME files for next natural orbital iteration
    if not task.get("natural_orbitals"):
        raise mcscript.exception.ScriptError("natural orbitals not enabled")

    work_dir = "work{:s}".format(postfix)
    obdme_info_filename = "mfdn.rppobdme.info"
    try:
        (J, g, n) = task["natorb_base_state"]
        obdme_filename = glob.glob(
            "{:s}/mfdn.statrobdme.seq*.2J{:02d}.n{:02d}.2T*".format(work_dir, 2*J, n)
            )
    except TypeError:
        obdme_filename = glob.glob(
            "{:s}/mfdn.statrobdme.seq{:03d}*".format(work_dir, task["natorb_base_state"])
            )

    print("Saving OBDME files for natural orbital generation...")
    mcscript.call(
        [
            "cp", "--verbose",
            os.path.join(work_dir, obdme_info_filename),
            environ.filenames.natorb_info_filename(postfix)
        ]
    )
    mcscript.call(
        [
            "cp", "--verbose",
            obdme_filename[0],
            environ.filenames.natorb_obdme_filename(postfix)
        ]
    )


def save_mfdn_output_out_only(task, postfix=""):
    """Collect and save MFDn output files only.

    Arguments:
        task (dict): as described in module docstring
        postfix (string, optional): identifier to add to generated files
    """
    # save quick inspection copies of mfdn.{res,out}
    descriptor = task["metadata"]["descriptor"]
    print("Saving basic output files...")
    work_dir = "work{:s}".format(postfix)
    filename_prefix = "{:s}-mfdn15-{:s}{:s}".format(mcscript.parameters.run.name, descriptor, postfix)
    res_filename = "{:s}.res".format(filename_prefix)
    mcscript.call(["cp", "--verbose", work_dir+"/mfdn.res", res_filename])
    out_filename = "{:s}.out".format(filename_prefix)
    mcscript.call(["cp", "--verbose", work_dir+"/mfdn.out", out_filename])

    # copy results out (if in multi-task run)
    if (mcscript.task.results_dir is not None):
        mcscript.call(
            [
                "cp",
                "--verbose",
                res_filename, out_filename,
                "--target-directory={}".format(mcscript.task.results_dir)
            ]
        )


def save_mfdn_output(task, postfix=""):
    """Collect and save MFDn output.

    Arguments:
        task (dict): as described in module docstring
        postfix (string, optional): identifier to add to generated files
    """
    # save quick inspection copies of mfdn.{res,out}
    descriptor = task["metadata"]["descriptor"]
    work_dir = "work{:s}".format(postfix)
    print("Saving basic output files...")
    filename_prefix = "{:s}-mfdn15-{:s}{:s}".format(mcscript.parameters.run.name, descriptor, postfix)
    res_filename = "{:s}.res".format(filename_prefix)
    mcscript.call(["cp", "--verbose", work_dir+"/mfdn.res", res_filename])
    out_filename = "{:s}.out".format(filename_prefix)
    mcscript.call(["cp", "--verbose", work_dir+"/mfdn.out", out_filename])

    # save full archive of input, log, and output files
    print("Saving full output files...")
    # logging
    archive_file_list = [
        environ.filenames.h2mixer_filename(postfix),
        "tbo_names{:s}.dat".format(postfix)
        ]
    # orbital information
    archive_file_list += [
        environ.filenames.orbitals_int_filename(postfix),
        environ.filenames.orbitals_filename(postfix),
        ]
    # transformation information
    archive_file_list += [
        environ.filenames.radial_xform_filename(postfix),
        # environ.filenames.radial_me_filename(postfix, operator_type, power),
        environ.filenames.radial_olap_int_filename(postfix),
        ]
    # Coulomb information:
    if task["use_coulomb"]:
        archive_file_list += [
            environ.filenames.orbitals_coul_filename(postfix),
            environ.filenames.radial_olap_coul_filename(postfix),
        ]
    # natural orbital information
    if task.get("natural_orbitals"):
        archive_file_list += [
            environ.filenames.natorb_info_filename(postfix),
            environ.filenames.natorb_obdme_filename(postfix),
            ]
        # glob for natural orbital xform
        archive_file_list += glob.glob(environ.filenames.natorb_xform_filename(postfix))
    # MFDn output
    archive_file_list += [
        work_dir+"/mfdn.input", work_dir+"/mfdn.out", work_dir+"/mfdn.res",
        work_dir+"/mfdn_partitioning.generated"
    ]
    if os.path.isfile(work_dir+"/mfdn_sp_orbitals.info"):
        archive_file_list += [work_dir+"/mfdn_sp_orbitals.info"]
    # renamed versions
    archive_file_list += [out_filename, res_filename]
    # MFDN obdme
    if (task["save_obdme"]):
        archive_file_list += glob.glob(work_dir+"/*obdme*")
    # observable output
    archive_file_list += glob.glob("em-gen.in")
    archive_file_list += glob.glob("obscalc-ob.*")
    # generate archive (outside work directory)
    archive_filename = "{:s}.tgz".format(filename_prefix)
    mcscript.call(
        [
            "tar", "zcvf", archive_filename,
            "--transform=s,{:s}/,,".format(work_dir),
            "--transform=s,^,{:s}/{:s}{:s}/,".format(mcscript.parameters.run.name, descriptor, postfix),
            "--show-transformed"
        ] + archive_file_list
    )

    # save wavefunctions (smwf files)
    if task.get("save_wavefunctions"):
        smwf_archive_file_list = glob.glob(work_dir+"/mfdn_smwf*")
        smwf_archive_file_list += glob.glob(work_dir+"/mfdn_MBgroups*")
        smwf_archive_filename = "{:s}-wf.tar".format(filename_prefix)
        mcscript.call(
            [
                "tar", "zcvf", smwf_archive_filename,
                "--transform=s,{:s}/,,".format(work_dir),
                "--transform=s,^,{:s}/{:s}{:s}/,".format(mcscript.parameters.run.name, descriptor, postfix),
                "--show-transformed"
            ] + smwf_archive_file_list
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
        if task.get("save_wavefunctions"):
            wavefunction_dir = os.path.join(mcscript.parameters.run.work_dir, "wavefunctions")
            mcscript.call(["mkdir", "-p", wavefunction_dir])
            mcscript.call(
                [
                    "mv",
                    "--verbose",
                    smwf_archive_filename,
                    "--target-directory={}".format(wavefunction_dir)
                ]
            )


def cleanup_mfdn_workdir(task, postfix=""):
    """Remove temporary MFDn work files.

    Arguments:
        task (dict): as described in module docstring
        postfix (string, optional): identifier to add to generated files
    """
    # cleanup of wave function files
    scratch_file_list = glob.glob("work{:s}/*".format(postfix))
    mcscript.call(["rm", "-vf"] + scratch_file_list)
