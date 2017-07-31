"""handlers.py -- task handlers for MFDn runs.

Patrick Fasano
University of Notre Dame

- 03/22/17 (pjf): Created, split from __init__.py.
- 04/07/17 (pjf): Update for mcscript namespace changes.
- 06/05/17 (pjf): Added basic handlers for oscillator and natural orbital runs.
- 06/07/17 (pjf): Clean up style.
- 06/22/17 (pjf): Update references to mcscript.exception.ScriptError.
- 07/31/17 (pjf): Move mfdn driver from handler argument to task dictionary.
"""
import os
import glob
import mcscript

from . import (
    config,
    radial,
    tbme,
    utils,
    mfdn_v14,
)


################################################################
# basic oscillator run
################################################################

def task_handler_oscillator(task, postfix=""):
    """Task handler for basic oscillator run.

    Arguments:
        task (dict): as described in module docstring
        postfix (string, optional): identifier to add to generated files
    """
    mfdn_driver = task.get("mfdn_driver")
    if mfdn_driver is None:
        mfdn_driver = mfdn_v14
    radial.set_up_orbitals_ho(task, postfix=postfix)
    radial.set_up_radial_analytic(task, postfix=postfix)
    tbme.generate_tbme(task, postfix=postfix)
    mfdn_driver.run_mfdn(task, postfix=postfix)
    mfdn_driver.save_mfdn_output(task, postfix=postfix)


################################################################
# basic natural orbital run
################################################################
def task_handler_natorb(task):
    """Task handler for basic oscillator+natural orbital run.

    Arguments:
        task (dict): as described in module docstring
    """
    mfdn_driver = task.get("mfdn_driver")
    if mfdn_driver is None:
        mfdn_driver = mfdn_v14

    # sanity checks
    if not task.get("natural_orbitals"):
        raise mcscript.exception.ScriptError("natural orbitals not enabled")

    natorb_base_state = task.get("natorb_base_state")
    if type(natorb_base_state) is not int:
        raise mcscript.exception.ScriptError("invalid natorb_base_state: {}".format(natorb_base_state))

    # first do base oscillator run
    task_handler_oscillator(task, postfix=utils.natural_orbital_indicator(0))

    # set correct basis mode
    task["basis_mode"] = config.BasisMode.kGeneric
    radial.set_up_orbitals_natorb(task=task, source_postfix=utils.natural_orbital_indicator(0), target_postfix=utils.natural_orbital_indicator(1))
    radial.set_up_radial_natorb(task=task, source_postfix=utils.natural_orbital_indicator(0), target_postfix=utils.natural_orbital_indicator(1))
    tbme.generate_tbme(task=task, postfix=utils.natural_orbital_indicator(1))
    mfdn_driver.run_mfdn(task=task, postfix=utils.natural_orbital_indicator(1))
    mfdn_driver.save_mfdn_output(task=task, postfix=utils.natural_orbital_indicator(1))


################################################################
# mfdn archiving
################################################################

def archive_handler_mfdn_res_only(task):
    """ Generate summary archive of MFDn results files.

    TODO: Update and test.
    """
    # write current toc
    toc_filename = mcscript.task.write_current_toc()

    # make archive -- results
    archive_filename = os.path.join(
        mcscript.task.archive_dir,
        "%s-results-%s.tgz" % (mcscript.parameters.run.name, mcscript.date_tag())
        )
    ## # store toc -- TODO once restructure subdirectories in tar file
    ## mcscript.call(
    ##     ["tar", "zcvf", archive_filename, toc_filename]
    ##     )
    os.chdir(mcscript.task.results_dir)
    result_files = glob.glob("*.res") + glob.glob("*.out") + glob.glob("*-emcalc-*.dat")
    mcscript.call(
        ["tar", "-zcvf", archive_filename] + result_files,
        cwd=mcscript.task.results_dir
    )

    # copy archive out to home results archive directory
    mcscript.call(
        ["cp", "-v", archive_filename, "-t", ncsm_config.data_dir_results_archive],
        cwd=mcscript.task.results_dir
    )
