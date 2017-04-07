import mcscript

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
        ##ncsm_config.data_dir_results_archive,
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
        ["tar", "-zcvf", archive_filename ] + result_files,
        cwd=mcscript.task.results_dir
    )

    # copy archive out to home results archive directory
    mcscript.call(
        ["cp","-v",archive_filename,"-t",ncsm_config.data_dir_results_archive],
        cwd=mcscript.task.results_dir
    )
