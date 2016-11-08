""" ncsm.py -- NCSM utility function for batch run control

  Expected dictionary keys for h2mixer + mfdn-h2 run:

        # nuclide parameters
        "nuclide" (tuple of int) : (Z,N) tuple
        "Nv" (int) : N of valence shell (for use in truncation)
        
        # Hamiltonian parameters
        "interaction" (string) : name for interaction (for use in descriptor and filenames)
        "coulomb" (bool) : logical value whether or not to include Coulomb
        "aLawson" (float) : aLawson
        "hwLawson" (float) : hwLawson (typically hw of basis)

        # basis parameters
        "hw_int" (float) : hw of source interaction
        "hw_int_coul" (float) : hw for Coulomb interaction
        "hw" (float) : hw of basis
        "basis" (tuple) : logical basis (radial_basis_p,1.,radial_basis_n,beta_n), to be scaled by hw
        "scaled_basis" (tuple) : deduced basis used in task descriptor, etc.
        "xform_cutoff" (tuple) : transform cutoff, as tuple ("ob"|"tb",N) 
        "xform_cutoff_coul" (tuple) : transform cutoff for Coulomb interaction, as tuple ("ob"|"tb",N) 
        "Nmax" (int) : Nmax
        "Nstep" (int) : Nstep (2 for single parity, 1 for mixed parity)
        "traditional_ho" (bool) : whether or not to assume traditional HO-basis run
            (a value of True switches to an abbreviated task descriptor and enable transition calculation in MFDn)  

        # diagonalization parameters
        "Mj" (float) : M-scheme Mj value
        "eigenvectors" (int) : number of eigenvectors to calculate
        "initial_vector" (int) : initial vector code for mfdn
        "lanczos" (int) : lanczos iterations
        "tolerance" (float) : diagonalization tolerance parameter

        # obdme parameters
        "obdme_multipolarity" (int) : maximum multipolarity for calculation of densities
        "obdme_reference_state_list" (list) : list of reference states (J,g,i) for density calculation, or "all2all"
        "keep_obdme" (bool) : optional, whether or not to save obdme files in archive

        # two-body observable parameters
        "obs-R20K20" (bool, optional) : whether or not to include center-of-mass diagnostic observables (default: False)
        "obs-am-sqr" (bool, optional) : whether or not to include squared angular momentum (default: False)
        "obs-H-components" (bool, optional) : whether or not to include Hamiltonian terms (default: False)
        
        # version parameters
        "mfdn_executable" (string) : mfdn executable name
        "mfdn_wrapper" (string) : wrapper code for given MFDn version

        # emcalc postprocessing parameters
        "em_multipolarity_list" (list, optional) : list of multipolarities ("E"|"M",lambda) 

  Expected dictionary keys for tbme to h2 conversion run:

        # interaction result naming
        "interaction" (string) : name for interaction (for use in descriptor and filenames)
        "truncation" (tuple) : truncation (for use in descriptor and filenames)
        "hw" (float) : hw of basis

        # interaction file parameters
        "int_dir" (string) : directory containing compressed interaction file
        "int_basename" (string) : base file name for interaction file


  Created by M. A. Caprio, University of Notre Dame.
  2/12/13 (mac): Originated.
  1/22/14 (mac): Python 3 update.
  8/3/14 (mac): Condense with mfdn_task and overhaul.  Upgrade for v14beta06.
  8/4/14 (mac): Add tbme-format conversion.
  12/9/14 (mac): Update h2gen and h2xform handlers.
  5/12/15 (mac): Update treatment of Coulomb in h2mixer runs (make Coulumb cutoff x
    form_cutoff_coul independent of from nuclear interaction cutoff, eliminate input of Coulomb file 
    if Coulomb off).
  5/13/15 (mac): Insert "future" statements for Python 2 legacy support.
  5/23/15 (mac): Add task parameter "keep_obdme".
  5/25/15 (mac): Add squared angular momentum two-body observables support (task parameter "obs-am-sqr").
  6/2/15 (mac): Update filename for squared angular momentum two-body observables.
  6/9/15 (mac): Sanity check on reference state parity.
  6/25/15 (mac): Add documentation. Replace all "twice angular momentum" parameters with true values.
  6/30/15 (mac): Save mfdn.out with results.
  7/20/15 (mac): Allow control over choice of partition file.
  8/4/15 (mac): Revise control over choice of partition file.
  9/22/16 (mac): Add Hamiltonian components as two-body observables.

"""

from __future__ import print_function, division
  
import os
import glob
import shutil
import datetime

import mcscript
import ncsm_config

################################################################
# interaction truncation utilities
################################################################

def interaction_truncation_for_Nmax(Nv,Nmax,standardize=True):
    """ interaction_truncation_for_NMax (Nv, Nmax) -> (N1b,N2b)

    Identifies (N1b,N2b) truncation needed for interaction to support given many-body run.

    Nv -- N for valence shell
    Nmax -- many-body Nmax
    standardize -- if True, returns a standardized (N1b,N2b) truncation (see below)

    The basic formulas are:
    
    N1b = Nv + Nmax
    N2b = 2*Nv + Nmax

    However, this function returns a standardized (N1b,N2b) truncation, from
    a more restricted set, to reduce the proliferation of stored interaction files:
    -- round odd (unnatural parity) Nmax to next higher even (natural parity) Nmax
    -- could also impose a sensible floor, such as Nmax = 6, to avoid unnecessary extra small interaction files, but perhaps this gets too complicated

    For instance, for p-shell nuclei (Nv=1)...

        Nmax  N1b-N2b
        0     1-2
        1     3-4
        2     3-4

    However, for the actual files used as input for MFDn, MFDn expects
    (as of version14_beta02) the exact minimum truncation needed, so
    standardization should be disabled for this truncation.

    """

    # standardize Nmax
    Nmax_base = 0  # value 0, so this restriction not currently active
    if (standardize):
        Nmax_used = max(Nmax,Nmax_base)
        Nmax_used = Nmax_used + (Nmax_used % 2)
    else:
        Nmax_used = Nmax

    # return truncation
    N1b = Nv + Nmax_used
    N2b = 2*Nv + Nmax_used
    return (N1b, N2b)
    

def regularize_truncation(truncation):
    """ regularize_truncation( (code, cutoff) ) -> (N1b,N2b)

    Converts a mnemonic description ("ob",N1b) or ("tb",N2b) to the
    appropriate (N1b,N2b) pair for a two-body interaction file.
    """

    (code, N) = truncation
    if (code == "ob"):
        return (N,2*N)
    elif (code == "tb"):
        return (N,N)

################################################################
# general basis naming conventions
################################################################

def basis_string_pn(basis):
    """ basis_string_pn(basis) -> name

    Takes basis specified as tuple (basis_name_p, beta_p,
    basis_name_n, beta_n) and converts to string for use in filename.
    """
    
    return "%s-%.3f-%s-%.3f" % basis

################################################################
# basis length parameter calculation
################################################################

def beta_for_hw(hw_int,hw):
    """ beta_for_hw(hw_int,hw) -> beta returns the length scale factor
    to scale from hw_int to hw
    """
    
    return (hw/hw_int)**(-0.5)

################################################################
# configuration
################################################################

def set_run_search_lists(radial=[],interaction=[],r2k2=[],xform=[]):
    """
    set_run_search_lists(radial=[],r2k2=[],xform=[]) sets lists of runs to search for data files
    """
    global radial_run_list, interaction_run_list, r2k2_run_list, xform_run_list

    radial_run_list = radial
    interaction_run_list = interaction
    r2k2_run_list = r2k2
    xform_run_list = xform

################################################################
# task descriptor for h2mixer + mfdn-h2 run
################################################################

def task_descriptor_format_2(current_task):
    """ format before scan by dilation """
    return "Z%d-N%d-%s-c%d-%s-hw%g-aL%g-hwL%g-Nmax%d" % (
        current_task["nuclide"][0], current_task["nuclide"][1],
        current_task["interaction"],
        current_task["coulomb"],
        mcscript.dashify((basis_string_pn(current_task["basis"]), mcscript.dashify(current_task["xform_cutoff"]))),
        current_task["hw"],
        current_task["aLawson"],
        current_task["hwLawson"],
        current_task["Nmax"]
        )

def task_descriptor_format_3(current_task):
    """ run0190 -- introducing scan by dilation from hw_int """
    return "Z%d-N%d-%s-c%d-%s-hw%.3f-aL%g-hwL%.3f-Nmax%d-%g-%s" % (
        current_task["nuclide"][0], current_task["nuclide"][1],
        current_task["interaction"],
        current_task["coulomb"],
        basis_string_pn(current_task["basis"]),
        current_task["hw"],
        current_task["aLawson"],
        current_task["hwLawson"],
        current_task["Nmax"],
        current_task["hw_int"], mcscript.dashify(current_task["xform_cutoff"])
        )

def task_descriptor_format_4(current_task):
    """ run0223 -- add twice M, Nstep flag, and Lanczos iterations to task descriptor """
    return "Z%d-N%d-%s-c%d-%s-hw%.3f-aL%g-hwL%.3f-Nmax%d%s-mmj%d-hwint%g-%s-lan%d" % (
        current_task["nuclide"][0], current_task["nuclide"][1],
        current_task["interaction"],
        current_task["coulomb"],
        basis_string_pn(current_task["basis"]),
        current_task["hw"],
        current_task["aLawson"],
        current_task["hwLawson"],
        current_task["Nmax"],
        mcscript.ifelse(current_task["Nstep"]==1,"x",""),
        current_task["2mj"],
        current_task["hw_int"],
        mcscript.dashify(current_task["xform_cutoff"]),
        current_task["lanczos"]
        )

def task_descriptor_format_5(current_task):
    """ run0301  -- rename twice M from mmj to MM, always pad Nmax to 2 digits

    Upgrade code to use new-style string format method.
    
    Abbreviate name for traditional ho run.

    Incorporate vc's FCI flag.
    """

    if(current_task["traditional_ho"]):
        template_string = (
            "Z{nuclide[0]}-N{nuclide[1]}-{interaction}-{coulomb:d}"
            "-hw{hw:.3f}"
            "-aL{aLawson:g}"
            "-Nmax{Nmax:02d}{mixed_parity_flag}{fci_flag}-MM{twice_mj:d}"
            "-lan{lanczos:d}"
            )
    else:
        template_string =(
            "Z{nuclide[0]}-N{nuclide[1]}-{interaction}-{coulomb:d}-{basis_string_pn}-hw{hw:.3f}"
            "-aL{aLawson:g}-hwL{hwLawson:.3f}"
            "-Nmax{Nmax:02d}{mixed_parity_flag}{fci_flag}-MM{twice_mj:d}"
            "-hwint{hw_int:g}-{xform_cutoff_name:s}-lan{lanczos:d}"
            )

    descriptor = template_string.format(
        basis_string_pn=basis_string_pn(current_task["basis"]),
        mixed_parity_flag=mcscript.ifelse(current_task["Nstep"]==1,"x",""),
        fci_flag = mcscript.ifelse(current_task.get("fci",False),"-fci",""),
        xform_cutoff_name=mcscript.dashify(current_task["xform_cutoff"]),
        twice_mj = int(2*current_task["Mj"]),
        **current_task
        )

    return descriptor

def task_descriptor_format_6(current_task):
    """ run0373  -- replace MM (integer) with Mj (float)
    """

    if(current_task["traditional_ho"]):
        template_string = (
            "Z{nuclide[0]}-N{nuclide[1]}-{interaction}-{coulomb:d}"
            "-hw{hw:.3f}"
            "-aL{aLawson:g}"
            "-Nmax{Nmax:02d}{mixed_parity_flag}{fci_flag}-Mj{Mj:03.1f}"
            "-lan{lanczos:d}"
            )
    else:
        template_string =(
            "Z{nuclide[0]}-N{nuclide[1]}-{interaction}-{coulomb:d}"
            "-{basis_string_pn}"  # general basis
            "-hw{hw:.3f}"
            "-aL{aLawson:g}"
            "-hwL{hwLawson:.3f}"  # general basis
            "-Nmax{Nmax:02d}{mixed_parity_flag}{fci_flag}-Mj{Mj:03.1f}"
            "-hwint{hw_int:g}-{xform_cutoff_name:s}" # general basis
            "-lan{lanczos:d}"
            )

    descriptor = template_string.format(
        basis_string_pn=basis_string_pn(current_task["basis"]),
        mixed_parity_flag=mcscript.ifelse(current_task["Nstep"]==1,"x",""),
        fci_flag = mcscript.ifelse(current_task.get("fci",False),"-fci",""),
        xform_cutoff_name=mcscript.dashify(current_task["xform_cutoff"]),
        ## Mj=current_task["2mj"]/2,
        **current_task
        )

    return descriptor

################################################################
# observable definitions
################################################################

angular_momentum_operator_list = [ 
##    "Lp", "Ln", "L", 
##    "Sp", "Sn", "S", 
##    "Jp", "Jn", "J"
    "L", "Sp", "Sn", "S", "J"
]


################################################################
# h2 production -- h2gen
################################################################

def task_pool_h2gen(current_task):
    """
    task_descriptor_h2gen(...) create task pool.
    """

    return "h2gen"

def task_descriptor_h2gen(current_task):
    """
    task_descriptor_h2gen(...) create task descriptor from basis description.
    """

    return basis_string_pn(current_task["basis"])


def task_handler_h2gen(current_task):
    """
    task_handler_h2gen(...) runs h2gen for the current task's basis.

    All Nmax values are iterated over in one task, for simplicity,
    since these r2k2 tasks are already generally short compared to the
    xform tasks.
    """

    # recover current basis
    basis = current_task["basis"]
    Nv = current_task["Nv"]
    Nmax_list = current_task["Nmax_list"]

    # iterate over Nmax
    for Nmax in Nmax_list:

        print ("h2gen:", basis_string_pn(basis), Nmax, "->", interaction_truncation_for_Nmax(Nv,Nmax))

        # set radial bases
        (basis_name_p, beta_p, basis_name_n, beta_n) = basis

 
        radial_basename_p_unresolved = "radial-me-%s" % (basis_name_p, )
        radial_basename_p = mcscript.subpath_search_basename (ncsm_config.data_dir_radial, radial_run_list, radial_basename_p_unresolved)
        if (radial_basename_p is None):
            print ("Cannot find radial file", radial_basename_p_unresolved, "in", ncsm_config.data_dir_radial, radial_run_list)
            sys.exit(1)
        radial_basename_n_unresolved = "radial-me-%s" % (basis_name_n, )
        radial_basename_n = mcscript.subpath_search_basename (ncsm_config.data_dir_radial, radial_run_list, radial_basename_n_unresolved)
        if (radial_basename_n is None):
            print ("Cannot find radial file", radial_basename_n_unresolved, "in", ncsm_config.data_dir_radial, radial_run_list)
            sys.exit(1)

        # set truncation
        (N1b, N2b) = interaction_truncation_for_Nmax(Nv,Nmax)

        # set output filname
        output_basename = mcscript.dashify(("r2k2",basis_string_pn(basis), N1b, N2b))

        # run h2gen
        h2gen_params = {
            "N1b" : N1b, "N2b" : N2b,
            "radial_basename_p" : radial_basename_p, "radial_basename_n" : radial_basename_n,
            "beta_p" : beta_p, "beta_n" : beta_n,
            "output_basename" : output_basename
            }
        call_h2gen (h2gen_params)

    # save results locally
    h2_filename_list = glob.glob("*.bin")
    mcscript.call([ "mv" , "--target-directory="+mcscript.task.results_dir ] + h2_filename_list)

    # copy files for immediate use
    #   normally would be done in archive phase 0, which is still available as fallback
    #   therefore, this is set as a "try but don't worry if you fail" action and will not result
    #   in the task being recorded as a failure -- archive phase 0 can/should be run as insurance
    ## target_dir = os.path.join(ncsm_config.data_dir_h2,mcscript.run.name)
    ## mcscript.call(["mkdir", "-p", target_dir],check_return=False)
    ## mcscript.call([ "cp", "-v", "--target-directory="+target_dir ] + h2_filename_list,cwd=mcscript.task.results_dir,check_return=False)

################################################################
# h2 production -- xform
################################################################

def task_pool_h2xform(current_task):
    """
    task_descriptor_h2xform(...) create task pool.
    """

    return "h2xform"

def task_pool_h2xform_Nmax(current_task):
    """
    task_descriptor_h2xform(...) create task pool showing last Nmax value.
    """

    return "Nmax{Nmax:02d}".format(Nmax=current_task["Nmax_list"][-1])


def task_descriptor_h2xform(current_task):
    """
    task_descriptor_h2xform(...) create task descriptor from basis description and scaling.
    """
    return "%s-%s-%g" % (
        current_task["interaction"],
        basis_string_pn(current_task["scaled_basis"]),
        current_task["hw_int"]
    )

def task_handler_h2xform(current_task):
    """
    task_handler_h2xform(...) runs h2xform for the current task's basis and interaction.
    """
    
    print ("h2xform", current_task["descriptor"])

    # recall parameters
    interaction = current_task["interaction"] 
    interaction_truncation = current_task["interaction_truncation"] 
    scaled_basis = current_task["scaled_basis"] 
    hw_int = current_task["hw_int"] 
    Nv = current_task["Nv"]
    xform_cutoff_list = current_task["xform_cutoff_list"]
    Nmax_list = current_task["Nmax_list"]

    # set radial bases
    (basis_name_p, beta_p, basis_name_n, beta_n) = scaled_basis
    radial_basename_p = mcscript.subpath_search_basename (ncsm_config.data_dir_radial, radial_run_list, "radial-xform-%s-%.3f" % (basis_name_p, beta_p))
    radial_basename_n = mcscript.subpath_search_basename (ncsm_config.data_dir_radial, radial_run_list, "radial-xform-%s-%.3f" % (basis_name_n, beta_n))

    # set input TBME file
    input_basename_unresolved = "%s-%s-%g" % (interaction, mcscript.dashify(interaction_truncation), hw_int)
    input_basename = mcscript.subpath_search_basename (ncsm_config.data_dir_h2, interaction_run_list, input_basename_unresolved)
    print ("Interaction TBME file:", input_basename_unresolved)
    if (input_basename is None):
        print ("Missing interaction...")
        exit(1)

    output_streams = []
    for cutoff in xform_cutoff_list:
        for Nmax in Nmax_list:

            # set output truncation
            output_truncation = interaction_truncation_for_Nmax(Nv,Nmax)

            # set output basename
            output_basename = "%s-%s-%s-%s-%g" % (interaction, basis_string_pn(scaled_basis), mcscript.dashify(cutoff), mcscript.dashify(output_truncation), hw_int)

            # record output parameters
            output_stream_params = {
                "basename" : output_basename,
                "radial_basename_p" : radial_basename_p, "radial_basename_n" : radial_basename_n,
                "cutoff" : cutoff,
                "truncation" : output_truncation
                }

            # append to list of output streams
            output_streams.append(output_stream_params)
            
    # run h2xform
    h2xform_params = { "input_basename" : input_basename, "output_streams" : output_streams }
    call_h2xform (h2xform_params)

    # save results locally
    h2_filename_list = glob.glob("*.bin")
    mcscript.call([ "mv" , "--target-directory="+mcscript.task.results_dir ] + h2_filename_list)

    # copy files for immediate use -- DISABLED
    #   normally would be done in archive phase 0, which is still available as fallback
    #   therefore, this is set as a "try but don't worry if you fail" action and will not result
    #   in the task being recorded as a failure -- archive phase 0 can/should be run as insurance
    # CAVEAT: this results in "double disk usage" if run directory and permanent h2 directory are on same filesystem
    if False:
        target_dir = os.path.join(ncsm_config.data_dir_h2,mcscript.run.name)
        mcscript.call(["mkdir", "-p", target_dir],check_return=False)
        mcscript.call(["cp", "-v", "--target-directory="+target_dir ] + h2_filename_list,cwd=mcscript.task.results_dir,check_return=False)

################################################################
# low-level h2 code invocation
################################################################

def call_h2gen (params):
    """ call_h2gen(params)

    Invokes h2gen, with input determined by the dictionary params, which should have the following key-value pairs:
    
    {
    "N1b" : N1b, "N2b" : N2b,
    "radial_basename_p" : radial_basename_p, "radial_basename_n" : radial_basename_n,
    "beta_p" : beta_p, "beta_n" : beta_n,
    "output_basename" : output_basename
    }
    """

    # generate input lines
    input_lines = [
        "%d %d" % ( params["N1b"], params["N2b"]),
        "%s %s" % ( params["radial_basename_p"], params["radial_basename_n"]),
        "%.3f %.3f" % ( params["beta_p"], params["beta_n"]),
        "%s" % ( params["output_basename"], )
        ]
    
    # run h2gen
    mcscript.call_serial([ ncsm_config.h2gen ], mcscript.run, input_lines=input_lines)

def call_h2xform (params):
    """ call_h2xform(params)

    Invokes h2xform, with input determined by the dictionary params,
    which should have the following key-value pairs:
    
    { "input_basename" : input_basename, "output_streams" : output_streams }

    In turn, output_streams should be a list of dictionaries, each of
    which should have the following key-value pairs:

    {
    "basename" : output_basename,
    "radial_basename_p" : radial_basename_p, "radial_basename_n" : radial_basename_n,
    "cutoff" : cutoff,
    "truncation" : output_truncation
    }

    The cutoff should be specified in the mnemonic form ("ob",N1b) or ("tb",N2b).

    """

    # generate input lines
    input_lines = [
        "%s" % ( params["input_basename"])
        ]

    for output_stream_params in params["output_streams"]:
        (N1b_cut, N2b_cut) = regularize_truncation(output_stream_params["cutoff"]);
        (N1b, N2b) = output_stream_params["truncation"];
        input_lines.append( 
            "%s %s %s %d %d %d %d" % (
                output_stream_params["basename"],
                output_stream_params["radial_basename_p"], output_stream_params["radial_basename_n"],
                N1b_cut, N2b_cut,
                N1b, N2b
                )
            )
    
    # run h2gen
    mcscript.call_serial([ ncsm_config.h2xform ], mcscript.run, input_lines=input_lines)

def call_h2mixer (params):
    """ call_h2mixer(params)

    Invokes h2xmixer, with input determined by the dictionary params,
    which should have the following keys:
    
    {
    "input_truncation" : (N1b_in,N2b_in),
    "output_truncation" : (N1b_out,N2b_out),
    "Z": Z, "N" : N, "hw" : hw,
    "r2k2_basename" : r2k2_basename,
    "input_basenames" : input_basenames
    "output_streams" : output_streams 
    }

    In turn, output_streams should be a list of dictionaries, each of
    which should have the following key-value pairs:

    {
    "basename" : output_basename,
    "components" : components
    }

    In turn, components should be a list of tuples, of arguments to "add",
    to be converted using str().

    """

    # generate input lines -- common
    input_lines = [
        "%d %d" % params["input_truncation"], 
        "%d %d" % params["output_truncation"], 
        "scale %d %d %g" % ( params["Z"], params["N"], params["hw"]),
        "in r2k2 %s" % ( params["r2k2_basename"], )
        ]

    # generate input lines -- additional input streams
    input_lines += [
        "in %s" % ( basename, ) 
        for basename in params["input_basenames"]
        ]

    # generate input lines -- output streams
    for output_stream in params["output_streams"]:
        input_lines += [
            "out %s" % ( output_stream["basename"], )
            ]
        input_lines += [
            "add %s" % ( mcscript.spacify(component), )
            for component in output_stream["components"]
        ]

    # run h2mixer
    mcscript.call_serial([ ncsm_config.h2mixer ], mcscript.run,
                         input_lines=input_lines)

################################################################
# h2 archiving
################################################################

def archive_handler_h2_cp ():

    # write current toc
    toc_filename = mcscript.task.write_current_toc()

    # copy files
    # may be absorbed on file-by-file basis into h2 xform handler
    # but this provides a fallback if original copy attempt is unsuccessful
    target_dir = os.path.join(ncsm_config.data_dir_h2,mcscript.run)
    mcscript.call(["mkdir", "-p", target_dir])
    mcscript.call(["cp", "-rv", mcscript.task.results_dir,target_dir,"--no-target-directory"])
    
def archive_handler_h2_archive ():

    # write current toc
    toc_filename = mcscript.task.write_current_toc()
    
    # make archive -- whole dir
    # note possible error return of tar if flags are modified, e.g., another archive phase running in parallel
    work_dir_parent = os.path.join(mcscript.task.task_root_dir,"..")
    archive_filename = os.path.join(
        ncsm_config.data_dir_h2_archive,
        "%s-archive-h2-%s.tgz" % (mcscript.run.name, mcscript.date_tag())
        )
    mcscript.call(
        ["tar", "zcvf", archive_filename,
         os.path.join(mcscript.run.name,toc_filename),
         os.path.join(mcscript.run.name,"flags"),
         os.path.join(mcscript.run.name,"output"),
         os.path.join(mcscript.run.name,"launch"),
         os.path.join(mcscript.run.name,"results")
         ],
        cwd=work_dir_parent,check_return=True
        )
    mcscript.call(["chown", ":m94", archive_filename])
    mcscript.call(["chmod", "g+r", archive_filename])

    # put to hsi
    hsi_subdir = "2014"
    hsi_arg = "lcd %s; cd %s; put %s" % (os.path.dirname(archive_filename), hsi_subdir, os.path.basename(archive_filename))
    mcscript.call(["hsi",hsi_arg])

##################################################################
# h2mixer 
##################################################################

def task_handler_h2mixer(current_task):
    """ Invoke mfdn in h2 run, given current task parameters.

    Expected dictionary keys: see initial module docstring
    """
    print ("h2mixer", current_task["descriptor"])

    # recall parameters
    nuclide = current_task["nuclide"] 
    interaction = current_task["interaction"] 
    coulomb = current_task["coulomb"] 
    hw_int = current_task["hw_int"] 
    aLawson = current_task["aLawson"] 
    hwLawson = current_task["hwLawson"]
    hw = current_task["hw"] 
    basis = current_task["basis"] 
    scaled_basis = current_task["scaled_basis"] 
    xform_cutoff = current_task["xform_cutoff"] 
    Nmax = current_task["Nmax"] 
    Nv = current_task["Nv"]

    # basic setup
    (Z, N) = nuclide
    input_truncation = interaction_truncation_for_Nmax(Nv,Nmax)
    r2k2_basename_unresolved = mcscript.dashify(("r2k2", basis_string_pn(basis), mcscript.dashify(input_truncation)))
    r2k2_basename = mcscript.subpath_search_basename (ncsm_config.data_dir_h2, r2k2_run_list, r2k2_basename_unresolved)
    output_truncation = interaction_truncation_for_Nmax(Nv,Nmax,standardize=False)
    h2mixer_params = {
        "input_truncation" : input_truncation,
        "output_truncation" : output_truncation,
        "Z": Z, "N" : N, "hw" : hw,
        "r2k2_basename" : r2k2_basename,
        "input_basenames" : [],
        "output_streams" : []
        }

    # input stream setup
    # input 5 -- built interaction matrix elements -- for scaled basis
    tbme_basename_VNN_unresolved = mcscript.dashify((interaction, basis_string_pn(scaled_basis), mcscript.dashify(xform_cutoff), mcscript.dashify(input_truncation), format(hw_int,"g")))
    h2mixer_params["input_basenames"].append(
        mcscript.subpath_search_basename (ncsm_config.data_dir_h2, xform_run_list, tbme_basename_VNN_unresolved)
        )
    # input 6 (OPTIONAL) -- built Coulomb matrix elements -- for unscaled basis
    if (coulomb):
        hw_int_coul = current_task["hw_int_coul"] 
        xform_cutoff_coul = current_task["xform_cutoff_coul"] 
        tbme_basename_VC_unresolved = mcscript.dashify(("VC", basis_string_pn(basis), mcscript.dashify(xform_cutoff_coul), mcscript.dashify(input_truncation), format(hw_int_coul,"g")))
        h2mixer_params["input_basenames"].append(
            mcscript.subpath_search_basename (ncsm_config.data_dir_h2, xform_run_list, tbme_basename_VC_unresolved)
        )
        
    # output stream -- Hamiltonian
    components = [
            ("Trel",), # Trel
            ("in", 5, 1.0), # VNN
            ("NCM", hwLawson, aLawson) # Lawson
            ]
    if (coulomb):
        components += [
            ("scaled", 6, hw_int_coul, -1, 1.0) # VC  -- add scaled <in#> <hw_ref> <degree> <multiplier> 
        ]
    h2mixer_params["output_streams"].append(
        {
            "basename" : "tbme-h2",
            "components" : components
        }
    )

    # output stream -- rrel2 (required observable)
    h2mixer_params["output_streams"].append(
        {
        "basename" : "tbme-rrel2",
        "components" : [ ("rrel2",) ]
        }
        )

    # output stream -- NCM(hwLawson)
    h2mixer_params["output_streams"].append(
        {
        "basename" : "tbme-NCM",
        "components" : [ ("NCM", hwLawson, 1.0) ]
        }
        )

    # output stream -- CM observables R20 & K20
    if (current_task.get("obs-R20K20",False)):
        h2mixer_params["output_streams"].append(
            {
                "basename" : "tbme-R20",
                "components" : [ ("R20",) ]
            }
        )
        h2mixer_params["output_streams"].append(
            {
                "basename" : "tbme-K20",
                "components" : [ ("K20",) ]
            }
        )

    # output stream -- squared angular momentum observables
    if (current_task.get("obs-am-sqr",False)):
        for op in angular_momentum_operator_list:
            h2mixer_params["output_streams"].append(
                {
                    "basename" : "tbme-{}2".format(op),
                    "components" : [ ("am-sqr", op) ]
                }
            )

    # output streams -- tbme for H components (debugging)
    if (current_task.get("obs-H-components",False)):
        h2mixer_params["output_streams"].append(
            {
                "basename" : "tbme-Trel",
                "components" : [ ("Trel",) ]
            }
        )
        h2mixer_params["output_streams"].append(
            {
                "basename" : "tbme-VNN",
                "components" : [ ("in", 5, 1.0) ]
            }
        )
        h2mixer_params["output_streams"].append(
            {
                "basename" : "tbme-coul",
                "components" : [ ("scaled", 6, hw_int_coul, -1, 1.0) ]
            }
        )


    # run h2mixer
    call_h2mixer (h2mixer_params)

##################################################################
# mfdn h2 run
##################################################################

def task_handler_mfdn_h2(current_task):
    """ Invoke mfdn in h2 run, given current task parameters.
    
    Expected dictionary keys: see initial module docstring
    """

    print ("MFDn", current_task["descriptor"])

    # optional dictionary keys
    # TODO: neaten up with setdefault method
    if ("obs-R20K20" not in current_task):
        current_task["obs-R20K20"] = False
    if ("obs-am-sqr" not in current_task):
        current_task["obs-am-sqr"] = False
    if ("obs-H-components" not in current_task):
        current_task["obs-H-components"] = False
    if ("em_multipolarity_list" not in current_task):
        current_task["em_multipolarity_list"] = []
    if ("keep_obdme" not in current_task):
        current_task["keep_obdme"] = True

    # set up MFDn basis truncation parameters
    Nv = current_task["Nv"]
    Nmin = current_task["Nmax"] % current_task["Nstep"]
    truncation = interaction_truncation_for_Nmax(Nv,current_task["Nmax"],standardize=False)
    ## BUG: through 130520: Nshell = truncation[1] + 1
    ## BUG: through 150805: Nshell given to MFDn was one step higher than needed for odd Nmax (default standardize="True")
    Nshell = truncation[0] + 1
    ## print("truncation",truncation,"Nshell",Nshell)


    # set up h2 filenames
    h2_basename = "tbme-h2"
    obs_basename_list = [ "tbme-rrel2", "tbme-NCM" ]
    if (current_task["obs-R20K20"]):
        obs_basename_list += [ "tbme-R20", "tbme-K20" ]
    if (current_task["obs-am-sqr"]):
        for op in angular_momentum_operator_list:
            obs_basename_list.append("tbme-{}2".format(op))
    if (current_task["obs-H-components"]):
        obs_basename_list += [ "tbme-Trel", "tbme-VNN", "tbme-coul"]

    # guard against pathetically common mistakes
    reference_state_list = current_task["obdme_reference_state_list"]
    for (J,g,i) in reference_state_list:
        twice_J = int(2*J)
        if ((twice_J%2) != (sum(current_task["nuclide"])%2)):
            raise ValueError("invalid angular momentum for reference state")
        if ((g != (current_task["Nmax"]%2)) and (current_task["Nstep"] != 1)):
            raise ValueError("invalid parity for reference state")

       
    # import partitioning file
    partitioning_suffix = os.environ.get("MFDN_PARTITIONING",None)
    if (partitioning_suffix is not None):
        partitioning_filename = os.path.join(
            ncsm_config.data_dir_partitioning,
            "mfdn_partitioning.info_Nsh{:02d}_{}".format(Nshell,partitioning_suffix)
        )
        print ("Checking for partition file %s..." % partitioning_filename)
        if (os.path.exists(partitioning_filename)):
            mcscript.call(["cp", partitioning_filename, "mfdn_partitioning.info"])
        else:
            raise mcscript.task.ScriptError("partition file not found")

    ## partitioning_filename = os.path.join(
    ##     ncsm_config.data_dir_partitioning,
    ##     "mfdn_partitioning.info_Nsh{}".format(Nshell)
    ## )

    # invoke mfdn
    current_task.update(
        {
            "Nmin" : Nmin,
            "Nshell" : Nshell,
            "h2_basename" : h2_basename,
            "obs_basename_list" : obs_basename_list
            }
        )
    mfdn_wrapper = current_task["mfdn_wrapper"]
    mfdn_wrapper(current_task)

    # process results file
    result_filename = "%s-mfdn-%s.res" % (mcscript.run.name, current_task["descriptor"])
    out_filename = "%s-mfdn-%s.out" % (mcscript.run.name, current_task["descriptor"])
    # remove detailed occupation listing
    # TO NEATEN: replace with shell-free subprocess call and pipes
    ## print ("Filtering mfdn.res...") 
    ## mcscript.call(["sed '/More/,$d' < mfdn.res > " + result_filename],shell=True) 
    ## mcscript.call(["cat < mfdn.res > " + result_filename],shell=True) 
    ## mcscript.call(["cat < mfdn.out > " + out_filename],shell=True) 
    mcscript.call(["cp","mfdn.res",result_filename]) 
    mcscript.call(["cp","mfdn.out",out_filename]) 

    # archive obdme and other data files
    # including unfiltered mfdn.res
    archive_filename = "%s-mfdn-%s.tgz" % (mcscript.run.name, current_task["descriptor"])
    print ("Archiving output files...")
    obdme_file_list = ["mfdn.out", "mfdn.dat", "mfdn.res", "mfdn_spstates.info"]
    if (current_task["keep_obdme"]):
        obdme_file_list += glob.glob("*obdme*") 
    print ("obdme files:", obdme_file_list)
    mcscript.call(["tar", "zcvf", archive_filename] + obdme_file_list)

    # move results out
    ##mcscript.call(["mv run* -t "+mcscript.task.results_dir],shell=True)
    mcscript.call(["mv", result_filename, out_filename, archive_filename, "--target-directory="+mcscript.task.results_dir])

    # cleanup of wave function files
    scratch_file_list = glob.glob("mfdn_smwf*") + glob.glob("mfdn_MBgroups0*")
    mcscript.call(["rm", "-vf"] + scratch_file_list)


##################################################################
# mfdn h2 run -- ad hoc salvage of full res files
# for runs 0313 to 0321
##################################################################

def task_handler_mfdn_h2_salvage(current_task):
    """ Salvage res file, given current task parameters.
    
    Regenerates obdme tar file, including raw mfdn.res.
    Also regenerates results file, without doing sed pruning.
     
    Expected dictionary keys: see initial module docstring
    """

    print ("MFDn", current_task["descriptor"])

    # process results file
    result_filename = "%s-mfdn-%s.res" % (mcscript.run.name, current_task["descriptor"])
    # remove detailed occupation listing
    # TO NEATEN: replace with shell-free subprocess call and pipes
    ## print ("Filtering mfdn.res...") 
    ## mcscript.call(["sed '/More/,$d' < mfdn.res > " + result_filename],shell=True) 
    mcscript.call(["cat < mfdn.res > " + result_filename],shell=True) 

    # archive obdme and other data files
    # including unfiltered mfdn.res
    archive_filename = "%s-mfdn-%s.tgz" % (mcscript.run.name, current_task["descriptor"])
    print ("Archiving output files...")
    obdme_file_list = ["mfdn.out", "mfdn.dat", "mfdn.res", "mfdn_spstates.info"] + glob.glob("*obdme*") 
    print ("obdme files:", obdme_file_list)
    mcscript.call(["tar", "zcvf", archive_filename] + obdme_file_list)

    # move results out
    ##mcscript.call(["mv run* -t "+mcscript.task.results_dir],shell=True)
    mcscript.call(["mv", result_filename, archive_filename, "--target-directory="+mcscript.task.results_dir])


##################################################################
# emcalc
##################################################################

def task_handler_emcalc(current_task):
    """
    TODO: replace 2mj parameter with Mj
    """
    # Performance note on NERSC: When run as a normal job, this
    # scripting involves many aprun call, each of which carries
    # considerable overhead (at least ~1 sec on hopper), so total execution
    # time is vastly slower than in epar mode, which has initial node setup
    # overhead but then runs the script locally on the compute node, so no aprun in needed.
    # Compare run0252 task 0007 (serial) at 909 sec with task 0009 (epar) at 37 sec.

    print ("emcalc", current_task["descriptor"])

    # set radial bases
    basis = current_task["basis"] 
    (basis_name_p, beta_p, basis_name_n, beta_n) = basis
    radial_basename_p_unresolved = "radial-me-%s" % (basis_name_p, )
    radial_basename_p = mcscript.subpath_search_basename (ncsm_config.data_dir_radial, radial_run_list, radial_basename_p_unresolved)
    if (radial_basename_p is None):
        print ("Cannot find radial file", radial_basename_p_unresolved, "in", ncsm_config.data_dir_radial, radial_run_list)
        sys.exit(1)
    radial_basename_n_unresolved = "radial-me-%s" % (basis_name_n, )
    radial_basename_n = mcscript.subpath_search_basename (ncsm_config.data_dir_radial, radial_run_list, radial_basename_n_unresolved)
    if (radial_basename_n is None):
        print ("Cannot find radial file", radial_basename_n_unresolved, "in", ncsm_config.data_dir_radial, radial_run_list)
        sys.exit(1)

    # basic setup
    emcalc_params = {
        "em_multipolarity_list" : current_task["em_multipolarity_list"],
        "2mj" : current_task["2mj"],
        "hw" : current_task["hw"],
        "beta_p" : beta_p, "beta_n" : beta_n,
        "radial_basename_p" : radial_basename_p, "radial_basename_n" : radial_basename_n
        }

    # recover density files from results archive (in case original run files lost/purged deleted)
    print ("Recovering obdme files...")
    if (glob.glob("*obdme*") == []):
        archive_filename = "%s-mfdn-%s.tgz" % (mcscript.run.name, current_task["descriptor"])
        mcscript.call(["tar", "xvf", os.path.join(mcscript.task.results_dir,archive_filename)])
    else:
        print ("   Found obdme files already in task directory...")

    # invoke emcalc
    call_emcalc (emcalc_params,current_task["descriptor"])

    # move results out
    main_result_filename = "%s-emcalc-%s.dat" % (mcscript.run.name, current_task["descriptor"])
    mcscript.call(["mv", "-v", main_result_filename, "--target-directory=" + mcscript.task.results_dir])

def call_emcalc (params,result_descriptor):
    """ call_emcalc(params)

    Input file format:
        em_sigma em_lambda
        hw
        beta_p beta_n
        <obdme orbital info file>
        <obdme data file>
        <proton radial r2 file> <neutron radial r2 file>

    Valentino's example:
        E 2                                                                                         
        25                                                                                          
        1 1                                                                                         
        /home/valentino/mfdn.obdme.info                                                             
        /home/valentino/mfdn.statobdme.seq005.2J08.2T00                                                
        /home/valentino/radial-me-HO /home/valentino/radial-me-HO                                    

    Invokes emcalc, with input determined by the dictionary params,
    which should have the following key-value pairs:
    
    {
        "em_multipolarity_list" : em_multipolarity_list,
        "hw" : current_task["hw"],
        "beta_p" : beta_p, "beta_n" : beta_n,
        "radial_basename_p" : radial_basename_p, "radial_basename_n" : radial_basename_n
    }

    This function scans all obdme files.

    """

    # find obdme files
    obdme_file_list = glob.glob("mfdn.obdme.seq*") + glob.glob("mfdn.statobdme.seq*")
    print ("obdme files:", obdme_file_list)
    
    # result file
    main_result_filename = "%s-emcalc-%s.dat" % (mcscript.run.name, result_descriptor)
    main_result_file = open(main_result_filename,"w")

    # loop over transitions
    for obdme_file in obdme_file_list:
        for em_multipolarity in params["em_multipolarity_list"]: 
            # generate input lines
            input_lines = [
                "%s %d" % em_multipolarity,
                "%g" % ( params["hw"]),
                "%.3f %.3f" % (params["beta_p"], params["beta_n"]),
                "mfdn.obdme.info",
                "%s" % (obdme_file),
                "%s %s" % (params["radial_basename_p"], params["radial_basename_n"])
                ]

            # run emcalc
            mcscript.call_serial([ncsm_config.emcalc], mcscript.run, input_lines=input_lines)

            # harvest results to main results file
            temporary_result_filename = "emcalc.out"
            temporary_result_file = open(temporary_result_filename,"r")
            shutil.copyfileobj(temporary_result_file,main_result_file)
            temporary_result_file.close()

    # finalize main result file
    main_result_file.close()

################################################################
# mfdn archiving
################################################################

def archive_handler_mfdn_res_only ():

    # write current toc
    toc_filename = mcscript.task.write_current_toc()

    # make archive -- results
    archive_filename = os.path.join(
        ##ncsm_config.data_dir_results_archive,
        mcscript.task.archive_dir,
        "%s-results-%s.tgz" % (mcscript.run.name, mcscript.date_tag())
        )
    ## # store toc -- TODO once restructure subdirectories in tar file
    ## mcscript.call(
    ##     ["tar", "zcvf", archive_filename, toc_filename]
    ##     )
    os.chdir (mcscript.task.results_dir)
    result_files = glob.glob("*.res") + glob.glob("*.out") + glob.glob("*-emcalc-*.dat")
    mcscript.call(["tar", "-zcvf", archive_filename ] + result_files, cwd=mcscript.task.results_dir)

    # copy archive out to home results archive directory
    mcscript.call(["cp","-v",archive_filename,"-t",ncsm_config.data_dir_results_archive], cwd=mcscript.task.results_dir)

def archive_handler_mfdn_archive ():

    # make archive -- whole dir
    archive_filename = mcscript.task.archive_handler_generic ()
   
    # put to hsi
    hsi_subdir = format(datetime.date.today().year,"04d")  # subdirectory named by year
    hsi_arg = "lcd %s; mkdir %s; cd %s; put %s" % (os.path.dirname(archive_filename), hsi_subdir, hsi_subdir, os.path.basename(archive_filename))
    mcscript.call(["hsi",hsi_arg])




##################################################################
# mfdn conversion
##################################################################

def task_descriptor_tbme_conv(current_task):
    descriptor = "{interaction:s}-{truncation_name:s}-{hw:g}".format(
        truncation_name=mcscript.dashify(current_task["truncation"]),
        **current_task
        )
    return descriptor

def task_handler_mfdn_tbme_conv(current_task):
    """ Invoke mfdn in tbme to h2 conversion, given current task parameters.

    With some inspiration from run0164 vxx conversion.
    
    Expected dictionary keys: see initial module docstring

    """

    print ("MFDn", current_task["descriptor"])

    # retrieve input file
    int_filename_compressed = os.path.join(
        current_task["int_dir"],
        "{int_basename}.gz".format(**current_task)
        )
    mcscript.call_serial(["getgz",int_filename_compressed],mcscript.run)
    mcscript.call(["ls","-Flah"])

    # augment parameters
    (N1bmax,N2bmax) = regularize_truncation(current_task["truncation"])
    current_task.update({"N1bmax" : N1bmax, "N2bmax" : N2bmax})
    
    # invoke mfdn
    mfdn_wrapper = current_task["mfdn_wrapper"]
    mfdn_wrapper(current_task)

    # move results out
    result_filename = "TBME_V2full.bin"
    saved_filename = os.path.join(mcscript.task.results_dir,"{descriptor}.bin".format(**current_task))
    mcscript.call(["mv", result_filename, saved_filename])

    # cleanup of h2 files
    scratch_file_list = glob.glob("*.bin")
    mcscript.call(["rm", "-vf"] + scratch_file_list)
