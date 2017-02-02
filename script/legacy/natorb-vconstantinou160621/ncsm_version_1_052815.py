""" ncsm.py -- NCSM utility function for batch run control

  Expected dictionary keys for h2mixer + mfdn-h2 run:

        # nuclide parameters
        "nuclide" : (Z,N) tuple
        "Nv" : N of valence shell (for use in truncation)
        
        # Hamiltonian parameters
        "interaction" : name for interaction (for use in descriptor and filenames)
        "coulomb" : logical value whether or not to include Coulomb
        "aLawson" : aLawson
        "hwLawson" : hwLawson (typically hw of basis)

        # basis parameters
        "hw_int" : hw of source interaction
        "hw_int_coul" : hw for Coulomb interaction
        "hw" : hw of basis
        "basis" : logical basis (radial_basis_p,1.,radial_basis_n,beta_n), to be scaled by hw
        "scaled_basis" : deduced basis used in task descriptor, etc.
        "xform_cutoff" : transform cutoff, as tuple ("ob"|"tb",N) 
        "xform_cutoff_coul" : transform cutoff for Coulomb interaction, as tuple ("ob"|"tb",N) 
        "Nmax" : Nmax
        "Nstep" : Nstep (2 for single parity, 1 for mixed parity)
        "traditional_ho" : whether or not to assume traditional HO-basis run
            (for abbreviated task descriptor and for MFDn transition calculation)  

        # diagonalization parameters
        "2mj" : twice the M-scheme M value
        "eigenvectors" : number of eigenvectors to calculate
        "initial_vector" : initial vector code for mfdn
        "lanczos" : lanczos iterations
        "tolerance" : diagonalization tolerance parameter

        # obdme parameters
        "obdme_twice_multipolarity" : twice maximum multipolarity for calculation of densities
        "obdme_reference_state_list" : list of reference states (J,g,i) for density calculation, or "all2all"
        "keep_obdme" : optional, whether or not to save obdme files in archive

        # two-body observable parameters
        "obs-R20K20" : optional, whether or not to include center-of-mass diagnostic observables (default: False)
        "obs-am-sqr" : optional, whether or not to include squared angular momentum (default: False)
        
        # version parameters
        "mfdn_executable" : mfdn executable name
        "mfdn_wrapper" : wrapper code for given MFDn version

        # emcalc postprocessing parameters
        "em_multipolarity_list" : optional, list of multipolarities ("E"|"M",lambda)

        # Full configuration run
        "fci" : optional, whether or not to perform an fci run (default: False) 

  Expected dictionary keys for tbme to h2 conversion run:

        # interaction result naming
        "interaction" : name for interaction (for use in descriptor and filenames)
        "truncation" : truncation (for use in descriptor and filenames)
        "hw" : hw of basis

        # interaction file parameters
        "int_dir" : directory containing compressed interaction file
        "int_basename" : base file name for interaction file

        # version parameters
        "mfdn_executable" : mfdn executable name
        "mfdn_wrapper" : wrapper code for given MFDn version



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
  5/28/15 (vc): Added fci functionality

  Last modified 5/25/15 (mac).

"""

from __future__ import print_function, division
  
import os
import glob
import shutil
import subprocess

# Added for natural orbitals
import math
#import numpy
import mcscript
import ncsm_config

################################################################
# interaction truncation utilities
################################################################

# V. Constantinou (05/28/15)
# Indroduce the boolean variable fci which defaults to False. If one wants to perform 
# an fci truncation then the variable must be set to true..

def interaction_truncation_for_Nmax(Nv, Nmax, standardize=True, fci=False):
    """ interaction_truncation_for_Nmax (Nv, Nmax) -> (N1b,N2b)

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

    In the case of a full configuration interaction run the basic formulas are:

    N1b = Nv + Nmax
    N2b = 2*N1b

    """

    # standardize Nmax
    Nmax_base = 0  # value 0, so this restriction not currently active
    if (standardize):
        Nmax_used = max(Nmax,Nmax_base)
        Nmax_used = Nmax_used + (Nmax_used % 2)
    else:
        Nmax_used = Nmax

    # return truncation
    
    if(fci):
        N1b = Nv + Nmax_used
        N2b = 2*N1b
    else:
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

def set_run_search_lists(radial=[],interaction=[],r2k2=[],xform=[], mfdn=[], nogen=[]):
    """
    set_run_search_lists(radial=[],r2k2=[],xform=[]) sets lists of runs to search for data files
    """
    # Adding mfdn_run_list to be used with nogen (06/01/15) V. Constantinou
    # Adding nogen_run_list to be used with noradial (06/01/15) V. Constantinou
    global radial_run_list, interaction_run_list, r2k2_run_list, xform_run_list, mfdn_run_list, nogen_run_list

    radial_run_list = radial
    interaction_run_list = interaction
    r2k2_run_list = r2k2
    xform_run_list = xform
    mfdn_run_list = mfdn
    nogen_run_list = nogen

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

# (05/25/15) V Constantinou. Adding a small flag to indicate that we are running in full configuration
# For now I wn't rename Nmax in the string to something else..

def task_descriptor_format_5(current_task):
    """ run0301  -- rename twice M from mmj to MM, always pad Nmax to 2 digits

    Upgrade code to use new-style string format method.
    
    Abbreviate name for traditional ho run.
    """

    if(current_task["traditional_ho"]):
        template_string = (
            "Z{nuclide[0]}-N{nuclide[1]}-{interaction}-{coulomb:d}"
            "-hw{hw:.3f}"
            "-aL{aLawson:g}"
            "-Nmax{Nmax:02d}{mixed_parity_flag}-MM{2mj:d}"
            "-lan{lanczos:d}"
            )
    else:
        template_string =(
            "Z{nuclide[0]}-N{nuclide[1]}-{interaction}-{coulomb:d}-{basis_string_pn}-hw{hw:.3f}"
            "-aL{aLawson:g}-hwL{hwLawson:.3f}"
            "-Nmax{Nmax:02d}{mixed_parity_flag}-MM{2mj:d}"
            "-hwint{hw_int:g}-{xform_cutoff_name:s}-lan{lanczos:d}"
            )

    descriptor = template_string.format(
        basis_string_pn=basis_string_pn(current_task["basis"]),
        mixed_parity_flag=mcscript.ifelse(current_task["Nstep"]==1,"x",""),
        xform_cutoff_name=mcscript.dashify(current_task["xform_cutoff"]),
        **current_task
        )

    if(current_task.get("fci",False)):
        descriptor = descriptor + "-fci"

    return descriptor

################################################################
# observable definitions
################################################################

angular_momentum_operator_list = [ 
    "Lp", "Ln", "L", 
    "Sp", "Sn", "S", 
    "Jp", "Jn", "J"
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

        if(current_task.get("fci",False)):
            print ("h2gen:", basis_string_pn(basis), Nmax, "->", interaction_truncation_for_Nmax(Nv,Nmax,standardize=True,fci=current_task["fci"]))
        else:
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
        if(current_task.get("fci",False)):
            (N1b, N2b) = interaction_truncation_for_Nmax(Nv,Nmax,standardize=True,fci=current_task["fci"])
        else:
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
            if(current_task.get("fci",False)):
                output_truncation = interaction_truncation_for_Nmax(Nv,Nmax,standardize=True,fci=current_task["fci"])
            else:
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
    if(current_task.get("fci",False)):
        input_truncation = interaction_truncation_for_Nmax(Nv,Nmax,standardize=True,fci=current_task["fci"])
    else:
        input_truncation = interaction_truncation_for_Nmax(Nv,Nmax)

    r2k2_basename_unresolved = mcscript.dashify(("r2k2", basis_string_pn(basis), mcscript.dashify(input_truncation)))
    r2k2_basename = mcscript.subpath_search_basename (ncsm_config.data_dir_h2, r2k2_run_list, r2k2_basename_unresolved)

    print(r2k2_basename_unresolved)
    print(r2k2_basename)
    
    if(current_task.get("fci",False)):
        output_truncation = interaction_truncation_for_Nmax(Nv,Nmax,standardize=False,fci=current_task["fci"])
    else:
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
    

    # New approach -- not really working.. August 22nd 2015. Add a different approach for the VC interaction. Uncomment in the future if you want to try again
    """
    hw_int_coul = current_task["hw_int_coul"] 
    xform_cutoff_coul = current_task["xform_cutoff_coul"] 
    scaled_basis_coul = current_task["scaled_basis_coul"]

    tbme_basename_VC_unresolved = mcscript.dashify(("VC", basis_string_pn(scaled_basis_coul), mcscript.dashify(xform_cutoff_coul), mcscript.dashify(input_truncation), format(hw_int_coul,"g")))
    h2mixer_params["input_basenames"].append(
        mcscript.subpath_search_basename (ncsm_config.data_dir_h2, xform_run_list, tbme_basename_VC_unresolved)
    )
    """

    # output stream -- Hamiltonian
    components = [
            ("Trel",), # Trel
            ("in", 5, 1.0), # VNN
            ("NCM", hwLawson, aLawson) # Lawson
            ]

    # Uncomment after testing..
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
                    "basename" : "tbme-{}".format(op),
                    "components" : [ ("am-sqr", op) ]
                }
            )

    # output streams -- tbme for H components (debugging)
    ## h2mixer_params["output_streams"].append(
    ##     {
    ##     "basename" : "tbme-h2-v2",
    ##     "components" :  [
    ##         ("in", 5, 1.0), 
    ##         ("scaled", 6, 20.0, -1,coulomb_factor)
    ##         ]
    ##     }
    ##     )
    ## h2mixer_params["output_streams"].append(
    ##     {
    ##     "basename" : "tbme-h2-Trel",
    ##     "components" :  [
    ##         ("Trel",)
    ##         ]
    ##     }
    ##     )
    ## h2mixer_params["output_streams"].append(
    ##     {
    ##     "basename" : "tbme-h2-Lawson",
    ##     "components" :  [
    ##         ("NCM", hwLawson, aLawson) # Lawson
    ##         ]
    ##     }
    ##     )


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
    if ("obs-R20K20" not in current_task):
        current_task["obs-R20K20"] = False
    if ("obs-am-sqr" not in current_task):
        current_task["obs-am-sqr"] = False
    if ("em_multipolarity_list" not in current_task):
        current_task["em_multipolarity_list"] = []
    if ("keep_obdme" not in current_task):
        current_task["keep_obdme"] = True


    # set up code parameters
    Nv = current_task["Nv"]
    Nmin = current_task["Nmax"] % current_task["Nstep"]
    truncation = interaction_truncation_for_Nmax(Nv,current_task["Nmax"])
    ## BUG: through 130520: Nshell = truncation[1] + 1
    Nshell = truncation[0] + 1

    # set up h2 filenames
    h2_basename = "tbme-h2"
    obs_basename_list = [ "tbme-rrel2", "tbme-NCM" ]
    if (current_task["obs-R20K20"]):
        obs_basename_list += [ "tbme-R20", "tbme-K20" ]
    if (current_task["obs-am-sqr"]):
        for op in angular_momentum_operator_list:
            obs_basename_list.append("tbme-{}".format(op))

    # import partitioning file
    partitioning_filename = os.path.join(ncsm_config.data_dir_partitioning,"mfdn_partitioning.info_Nsh"+str(Nshell))
    print ("Checking for partition file %s..." % partitioning_filename)
    if (os.path.exists(partitioning_filename)):
        mcscript.call(["cp", partitioning_filename, "mfdn_partitioning.info"])
    else:
        print ("Not found.")
        
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
    # remove detailed occupation listing
    # TO NEATEN: replace with shell-free subprocess call and pipes
    ## print ("Filtering mfdn.res...") 
    ## mcscript.call(["sed '/More/,$d' < mfdn.res > " + result_filename],shell=True) 
    mcscript.call(["cat < mfdn.res > " + result_filename],shell=True) 

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
    mcscript.call(["mv", result_filename, archive_filename, "--target-directory="+mcscript.task.results_dir])

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

# September 20 2015 fixing the script to work with new reduced obdme


def task_handler_emcalc(current_task):

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
        #"2mj" : current_task["2mj"],
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
        /home/valentino/mfdn.rppobdme.info                                                             
        /home/valentino/mfdn.statrobdme.seq004.2J04.p0.n01.2T02                                                
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
    obdme_file_list = glob.glob("*seq*")
    print ("obdme files:", obdme_file_list)
    
    # result file
    main_result_filename = "%s-emcalc-%s.dat" % (mcscript.run.name, result_descriptor)
    main_result_file = open(main_result_filename,"w")

    #print(params["em_multipolarity_list"])


    # loop over transitions
    for obdme_file in obdme_file_list:
        for em_multipolarity in params["em_multipolarity_list"]: 
            # generate input lines
            #print(em_multipolarity)

            input_lines = [
                "%s %d" % em_multipolarity,
                "%g" % ( params["hw"]),
                "%.3f %.3f" % (params["beta_p"], params["beta_n"]),
                "mfdn.rppobdme.info",
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
    result_files = glob.glob("*.res") + glob.glob("*-emcalc-*.dat")
    mcscript.call(["tar", "-zcvf", archive_filename ] + result_files, cwd=mcscript.task.results_dir)

    # copy archive out to home results archive directory
    mcscript.call(["cp","-v",archive_filename,"-t",ncsm_config.data_dir_results_archive], cwd=mcscript.task.results_dir)

def archive_handler_mfdn_archive ():

    # make archive -- whole dir
    archive_filename = mcscript.task.archive_handler_generic ()
   
    # put to hsi
    hsi_subdir = "2016"
    hsi_arg = "lcd %s; cd %s; put %s" % (os.path.dirname(archive_filename), hsi_subdir, os.path.basename(archive_filename))
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

################################################################
# natural orbitals production -- nogen
################################################################

def task_pool_nogen(current_task):
    """
    task_pool_nogen(current_task): Returns the task name
    """

    return "nogen"

def task_handler_nogen(current_task):
    """ task_handler_nogen(current_task): Reads the dictionary keys found in current_task,
    it locates the zip file containing the obdme files, it unzips them and it creates
    the natural orbitals, and the occupation probabilities..
    """
    # Quick and dirty way to avoid searching for density files with the fci flag attached to them when we perform FCI
    # runs with the natural orbital basis..

    #for kk in range(len(mfdn_run_list)):

    current_task["fci"] = False
    #current_task["hw_int"] = 40.
    #current_task["xform_cutoff"] = ("ob",13)

    print ("nogen", task_descriptor_format_5(current_task))  

    # Find the *.tgz files from an mfdn run..
    # Append the subdirectory "results" in the search path. The search path should become an 
    # environmental variable declared in your .cshrc.ext file and read in by ncsm_config.py

    # subdirectory = mfdn_run_list[kk] + "/results"
    archive_filename = "%s-mfdn-%s.tgz" % (mfdn_run_list[0], task_descriptor_format_5(current_task))
    obdme_archive = os.path.join("/scratch1/scratchdirs/cconsta1/runs",mfdn_run_list[0],"results",archive_filename)  

    # current_task["hw_int"] = 10.
    current_task["fci"] = True

    print(obdme_archive)

    print("The hw value is: ", current_task["hw"])
 
    # TODO: Replace subprocess.call with a safer option
    # The reason for using subprocess.call invoking the shell is to be able to use a wildcard and extract only the relevant statobdme files from
    # the archive

    subprocess.call("tar -xzf" + obdme_archive + " --wildcards --no-anchored '*statrobdme*'",shell=True)
    subprocess.call("tar -xzf" + obdme_archive + " --wildcards --no-anchored '*mfdn.rppobdme.info*'",shell=True)
    # mcscript.call(["tar", "xvf", obdme_archive])

    # The natural orbitals are different for protons and neutrons 

    species_list = ["p","n"] 

    # Loop over the two species and declare the input parameters for nogen

    for species in species_list:
        nogen_params = {
            "Nmax" : current_task["Nmax"], "Nv" : current_task["Nv"], "species": species
            }

        # invoke nogen

        print(nogen_params, "\n", task_descriptor_format_5(current_task))

        call_nogen (nogen_params,task_descriptor_format_5(current_task))

        # save results locally

        natural_orbitals_filename_list = glob.glob("natural-orbital*"+species+"*")
        occupations_orbitals_filename_list = glob.glob("occupations-orbital*"+species+"*")
   
        # Changing the filename
  
        for natural_orbitals_file in natural_orbitals_filename_list:
            new_no_filename = natural_orbitals_file[:-4] + "-Nmax" + "{:02d}".format(current_task["Nmax"]) + "-hw" + "{:.3f}".format(current_task["hw"])+"-iter-0.dat"
            subprocess.call("cp " + natural_orbitals_file + " " + new_no_filename, shell=True)
            subprocess.call("mv " +  new_no_filename + " ../results/", shell=True)
        
        for occupations_orbitals_file in occupations_orbitals_filename_list:
            new_occupations_filename = occupations_orbitals_file[:-4] + "-Nmax" + "{:02d}".format(current_task["Nmax"]) + "-hw" + "{:.3f}".format(current_task["hw"])+"-iter-0.dat"
            subprocess.call("cp " + occupations_orbitals_file + " " + new_occupations_filename, shell=True)  
            subprocess.call("mv " +  new_occupations_filename + " ../results/", shell=True) 

    # delete all old densities after a successful run

    subprocess.call("rm mfdn.*", shell=True)
          
           
def call_nogen (params,result_descriptor):
    """ call_nogen(params,descriptor): Caller function that calls nogen
    """
    # find obdme files
    # TODO: Make the caller function flexible to use the obdme from states other than the ground..

    obdme_file = glob.glob("mfdn.statrobdme.seq001*") 

    print(obdme_file)
    
    input_lines = [
        "%d %d %s" % (params["Nmax"],params["Nv"], params["species"]),
        "mfdn.rppobdme.info",
        "%s" % obdme_file[0]
        ]

    print(input_lines)
    print("Run is: ", mcscript.run.name)

    mcscript.call_serial(["/global/u2/c/cconsta1/noutils/nogen"], mcscript.run, input_lines=input_lines)



################################################################
# producing the kinematic files in the natural orbital basis -- noradial
################################################################

def task_pool_noradial(current_task):

    return "noradial"

def task_handler_noradial(current_task):
    """ task_handler_noradial(): Handles the noradial execution providing all the necessary parameters
    """
    print ("noradial", current_task["descriptor"])   


    # TODO: Replace with environmental variable
    radial_files_path = "/scratch1/scratchdirs/cconsta1/data/radial-basis/130704/"

    # noradial is designed to work with one radial file at a time. Therefore here we provide a list
    # of all the necessary kinematic files that need to be converted over to the natural orbital basis 

    # TODO: Make the radial files more flexible. i.e, if one wants to use the natural orbitals for the CS basis
    # he should be able to do so..
    
    basis = current_task["basis"][0]

    radial_files_list = ['radial-me-' + basis + '-r0.dat', 'radial-me-' + basis + '-r1.dat', 'radial-me-'+ basis +'-r2.dat', 'radial-me-'+ basis +'-k1.dat','radial-me-'+ basis +'-k2.dat']

    species_list = ["p","n"]    

    for radial_file in radial_files_list:
        for species in species_list:
            input_radial_basename = radial_files_path + radial_file

            print(radial_file)

            noradial_params = {
                "Nmax" : current_task["Nmax"], "Nv" : current_task["Nv"], "species": species,
                "input_radial_basename": input_radial_basename
                }

            # invoke noradial

            print(noradial_params, "\n", current_task["descriptor"])

            call_noradial (noradial_params,current_task["descriptor"])

def call_noradial (params,result_descriptor):
    
    input_lines = [
        "%d %d %s" % (params["Nmax"],params["Nv"], params["species"]),
        "%s" % params["input_radial_basename"]
        ]

    print(input_lines)

    mcscript.call_serial(["/global/u2/c/cconsta1/noutils/noradial"], mcscript.run, input_lines=input_lines)  		   

################################################################
# natural orbitals xform files preparation -- noxform
################################################################

def task_pool_noxform(current_task):

    return "noxform"

def task_handler_noxform(current_task):

    print ("noxform", current_task["descriptor"])   

    radial_files_path = "/scratch1/scratchdirs/cconsta1/data/radial-basis/130704/"

    beta_interaction = "%.3f" % math.sqrt(current_task["hw_int"]/current_task["hw"])

    beta_interaction_coul = "%.3f" % math.sqrt(current_task["hw_int_coul"]/20.)

    print(beta_interaction," ",beta_interaction_coul)

    basis = current_task["basis"][0]

    #radial_xform_files_list = ["radial-xform-"+basis+"-" + beta_interaction + ".dat", "radial-xform-"+basis+"-" + beta_interaction_coul + ".dat"]  

    species_list = ["p","n"]  

    # JISP16 interaction transformation...

    #for radial_xform_file in radial_xform_files_list:
    for species in species_list:
        input_radial_basename = radial_files_path + "radial-xform-"+basis+"-" + beta_interaction + ".dat"

        print(input_radial_basename)

        noxform_params = {
            "Nmax" : current_task["Nmax"], "Nv" : current_task["Nv"], "species": species,
            "input_radial_basename": input_radial_basename
            }

        # invoke noxform

        print(noxform_params, "\n", current_task["descriptor"])

        call_noxform (noxform_params,current_task["descriptor"])

        subprocess.call("mv " +  "radial-xform-NO-"+species+"-" + beta_interaction + ".dat" + " radial-xform-NO-"+species+"-JISP16-" + beta_interaction + ".dat", shell=True)

    # VC interaction transformation...
    
    #print("Old hw: ", current_task["hw"], "\n")
    #subprocess.call("rm *orbital*", shell=True)
    #current_task["hw"] = 20.
    #current_task["hwLawson"] = 20.
    #print("New hw: ", current_task["hw"], "\n")
    
    #print("Running for hw = 20 MeV NOs..")
    #task_handler_nogen(current_task)


    #for radial_xform_file in radial_xform_files_list:
    for species in species_list:
        input_radial_basename = radial_files_path + "radial-xform-"+basis+"-" + beta_interaction_coul + ".dat"

        print(input_radial_basename)

        noxform_params = {
            "Nmax" : current_task["Nmax"], "Nv" : current_task["Nv"], "species": species,
            "input_radial_basename": input_radial_basename
            }

        # invoke noxform

        print(noxform_params, "\n", current_task["descriptor"])

        call_noxform (noxform_params,current_task["descriptor"])

        subprocess.call("mv " +  "radial-xform-NO-"+species+"-" + beta_interaction_coul + ".dat" + " radial-xform-NO-"+species+"-VC-" + beta_interaction_coul + ".dat", shell=True)

def call_noxform (params,result_descriptor):
    
    input_lines = [
        "%d %d %s" % (params["Nmax"],params["Nv"], params["species"]),
        "%s" % params["input_radial_basename"]
        ]

    print(input_lines)

    mcscript.call_serial(["/global/u2/c/cconsta1/noutils/noxform"], mcscript.run, input_lines=input_lines) 

################################################################
# natural orbitals -- wrapper handler
################################################################

def task_handler_no_wrapper(current_task):
    """ task_handler_no_wrapper() calls the task handlers for nogen,
    noradial and noxform so that all calculations can be performed in one 
    phase. Creating the natural orbitals, the radial files and the xform
    files is a quick process..
    """

    task_handler_nogen(current_task)
    task_handler_noradial(current_task)
    task_handler_noxform(current_task)

		# remove the natural orbital files from task directory to allow for iterations..

    subprocess.call("rm natural-orbital*", shell=True)
    subprocess.call("rm occupations-orbital*", shell=True)

		# however zip them into result directory..

    # archive_handler_natural_orbitals()

    # delete all densities

    #subprocess.call("rm mfdn.*", shell=True)

# Iterate natural orbitals new approach.. Just multiply the NOs using numpy..

"""
def task_handler_noiter(current_task):

    # obdme_archive = os.path.join("/scratch/scratchdirs/cconsta1/runs",mfdn_run_list[kk],"results",archive_filename)

    path_to_no_orbitals = os.path.join("/scratch/scratchdirs/cconsta1/runs", mcscript.run.name,"results")

    lmax = current_task["Nmax"] + current_task["Nv"]

    species_list = ["p", "n"]

    for l in range(lmax+1):
        for species in species_list:

            orbital_list_jmin = []
            orbital_list_jmax = []
            
            for kk in range(len(mfdn_run_list)):
                
                # Create input file
                # natural-orbital-2l8-2j9-p-Nmax04-hw40.000-iter-0.dat
                
                if(l==0):
                    natural_orbital_file = "natural-orbital-2l"+str(2*l)+"-2j-"+str(2*l+1)+"-"+species+"-Nmax" + "{:02d}".format(current_task["Nmax"]) + "-hw" + "{:.3f}".format(current_task["hw"])+"-iter-"+str(kk)+".dat"
                    orbital_list_jmin.append(natural_orbital_file)                
                else:
                    natural_orbital_file_jmin = "natural-orbital-2l"+str(2*l)+"-2j-"+str(2*l-1)+"-"+species+"-Nmax" + "{:02d}".format(current_task["Nmax"]) + "-hw" + "{:.3f}".format(current_task["hw"])+"-iter-"+str(kk)+".dat"
                    natural_orbital_file_jmax = "natural-orbital-2l"+str(2*l)+"-2j-"+str(2*l+1)+"-"+species+"-Nmax" + "{:02d}".format(current_task["Nmax"]) + "-hw" + "{:.3f}".format(current_task["hw"])+"-iter-"+str(kk)+".dat"
              
                    orbital_list_jmin.append(natural_orbital_file_jmin)
                    orbital_list_jmax.append(natural_orbital_file_jmax)

            if(l==0):
                if(len(orbital_list_jmin)==1):
                    no_orbital = numpy.loadtxt(path_to_no_orbitals+orbital_list_jmin[0])
                    no_orbital_filename = "natural-orbital-2l"+str(2*l)+"-2j-"+str(2*l-1)+"-"+species+".dat"
                    numpy.savetxt(no_orbital_filename,no_orbital)
                else:
                    no_orbital_1 = 

"""

# Iterate the natural orbitals -- Quick and dirty approach but its code duplication... 09/15/15

################################################################
# natural orbitals production -- nogen -- iteration phase
################################################################

def task_handler_nogen_iter(current_task):
    """ task_handler_nogen(current_task): Reads the dictionary keys found in current_task,
    it locates the zip file containing the obdme files, it unzips them and it creates
    the natural orbitals, and the occupation probabilities..
    """
    # Quick and dirty way to avoid searching for density files with the fci flag attached to them when we perform FCI
    # runs with the natural orbital basis..

    # current_task["descriptor"]

    print ("nogen_iter", task_descriptor_format_5(current_task))  

    # Find the *.tgz files from an mfdn run..
    # Append the subdirectory "results" in the search path. The search path should become an 
    # environmental variable declared in your .cshrc.ext file and read in by ncsm_config.py

    # subdirectory = mfdn_run_list[0] + "/results"

		# mfdn_run_list[1] is the second iteration

    # task_descriptor_format_5(current_task)

    archive_filename = "%s-mfdn-%s.tgz" % (mfdn_run_list[1], task_descriptor_format_5(current_task))

    print(archive_filename)

    obdme_archive = os.path.join("/scratch1/scratchdirs/cconsta1/runs",mfdn_run_list[1],"results",archive_filename)  

    print(obdme_archive)
 
    # TODO: Replace subprocess.call with a safer option
    # The reason for using subprocess.call invoking the shell is to be able to use a wildcard and extract only the relevant statobdme files from
    # the archive

    subprocess.call("tar -xzf" + obdme_archive + " --wildcards --no-anchored '*statrobdme*'",shell=True)
    subprocess.call("tar -xzf" + obdme_archive + " --wildcards --no-anchored '*mfdn.rppobdme.info*'",shell=True)
    # mcscript.call(["tar", "xvf", obdme_archive])

    # The natural orbitals are different for protons and neutrons 

    species_list = ["p","n"] 

    # Loop over the two species and declare the input parameters for nogen

    for species in species_list:
        nogen_params = {
            "Nmax" : current_task["Nmax"], "Nv" : current_task["Nv"], "species": species
            }

        # invoke nogen

        print(nogen_params, "\n", task_descriptor_format_5(current_task))

        call_nogen (nogen_params, task_descriptor_format_5(current_task))

        # save results locally

        natural_orbitals_filename_list = glob.glob("natural-orbital*"+species+"*")
        occupations_orbitals_filename_list = glob.glob("occupations-orbital*"+species+"*")
   
        # Changing the filename
  
        for natural_orbitals_file in natural_orbitals_filename_list:
            new_no_filename = natural_orbitals_file[:-4] + "-Nmax" + "{:02d}".format(current_task["Nmax"]) + "-hw" + "{:.3f}".format(current_task["hw"])+"-iter1.dat"
            subprocess.call("cp " + natural_orbitals_file + " " + new_no_filename, shell=True)
            subprocess.call("mv " +  new_no_filename + " ../results/", shell=True)
        
        for occupations_orbitals_file in occupations_orbitals_filename_list:
            new_occupations_filename = occupations_orbitals_file[:-4] + "-Nmax" + "{:02d}".format(current_task["Nmax"]) + "-hw" + "{:.3f}".format(current_task["hw"])+"-iter1.dat"
            subprocess.call("cp " + occupations_orbitals_file + " " + new_occupations_filename, shell=True)  
            subprocess.call("mv " +  new_occupations_filename + " ../results/", shell=True)      

    # delete all old densities after a successful run

    subprocess.call("rm mfdn.*", shell=True)      

################################################################
# producing the kinematic files in the natural orbital basis -- noradial -- iteration
################################################################

def task_handler_noradial_iter(current_task):
    """ task_handler_noradial(): Handles the noradial execution providing all the necessary parameters
    """
    print ("noradial", current_task["descriptor"]) 

    # noradial is designed to work with one radial file at a time. Therefore here we provide a list
    # of all the necessary kinematic files that need to be converted over to the natural orbital basis 

    radial_files_list = ["r1.dat", "r2.dat", "k1.dat","k2.dat"]

    species_list = ["p","n"]    

    for radial_file in radial_files_list:
        for species in species_list:

            print(radial_file)

            noradial_params = {
                "Nmax" : current_task["Nmax"], "Nv" : current_task["Nv"], "species": species,
                "input_radial_basename": "radial-me-NO-"+ species + "-" + radial_file
                }

            # invoke noradial

            print(noradial_params, "\n", current_task["descriptor"])

            call_noradial_iter (noradial_params,current_task["descriptor"])  		   


def call_noradial_iter (params,result_descriptor):
    
    input_lines = [
        "%d %d %s" % (params["Nmax"],params["Nv"], params["species"]),
        "%s" % params["input_radial_basename"]
        ]

    print(input_lines)

    mcscript.call_serial(["/global/u2/c/cconsta1/noutils/noradial_iter"], mcscript.run, input_lines=input_lines)


################################################################
# natural orbitals xform files preparation -- noxform -- iteration
################################################################

def task_handler_noxform_iter(current_task):

    print ("noxform", current_task["descriptor"])   

    beta_interaction = "%.3f" % math.sqrt(current_task["hw_int"]/current_task["hw"])

    beta_interaction_coul = "%.3f" % math.sqrt(current_task["hw_int_coul"]/20.)

    print(beta_interaction," ",beta_interaction_coul)

    # radial_xform_files_list = [ beta_interaction + ".dat", beta_interaction_coul + ".dat"]  

    species_list = ["p","n"]  

    # JISP16 transformation..

    #for radial_xform_file in radial_xform_files_list:
    for species in species_list:
        noxform_iter_params = {
            "Nmax" : current_task["Nmax"], "Nv" : current_task["Nv"], "species": species,
            "input_radial_basename": "radial-xform-NO-" + species + "-JISP16-" + beta_interaction + ".dat"
            }

        # invoke noxform_iter

        print(noxform_iter_params, "\n", current_task["descriptor"])

        call_noxform_iter (noxform_iter_params,current_task["descriptor"])

        subprocess.call("mv " +  "radial-xform-NO-"+species+"-" + beta_interaction + ".dat" + " radial-xform-NO-"+species+"-JISP16-" + beta_interaction + ".dat", shell=True)

    # VC interaction transformation...
    
    #print("Old hw: ", current_task["hw"], "\n")
    #subprocess.call("rm *orbital*", shell=True)
    #current_task["hw"] = 20.
    #current_task["hwLawson"] = 20.
    #print("New hw: ", current_task["hw"], "\n")
    
    #print("Running for hw = 20 MeV NOs..")
    #task_handler_nogen_iter(current_task)


    #for radial_xform_file in radial_xform_files_list:
    for species in species_list:
        input_radial_basename = "radial-xform-NO-" + species + "-VC-" + beta_interaction_coul + ".dat"

        print(input_radial_basename)

        noxform_iter_params = {
            "Nmax" : current_task["Nmax"], "Nv" : current_task["Nv"], "species": species,
            "input_radial_basename": input_radial_basename
            }

        # invoke noxform

        print(noxform_iter_params, "\n", current_task["descriptor"])

        call_noxform_iter (noxform_iter_params,current_task["descriptor"])

        subprocess.call("mv " +  "radial-xform-NO-"+species+"-" + beta_interaction_coul + ".dat" + " radial-xform-NO-"+species+"-VC-" + beta_interaction_coul + ".dat", shell=True)

def call_noxform_iter (params,result_descriptor):
    
    input_lines = [
        "%d %d %s" % (params["Nmax"],params["Nv"], params["species"]),
        "%s" % params["input_radial_basename"]
        ]

    print(input_lines)

    mcscript.call_serial(["/global/u2/c/cconsta1/noutils/noxform_iter"], mcscript.run, input_lines=input_lines)

################################################################
# natural orbitals -- wrapper handler -- iteration
################################################################

def task_handler_no_wrapper_iter(current_task):
    """ task_handler_no_wrapper() calls the task handlers for nogen,
    noradial and noxform so that all calculations can be performed in one 
    phase. Creating the natural orbitals, the radial files and the xform
    files is a quick process..
    """

    task_handler_nogen_iter(current_task)
    task_handler_noradial_iter(current_task)
    task_handler_noxform_iter(current_task)

		# remove the natural orbital files from task directory to allow for iterations..

    # subprocess.call("rm *orbital*", shell=True)

		# however zip them into result directory..

    # archive_handler_natural_orbitals()

# end iteration code

################################################################
# h2 production -- h2gen-NO
################################################################

def task_pool_h2gen_no(current_task):
    """
    task_descriptor_h2gen_no(...) create task pool.
    """

    return "h2gen-NO"

##def task_descriptor_h2gen_no(current_task):
#    """
#    task_descriptor_h2gen_no(...) create task descriptor from basis description.
#    """

#    return basis_string_pn(current_task["basis"])


def task_handler_h2gen_no(current_task):
    """
    task_handler_h2gen_no(...) runs h2gen for the current task parameters and 
    the NO basis.    
    """
    # recover current basis
    # NO update. Must add a new dictionary key called no_basis
    basis = current_task["basis"]
    no_basis = current_task["no_basis"]
    Nv = current_task["Nv"]
    Nmax = current_task["Nmax"]

    if(current_task.get("fci",False)):
        print ("h2gen:", basis_string_pn(basis), Nmax, "->", interaction_truncation_for_Nmax(Nv,Nmax,standardize=True,fci=current_task["fci"]))
    else:
        print ("h2gen:", basis_string_pn(basis), Nmax, "->", interaction_truncation_for_Nmax(Nv,Nmax))
        # set radial bases

    (no_basis_name_p, beta_p, no_basis_name_n, beta_n) = no_basis

    radial_basename_p = "radial-me-NO-p"  
    
    radial_basename_n = "radial-me-NO-n" 
    # radial_basename_n = mcscript.subpath_search_basename (ncsm_config.data_dir_radial, radial_run_list, radial_basename_n_unresolved)
    # if (radial_basename_n is None):
    #    print ("Cannot find radial file", radial_basename_n_unresolved, "in", ncsm_config.data_dir_radial, radial_run_list)
    #    sys.exit(1)

    # set truncation
    if(current_task.get("fci",False)):
        (N1b, N2b) = interaction_truncation_for_Nmax(Nv,Nmax,standardize=True,fci=current_task["fci"])
    else:
        (N1b, N2b) = interaction_truncation_for_Nmax(Nv,Nmax)

    # set output filname
    output_basename = mcscript.dashify(("r2k2",basis_string_pn(no_basis), N1b, N2b))

    # run h2gen
    h2gen_params = {
        "N1b" : N1b, "N2b" : N2b,
        "radial_basename_p" : radial_basename_p, "radial_basename_n" : radial_basename_n,
        "beta_p" : beta_p, "beta_n" : beta_n,
        "output_basename" : output_basename
        }
    call_h2gen (h2gen_params)


################################################################
# h2 production -- xform -- natural orbital basis
################################################################

def task_pool_h2xform_no(current_task):
    """
    task_pool_h2xform_no(...) create task pool.
    """

    return "h2xform-NO"

def task_pool_h2xform_no_Nmax(current_task):
    """
    task_pool_h2xform_no_Nmax(...) create task pool showing last Nmax value.
    """

    return "Nmax{Nmax:02d}".format(Nmax=current_task["Nmax_list"][-1])


def task_descriptor_h2xform_no(current_task):
    """
    task_descriptor_h2xform(...) create task descriptor from basis description and scaling.
    """
    return "%s-%s-%g" % (
        current_task["interaction"],
        basis_string_pn(current_task["scaled_basis"]),
        current_task["hw_int"]
    )

def task_handler_h2xform_no(current_task):
    """
    task_handler_h2xform_no(...) runs h2xform for the current task's basis and interaction.
    """
    
    print ("h2xform_no", current_task["descriptor"])

    # recall parameters

    interaction = current_task["interaction"]
    interaction_truncation = current_task["interaction_truncation"] 
    no_scaled_basis = current_task["no_scaled_basis"] 
    no_basis = current_task["no_basis"]
    hw_int = current_task["hw_int"] 
    Nv = current_task["Nv"]
    xform_cutoff_list = current_task["xform_cutoff_list"]
    Nmax = current_task["Nmax"]   
    cutoff = current_task["xform_cutoff"]

    # set radial bases
    (basis_name_p, beta_p, basis_name_n, beta_n) = no_scaled_basis

    # radial basename for the strong interaction
    radial_basename_p = "radial-xform-%s-%s-%.3f" % (basis_name_p, "JISP16" ,beta_p)
    radial_basename_n = "radial-xform-%s-%s-%.3f" % (basis_name_n, "JISP16" ,beta_n)

    print(radial_basename_p)

    # set input TBME file
    input_basename_unresolved = "%s-%s-%g" % (interaction, mcscript.dashify(interaction_truncation), hw_int)
    input_basename = mcscript.subpath_search_basename (ncsm_config.data_dir_h2, interaction_run_list, input_basename_unresolved)
    print ("Interaction TBME file:", input_basename_unresolved)
    if (input_basename is None):
        print ("Missing interaction...")
        exit(1)

    output_streams = []

    #for cutoff in xform_cutoff_list:
    #for Nmax in Nmax_list:

    # set output truncation
    if(current_task.get("fci",False)):
            output_truncation = interaction_truncation_for_Nmax(Nv,Nmax,standardize=True,fci=current_task["fci"])
    else:
            output_truncation = interaction_truncation_for_Nmax(Nv,Nmax)

    # set output basename
    output_basename = "%s-%s-%s-%s-%g" % (interaction, basis_string_pn(no_scaled_basis), mcscript.dashify(cutoff), mcscript.dashify(output_truncation), hw_int)

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

    #TODO: Figure out a better way to implement converting the two interaction files at the same time because this is
    # pure code duplication..

    # Coulomb interaction
    
    interaction = "VC"
    hw_int_coul = current_task["hw_int_coul"]
    (basis_name_p, beta_p, basis_name_n, beta_n) = no_basis
    
    radial_basename_p = "radial-xform-%s-%s-%.3f" % (basis_name_p, "VC", beta_p)
    radial_basename_n = "radial-xform-%s-%s-%.3f" % (basis_name_n, "VC", beta_n)

    # set input TBME file
    input_basename_unresolved = "%s-%s-%g" % (interaction, mcscript.dashify(interaction_truncation), hw_int_coul)
    input_basename = mcscript.subpath_search_basename (ncsm_config.data_dir_h2, interaction_run_list, input_basename_unresolved)
    print ("Interaction TBME file:", input_basename_unresolved)
    if (input_basename is None):
        print ("Missing interaction...")
        exit(1)

    output_streams = []

    #for cutoff in xform_cutoff_list:
    #for Nmax in Nmax_list:

    # set output truncation
    if(current_task.get("fci",False)):
            output_truncation = interaction_truncation_for_Nmax(Nv,Nmax,standardize=True,fci=current_task["fci"])
    else:
            output_truncation = interaction_truncation_for_Nmax(Nv,Nmax)

    # set output basename
    output_basename = "%s-%s-%s-%s-%g" % (interaction, basis_string_pn(no_basis) , mcscript.dashify(cutoff), mcscript.dashify(output_truncation), hw_int_coul)

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
    
      
################################################################
# h2 production -- h2mixer -- natural orbital basis
################################################################
  

def task_handler_h2mixer_no(current_task):
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
    no_basis = current_task["no_basis"] 
    no_scaled_basis = current_task["no_scaled_basis"] 
    xform_cutoff = current_task["xform_cutoff"] 
    Nmax = current_task["Nmax"] 
    Nv = current_task["Nv"]

    # basic setup

    (Z, N) = nuclide
    if(current_task.get("fci",False)):
        input_truncation = interaction_truncation_for_Nmax(Nv,Nmax,standardize=True,fci=current_task["fci"])
    else:
        input_truncation = interaction_truncation_for_Nmax(Nv,Nmax)

    r2k2_basename_unresolved = mcscript.dashify(("r2k2", basis_string_pn(no_basis), mcscript.dashify(input_truncation)))
    #r2k2_basename = mcscript.subpath_search_basename (ncsm_config.data_dir_h2, r2k2_run_list, r2k2_basename_unresolved)

    print(r2k2_basename_unresolved)
    #print(r2k2_basename)
    
    if(current_task.get("fci",False)):
        output_truncation = interaction_truncation_for_Nmax(Nv,Nmax,standardize=False,fci=current_task["fci"])
    else:
        output_truncation = interaction_truncation_for_Nmax(Nv,Nmax,standardize=False)

    h2mixer_params = {
        "input_truncation" : input_truncation,
        "output_truncation" : output_truncation,
        "Z": Z, "N" : N, "hw" : hw,
        "r2k2_basename" : r2k2_basename_unresolved,
        "input_basenames" : [],
        "output_streams" : []
        }

    # input stream setup
    # input 5 -- built interaction matrix elements -- for scaled basis
    tbme_basename_VNN_unresolved = mcscript.dashify((interaction, basis_string_pn(no_scaled_basis), mcscript.dashify(xform_cutoff), mcscript.dashify(input_truncation), format(hw_int,"g")))
    h2mixer_params["input_basenames"].append(
        #mcscript.subpath_search_basename (ncsm_config.data_dir_h2, xform_run_list, tbme_basename_VNN_unresolved)
        tbme_basename_VNN_unresolved
        )
    # input 6 (OPTIONAL) -- built Coulomb matrix elements -- for unscaled basis
    if (coulomb):
        hw_int_coul = current_task["hw_int_coul"] 
        xform_cutoff_coul = current_task["xform_cutoff_coul"] 
        tbme_basename_VC_unresolved = mcscript.dashify(("VC", basis_string_pn(no_basis), mcscript.dashify(xform_cutoff_coul), mcscript.dashify(input_truncation), format(hw_int_coul,"g")))
        h2mixer_params["input_basenames"].append(
            #mcscript.subpath_search_basename (ncsm_config.data_dir_h2, xform_run_list, tbme_basename_VC_unresolved)
            tbme_basename_VC_unresolved
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
                    "basename" : "tbme-{}".format(op),
                    "components" : [ ("am-sqr", op) ]
                }
            )

    # output streams -- tbme for H components (debugging)
    ## h2mixer_params["output_streams"].append(
    ##     {
    ##     "basename" : "tbme-h2-v2",
    ##     "components" :  [
    ##         ("in", 5, 1.0), 
    ##         ("scaled", 6, 20.0, -1,coulomb_factor)
    ##         ]
    ##     }
    ##     )
    ## h2mixer_params["output_streams"].append(
    ##     {
    ##     "basename" : "tbme-h2-Trel",
    ##     "components" :  [
    ##         ("Trel",)
    ##         ]
    ##     }
    ##     )
    ## h2mixer_params["output_streams"].append(
    ##     {
    ##     "basename" : "tbme-h2-Lawson",
    ##     "components" :  [
    ##         ("NCM", hwLawson, aLawson) # Lawson
    ##         ]
    ##     }
    ##     )


    # run h2mixer
    call_h2mixer (h2mixer_params)

##################################################################
# emcalc  -- natural orbitals version
##################################################################

# For transitions in the NO basis we need to use the newly generated 
# radial files and of course the densities from the NO runs
# Also the executable should be emcalc_lj

def task_handler_emcalc_no(current_task):

    # Performance note on NERSC: When run as a normal job, this
    # scripting involves many aprun call, each of which carries
    # considerable overhead (at least ~1 sec on hopper), so total execution
    # time is vastly slower than in epar mode, which has initial node setup
    # overhead but then runs the script locally on the compute node, so no aprun in needed.
    # Compare run0252 task 0007 (serial) at 909 sec with task 0009 (epar) at 37 sec.

    print ("emcalc", current_task["descriptor"])

    # set radial bases
    basis = current_task["no_basis"] 
    (basis_name_p, beta_p, basis_name_n, beta_n) = basis
    radial_basename_p = "radial-me-%s" % (basis_name_p, )
    #radial_basename_p = mcscript.subpath_search_basename (ncsm_config.data_dir_radial, radial_run_list, radial_basename_p_unresolved)
    if (radial_basename_p is None):
        print ("Cannot find radial file", radial_basename_p)
        sys.exit(1)
    radial_basename_n = "radial-me-%s" % (basis_name_n, )
    #radial_basename_n = mcscript.subpath_search_basename (ncsm_config.data_dir_radial, radial_run_list, radial_basename_n_unresolved)
    if (radial_basename_n is None):
        print ("Cannot find radial file", radial_basename_n)
        sys.exit(1)

    # basic setup
    emcalc_params = {
        "em_multipolarity_list" : current_task["em_multipolarity_list"],
        #"2mj" : current_task["2mj"],
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
    call_emcalc_no (emcalc_params,current_task["descriptor"])

    # move results out
    main_result_filename = "%s-emcalc-%s.dat" % (mcscript.run.name, current_task["descriptor"])
    mcscript.call(["mv", "-v", main_result_filename, "--target-directory=" + mcscript.task.results_dir])

def call_emcalc_no (params,result_descriptor):
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
        /home/valentino/mfdn.rppobdme.info                                                             
        /home/valentino/mfdn.statrobdme.seq004.2J04.p0.n01.2T02                                                
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
    obdme_file_list = glob.glob("*seq*")
    print ("obdme files:", obdme_file_list)
    
    # result file
    main_result_filename = "%s-emcalc-%s.dat" % (mcscript.run.name, result_descriptor)
    main_result_file = open(main_result_filename,"w")

    #print(params["em_multipolarity_list"])


    # loop over transitions
    for obdme_file in obdme_file_list:
        for em_multipolarity in params["em_multipolarity_list"]: 
            # generate input lines
            #print(em_multipolarity)

            input_lines = [
                "%s %d" % em_multipolarity,
                "%g" % ( params["hw"]),
                "%.3f %.3f" % (params["beta_p"], params["beta_n"]),
                "mfdn.rppobdme.info",
                "%s" % (obdme_file),
                "%s %s" % (params["radial_basename_p"], params["radial_basename_n"])
                ]

            # run emcalc
            mcscript.call_serial([ncsm_config.emcalc_lj], mcscript.run, input_lines=input_lines)

            # harvest results to main results file
            temporary_result_filename = "emcalc.out"
            temporary_result_file = open(temporary_result_filename,"r")
            shutil.copyfileobj(temporary_result_file,main_result_file)
            temporary_result_file.close()

    # finalize main result file
    main_result_file.close()





def archive_handler_natural_orbitals ():    

    # make archive -- results
    archive_filename = os.path.join(
        ##ncsm_config.data_dir_results_archive,
        mcscript.task.archive_dir,
        "%s-natural-orbitals-%s.tgz" % (mcscript.run.name, mcscript.date_tag())
        )
    ## # store toc -- TODO once restructure subdirectories in tar file
    ## mcscript.call(
    ##     ["tar", "zcvf", archive_filename, toc_filename]
    ##     )
    os.chdir (mcscript.task.results_dir)
    result_files = glob.glob("*orbital*") 
    mcscript.call(["tar", "-zcvf", archive_filename ] + result_files, cwd=mcscript.task.results_dir)

    # copy archive out to home results archive directory
    mcscript.call(["cp","-v",archive_filename,"-t",ncsm_config.data_dir_results_archive], cwd=mcscript.task.results_dir)
		
		# delete files from results directory..
    # subprocess.call("rm *orbital*", shell=True)


