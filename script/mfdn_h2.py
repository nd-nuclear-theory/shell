"""mfdn_h2.py -- define scripting for h2mixer+MFDn runs

  Expected dictionary keys for h2mixer + mfdn run:

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


  M. A. Caprio
  University of Notre Dame

  - 12/14/16 (mac): Created, drawing on ncsm.py (created 2/12/13) and
    shell package example code generate_input.py (created 11/7/16).

"""
  
import os
import glob
import shutil
import datetime

import mcscript
import ncsm_config


################################################################
# utility calculations
################################################################

def weight_max_string(truncation):
    """ Convert (rank,cutoff) to "wp wn wpp wnn wpn" string.

    >>> weight_max_string(("ob",4))
        "4 4 8 8 8"
    >>> weight_max_string(("tb",4))
        "4 4 4 4 4"
    """

    (code, N) = truncation
    if (code == "ob"):
        cutoffs = (N,2*N)
    elif (code == "tb"):
        cutoffs = (N,N)

    return "{0[0]} {0[0]} {0[1]} {0[1]} {0[1]}".format(cutoffs)

def oscillator_length(hw):
    """ Calculate oscillator length for given oscillator frequency.

    b(hw) = (hbar c)/[(m_N c^2) (hbar omega)]^(1/2)

    Arguments:
        hw (numeric): hbar omega in MeV

    Returns:
        (float): b in fm
    """
    
    return k_hbar_c/math.sqrt(k_mN_csqr*hw)

################################################################
# configuration
################################################################

class Configuration(object): 
    """Object to MFDn environment configuration parameters into common
    name space.

    Defines several variables based on information provided through
    the environment:
    
        run: 

    """
    def __init__(self):

        # environment
        data_dir_h2 = os.environ.get("SHELL_DATA_DIR_H2")
        print("SHELL_DATA_DIR_H2",data_dir_h2)
        install_dir = os.environ.get("SHELL_INSTALL_DIR")
        print("SHELL_INSTALL_DIR",install_dir)
        interaction_run_list = []

    def shell_filename(self,code_name):
        """Construct filename for shell package executable."""
        return os.path.join(install_dir,"bin",code_name)


configuration = Configuration()

################################################################
# task descriptor for h2mixer + mfdn run
################################################################


def task_descriptor_format_7(current_task):
    """ runxxxx -- (in progress)
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
        **current_task
        )

    return descriptor

################################################################
# observable definitions
################################################################

angular_momentum_operator_list = [ 
    "L", "Sp", "Sn", "S", "J"
]


def task_handler_ho(current_task):
    """ Invoke mfdn in h2 run, given current task parameters.

    Expected dictionary keys: see initial module docstring
    """

    # temporary filenames
    orbitals_int_filename = "orbitals-int.dat"  # for overlaps from interaction tbmes
    orbitals_filename = "orbitals.dat"  # for target space

    # orbital generation -- interaction files
    mcscript.utils.call(
        [
            configuration.shell_filename("orbital-gen"),
            "--oscillator",
            "{Nmax_orb_int:d}".format(current_task),
            "{:s}".format(orbitals_int_filename)
        ]

    # orbital generation -- target space
    mcscript.utils.call(
        [
            configuration.shell_filename("orbital-gen"),
            "--oscillator",
            "{Nmax_orb:d}".format(current_task),
            "{:s}".format(orbitals_filename)
        ]

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
           
    
if (___name___ == "___MAIN___"):
    pass
