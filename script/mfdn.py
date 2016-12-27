"""mfdn.py -- define scripting for h2mixer+MFDn runs

    Expected dictionary keys for h2mixer + mfdn run:

        # nuclide parameters
        "nuclide" (tuple of int): (Z,N) tuple
        
        # Hamiltonian parameters
        "interaction" (string): name for interaction (for use in descriptor and filenames)
        "coulomb" (bool): whether or not to include Coulomb
        "aLawson" (float): aLawson
        "hwLawson" (float): hwLawson (typically hw of basis)

        # basis parameters
        "basis_mode" (int): enumerated value indicating direct oscillator, dilated oscillator,
            or generic run mode
        "hw" (float): hw of basis
        "hw_int" (float): hw of basis for source interaction TBMEs
        "truncation_int" (tuple): input interaction TBME cutoff, as tuple ("ob"|"tb",N) 
        "hw_coul" (float): hw of basis for source Coulomb TBMEs
        "hw_coul_prime" (float): hw to which to rescale Coulomb TBMEs before two-body transformation
        "truncation_coul" (tuple): input Coulomb tbme cutoff, as tuple ("ob"|"tb",N) 

        # transformation parameters
        "xform_truncation_int" (tuple): transform cutoff for interaction, as tuple ("ob"|"tb",N) 
        "xform_truncation_coul" (tuple): transform cutoff for Coulomb, as tuple ("ob"|"tb",N) 

        # legacy -- to adapt or remove
        "basis" (tuple): logical basis (radial_basis_p,1.,radial_basis_n,beta_n), to be scaled by hw
        "scaled_basis" (tuple): deduced basis used in task descriptor, etc.

        # traditional oscillator many-body truncation
        "ho_truncation" (bool): whether or not to assume traditional HO-style many-body truncation
        "Nv" (int): N of valence shell (for use in truncation)
        "Nmax" (int): Nmax
        "Nstep" (int): Nstep (2 for single parity, 1 for mixed parity)

        # diagonalization parameters
        "Mj" (float): M-scheme Mj value
        "eigenvectors" (int): number of eigenvectors to calculate
        "initial_vector" (int): initial vector code for mfdn
        "lanczos" (int): lanczos iterations
        "tolerance" (float): diagonalization tolerance parameter

        # obdme parameters
        "obdme_multipolarity" (int): maximum multipolarity for calculation of densities
        "obdme_reference_state_list" (list): list of reference states (J,g,i) for density calculation, or "all2all"
        "save_obdme" (bool): whether or not to save obdme files in archive

        # two-body observable parameters
        "obs-R20K20" (bool, optional): whether or not to include center-of-mass diagnostic observables (default: False)
        "obs-am-sqr" (bool, optional): whether or not to include squared angular momentum (default: False)
        "obs-H-components" (bool, optional): whether or not to include Hamiltonian terms (default: False)
        
        # version parameters
        "mfdn_executable" (string): mfdn executable name
        "mfdn_wrapper" (string): wrapper code for given MFDn version

        # emcalc postprocessing parameters
        "em_multipolarity_list" (list, optional): list of multipolarities ("E"|"M",lambda) 


  Mark A. Caprio
  University of Notre Dame

  - 12/14/16 (mac): Created, drawing on ncsm.py (created 2/12/13) and
    shell package example code generate_input.py (created 11/7/16).

"""
  
import datetime
import glob
import os
import shutil

import mcscript

################################################################
# radial basis modes
################################################################

# General modes of operation for radial basis
#
# k_radial_direct:
#   - final basis is oscillator basis (hw)
#   - source basis for VNN is oscillator basis of same
#     oscillator length (hw_int=hw); therefore no transformation
#     needed on VNN TBMEs
#   - Coulomb TBMEs need only scaling for dilation (hw_c -> hw)
#   - MFDn can use built-in oscillator OBMEs for observables
#
# k_radial_dilated:
#   - final basis is oscillator basis (hw)
#   - source basis for VNN is oscillator basis of different
#     oscillator length; therefore transformation
#     needed on VNN TBMEs (hw_int -> hw)
#   - Coulomb TBMEs need only scaling for dilation (hw_c -> hw)
#   - MFDn can use built-in oscillator OBMEs for observables
#
# k_radial_generic:
#   - final basis is not assumed to be oscillator basis
#     (still has nominal hw to define a length parameter)
#   - transformation needed on VNN TBMEs (hw_int HO -> hw generic)
#   - Coulomb TBMEs may be rescaled (hw_c -> hw_cp) but then need
#     transformation (hw_cp HO -> hw generic)
#   - MFDn *cannot* use built-in oscillator OBMEs for observables

k_radial_direct = 0
k_radial_dilated = 1
k_radial_generic = 2

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
    """Object to collect MFDn environment configuration parameters into
    common name space.

    Defines several variables based on information provided through
    the environment:
    
        run:

    Methods:

        shell_filename(name): Generate full qualified filename for
        shell utility code, given root name.

        interaction_filename(name): Generate full qualified filename
        for interaction h2 file, given root name.

    """
    def __init__(self):

        # environment
        self.data_dir_h2 = os.environ.get("SHELL_DATA_DIR_H2")
        self.install_dir = os.environ.get("SHELL_INSTALL_DIR")
        self.interaction_run_list = []

    def shell_filename(self,name):
        """Construct filename for shell package executable."""
        return os.path.join(self.install_dir,"bin",name)

    def interaction_filename(self,name):
        """Construct filename for interaction h2 file."""
        return os.path.join(self.data_dir_h2,name)


configuration = Configuration()

################################################################
# task descriptor for h2mixer + mfdn run
################################################################


def task_descriptor_format_7(task):
    """ runxxxx -- (in progress)
    """

    if(task["traditional_ho"]):
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
        basis_string_pn=basis_string_pn(task["basis"]),
        mixed_parity_flag=mcscript.ifelse(task["Nstep"]==1,"x",""),
        fci_flag = mcscript.ifelse(task.get("fci",False),"-fci",""),
        xform_cutoff_name=mcscript.dashify(task["xform_cutoff"]),
        **task
        )

    return descriptor

################################################################
# observable definitions
################################################################

angular_momentum_operator_list = [ 
    "L", "Sp", "Sn", "S", "J"
]


##################################################################
# traditional ho run
##################################################################

def task_handler_mfdn_analytic(task):
    """Carry out MFDn run in analytic basis: direct oscillator, dilated
    oscillator, or generic (TODO).

    Limitation: Currently only supports harmonic oscillator style
    truncation.

    Expected dictionary keys: see initial module docstring

    """

    print("task:",task)

    # orbital filenames
    orbitals_int_filename = "orbitals-int.dat"  # orbitals for interaction tbme basis
    orbitals_coul_filename = "orbitals-coul.dat"  # orbitals for Coulomb tbme basis
    orbitals_filename = "orbitals.dat"  # orbitals for target basis

    # radial matrix element filenames (for kinematic radial matrix elements)
    #
    # "{}{}" will be replaced by {"r1","r2","k1","k2"}
    radial_me_filename_template = "radial-me-{}{}.dat"

    # radial overlap filenames
    radial_olap_int_filename = "radial-olap-int.dat"  # overlaps from interaction tbme basis
    radial_olap_coul_filename = "radial-olap-coul.dat"  # overlaps from Coulomb tbme basis

    # validate basis mode
    if (task["basis_mode"] not in {k_radial_direct,k_radial_dilated}):  # no k_radial_generic yet
        raise ValueError("invalid basis mode {basis_mode}".format(**task))
    if (not task["ho_truncation"]):
        raise ValueError("expecting ho_truncation to be True but found {ho_truncation}".format(**task))
    
    # basis radial code -- expected by radial_utils codes
    basis_radial_code = "oscillator"  # TO GENERALIZE: if not oscillator basis

    # generate orbitals -- interaction bases
    mcscript.call(
        [
            configuration.shell_filename("orbital-gen"),
            "--oscillator",
            "{truncation_int[1]:d}".format(**task),
            "{:s}".format(orbitals_int_filename)
        ]
    )
    if (task["coulomb"]):
        mcscript.call(
            [
                configuration.shell_filename("orbital-gen"),
                "--oscillator",
                "{truncation_coul[1]:d}".format(**task),
                "{:s}".format(orbitals_coul_filename)
            ]
        )

    # generate orbitals -- target basis
    if (task["ho_truncation"]):
        Nmax_orb = task["Nmax"] + task["Nv"]
        mcscript.call(
            [
                configuration.shell_filename("orbital-gen"),
                "--oscillator",
                "{Nmax_orb:d}".format(Nmax_orb=Nmax_orb),
                "{:s}".format(orbitals_filename)
            ]
        )

    # generate radial integrals
    for operator_type in ["r","k"]:
        for power in [1,2]:
            radial_me_filename = radial_me_filename_template.format(operator_type,power)
            mcscript.call(
                [
                    configuration.shell_filename("radial-gen"),
                    "--kinematic",
                    "{:s}".format(operator_type),
                    "{:d}".format(power),
                    basis_radial_code,
                    orbitals_filename,
                    radial_me_filename
                ]
            )

    # generate radial overlaps
    if (task["basis_mode"] not in {k_radial_direct}):
        b_ratio = math.sqrt(task["hw_int"]/task["hw"])
        mcscript.call(
            [
                configuration.shell_filename("radial-gen"),
                "--overlaps",
                "{:g}".format(b_ratio),
                basis_radial_code,
                orbitals_int_filename,
                radial_olap_int_filename
            ]
        )
    if (task["coulomb"] and (task["basis_mode"] not in {k_radial_direct,k_radial_dilated})):
        b_ratio = math.sqrt(task["hw_coul_prime"]/task["hw"])
        mcscript.call(
            [
                configuration.shell_filename("radial-gen"),
                "--overlaps",
                "{:g}".format(b_ratio),
                basis_radial_code,
                orbitals_coul_filename,
                radial_olap_coul_filename
            ]
        )

## ##################################################################
## # h2mixer 
## ##################################################################
## 
## def task_handler_h2mixer(task):
##     """ Invoke mfdn in h2 run, given current task parameters.
## 
##     Expected dictionary keys: see initial module docstring
##     """
##     print ("h2mixer", task["descriptor"])
## 
##     # recall parameters
##     nuclide = task["nuclide"] 
##     interaction = task["interaction"] 
##     coulomb = task["coulomb"] 
##     hw_int = task["hw_int"] 
##     aLawson = task["aLawson"] 
##     hwLawson = task["hwLawson"]
##     hw = task["hw"] 
##     basis = task["basis"] 
##     scaled_basis = task["scaled_basis"] 
##     xform_cutoff = task["xform_cutoff"] 
##     Nmax = task["Nmax"] 
##     Nv = task["Nv"]
## 
##     # basic setup
##     (Z, N) = nuclide
##     input_truncation = interaction_truncation_for_Nmax(Nv,Nmax)
##     r2k2_basename_unresolved = mcscript.dashify(("r2k2", basis_string_pn(basis), mcscript.dashify(input_truncation)))
##     r2k2_basename = mcscript.subpath_search_basename (ncsm_config.data_dir_h2, r2k2_run_list, r2k2_basename_unresolved)
##     output_truncation = interaction_truncation_for_Nmax(Nv,Nmax,standardize=False)
##     h2mixer_params = {
##         "input_truncation" : input_truncation,
##         "output_truncation" : output_truncation,
##         "Z": Z, "N" : N, "hw" : hw,
##         "r2k2_basename" : r2k2_basename,
##         "input_basenames" : [],
##         "output_streams" : []
##         }
## 
##     # input stream setup
##     # input 5 -- built interaction matrix elements -- for scaled basis
##     tbme_basename_VNN_unresolved = mcscript.dashify((interaction, basis_string_pn(scaled_basis), mcscript.dashify(xform_cutoff), mcscript.dashify(input_truncation), format(hw_int,"g")))
##     h2mixer_params["input_basenames"].append(
##         mcscript.subpath_search_basename (ncsm_config.data_dir_h2, xform_run_list, tbme_basename_VNN_unresolved)
##         )
##     # input 6 (OPTIONAL) -- built Coulomb matrix elements -- for unscaled basis
##     if (coulomb):
##         hw_int_coul = task["hw_int_coul"] 
##         xform_cutoff_coul = task["xform_cutoff_coul"] 
##         tbme_basename_VC_unresolved = mcscript.dashify(("VC", basis_string_pn(basis), mcscript.dashify(xform_cutoff_coul), mcscript.dashify(input_truncation), format(hw_int_coul,"g")))
##         h2mixer_params["input_basenames"].append(
##             mcscript.subpath_search_basename (ncsm_config.data_dir_h2, xform_run_list, tbme_basename_VC_unresolved)
##         )
##         
##     # output stream -- Hamiltonian
##     components = [
##             ("Trel",), # Trel
##             ("in", 5, 1.0), # VNN
##             ("NCM", hwLawson, aLawson) # Lawson
##             ]
##     if (coulomb):
##         components += [
##             ("scaled", 6, hw_int_coul, -1, 1.0) # VC  -- add scaled <in#> <hw_ref> <degree> <multiplier> 
##         ]
##     h2mixer_params["output_streams"].append(
##         {
##             "basename" : "tbme-h2",
##             "components" : components
##         }
##     )
## 
##     # output stream -- rrel2 (required observable)
##     h2mixer_params["output_streams"].append(
##         {
##         "basename" : "tbme-rrel2",
##         "components" : [ ("rrel2",) ]
##         }
##         )
## 
##     # output stream -- NCM(hwLawson)
##     h2mixer_params["output_streams"].append(
##         {
##         "basename" : "tbme-NCM",
##         "components" : [ ("NCM", hwLawson, 1.0) ]
##         }
##         )
## 
##     # output stream -- CM observables R20 & K20
##     if (task.get("obs-R20K20",False)):
##         h2mixer_params["output_streams"].append(
##             {
##                 "basename" : "tbme-R20",
##                 "components" : [ ("R20",) ]
##             }
##         )
##         h2mixer_params["output_streams"].append(
##             {
##                 "basename" : "tbme-K20",
##                 "components" : [ ("K20",) ]
##             }
##         )
## 
##     # output stream -- squared angular momentum observables
##     if (task.get("obs-am-sqr",False)):
##         for op in angular_momentum_operator_list:
##             h2mixer_params["output_streams"].append(
##                 {
##                     "basename" : "tbme-{}2".format(op),
##                     "components" : [ ("am-sqr", op) ]
##                 }
##             )
## 
##     # output streams -- tbme for H components (debugging)
##     if (task.get("obs-H-components",False)):
##         h2mixer_params["output_streams"].append(
##             {
##                 "basename" : "tbme-Trel",
##                 "components" : [ ("Trel",) ]
##             }
##         )
##         h2mixer_params["output_streams"].append(
##             {
##                 "basename" : "tbme-VNN",
##                 "components" : [ ("in", 5, 1.0) ]
##             }
##         )
##         h2mixer_params["output_streams"].append(
##             {
##                 "basename" : "tbme-coul",
##                 "components" : [ ("scaled", 6, hw_int_coul, -1, 1.0) ]
##             }
##         )
## 
## 
##     # run h2mixer
##     call_h2mixer (h2mixer_params)
## 
## ##################################################################
## # mfdn h2 run
## ##################################################################
## 
## def task_handler_mfdn_h2(task):
##     """ Invoke mfdn in h2 run, given current task parameters.
##     
##     Expected dictionary keys: see initial module docstring
##     """
## 
##     print ("MFDn", task["descriptor"])
## 
##     # optional dictionary keys
##     # TODO: neaten up with setdefault method
##     if ("obs-R20K20" not in task):
##         task["obs-R20K20"] = False
##     if ("obs-am-sqr" not in task):
##         task["obs-am-sqr"] = False
##     if ("obs-H-components" not in task):
##         task["obs-H-components"] = False
##     if ("em_multipolarity_list" not in task):
##         task["em_multipolarity_list"] = []
##     if ("keep_obdme" not in task):
##         task["keep_obdme"] = True
## 
##     # set up MFDn basis truncation parameters
##     Nv = task["Nv"]
##     Nmin = task["Nmax"] % task["Nstep"]
##     truncation = interaction_truncation_for_Nmax(Nv,task["Nmax"],standardize=False)
##     ## BUG: through 130520: Nshell = truncation[1] + 1
##     ## BUG: through 150805: Nshell given to MFDn was one step higher than needed for odd Nmax (default standardize="True")
##     Nshell = truncation[0] + 1
##     ## print("truncation",truncation,"Nshell",Nshell)
## 
## 
##     # set up h2 filenames
##     h2_basename = "tbme-h2"
##     obs_basename_list = [ "tbme-rrel2", "tbme-NCM" ]
##     if (task["obs-R20K20"]):
##         obs_basename_list += [ "tbme-R20", "tbme-K20" ]
##     if (task["obs-am-sqr"]):
##         for op in angular_momentum_operator_list:
##             obs_basename_list.append("tbme-{}2".format(op))
##     if (task["obs-H-components"]):
##         obs_basename_list += [ "tbme-Trel", "tbme-VNN", "tbme-coul"]
## 
##     # guard against pathetically common mistakes
##     reference_state_list = task["obdme_reference_state_list"]
##     for (J,g,i) in reference_state_list:
##         twice_J = int(2*J)
##         if ((twice_J%2) != (sum(task["nuclide"])%2)):
##             raise ValueError("invalid angular momentum for reference state")
##         if ((g != (task["Nmax"]%2)) and (task["Nstep"] != 1)):
##             raise ValueError("invalid parity for reference state")
## 
##        
##     # import partitioning file
##     partitioning_suffix = os.environ.get("MFDN_PARTITIONING",None)
##     if (partitioning_suffix is not None):
##         partitioning_filename = os.path.join(
##             ncsm_config.data_dir_partitioning,
##             "mfdn_partitioning.info_Nsh{:02d}_{}".format(Nshell,partitioning_suffix)
##         )
##         print ("Checking for partition file %s..." % partitioning_filename)
##         if (os.path.exists(partitioning_filename)):
##             mcscript.call(["cp", partitioning_filename, "mfdn_partitioning.info"])
##         else:
##             raise mcscript.task.ScriptError("partition file not found")
## 
##     ## partitioning_filename = os.path.join(
##     ##     ncsm_config.data_dir_partitioning,
##     ##     "mfdn_partitioning.info_Nsh{}".format(Nshell)
##     ## )
## 
##     # invoke mfdn
##     task.update(
##         {
##             "Nmin" : Nmin,
##             "Nshell" : Nshell,
##             "h2_basename" : h2_basename,
##             "obs_basename_list" : obs_basename_list
##             }
##         )
##     mfdn_wrapper = task["mfdn_wrapper"]
##     mfdn_wrapper(task)
## 
##     # process results file
##     result_filename = "%s-mfdn-%s.res" % (mcscript.run.name, task["descriptor"])
##     out_filename = "%s-mfdn-%s.out" % (mcscript.run.name, task["descriptor"])
##     # remove detailed occupation listing
##     # TO NEATEN: replace with shell-free subprocess call and pipes
##     ## print ("Filtering mfdn.res...") 
##     ## mcscript.call(["sed '/More/,$d' < mfdn.res > " + result_filename],shell=True) 
##     ## mcscript.call(["cat < mfdn.res > " + result_filename],shell=True) 
##     ## mcscript.call(["cat < mfdn.out > " + out_filename],shell=True) 
##     mcscript.call(["cp","mfdn.res",result_filename]) 
##     mcscript.call(["cp","mfdn.out",out_filename]) 
## 
##     # archive obdme and other data files
##     # including unfiltered mfdn.res
##     archive_filename = "%s-mfdn-%s.tgz" % (mcscript.run.name, task["descriptor"])
##     print ("Archiving output files...")
##     obdme_file_list = ["mfdn.out", "mfdn.dat", "mfdn.res", "mfdn_spstates.info"]
##     if (task["keep_obdme"]):
##         obdme_file_list += glob.glob("*obdme*") 
##     print ("obdme files:", obdme_file_list)
##     mcscript.call(["tar", "zcvf", archive_filename] + obdme_file_list)
## 
##     # move results out
##     ##mcscript.call(["mv run* -t "+mcscript.task.results_dir],shell=True)
##     mcscript.call(["mv", result_filename, out_filename, archive_filename, "--target-directory="+mcscript.task.results_dir])
## 
##     # cleanup of wave function files
##     scratch_file_list = glob.glob("mfdn_smwf*") + glob.glob("mfdn_MBgroups0*")
##     mcscript.call(["rm", "-vf"] + scratch_file_list)
           
    
if (__name__ == "__MAIN__"):
    pass
