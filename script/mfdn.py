"""mfdn.py -- define scripting for h2mixer+MFDn runs

    Expected dictionary keys for h2mixer + mfdn run:

        # nuclide parameters
        "nuclide" (tuple of int): (Z,N) tuple
        
        # Hamiltonian parameters
        "interaction" (string): name for interaction (for use in descriptor and filenames)
        "use_coulomb" (bool): whether or not to include Coulomb
        "a_cm" (float): coefficient of N_cm for Lawson term
        "hw_cm" (float): hw of N_cm for Lawson term (None: use hw of basis)

        # input TBME parameters
        "truncation_int" (tuple): input interaction TBME cutoff, as tuple ("ob"|"tb",N) 
        "hw_int" (float): hw of basis for source interaction TBMEs
        "truncation_coul" (tuple): input Coulomb TBME cutoff, as tuple ("ob"|"tb",N) 
        "hw_coul" (float): hw of basis for source Coulomb TBMEs

        # basis parameters
        "basis_mode" (int): enumerated value indicating direct
            oscillator (k_basis_mode_direct), dilated oscillator
            (k_basis_mode_dilated), or generic run mode
            (k_basis_mode_generic), as explained further where these
            constants are defined below; beware that natural orbital
            runs must be treated as generic even if the initial
            (starting basis) run is of oscillator type
        "hw" (float): hw of basis

        # transformation parameters
        "xform_truncation_int" (tuple): transform cutoff for interaction, as tuple ("ob"|"tb",N) 
        "xform_truncation_coul" (tuple): transform cutoff for Coulomb, as tuple ("ob"|"tb",N)
        "hw_coul_rescaled" (float): hw to which to rescale Coulomb TBMEs before two-body
            transformation (None: use hw of basis)
             - direct oscillator run: naturally hw to avoid any two-body transformation
             - dilated oscillator run: naturally hw to avoid any two-body transformation, but one could
               also prefer hw_int for uniformity in the two-body transformation (at the expense of
               introducing some transformation error in the Coulomb interaction)
             - generic run: naturally hw_int for uniformity in the two-body transformation
        "target_truncation" (tuple): truncation of target TBMEs, as weight_max tuple
             (None: deduce automatically from valence shell and many-body truncation information)

        # legacy -- to adapt or remove
        "basis" (tuple): logical basis (radial_basis_p,1.,radial_basis_n,beta_n), to be scaled by hw
        "scaled_basis" (tuple): deduced basis used in task descriptor, etc.

        # traditional oscillator truncation
        "ho_truncation" (bool): whether or not to assume traditional HO-style truncation
            - weighting is by N=2n+l
            - two-body truncation is by a one-body or two-body N cutoff (compatible with MFDn
              version 14 and h2 format 0), rather than any more general (Wp,Wn,Wpp,Wnn,Wpn)
              cutoff (which requires MFDn version 15 and h2 format 15099)
            - many-body truncation is by a one-body (FCI) or many-body (Nmax) N cutoff
        "many_body_truncation" (str): many-body truncation rank ("FCI" for one-body, or
            "Nmax" for A-body)
        "Nv" (int): N of valence shell (for use in truncation)
        "Nmax" (int): many-body oscillator cutoff (interpreted as one-body Nmax_orb for "FCI"
             truncation, or many-body excitation cutoff Nmax for "Nmax" truncation)
        "Nstep" (int): Nstep (2 for single parity, 1 for mixed parity)
 
        # diagonalization parameters
        "Mj" (float): M-scheme angular momentum projection value (as true value,
            possibly half-integral, *not* doubled!)
        "eigenvectors" (int): number of eigenvectors to calculate
        "initial_vector" (int): initial vector code for mfdn
        "lanczos" (int): lanczos iterations
        "tolerance" (float): diagonalization tolerance parameter
        "partition_filename" (str): filename for partition file to use with MFDn (None: no
            partition file); for now absolute path is required, but path search protocol may
            be restored in future

        # obdme parameters
        "obdme_multipolarity" (int): maximum multipolarity for calculation of densities
        "obdme_reference_state_list" (list): list of reference states (J,g,i) for density calculation, or "all2all"
        "save_obdme" (bool): whether or not to save obdme files in archive

        # two-body observables
        "observable_sets" (list of str): codes for predefined observable sets to include
            "H-components": Hamiltonian terms
            "am-sqr": squared angular momenta
            "R20K20": center-of-mass diagnostic observables (TODO)
        
        # version parameters
        "h2_format" (int): h2 file format to use (values include: 0, 15099)
        "mfdn_executable" (string): mfdn executable name

        # emcalc postprocessing parameters
        "em_multipolarity_list" (list, optional): list of multipolarities ("E"|"M",lambda) 

        # natural orbital postprocessing
        "enable_natural_orbitals" (bool): whether or not to generate natural orbital auxiliary files
          
  Mark A. Caprio
  University of Notre Dame

  - 12/14/16 (mac): Created, drawing on ncsm.py (created 2/12/13) and
    shell package example code generate_input.py (created 11/7/16).
  - 12/27/16 (mac): Rough in scripting of MFDn v14 run.
  - 12/29/16 (mac): Complete basic scripting for MFDn v14 oscillator-type runs.

"""
  
import datetime
import glob
import math
import os
import shutil

import mcscript

################################################################
# radial basis modes
################################################################

# General modes of operation for radial basis
#
# k_basis_mode_direct:
#   - final basis is oscillator basis (hw)
#   - source basis for VNN is oscillator basis of same
#     oscillator length (hw_int=hw); therefore no transformation
#     needed on VNN TBMEs
#   - Coulomb TBMEs need only scaling for dilation (hw_c -> hw)
#   - MFDn can use built-in oscillator OBMEs for observables
#
# k_basis_mode_dilated:
#   - final basis is oscillator basis (hw)
#   - source basis for VNN is oscillator basis of different
#     oscillator length; therefore transformation
#     needed on VNN TBMEs (hw_int -> hw)
#   - Coulomb TBMEs need only scaling for dilation (hw_c -> hw)
#   - MFDn can use built-in oscillator OBMEs for observables
#
# k_basis_mode_generic:
#   - final basis is not assumed to be oscillator basis
#     (still has nominal hw to define a length parameter)
#   - transformation needed on VNN TBMEs (hw_int HO -> hw generic)
#   - Coulomb TBMEs may be rescaled (hw_c -> hw_cp) but then need
#     transformation (hw_cp HO -> hw generic)
#   - MFDn *cannot* use built-in oscillator OBMEs for observables

k_basis_mode_direct = 0
k_basis_mode_dilated = 1
k_basis_mode_generic = 2

################################################################
# physical constants
################################################################

k_mN_csqr = 938.92  # (m_N c^2)~938.92 MeV
k_hbar_c = 197.327  # (hbar c)~197.327 Mev fm
    

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

    Data members:

        install_dir (str): installation directory for shell project ("SHELL_INSTALL_DIR")
        mfdn_dir (str): base directory for finding MFDn executables ("SHELL_MFDN_DIR")
            EX: base directory "${HOME}/projects/mfdn" may contain "mfdn-v14-beta06-newmake/xmfdn-h2-lan"
            and "beta00/xmfdn-h2-lan"
        data_dir_h2_list (str): base directories for interaction tbme files ("SHELL_DATA_DIR_H2")
            Environment variable is interpreted as a PATH-style colon-delimited list.
        interaction_run_list (list of str): subdirectories for interaction tbme files
            (to be set by calling run script)

    Methods:

        shell_filename(name): Generate full qualified filename for
            shell utility executable, given root name.

        mfdn_filename(name): Generate full qualified filename for
            MFDn executable, given root name.

        interaction_filename(name): Generate full qualified filename
            for interaction h2 file, given root name.

    """
    def __init__(self):

        # environment
        self.install_dir = os.environ.get("SHELL_INSTALL_DIR")
        self.mfdn_dir = os.environ.get("SHELL_MFDN_DIR")
        self.data_dir_h2_list = os.environ.get("SHELL_DATA_DIR_H2").split(":")
        self.interaction_run_list = []

    def shell_filename(self,name):
        """Construct filename for shell package executable."""
        return os.path.join(self.install_dir,"bin",name)

    def mfdn_filename(self,name):
        """Construct filename for MFDn executable."""
        return os.path.join(self.mfdn_dir,name)

    def interaction_filename(self,name):
        """Construct filename for interaction h2 file."""
        ##return os.path.join(self.data_dir_h2,name)
        return mcscript.utils.search_in_subdirectories(self.data_dir_h2_list,self.interaction_run_list,name)


configuration = Configuration()

################################################################
# task descriptor for h2mixer + mfdn run
################################################################


def task_descriptor_7(task):
    """ Task descriptor format 7

        Overhaul for new h2utils scripting:
        - Strip back down to basic form for oscillator-like runs only.
        - Adjust some field labels.
        - Add tolerance.
    """

    if(
        task["ho_truncation"]
        and
        task["basis_mode"] in {k_basis_mode_direct,k_basis_mode_dilated}
    ):
        # traditional oscillator run
        template_string = (
            "Z{nuclide[0]}-N{nuclide[1]}-{interaction}-coul{coulomb_flag:d}"
            "-hw{hw:.3f}"
            "-a_cm{a_cm:g}"
            "-Nmax{Nmax:02d}{mixed_parity_indicator}{fci_indicator}-Mj{Mj:03.1f}"
            "-lan{lanczos:d}-tol{tolerance:.1e}"
            "{natural_orbital_indicator}"
            )
    else:
        raise mcscript.ScriptError("mode not supported by task descriptor")

    if (task["many_body_truncation"]=="fci"):
        fci_indicator = "fci"
    else:
        fci_indicator = ""
    mixed_parity_indicator = mcscript.utils.ifelse(task["Nstep"]==1,"x","")
    coulomb_flag = int(task["use_coulomb"])

    natural_orbital_iteration = task.get("natural_orbital_iteration")
    if (natural_orbital_iteration is None):
        natural_orbital_indicator = ""
    else:
        natural_orbital_indicator = "-no{:1d}".format(natural_orbital_iteration)

    descriptor = template_string.format(
        coulomb_flag=coulomb_flag,
        mixed_parity_indicator=mixed_parity_indicator,
        fci_indicator=fci_indicator,
        natural_orbital_indicator=natural_orbital_indicator,
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

def set_up_orbitals(task):
    """Set up source and target orbitals for MFDn run.

    Arguments:
        task (dict): as described in module docstring

    Limitation: Currently only supports harmonic oscillator style
    truncation.
    """

    # orbital filenames
    orbitals_int_filename = "orbitals-int.dat"  # orbitals for interaction tbme basis
    orbitals_coul_filename = "orbitals-coul.dat"  # orbitals for Coulomb tbme basis
    orbitals_filename = "orbitals.dat"  # orbitals for target basis

    # validate basis mode
    if (not task["ho_truncation"]):
        raise ValueError("expecting ho_truncation to be True but found {ho_truncation}".format(**task))

    # generate orbitals -- interaction bases
    mcscript.call(
        [
            configuration.shell_filename("orbital-gen"),
            "--oscillator",
            "{truncation_int[1]:d}".format(**task),
            "{:s}".format(orbitals_int_filename)
        ]
    )
    if (task["use_coulomb"]):
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
        if (task["many_body_truncation"]=="Nmax"):
            Nmax_orb = task["Nmax"] + task["Nv"]
        elif (task["many_body_truncation"]=="FCI"):
            Nmax_orb = task["Nmax"]
        mcscript.call(
            [
                configuration.shell_filename("orbital-gen"),
                "--oscillator",
                "{Nmax_orb:d}".format(Nmax_orb=Nmax_orb),
                "{:s}".format(orbitals_filename)
            ]
        )

def set_up_radial_analytic(task):
    """Generate radial integrals and overlaps by integration for MFDn run
    in analytic basis.

    Operation mode may in general be direct oscillator, dilated
    oscillator, or generic (TODO).

    Arguments:
        task (dict): as described in module docstring

    """

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
    if (task["basis_mode"] not in {k_basis_mode_direct,k_basis_mode_dilated}):  # no k_basis_mode_generic yet
        raise ValueError("invalid basis mode {basis_mode}".format(**task))
    
    # basis radial code -- expected by radial_utils codes
    basis_radial_code = "oscillator"  # TO GENERALIZE: if not oscillator basis

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
                ],
                mode = mcscript.call.serial  # TODO: upgrade to SMP
            )

    # generate radial overlaps
    if (task["basis_mode"] not in {k_basis_mode_direct}):
        b_ratio = math.sqrt(task["hw_int"]/task["hw"])
        mcscript.call(
            [
                configuration.shell_filename("radial-gen"),
                "--overlaps",
                "{:g}".format(b_ratio),
                basis_radial_code,
                orbitals_int_filename,
                radial_olap_int_filename
            ],
            mode = mcscript.call.serial  # TODO: upgrade to SMP
        )
    if (task["use_coulomb"] and (task["basis_mode"] not in {k_basis_mode_direct,k_basis_mode_dilated})):
        b_ratio = math.sqrt(task["hw_coul_rescaled"]/task["hw"])
        mcscript.call(
            [
                configuration.shell_filename("radial-gen"),
                "--overlaps",
                "{:g}".format(b_ratio),
                basis_radial_code,
                orbitals_coul_filename,
                radial_olap_coul_filename
            ],
            mode = mcscript.call.serial  # TODO: upgrade to SMP
        )

def generate_tbme(task):
    """Generate TBMEs for MFDn run.

    Operation mode may in general be direct oscillator, dilated
    oscillator, or generic (TODO).

    Arguments:
        task (dict): as described in module docstring

    """

    # validate basis mode
    if (not task["ho_truncation"]):
        raise ValueError("expecting ho_truncation to be True but found {ho_truncation}".format(**task))

    # extract parameters for convenience
    A = sum(task["nuclide"])
    a_cm = task["a_cm"]
    hw = task["hw"]
    hw_cm = task["hw_cm"]
    if (hw_cm is None):
        hw_cm = hw
    hw_coul = task["hw_coul"]
    hw_coul_rescaled = task["hw_coul_rescaled"]
    if (hw_coul_rescaled is None):
        hw_coul_rescaled = hw

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

    # accumulate h2mixer input lines
    lines = []

    # initial comment
    lines.append("# task: {}".format(task))
    lines.append("")

    # global mode definitions
    target_truncation = task["target_truncation"]
    if (target_truncation is None):
        # automatic derivation
        if (task["ho_truncation"]):
            if (task["many_body_truncation"]=="Nmax"):
                N2_max = 2*task["Nv"]+task["Nmax"]
                target_weight_max = weight_max_string(("tb",N2_max))
            elif (task["many_body_truncation"]=="FCI"):
                N1_max = task["Nmax"]
                target_weight_max = weight_max_string(("ob",N1_max))
        else:
            # calculation of required weight_max will require external program for occupation counting
            pass
    else:
        # given value
        target_weight_max = target_truncation
    lines.append("set-target-indexing {orbitals_filename} {target_weight_max}".format(
        orbitals_filename=orbitals_filename,
        target_weight_max=target_weight_max,
        **task
    ))
    lines.append("set-target-multipolarity 0 0 0")
    lines.append("set-output-format {h2_format}".format(**task))
    lines.append("set-mass {A}".format(A=A,**task))
    lines.append("")

    # radial operator inputs
    for operator_type in ["r","k"]:
        for power in [1,2]:
            radial_me_filename = radial_me_filename_template.format(operator_type,power)
            lines.append("define-radial-operator {} {} {}".format(operator_type,power,radial_me_filename))
    lines.append("")
    
    # sources: identity and kinematic operators
    lines.append("define-source operator identity")
    lines.append("define-source operator Ursqr")
    lines.append("define-source operator Vr1r2")
    lines.append("define-source operator Uksqr")
    lines.append("define-source operator Vk1k2")
    lines.append("")

    # sources: VNN
    VNN_filename = configuration.interaction_filename(
        "{}-{}-{:g}.bin".format(
            task["interaction"],
            mcscript.utils.dashify(task["truncation_int"]),
            task["hw_int"]
        )
    )
    if (task["basis_mode"]==k_basis_mode_direct):
        lines.append("define-source input VNN {VNN_filename}".format(VNN_filename=VNN_filename,**task))
    else:
        xform_weight_max_int = weight_max_string(task["xform_truncation_int"])
        lines.append("define-source xform VNN {VNN_filename} {xform_weight_max_int} {radial_olap_int_filename}".format(
            VNN_filename=VNN_filename,
            xform_weight_max_int=xform_weight_max_int,
            radial_olap_int_filename=radial_olap_int_filename,
            **task
        ))

    # sources: Coulomb
    #
    # Note: This is the "unscaled" Coulomb, still awaiting the scaling
    # factor from dilation.
    if (task["use_coulomb"]):
        VC_filename = configuration.interaction_filename(
            "{}-{}-{:g}.bin".format(
                "VC",
                mcscript.utils.dashify(task["truncation_coul"]),
                task["hw_coul"]
            )
        )
        if (task["basis_mode"] in {k_basis_mode_direct,k_basis_mode_dilated}):
            lines.append("define-source input VC_unscaled {VC_filename}".format(VC_filename=VC_filename,**task))
        else:
            xform_weight_max_coul = weight_max_string(task["xform_truncation_coul"])
            lines.append("define-source xform VC_unscaled {VC_filename} {xform_weight_max_coul} {radial_olap_coul_filename}".format(
                VC_filename=VC_filename,
                xform_weight_max_coul=xform_weight_max_coul,
                radial_olap_coul_filename=radial_olap_coul_filename,
                **task
            ))

    lines.append("")


    # target: Hamiltonian
    coef_Uksqr = (A-1)/(2*A)*hw+1/(2*A)*a_cm*(hw/hw_cm)
    coef_Vk1k2 = -1/A*hw+1/A*a_cm*(hw/hw_cm)
    coef_Ursqr = 1/(2*A)*a_cm*(hw_cm/hw)
    coef_Vr1r2 = 1/A*a_cm*(hw_cm/hw)
    coef_identity = -3/2*a_cm
    lines.append("define-target tbme-H.bin")
    lines.append("  add-source Ursqr {:e}".format(coef_Ursqr))
    lines.append("  add-source Vr1r2 {:e}".format(coef_Vr1r2))
    lines.append("  add-source Uksqr {:e}".format(coef_Uksqr))
    lines.append("  add-source Vk1k2 {:e}".format(coef_Vk1k2))
    lines.append("  add-source identity {:e}".format(coef_identity))
    lines.append("  add-source VNN {}".format(1.))
    if (task["use_coulomb"]):
        coef_VC = math.sqrt(hw_coul_rescaled/hw_coul)
        lines.append("  add-source VC_unscaled {:e}".format(coef_VC))

    # target: radius squared
    coef_Ursqr = (A-1)*(oscillator_length(hw)/A)**2
    coef_Vr1r2 = -2*(oscillator_length(hw)/A)**2
    lines.append("define-target tbme-rrel2.bin")
    lines.append("  add-source Ursqr {:e}".format(coef_Ursqr))
    lines.append("  add-source Vr1r2 {:e}".format(coef_Vr1r2))

    # target: NCM
    #
    # These are the "a_cm" terms from the Hamiltonian (as above), but
    # without including the factor of "a_cm".
    coef_Uksqr = 1/(2*A)*(hw/hw_cm)
    coef_Vk1k2 = 1/A*(hw/hw_cm)
    coef_Ursqr = 1/(2*A)*(hw_cm/hw)
    coef_Vr1r2 = 1/A*(hw_cm/hw)
    coef_identity = -3/2
    lines.append("define-target tbme-Ncm.bin")
    lines.append("  add-source Ursqr {:e}".format(coef_Ursqr))
    lines.append("  add-source Vr1r2 {:e}".format(coef_Vr1r2))
    lines.append("  add-source Uksqr {:e}".format(coef_Uksqr))
    lines.append("  add-source Vk1k2 {:e}".format(coef_Vk1k2))
    lines.append("  add-source identity {:e}".format(coef_identity))

    # optional observable sets

    # Hamiltonian components
    if ("H-components" in task["observable_sets"]):

        # target: Trel (diagnostic)
        #
        # These are the "non-a_cm" terms from the Hamiltonian (as above).
        coef_Uksqr = (A-1)/(2*A)*hw
        coef_Vk1k2 = -1/A*hw
        lines.append("define-target tbme-Trel.bin")
        lines.append("  add-source Uksqr {:e}".format(coef_Uksqr))
        lines.append("  add-source Vk1k2 {:e}".format(coef_Vk1k2))

        # target: VNN (diagnostic)
        lines.append("define-target tbme-VNN.bin")
        lines.append("  add-source VNN {}".format(1.))

        # target: VC (diagnostic)
        if (task["use_coulomb"]):
            lines.append("define-target tbme-VC.bin")
            coef_VC = math.sqrt(hw_coul_rescaled/hw_coul)
            lines.append("  add-source VC_unscaled {:e}".format(coef_VC))

    # squared angular momenta
    if ("am-sqr" in task["observable_sets"]):

        # sources: squared angular momentum operators
        lines.append("define-source operator L")
        lines.append("define-source operator Sp")
        lines.append("define-source operator Sn")
        lines.append("define-source operator S")
        lines.append("define-source operator J")

        # targets: squared angular momentum operators
        lines.append("define-target tbme-L.bin")
        lines.append("  add-source L {:e}".format(1))
        lines.append("define-target tbme-Sp.bin")
        lines.append("  add-source Sp {:e}".format(1))
        lines.append("define-target tbme-Sn.bin")
        lines.append("  add-source Sn {:e}".format(1))
        lines.append("define-target tbme-S.bin")
        lines.append("  add-source S {:e}".format(1))
        lines.append("define-target tbme-J.bin")
        lines.append("  add-source J {:e}".format(1))

    # ensure terminal line
    lines.append("")
    
    # diagnostic: log input lines to file
    #
    # This is purely for easy diagnostic purposes, since lines will be
    # fed directly to h2mixer as stdin below.
    mcscript.utils.write_input("h2mixer.in",input_lines=lines,verbose=False)

    # invoke h2mixer
    mcscript.call(
        [
            configuration.shell_filename("h2mixer")
        ],
        input_lines=lines,        
        mode = "serial"  # SMP only
    )

def run_mfdn_v14_b06(task):
    """Generate input file and execute MFDn version 14 beta 06.

    Arguments:
        task (dict): as described in module docstring
    
    Raises:
        mcscript.ScriptError: if MFDn output not found

    """

    # validate basis mode
    if (not task["ho_truncation"]):
        raise ValueError("expecting ho_truncation to be True but found {ho_truncation}".format(**task))

    # accumulate MFDn input lines
    lines = []

    # base parameters
    twice_Mj = int(2*task["Mj"])
    if (task["many_body_truncation"]=="Nmax"):
        Nmax_orb = task["Nmax"] + task["Nv"]
    elif (task["many_body_truncation"]=="FCI"):
        Nmax_orb = task["Nmax"]
    Nshell = Nmax_orb+1
    if (task["basis_mode"] in {k_basis_mode_direct,k_basis_mode_dilated}):  
        hw_for_trans = task["hw"]
    else:
        hw_for_trans = 0  # disable MFDn hard-coded oscillator one-body observables
    ## ndiag = int(os.environ.get("MFDN_NDIAG",0))  # allows override for spares, not so elegant
    ndiag = 0
    if (task["Nstep"]==2):
        Nmin = task["Nmax"]%2
    else:
        Nmin = 1

    lines.append("{:d}  # IFLAGMBSI".format(0))
    lines.append("{ndiag:d}  # ndiag (0: no spares, automatic ndiag)".format(ndiag=ndiag,**task))
    lines.append("{:d}  # number of classes".format(2))
    lines.append("{nuclide[0]:d}  # protons (class 1 particles)".format(**task))
    lines.append("{nuclide[1]:d}  # protons (class 2 particles)".format(**task))
    lines.append("1 {Nshell:d}  # min, max # S.P. shells for class 1 particles".format(Nshell=Nshell,**task))
    lines.append("1 {Nshell:d}  # min, max # S.P. shells for class 2 particles".format(Nshell=Nshell,**task))
    lines.append("{Nmin:d} {Nmax:d} {Nstep:d}  # N_min, N_max, delta_N".format(Nmin=Nmin,**task))
    lines.append("{:d}   # Total 2 M_j".format(twice_Mj))
    lines.append("{eigenvectors:d} {lanczos:d} {initial_vector:d} {tolerance:e}  # number of eigenvalues/vectors, max number of its, ...)".format(**task))
    lines.append("{:d} {:d}  # rank of input Hamiltonian/interaction".format(2,2))
    lines.append("{hw_for_trans:g} {k_mN_csqr:g}  # h-bar*omega, Nucleon mass (MeV) ".format(
        hw_for_trans=hw_for_trans,k_mN_csqr=k_mN_csqr,**task
    ))

    # tbo: collect tbo names
    obs_basename_list = ["tbme-rrel2","tbme-Ncm"]
    if ("H-components" in task["observable_sets"]):
        obs_basename_list += ["tbme-Trel","tbme-VNN"]
        if (task["use_coulomb"]):
            obs_basename_list += ["tbme-VC"]
    if ("am-sqr" in task["observable_sets"]):
        obs_basename_list += ["tbme-L","tbme-Sp","tbme-Sn","tbme-S","tbme-J"]

    # tbo: log tbo names in separate file to aid future data analysis
    mcscript.utils.write_input("tbo_names.dat",input_lines=obs_basename_list)

    # tbo: write list of operators
    lines.append("tbme-H")  # Hamiltonian basename
    num_obs = 2 + len(obs_basename_list)
    lines.append("{:d}   # number of observables (J, T, R2, ...)".format(num_obs))
    lines += obs_basename_list
    
    # obdme: parameters
    lines.append("{enable_obd:d} {twice_multipolarity:d} # static one-body density matrix elements (0: no one-body densities), twice multipolarity".format(
        enable_obd=1,twice_multipolarity=2*task["obdme_multipolarity"]
    ))
    lines.append("{num_reference_states:d} {max_delta_J:d} # number of reference states for transitions (0: no transitions, -1: all2all), max delta2J (?)".format(
        num_reference_states=len(task["obdme_reference_state_list"]),
        max_delta_J=2*task["obdme_multipolarity"]
    ))

    # obdme: validate reference state list
    #
    # guard against pathetically common mistakes
    for (J,g_rel,i) in task["obdme_reference_state_list"]:
        # validate integer/half-integer character of angular momentum
        twice_J = int(2*J)
        if ((twice_J%2) != (sum(task["nuclide"])%2)):
            raise ValueError("invalid angular momentum for reference state")
        # validate grade (here taken relative to natural grade)
        if ((g_rel != (task["Nmax"]%2)) and (task["Nstep"] != 1)):
            raise ValueError("invalid parity for reference state")

    # obdme: write reference state list
    for reference_state in task["obdme_reference_state_list"]:
        lines.append("{:d} {:d} {:d}".format(2*reference_state[0],reference_state[1],reference_state[2]))

    # ensure terminal line
    lines.append("")

    # generate MFDn input file
    mcscript.utils.write_input("mfdn.dat",input_lines=lines)
    
    # import partitioning file
    if (task["partition_filename"] is not None):
        if (not os.path.exists(task["partition_filename"])):
            raise mcscript.ScriptError("partition file not found")
        mcscript.call(["cp","--verbose",task["partition_filename"],"mfdn_partitioning.info"])

    # invoke MFDn
    mcscript.call(
        [
            configuration.mfdn_filename(task["mfdn_executable"])
        ],
        mode = mcscript.call.hybrid,
        check_return=True
    )

    # test for basic indications of success
    if (not os.path.exists("mfdn.out")):
        raise mcscript.ScriptError("mfdn.out not found")
    if (not os.path.exists("mfdn.res")):
        raise mcscript.ScriptError("mfdn.res not found")

def save_mfdn_output(task):
    """Generate input file and execute MFDn version 14 beta 06.

    Arguments:
        task (dict): as described in module docstring
    
    Raises:
        mcscript.ScriptError: if MFDn output not found

    """

    # save quick inspection copies of mfdn.{res,out}
    print("Saving basic output files...")
    res_filename = "{:s}-mfdn-{:s}.res".format(mcscript.run.name,task["descriptor"])
    mcscript.call(["cp","--verbose","mfdn.res",res_filename]) 
    out_filename = "{:s}-mfdn-{:s}.out".format(mcscript.run.name,task["descriptor"])
    mcscript.call(["cp","--verbose","mfdn.out",out_filename]) 

    # save full archive of input, log, and output files
    print("Saving full output files...")
    # logging
    archive_file_list = [
        "h2mixer.in","tbo_names.dat"
        ]
    # MFDn output
    archive_file_list += [
        "mfdn.dat","mfdn.out","mfdn.res","mfdn_partitioning.generated","mfdn_spstates.info"
    ]
    # renamed versions
    archive_file_list += [out_filename,res_filename]
    # MFDN obdme
    if (task["save_obdme"]):
        archive_file_list += glob.glob("*obdme*")
    # generate archive
    archive_filename = "{:s}-mfdn-{:s}.tgz".format(mcscript.run.name,task["descriptor"])
    mcscript.call(["tar", "zcvf", archive_filename] + archive_file_list)

    # copy results out (if in multi-task run)
    if (mcscript.task.results_dir is not None):
        mcscript.call(
            [
                "cp",
                "--verbose",
                res_filename,out_filename,archive_filename,
                "--target-directory={}".format(mcscript.task.results_dir)
            ]
        )

    # cleanup of wave function files
    scratch_file_list = glob.glob("mfdn_smwf*") + glob.glob("mfdn_MBgroups0*")
    mcscript.call(["rm", "-vf"] + scratch_file_list)
           
    
if (__name__ == "__MAIN__"):
    pass
