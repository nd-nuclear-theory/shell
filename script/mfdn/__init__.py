"""mfdn -- define scripting for h2mixer+MFDn runs

    Expected dictionary keys for h2mixer + mfdn run:

        # nuclide parameters
        "nuclide" (tuple of int): (Z,N) tuple

        # Hamiltonian parameters
        "interaction" (string): name for interaction (for use in descriptor and filenames)
        "use_coulomb" (bool): whether or not to include Coulomb
        "a_cm" (float): coefficient of N_cm for Lawson term
        "hw_cm" (float): hw of N_cm for Lawson term (None: use hw of basis)
        "hamiltonian" (CoefficientDict, optional): specification of Hamiltonian as a
            CoefficientDict of two-body operators passed as sources to h2mixer (see mfdn.operators)

        # input TBME parameters
        "truncation_int" (tuple): input interaction TBME cutoff, as tuple ("ob"|"tb",N)
        "hw_int" (float): hw of basis for source interaction TBMEs
        "truncation_coul" (tuple): input Coulomb TBME cutoff, as tuple ("ob"|"tb",N)
        "hw_coul" (float): hw of basis for source Coulomb TBMEs

        # basis parameters
        "basis_mode" (BasisMode): enumerated value indicating direct
            oscillator (k_basis_mode_direct), dilated oscillator
            (k_basis_mode_dilated), or generic run mode
            (k_basis_mode_generic), as explained further where this
            enum is defined below.
        "hw" (float): hw of basis

        # transformation parameters
        "xform_truncation_int" (tuple, optional): transform cutoff for interaction, as tuple ("ob"|"tb",N)
        "xform_truncation_coul" (tuple, optional): transform cutoff for Coulomb, as tuple ("ob"|"tb",N)
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

        # truncation parameters
        "truncation_mode" (TruncationMode): enumerated value indicating oscillator-like
            truncation or triangular truncation; see docstring of TruncationMode for
            information.
        "truncation_parameters" (dict): truncation parameters, specific to each enumerated
            truncation type; see docstring of TruncationMode for full documentation

        # diagonalization parameters
        "Mj" (float): M-scheme angular momentum projection value (as true value,
            possibly half-integral, *not* doubled!)
        "eigenvectors" (int): number of eigenvectors to calculate
        "initial_vector" (int): initial vector code for mfdn
        "lanczos" (int): lanczos iterations
        "tolerance" (float): diagonalization tolerance parameter
        "ndiag" (int): number of spare diagonal nodes
        "partition_filename" (str): filename for partition file to use with MFDn (None: no
            partition file); for now absolute path is required, but path search protocol may
            be restored in future

        # obdme parameters
        "obdme_multipolarity" (int): maximum multipolarity for calculation of densities
        "obdme_reference_state_list" (list): list of reference states (J,g,i) for density calculation, or "all2all"
        "save_obdme" (bool): whether or not to save obdme files in archive

        # two-body observables
        "observables" (dict of {"basename": CoefficientDict}): additional observable definitions
            (see mfdn.operators)
        "observable_sets" (list of str): codes for predefined observable sets to include
            "H-components": Hamiltonian terms
            "am-sqr": squared angular momenta
            "R20K20": center-of-mass diagnostic observables (TODO)

        # version parameters
        "h2_format" (int): h2 file format to use (values include: 0, 15099)
        "mfdn_executable" (string): mfdn executable name

        # natural orbital parameters
        "natorb_iteration" (int): 0-based iteration number
        "natorb_base_state" (int): MFDn sequence number of state off which to
            build natural orbitals

        # emcalc postprocessing parameters
        "em_multipolarity_list" (list, optional): list of multipolarities ("E"|"M",lambda)

  Mark A. Caprio
  University of Notre Dame

  - 12/14/16 (mac): Created, drawing on ncsm.py (created 2/12/13) and
    shell package example code generate_input.py (created 11/7/16).
  - 12/27/16 (mac): Rough in scripting of MFDn v14 run.
  - 12/29/16 (mac): Complete basic scripting for MFDn v14 oscillator-type runs.
  - 1/22-28/17 (pjf): Implement iterated natural orbitals.
    + Turn k_basis_mode_* into Python enum BasisMode
    + Add TruncationMode enumeration
    + Add FilenameConfiguration to hold filename patterns
    + Updated docs to more clearly differentiate between basis modes and
      truncation modes. Basis modes only control the 0th natural orbital run,
      since the nth (n>0) run is always "generic". Truncation modes control
      which states are included in the single and two-body bases.
    + MFDn is now run in a temporary work/ subdirectory. This ensures that MFDn
      can't get confused by earlier natural orbital runs.
    + Rename save_mfdn_output -> save_mfdn_v14_output.
  - 1/30/17 (mac): Downgrade syntax to python 3.4.
  - 1/31/17 (mac): Fix one-body truncation on Hamiltonian tbme files.
  - 2/3/17 (pjf):
    + Make "xform_truncation_int" and "xform_truncation_coul" optional.
    + Fix save_mfdn_v14_output() when Coulomb is turned off.
    + Fix natural orbital indicator in task_descriptor_7.
  - 2/10/17 (pjf):
    + Fix "fci" -> "-fci" flag
    + Add ndiag parameter to task dictionary.
  - 2/20/17 (pjf): Rename mfdn.py -> mfdn/__init__.py
"""

import datetime
import glob
import math
import os
import shutil
import enum

import mcscript

from .utils import *
import mfdn.operators

################################################################
# radial basis modes
################################################################


@enum.unique
class BasisMode(enum.Enum):
    """General modes of operation for radial basis

    kDirect (alias k_basis_mode_direct):
      - no0 basis is oscillator basis (hw)
      - source basis for VNN is oscillator basis of same
        oscillator length (hw_int=hw); therefore no transformation
        needed on VNN TBMEs
      - Coulomb TBMEs need only scaling for dilation (hw_c -> hw)
      - MFDn can use built-in oscillator OBMEs for observables

    kDilated (alias k_basis_mode_dilated):
      - no0 basis is oscillator basis (hw)
      - source basis for VNN is oscillator basis of different
        oscillator length; therefore transformation
        needed on VNN TBMEs (hw_int -> hw)
      - Coulomb TBMEs need only scaling for dilation (hw_c -> hw)
      - MFDn can use built-in oscillator OBMEs for observables

    kGeneric (alias k_basis_mode_generic):
      - no0 basis is not assumed to be oscillator basis
        (still has nominal hw to define a length parameter)
      - transformation needed on VNN TBMEs (hw_int HO -> hw generic)
      - Coulomb TBMEs may be rescaled (hw_c -> hw_cp) but then need
        transformation (hw_cp HO -> hw generic)
      - MFDn *cannot* use built-in oscillator OBMEs for observables

    """
    kDirect = 0
    kDilated = 1
    kGeneric = 2

k_basis_mode_direct = BasisMode.kDirect
k_basis_mode_dilated = BasisMode.kDilated
k_basis_mode_generic = BasisMode.kGeneric

################################################################
# truncation modes
################################################################

@enum.unique
class TruncationMode(enum.Enum):
    """General truncation modes for radial basis

    k_truncation_mode_ho:
        - traditional Nmax truncation; weight is (2*n + l)
        - compatible with MFDn v14+
        - "truncation_parameters" (dict):
            - "many_body_truncation" (str): many-body truncation rank ("FCI" for one-body, or
                "Nmax" for A-body)
            - "Nv" (int): N of valence shell (for use in truncation)
            - "Nmax" (int): many-body oscillator cutoff (interpreted as one-body Nmax_orb for "FCI"
                 truncation, or many-body excitation cutoff Nmax for "Nmax" truncation)
            - "Nstep" (int): Nstep (2 for single parity, 1 for mixed parity)


    k_truncation_mode_triangular:
        - weight is (n_coeff*n + l_coeff*l)
        - compatible with MFDN v15+
        - "truncation_parameters" (dict):
            - "n_coeff": coefficient in front of n
            - "l_coeff": coefficient in front of l
            - "WTmax": maximum weight

    """

    kHO = 0
    kTriangular = 1

k_truncation_mode_ho = TruncationMode.kHO
k_truncation_mode_triangular = TruncationMode.kTriangular


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

class FilenameConfiguration(object):
    """Object to collect filename configuration into common namespace.

    Data members:

        orbitals_int_filename_template (str): filename template for interaction
            tbme basis orbitals
        orbitals_coul_filename_template (str): filename template for Coulomb tbme
            basis orbitals
        orbitals_filename_template (str): filename template for target basis orbitals
        radial_xform_filename_template (str): filename template for change of basis xform
        radial_me_filename_template (str): filename template for radial matrix elements
        radial_olap_int_filename_template (str): filename template for overlaps from interaction tbme basis
        radial_olap_coul_filename_template (str): filename template for overlaps from Coulomb tbme basis
        h2mixer_filename_template (str): filename template for h2mixer input
        natorb_info_filename_template (str): filename template for OBDME info for building natural orbitals
        natorb_obdme_filename_template (str): filename template for static OBDME file for building natural orbitals
        natorb_xform_filename_template (str): filename template for natural orbital xform from previous basis

    Methods:

        orbitals_int_filename(natural_orbital_iteration): Generate filename for
            interaction tbme basisvorbitals, given natural orbital iteration.

        orbital_coul_filename(natural_orbital_iteration): Generate filename for
            Coulomb tbme basis orbitals, given natural orbital iteration.

        orbital_filename(natural_orbital_iteration): Generate filename for target
            basis orbitals, given natural orbital iteration.

        radial_xform_filename(natural_orbital_iteration): Generate filename for
            change of basis xform, given natural orbital iteration.

        radial_me_filename(natural_orbital_iteration, operator_type, power):
            Generate filename for radial matrix elements, given natural orbital
            iteration, operator type, and power.

        radial_olap_int_filename(natural_orbital_iteration): Generate filename
            for overlaps from interaction tbme basis, given natural orbital iteration.

        radial_olap_coul_filename(natural_orbital_iteration): Generate filename
            for overlaps from Coulomb tbme basis, given natural orbital iteration.

        h2mixer_filename(natural_orbital_iteration): Generate filename h2mixer
            input, given natural_orbital_iteration.

        def natorb_info_filename(natural_orbital_iteration): Generate filename for
            MFDn OBDME info output (from which to build natural orbitals), given
            natural_orbital_iteration.

        def natorb_obdme_filename(natural_orbital_iteration): Generate filename
            for MFDn static OBDME output (from which to build natural orbitals),
            given natural_orbital_iteration.

        def natorb_xform_filename(natural_orbital_iteration): Generate filename
            for xform from previous natural orbit to current natural orbit, given
            natural_orbital_iteration.

    """

    # orbital filename templates
    orbitals_int_filename_template = "orbitals-int{:s}.dat"
    orbitals_coul_filename_template = "orbitals-coul{:s}.dat"
    orbitals_filename_template = "orbitals{:s}.dat"
    radial_xform_filename_template = "radial-xform{:s}.dat"
    radial_me_filename_template = "radial-me-{}{}{:s}.dat" # "{}{}" will be replaced by {"r1","r2","k1","k2"}
    radial_olap_int_filename_template = "radial-olap-int{:s}.dat"
    radial_olap_coul_filename_template = "radial-olap-coul{:s}.dat"
    h2mixer_filename_template = "h2mixer{:s}.in"
    natorb_info_filename_template = "natorb-obdme{:s}.info"
    natorb_obdme_filename_template = "natorb-obdme{:s}.dat"
    natorb_xform_filename_template = "natorb-xform{:s}.dat"

    def orbitals_int_filename(self, natural_orbital_iteration):
        """Construct filename for interaction tbme basis orbitals."""
        # don't make the interaction orbital filename dependent on iteration
        return self.orbitals_int_filename_template.format("")

    def orbitals_coul_filename(self, natural_orbital_iteration):
        """Construct filename for Coulomb tbme basis orbitals."""
        # don't make the Coulomb orbital filename dependent on iteration
        return self.orbitals_coul_filename_template.format("")

    def orbitals_filename(self, natural_orbital_iteration):
        """Construct filename for target basis orbitals."""
        return self.orbitals_filename_template.format(natural_orbital_indicator(natural_orbital_iteration))

    def radial_xform_filename(self, natural_orbital_iteration):
        """Construct filename for change of basis xform."""
        return self.radial_xform_filename_template.format(natural_orbital_indicator(natural_orbital_iteration))

    def radial_me_filename(self, natural_orbital_iteration, operator_type, power):
        """Construct filename for radial matrix elements."""
        return self.radial_me_filename_template.format(operator_type,power,natural_orbital_indicator(natural_orbital_iteration))

    def radial_olap_int_filename(self, natural_orbital_iteration):
        """Construct filename for overlaps from interaction tbme basis."""
        return self.radial_olap_int_filename_template.format(natural_orbital_indicator(natural_orbital_iteration))

    def radial_olap_coul_filename(self, natural_orbital_iteration):
        """Construct filename for overlaps from Coulomb tbme basis."""
        return self.radial_olap_coul_filename_template.format(natural_orbital_indicator(natural_orbital_iteration))

    def h2mixer_filename(self, natural_orbital_iteration):
        """Construct filename for h2mixer input."""
        return self.h2mixer_filename_template.format(natural_orbital_indicator(natural_orbital_iteration))

    def natorb_info_filename(self, natural_orbital_iteration):
        """Construct filename for MFDn OBDME info output."""
        return self.natorb_info_filename_template.format(natural_orbital_indicator(natural_orbital_iteration))

    def natorb_obdme_filename(self, natural_orbital_iteration):
        """Construct filename for MFDn OBDME output."""
        return self.natorb_obdme_filename_template.format(natural_orbital_indicator(natural_orbital_iteration))

    def natorb_xform_filename(self, natural_orbital_iteration):
        """Construct filename for natural orbital xform."""
        return self.natorb_xform_filename_template.format(natural_orbital_indicator(natural_orbital_iteration))

filenames = FilenameConfiguration()

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
        task["truncation_mode"] is k_truncation_mode_ho
        and
        task["basis_mode"] in {k_basis_mode_direct,k_basis_mode_dilated}
    ):
        # traditional oscillator run
        template_string = (
            "Z{nuclide[0]}-N{nuclide[1]}-{interaction}-coul{coulomb_flag:d}"
            "-hw{hw:06.3f}"
            "-a_cm{a_cm:g}"
            "-Nmax{Nmax:02d}{mixed_parity_indicator}{fci_indicator}-Mj{Mj:03.1f}"
            "-lan{lanczos:d}-tol{tolerance:.1e}"
            "{natural_orbital_indicator}"
            )
    else:
        raise mcscript.ScriptError("mode not supported by task descriptor")

    truncation_parameters = task["truncation_parameters"]
    if (truncation_parameters["many_body_truncation"]=="fci"):
        fci_indicator = "-fci"
    else:
        fci_indicator = ""
    mixed_parity_indicator = mcscript.utils.ifelse(truncation_parameters["Nstep"]==1,"x","")
    coulomb_flag = int(task["use_coulomb"])
    natural_orbital_iteration = task.get("natorb_iteration")
    natural_orbital_str = ("-natorb" if (natural_orbital_iteration is not None) else "")

    descriptor = template_string.format(
        coulomb_flag=coulomb_flag,
        mixed_parity_indicator=mixed_parity_indicator,
        fci_indicator=fci_indicator,
        natural_orbital_indicator=natural_orbital_str,
        **mcscript.utils.dict_union(task,truncation_parameters)
        )

    return descriptor


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

    # validate truncation mode
    if (task["truncation_mode"] is not k_truncation_mode_ho):
        raise ValueError("expecting truncation_mode to be {} but found {ho_truncation}".format(k_truncation_mode_ho,**task))

    natural_orbital_iteration = task.get("natorb_iteration")

    # generating orbitals for a non-natural orbitals run is the same
    # as generating orbitals for the 0th natural orbital run
    # generate orbitals and basis xform -- target basis
    if (natural_orbital_iteration in {None,0}):
        # generate orbitals -- interaction bases
        mcscript.call(
            [
                configuration.shell_filename("orbital-gen"),
                "--oscillator",
                "{truncation_int[1]:d}".format(**task),
                "{:s}".format(filenames.orbitals_int_filename(natural_orbital_iteration))
            ]
        )
        if (task["use_coulomb"]):
            mcscript.call(
                [
                    configuration.shell_filename("orbital-gen"),
                    "--oscillator",
                    "{truncation_coul[1]:d}".format(**task),
                    "{:s}".format(filenames.orbitals_coul_filename(natural_orbital_iteration))
                ]
            )

        # generate orbitals -- target basis
        if (task["truncation_mode"] is k_truncation_mode_ho):
            truncation_parameters = task["truncation_parameters"]
            if (truncation_parameters["many_body_truncation"]=="Nmax"):
                Nmax_orb = truncation_parameters["Nmax"] + truncation_parameters["Nv"]
            elif (truncation_parameters["many_body_truncation"]=="FCI"):
                Nmax_orb = truncation_parameters["Nmax"]
            mcscript.call(
                [
                    configuration.shell_filename("orbital-gen"),
                    "--oscillator",
                    "{Nmax_orb:d}".format(Nmax_orb=Nmax_orb),
                    "{:s}".format(filenames.orbitals_filename(natural_orbital_iteration))
                ]
            )
            mcscript.call(
                [
                    configuration.shell_filename("radial-gen"),
                    "--identity",
                    filenames.orbitals_filename(natural_orbital_iteration),
                    filenames.radial_xform_filename(natural_orbital_iteration)
                ]
            )
        else:
            # TODO implement non-HO truncation
            raise mcscript.ScriptError("expecting ho_truncation to be True but found {ho_truncation}".format(**task))
    else:
        mcscript.call(
            [
                configuration.shell_filename("natorb-gen"),
                filenames.orbitals_filename(natural_orbital_iteration-1),
                filenames.natorb_info_filename(natural_orbital_iteration-1),
                filenames.natorb_obdme_filename(natural_orbital_iteration-1),
                filenames.natorb_xform_filename(natural_orbital_iteration),
                filenames.orbitals_filename(natural_orbital_iteration)
            ]
        )

def set_up_radial(task):
    """Generate radial integrals and overlaps for MFDn run.

    Operation mode may in general be direct oscillator, dilated
    oscillator, or generic.

    Arguments:
        task (dict): as described in module_docstring

    """

    # get natural orbital iteration; None should be treated the same as 0
    natural_orbital_iteration = task.get("natorb_iteration")
    if (natural_orbital_iteration in {None,0}):
        return set_up_radial_analytic(task)
    elif (natural_orbital_iteration > 0):
        return set_up_radial_natorb(task)

def set_up_radial_analytic(task):
    """Generate radial integrals and overlaps by integration for MFDn run
    in analytic basis.

    Operation mode may in general be direct oscillator, dilated
    oscillator, or generic (TODO).

    Arguments:
        task (dict): as described in module docstring

    """

    # validate basis mode
    if (task["basis_mode"] not in {k_basis_mode_direct,k_basis_mode_dilated}):  # no k_basis_mode_generic yet
        raise ValueError("invalid basis mode {basis_mode}".format(**task))

    # get natural orbital iteration
    natural_orbital_iteration = task.get("natorb_iteration")

    # basis radial code -- expected by radial_utils codes
    basis_radial_code = "oscillator"  # TO GENERALIZE: if not oscillator basis

    # generate radial integrals
    for operator_type in ["r","k"]:
        for power in [1,2]:
            mcscript.call(
                [
                    configuration.shell_filename("radial-gen"),
                    "--kinematic",
                    "{:s}".format(operator_type),
                    "{:d}".format(power),
                    basis_radial_code,
                    filenames.orbitals_filename(natural_orbital_iteration),
                    filenames.radial_me_filename(natural_orbital_iteration, operator_type, power)
                ],
                mode = mcscript.call.serial
            )

    # generate radial overlaps -- generate trivial identities if applicable
    if (task["basis_mode"] in {k_basis_mode_direct}):
        mcscript.call(
            [
                configuration.shell_filename("radial-gen"),
                "--identity",
                filenames.orbitals_int_filename(natural_orbital_iteration),
                filenames.orbitals_filename(natural_orbital_iteration),
                filenames.radial_olap_int_filename(natural_orbital_iteration)
            ],
            mode = mcscript.call.serial
        )
    else:
        b_ratio = math.sqrt(task["hw_int"]/task["hw"])
        mcscript.call(
            [
                configuration.shell_filename("radial-gen"),
                "--overlaps",
                "{:g}".format(b_ratio),
                basis_radial_code,
                filenames.orbitals_int_filename(natural_orbital_iteration),
                filenames.orbitals_filename(natural_orbital_iteration),
                filenames.radial_olap_int_filename(natural_orbital_iteration)
            ],
            mode = mcscript.call.serial
        )
    if (task["use_coulomb"]):
        if (task["basis_mode"] in {k_basis_mode_direct,k_basis_mode_dilated}):
            mcscript.call(
                [
                    configuration.shell_filename("radial-gen"),
                    "--identity",
                    filenames.orbitals_coul_filename(natural_orbital_iteration),
                    filenames.orbitals_filename(natural_orbital_iteration),
                    filenames.radial_olap_coul_filename(natural_orbital_iteration)
                ],
                mode = mcscript.call.serial
            )
        else:
            if task.get("hw_coul_rescaled") is None:
                b_ratio = 1
            else:
                b_ratio = math.sqrt(task["hw_coul_rescaled"]/task["hw"])
            mcscript.call(
                [
                    configuration.shell_filename("radial-gen"),
                    "--overlaps",
                    "{:g}".format(b_ratio),
                    basis_radial_code,
                    filenames.orbitals_coul_filename(natural_orbital_iteration),
                    filenames.orbitals_filename(natural_orbital_iteration),
                    filenames.radial_olap_coul_filename(natural_orbital_iteration)
                ],
                mode = mcscript.call.serial
            )

def set_up_radial_natorb(task):
    """Generate radial integrals and overlaps by transformation for MFDn run
    in natural orbital basis.

    Operation mode must be generic.

    Arguments:
        task (dict): as described in module docstring

    """

    # get natural orbital iteration
    natural_orbital_iteration = task.get("natorb_iteration")

    if (natural_orbital_iteration in {None,0}):
        raise mcscript.ScriptError(
            "cannot set up natural orbital radial files in natorb_iteration {}".format(natural_orbital_iteration)
        )

    # compose radial transform
    mcscript.call(
        [
            configuration.shell_filename("radial-compose"),
            filenames.radial_xform_filename(natural_orbital_iteration-1),
            filenames.natorb_xform_filename(natural_orbital_iteration),
            filenames.radial_xform_filename(natural_orbital_iteration)
        ]
    )

    # compose interaction transform
    mcscript.call(
        [
            configuration.shell_filename("radial-compose"),
            filenames.radial_olap_int_filename(natural_orbital_iteration-1),
            filenames.natorb_xform_filename(natural_orbital_iteration),
            filenames.radial_olap_int_filename(natural_orbital_iteration)
        ]
    )

    # compose Coulomb transform
    mcscript.call(
        [
            configuration.shell_filename("radial-compose"),
            filenames.radial_olap_coul_filename(natural_orbital_iteration-1),
            filenames.natorb_xform_filename(natural_orbital_iteration),
            filenames.radial_olap_coul_filename(natural_orbital_iteration)
        ]
    )

    # transform radial integrals
    for operator_type in ["r","k"]:
        for power in [1,2]:
            mcscript.call(
                [
                    configuration.shell_filename("radial-xform"),
                    filenames.orbitals_filename(natural_orbital_iteration),
                    filenames.radial_xform_filename(natural_orbital_iteration),
                    filenames.radial_me_filename(0, operator_type, power),
                    filenames.radial_me_filename(natural_orbital_iteration, operator_type, power)
                ],
                mode = mcscript.call.serial
            )

def generate_tbme(task):
    """Generate TBMEs for MFDn run.

    Operation mode may in general be direct oscillator, dilated
    oscillator, or generic (TODO).

    Arguments:
        task (dict): as described in module docstring

    """

    # validate basis mode
    if (task["truncation_mode"] is not k_truncation_mode_ho):
        raise ValueError("expecting truncation_mode to be {} but found {ho_truncation}".format(k_truncation_mode_ho,**task))

    # extract parameters for convenience
    natural_orbital_iteration = task.get("natorb_iteration")
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
    xform_truncation_int = task.get("xform_truncation_int")
    if (xform_truncation_int is None):
        xform_truncation_int = task["truncation_int"]
    xform_truncation_coul = task.get("xform_truncation_coul")
    if (xform_truncation_coul is None):
        xform_truncation_coul = task["truncation_coul"]

    # accumulate h2mixer targets
    targets = {}

    # target: Hamiltonian
    if (task.get("hamiltonian")):
        targets["tbme-H"] = task["hamiltonian"]
    else:
        targets["tbme-H"] = mfdn.operators.Hamiltonian(
            A=A, hw=hw, a_cm=a_cm, bsqr_intr=hw_cm/hw,
            use_coulomb=task["use_coulomb"], bsqr_coul=hw_coul_rescaled/hw_coul
        )

    # accumulate observables
    if (task.get("observables")):
        targets.update(task.get("observables"))

    # target: radius squared
    if ("tbme-rrel2" not in targets.keys()):
        targets["tbme-rrel2"] = mfdn.operators.rrel2(A, hw)
    # target: Ncm
    if ("tbme-Ncm" not in targets.keys()):
        targets["tbme-Ncm"] = mfdn.operators.Ncm(A, hw/hw_cm)

    # optional observable sets
    # Hamiltonian components
    if ("H-components" in task["observable_sets"]):
        # target: Trel (diagnostic)
        targets["tbme-Trel"] = mfdn.operators.Trel(A, hw)
        # target: VNN (diagnostic)
        targets["tbme-VNN"] = mfdn.operators.VNN()
        # target: VC (diagnostic)
        if (task["use_coulomb"]):
            targets["tbme-VC"] = mfdn.operators.VC(hw_coul_rescaled/hw_coul)
    # squared angular momenta
    if ("am-sqr" in task["observable_sets"]):
        targets["tbme-L"] = mfdn.operators.L()
        targets["tbme-Sp"] = mfdn.operators.Sp()
        targets["tbme-Sn"] = mfdn.operators.Sn()
        targets["tbme-S"] = mfdn.operators.S()
        targets["tbme-J"] = mfdn.operators.J()

    # get set of required sources
    required_sources = set()
    required_sources.update(*[op.keys() for op in targets.values()])

    # accumulate h2mixer input lines
    lines = []

    # initial comment
    lines.append("# task: {}".format(task))
    lines.append("")

    # global mode definitions
    target_truncation = task["target_truncation"]
    if (target_truncation is None):
        # automatic derivation
        if (task["truncation_mode"] is k_truncation_mode_ho):
            truncation_parameters = task["truncation_parameters"]
            if (truncation_parameters["many_body_truncation"]=="Nmax"):
                # important: truncation of orbitals file, one-body
                # truncation of interaction file, and MFDn
                # single-particle shells (beware 1-based) must agree
                N1_max = truncation_parameters["Nv"]+truncation_parameters["Nmax"]
                N2_max = 2*truncation_parameters["Nv"]+truncation_parameters["Nmax"]
                target_weight_max = weight_max_string((N1_max,N2_max))
            elif (truncation_parameters["many_body_truncation"]=="FCI"):
                N1_max = truncation_parameters["Nmax"]
                target_weight_max = weight_max_string(("ob",N1_max))
        else:
            # calculation of required weight_max will require external program for occupation counting
            pass
    else:
        # given value
        target_weight_max = target_truncation
    lines.append("set-target-indexing {orbitals_filename} {target_weight_max}".format(
        orbitals_filename=filenames.orbitals_filename(natural_orbital_iteration),
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
            radial_me_filename = filenames.radial_me_filename(natural_orbital_iteration, operator_type, power)
            lines.append("define-radial-operator {} {} {}".format(operator_type,power,radial_me_filename))
    lines.append("")

    # sources: h2mixer built-ins
    builtin_sources = (mfdn.operators.kinematic_operator_set | mfdn.operators.angular_momentum_operator_set)
    for source in sorted(builtin_sources & required_sources):
        lines.append("define-source operator " + source)
    lines.append("")

    # sources: VNN
    if ("VNN" in required_sources):
        VNN_filename = configuration.interaction_filename(
            "{}-{}-{:g}.bin".format(
                task["interaction"],
                mcscript.utils.dashify(task["truncation_int"]),
                task["hw_int"]
            )
        )
        if (task["basis_mode"]==k_basis_mode_direct and natural_orbital_iteration in {None,0}):
            lines.append("define-source input VNN {VNN_filename}".format(VNN_filename=VNN_filename,**task))
        else:
            xform_weight_max_int = weight_max_string(xform_truncation_int)
            lines.append("define-source xform VNN {VNN_filename} {xform_weight_max_int} {radial_olap_int_filename}".format(
                VNN_filename=VNN_filename,
                xform_weight_max_int=xform_weight_max_int,
                radial_olap_int_filename=filenames.radial_olap_int_filename(natural_orbital_iteration),
                **task
            ))

    # sources: Coulomb
    #
    # Note: This is the "unscaled" Coulomb, still awaiting the scaling
    # factor from dilation.
    if ("VC_unscaled" in required_sources):
        VC_filename = configuration.interaction_filename(
            "{}-{}-{:g}.bin".format(
                "VC",
                mcscript.utils.dashify(task["truncation_coul"]),
                task["hw_coul"]
            )
        )
        if (task["basis_mode"] in {k_basis_mode_direct,k_basis_mode_dilated} and natural_orbital_iteration in {None,0}):
            lines.append("define-source input VC_unscaled {VC_filename}".format(VC_filename=VC_filename,**task))
        else:
            xform_weight_max_coul = weight_max_string(xform_truncation_coul)
            lines.append("define-source xform VC_unscaled {VC_filename} {xform_weight_max_coul} {radial_olap_coul_filename}".format(
                VC_filename=VC_filename,
                xform_weight_max_coul=xform_weight_max_coul,
                radial_olap_coul_filename=filenames.radial_olap_coul_filename(natural_orbital_iteration),
                **task
            ))

    lines.append("")


    # targets: generate h2mixer input
    for (basename, operator) in targets.items():
        lines.append("define-target work/"+basename+".bin")
        for (source, coefficient) in operator.items():
            lines.append("  add-source {:s} {:e}".format(source, coefficient))
        lines.append("")

    # ensure terminal line
    lines.append("")

    # diagnostic: log input lines to file
    #
    # This is purely for easy diagnostic purposes, since lines will be
    # fed directly to h2mixer as stdin below.
    mcscript.utils.write_input(filenames.h2mixer_filename(natural_orbital_iteration),input_lines=lines,verbose=False)

    # create work directory if it doesn't exist yet (-p)
    mcscript.call(["mkdir","-p","work"])

    # invoke h2mixer
    mcscript.call(
        [
            configuration.shell_filename("h2mixer")
        ],
        input_lines=lines,
        mode = mcscript.call.serial
    )

def run_mfdn_v14_b06(task):
    """Generate input file and execute MFDn version 14 beta 06.

    Arguments:
        task (dict): as described in module docstring

    Raises:
        mcscript.ScriptError: if MFDn output not found

    """

    # validate truncation mode
    if (task["truncation_mode"] is not k_truncation_mode_ho):
        raise ValueError("expecting truncation_mode to be {} but found {ho_truncation}".format(k_truncation_mode_ho,**task))

    # accumulate MFDn input lines
    lines = []

    # base parameters
    natural_orbital_iteration = task.get("natorb_iteration")
    truncation_parameters = task["truncation_parameters"]
    twice_Mj = int(2*task["Mj"])
    if (truncation_parameters["many_body_truncation"]=="Nmax"):
        Nmax_orb = truncation_parameters["Nmax"] + truncation_parameters["Nv"]
    elif (truncation_parameters["many_body_truncation"]=="FCI"):
        Nmax_orb = truncation_parameters["Nmax"]
    Nshell = Nmax_orb+1
    if (task["basis_mode"] in {k_basis_mode_direct,k_basis_mode_dilated} and natural_orbital_iteration in {None,0}):
        hw_for_trans = task["hw"]
    else:
        hw_for_trans = 0  # disable MFDn hard-coded oscillator one-body observables
    ## ndiag = int(os.environ.get("MFDN_NDIAG",0))  # allows override for spares, not so elegant
    ndiag = task.get("ndiag")
    if ndiag is None:
        ndiag = 0
    if (truncation_parameters["Nstep"]==2):
        Nmin = truncation_parameters["Nmax"]%2
    else:
        Nmin = 1

    lines.append("{:d}  # IFLAGMBSI".format(0))
    lines.append("{ndiag:d}  # ndiag (0: no spares, automatic ndiag)".format(ndiag=ndiag,**task))
    lines.append("{:d}  # number of classes".format(2))
    lines.append("{nuclide[0]:d}  # protons (class 1 particles)".format(**task))
    lines.append("{nuclide[1]:d}  # protons (class 2 particles)".format(**task))
    lines.append("1 {Nshell:d}  # min, max # S.P. shells for class 1 particles".format(Nshell=Nshell,**task))
    lines.append("1 {Nshell:d}  # min, max # S.P. shells for class 2 particles".format(Nshell=Nshell,**task))
    lines.append("{Nmin:d} {Nmax:d} {Nstep:d}  # N_min, N_max, delta_N".format(Nmin=Nmin,**mcscript.utils.dict_union(task,truncation_parameters)))
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
    obs_basename_list += list(task["observables"].keys())

    # tbo: log tbo names in separate file to aid future data analysis
    mcscript.utils.write_input("tbo_names.dat",input_lines=obs_basename_list)

    # tbo: write list of operators
    lines.append("tbme-H")  # Hamiltonian basename
    num_obs = 2 + len(obs_basename_list)
    if (num_obs > 8):
        raise mcscript.ScriptError("Too many observables for MFDn v14")
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
        if ((g_rel != (truncation_parameters["Nmax"]%2)) and (truncation_parameters["Nstep"] != 1)):
            raise ValueError("invalid parity for reference state")

    # obdme: write reference state list
    for reference_state in task["obdme_reference_state_list"]:
        lines.append("{:d} {:d} {:d}".format(2*reference_state[0],reference_state[1],reference_state[2]))

    # ensure terminal line
    lines.append("")

    # generate MFDn input file
    mcscript.utils.write_input("work/mfdn.dat",input_lines=lines)

    # import partitioning file
    if (task["partition_filename"] is not None):
        if (not os.path.exists(task["partition_filename"])):
            raise mcscript.ScriptError("partition file not found")
        mcscript.call(["cp","--verbose",task["partition_filename"],"work/mfdn_partitioning.info"])

    # enter work directory
    os.chdir("work")

    # invoke MFDn
    mcscript.call(
        [
            configuration.mfdn_filename(task["mfdn_executable"])
        ],
        # mode = mcscript.call.hybrid,
        check_return=True
    )

    # test for basic indications of success
    if (not os.path.exists("mfdn.out")):
        raise mcscript.ScriptError("mfdn.out not found")
    if (not os.path.exists("mfdn.res")):
        raise mcscript.ScriptError("mfdn.res not found")

    # leave work directory
    os.chdir("..")

def save_mfdn_v14_output(task):
    """Generate input file and execute MFDn version 14 beta 06.

    Arguments:
        task (dict): as described in module docstring

    Raises:
        mcscript.ScriptError: if MFDn output not found

    """

    # save quick inspection copies of mfdn.{res,out}
    natural_orbital_iteration = task.get("natorb_iteration")
    descriptor = task["descriptor"] + natural_orbital_indicator(natural_orbital_iteration)
    print("Saving basic output files...")
    res_filename = "{:s}-mfdn-{:s}.res".format(mcscript.run.name,descriptor)
    mcscript.call(["cp","--verbose","work/mfdn.res",res_filename])
    out_filename = "{:s}-mfdn-{:s}.out".format(mcscript.run.name,descriptor)
    mcscript.call(["cp","--verbose","work/mfdn.out",out_filename])

    # save OBDME files for next natural orbital iteration
    if (natural_orbital_iteration is not None):
        print("Saving OBDME files for natural orbital generation...")
        obdme_info_filename = "mfdn.rppobdme.info"
        mcscript.call(
            [
                "cp","--verbose",
                "work/{}".format(obdme_info_filename),
                filenames.natorb_info_filename(natural_orbital_iteration)
            ]
        )
        obdme_filename = glob.glob("work/mfdn.statrobdme.seq{:03d}*".format(task["natorb_base_state"]))
        mcscript.call(
            [
                "cp","--verbose",
                obdme_filename[0],
                filenames.natorb_obdme_filename(natural_orbital_iteration)
            ]
        )

    # save full archive of input, log, and output files
    print("Saving full output files...")
    # logging
    archive_file_list = [
        filenames.h2mixer_filename(natural_orbital_iteration),
        "tbo_names.dat"
        ]
    # orbital information
    archive_file_list += [
        filenames.orbitals_int_filename(natural_orbital_iteration),
        filenames.orbitals_filename(natural_orbital_iteration),
        ]
    # transformation information
    archive_file_list += [
        filenames.radial_xform_filename(natural_orbital_iteration),
        # filenames.radial_me_filename(natural_orbital_iteration, operator_type, power),
        filenames.radial_olap_int_filename(natural_orbital_iteration),
        ]
    # Coulomb information:
    if task["use_coulomb"]:
        archive_file_list += [
            filenames.orbitals_coul_filename(natural_orbital_iteration),
            filenames.radial_olap_coul_filename(natural_orbital_iteration),
        ]
    # natural orbital information
    if natural_orbital_iteration not in {None}:
        archive_file_list += [
            filenames.natorb_info_filename(natural_orbital_iteration),
            filenames.natorb_obdme_filename(natural_orbital_iteration),
            ]
    if (natural_orbital_iteration not in {None,0}):
        archive_file_list += [
            filenames.natorb_xform_filename(natural_orbital_iteration),
            ]
    # MFDn output
    archive_file_list += [
        "work/mfdn.dat","work/mfdn.out","work/mfdn.res","work/mfdn_partitioning.generated","work/mfdn_spstates.info"
    ]
    # renamed versions
    archive_file_list += [out_filename,res_filename]
    # MFDN obdme
    if (task["save_obdme"]):
        archive_file_list += glob.glob("work/*obdme*")
    # generate archive (outside work directory)
    archive_filename = "{:s}-mfdn-{:s}.tgz".format(mcscript.run.name,descriptor)
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
                res_filename,out_filename,archive_filename,
                "--target-directory={}".format(mcscript.task.results_dir)
            ]
        )

    # cleanup of wave function files
    scratch_file_list = glob.glob("work/*")
    mcscript.call(["rm", "-vf"] + scratch_file_list)


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
        "%s-results-%s.tgz" % (mcscript.run.name, mcscript.date_tag())
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



if (__name__ == "__MAIN__"):
    pass
