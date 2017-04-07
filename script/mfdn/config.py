import enum
import os

import mcscript.utils

from . import utils

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

    kHO (alias k_truncation_mode_ho):
        - traditional Nmax truncation; weight is (2*n + l)
        - compatible with MFDn v14+
        - "truncation_parameters" (dict):
            - "many_body_truncation" (str): many-body truncation rank ("FCI" for one-body, or
                "Nmax" for A-body)
            - "Nv" (int): N of valence shell (for use in truncation)
            - "Nmax" (int): many-body oscillator cutoff (interpreted as one-body Nmax_orb for "FCI"
                 truncation, or many-body excitation cutoff Nmax for "Nmax" truncation)
            - "Nstep" (int): Nstep (2 for single parity, 1 for mixed parity)


    kTriangular (alias k_truncation_mode_triangular):
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

class Environment(object):
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


environ = Environment()

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
        return self.orbitals_filename_template.format(utils.natural_orbital_indicator(natural_orbital_iteration))

    def radial_xform_filename(self, natural_orbital_iteration):
        """Construct filename for change of basis xform."""
        return self.radial_xform_filename_template.format(utils.natural_orbital_indicator(natural_orbital_iteration))

    def radial_me_filename(self, natural_orbital_iteration, operator_type, power):
        """Construct filename for radial matrix elements."""
        return self.radial_me_filename_template.format(operator_type,power,utils.natural_orbital_indicator(natural_orbital_iteration))

    def radial_olap_int_filename(self, natural_orbital_iteration):
        """Construct filename for overlaps from interaction tbme basis."""
        return self.radial_olap_int_filename_template.format(utils.natural_orbital_indicator(natural_orbital_iteration))

    def radial_olap_coul_filename(self, natural_orbital_iteration):
        """Construct filename for overlaps from Coulomb tbme basis."""
        return self.radial_olap_coul_filename_template.format(utils.natural_orbital_indicator(natural_orbital_iteration))

    def h2mixer_filename(self, natural_orbital_iteration):
        """Construct filename for h2mixer input."""
        return self.h2mixer_filename_template.format(utils.natural_orbital_indicator(natural_orbital_iteration))

    def natorb_info_filename(self, natural_orbital_iteration):
        """Construct filename for MFDn OBDME info output."""
        return self.natorb_info_filename_template.format(utils.natural_orbital_indicator(natural_orbital_iteration))

    def natorb_obdme_filename(self, natural_orbital_iteration):
        """Construct filename for MFDn OBDME output."""
        return self.natorb_obdme_filename_template.format(utils.natural_orbital_indicator(natural_orbital_iteration))

    def natorb_xform_filename(self, natural_orbital_iteration):
        """Construct filename for natural orbital xform."""
        return self.natorb_xform_filename_template.format(utils.natural_orbital_indicator(natural_orbital_iteration))

filenames = FilenameConfiguration()
