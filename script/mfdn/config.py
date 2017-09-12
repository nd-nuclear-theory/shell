"""config.py -- define configuration options for MFDn scripting

Patrick Fasano
University of Notre Dame

- 3/22/17 (pjf): Created, split from __init__.py.
- 4/7/17 (pjf): Rename Configuration -> Environment.
- 6/3/17 (pjf): Remove dependence of filenames on natural orbital iteration.
- 6/5/17 (pjf): Clean up formatting.
- 8/11/17 (pjf): Split TruncationMode into SingleParticleTruncationMode and
    ManyBodyTruncationMode.
- 08/26/17 (pjf): Add parity flag for WeightMax many-body truncation mode.
"""

import os
import enum

import mcscript.parameters
import mcscript.utils


################################################################
# radial basis modes
################################################################

@enum.unique
class BasisMode(enum.Enum):
    """General modes of operation for radial basis

    kDirect:
      - no0 basis is oscillator basis (hw)
      - source basis for VNN is oscillator basis of same
        oscillator length (hw_int=hw); therefore no transformation
        needed on VNN TBMEs
      - Coulomb TBMEs need only scaling for dilation (hw_c -> hw)
      - MFDn can use built-in oscillator OBMEs for observables
    kDilated:
      - no0 basis is oscillator basis (hw)
      - source basis for VNN is oscillator basis of different
        oscillator length; therefore transformation
        needed on VNN TBMEs (hw_int -> hw)
      - Coulomb TBMEs need only scaling for dilation (hw_c -> hw)
      - MFDn can use built-in oscillator OBMEs for observables

    kGeneric:
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


################################################################
# truncation modes
################################################################

@enum.unique
class SingleParticleTruncationMode(enum.Enum):
    """General truncation modes for single-particle basis

    kNmax:
        - traditional Nmax truncation; weight is (2*n + l)
        - compatible with MFDn v14+
        - "truncation_parameters" (dict) must contain:
            - "Nv" (int): N of valence shell (for use in truncation)
            - "Nmax" (int): single-particle excitation oscillator cutoff
                (interpreted as one-body Nmax_orb for "FCI" truncation, or
                many-body excitation cutoff Nmax for "Nmax" truncation)
            - "Nstep" (int): Nstep (2 for single parity, 1 for mixed parity)


    kTriangular:
        - weight is (n_coeff*n + l_coeff*l)
        - compatible with MFDn v15+
        - "truncation_parameters" (dict) must contain:
            - "n_coeff" (float): coefficient in front of n
            - "l_coeff" (float): coefficient in front of l
            - "sp_weight_max" (float): maximum weight for single-particle orbitals

    """

    kNmax = 1
    kTriangular = 2


@enum.unique
class ManyBodyTruncationMode(enum.Enum):
    """General truncation modes for many-body basis

    kNmax:
        - traditional Nmax truncation; weight is (2*n + l)
        - compatible with MFDn v14+
        - must be used with SingleParticleTruncationMode.kNmax
        - "truncation_parameters" (dict) must contain:
            - "Nv" (int): N of valence shell (for use in truncation)
            - "Nmax" (int): many-body excitation cutoff
            - "Nstep" (int): Nstep (2 for single parity, 1 for mixed parity)


    kWeightMax:
        - compatible with MFDn v15+
        - "truncation_parameters" (dict) must contain:
            - "mb_weight_max" (float): maximum weight for many-body states
            - "parity" (int): absolute parity for run (+1, 0, or -1)


    kFCI:
        - compatible with MFDn v14+
        - many-body basis constrained only by single-particle basis
        - "truncation_parameters" (dict) must contain:
            - "parity" (int): absolute parity for run (+1, 0, or -1)
    """

    kNmax = 1
    kWeightMax = 2
    kFCI = 3


################################################################
# configuration
################################################################

class Environment(object):
    """Object to collect MFDn environment configuration parameters into common name space.

    Data members:

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
        self.mfdn_dir = os.environ.get("SHELL_MFDN_DIR")
        self.data_dir_h2_list = os.environ.get("SHELL_DATA_DIR_H2").split(":")
        self.interaction_run_list = []

    def shell_filename(self, name):
        """Construct filename for shell package executable."""
        return os.path.join(mcscript.parameters.run.install_dir, "shell", "bin", name)

    def mfdn_filename(self, name):
        """Construct filename for MFDn executable."""
        return os.path.join(self.mfdn_dir, name)

    def interaction_filename(self, name):
        """Construct filename for interaction h2 file."""
        return mcscript.utils.search_in_subdirectories(self.data_dir_h2_list, self.interaction_run_list, name)


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

        orbitals_int_filename(postfix): Generate filename for interaction tbme
            basis orbitals, given postfix.

        orbital_coul_filename(postfix): Generate filename for Coulomb tbme basis
            orbitals, given postfix.

        orbital_filename(postfix): Generate filename for target basis orbitals,
            given postfix.

        radial_xform_filename(postfix): Generate filename for change of basis
            xform, given postfix.

        radial_me_filename(postfix, operator_type, power): Generate filename for
            radial matrix elements, given natural orbital iteration, operator
            type, and power.

        radial_olap_int_filename(postfix): Generate filename for overlaps from
            interaction tbme basis, given postfix.

        radial_olap_coul_filename(postfix): Generate filename for overlaps from
            Coulomb tbme basis, given postfix.

        h2mixer_filename(postfix): Generate filename h2mixer input, given
            postfix.

        def natorb_info_filename(postfix): Generate filename for MFDn OBDME info
            output (from which to build natural orbitals), given postfix.

        def natorb_obdme_filename(postfix): Generate filename for MFDn static
            OBDME output (from which to build natural orbitals), given postfix.

        def natorb_xform_filename(postfix): Generate filename for xform from
            previous natural orbit to current natural orbit, given postfix.
    """

    # orbital filename templates
    orbitals_int_filename_template = "orbitals-int{:s}.dat"
    orbitals_coul_filename_template = "orbitals-coul{:s}.dat"
    orbitals_filename_template = "orbitals{:s}.dat"
    radial_xform_filename_template = "radial-xform{:s}.dat"
    radial_me_filename_template = "radial-me-{}{}{:s}.dat"  # "{}{}" will be replaced by {"r1","r2","k1","k2"}
    radial_olap_int_filename_template = "radial-olap-int{:s}.dat"
    radial_olap_coul_filename_template = "radial-olap-coul{:s}.dat"
    h2mixer_filename_template = "h2mixer{:s}.in"
    natorb_info_filename_template = "natorb-obdme{:s}.info"
    natorb_obdme_filename_template = "natorb-obdme{:s}.dat"
    natorb_xform_filename_template = "natorb-xform{:s}.dat"

    def orbitals_int_filename(self, postfix):
        """Construct filename for interaction tbme basis orbitals."""
        # don't make the interaction orbital filename dependent on iteration
        return self.orbitals_int_filename_template.format("")

    def orbitals_coul_filename(self, postfix):
        """Construct filename for Coulomb tbme basis orbitals."""
        # don't make the Coulomb orbital filename dependent on iteration
        return self.orbitals_coul_filename_template.format("")

    def orbitals_filename(self, postfix):
        """Construct filename for target basis orbitals."""
        return self.orbitals_filename_template.format(postfix)

    def radial_xform_filename(self, postfix):
        """Construct filename for change of basis xform."""
        return self.radial_xform_filename_template.format(postfix)

    def radial_me_filename(self, postfix, operator_type, power):
        """Construct filename for radial matrix elements."""
        return self.radial_me_filename_template.format(operator_type, power, postfix)

    def radial_olap_int_filename(self, postfix):
        """Construct filename for overlaps from interaction tbme basis."""
        return self.radial_olap_int_filename_template.format(postfix)

    def radial_olap_coul_filename(self, postfix):
        """Construct filename for overlaps from Coulomb tbme basis."""
        return self.radial_olap_coul_filename_template.format(postfix)

    def h2mixer_filename(self, postfix):
        """Construct filename for h2mixer input."""
        return self.h2mixer_filename_template.format(postfix)

    def natorb_info_filename(self, postfix):
        """Construct filename for MFDn OBDME info output."""
        return self.natorb_info_filename_template.format(postfix)

    def natorb_obdme_filename(self, postfix):
        """Construct filename for MFDn OBDME output."""
        return self.natorb_obdme_filename_template.format(postfix)

    def natorb_xform_filename(self, postfix):
        """Construct filename for natural orbital xform."""
        return self.natorb_xform_filename_template.format(postfix)


filenames = FilenameConfiguration()