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
            oscillator (modes.BasisMode.kDirect), dilated oscillator
            (modes.BasisMode.kDilated), or generic run mode
            (modes.BasisMode.kGeneric), as explained further where this
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

        # truncation parameters
        "sp_truncation_mode" (modes.SingleParticleTruncationMode): enumerated value
            indicating single-particle basis truncation; see docstring of
            SingleParticleTruncationMode for information.
        "mb_truncation_mode" (modes.ManyBodyTruncationMode): enumerated value
            indicating many-body basis truncation; see docstring of
            ManyBodyTruncationMode for information.
        "truncation_parameters" (dict): truncation parameters, specific to each enumerated
            truncation type; see docstrings of SingleParticleTruncationMode and
            ManyBodyTruncationMode for full documentation

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
        "observables" (list of ("basename", CoefficientDict) tuples): additional observable definitions
            (see mfdn.operators)
        "observable_sets" (list of str): codes for predefined observable sets to include
            "H-components": Hamiltonian terms
            "am-sqr": squared angular momenta
            "isospin": isospin observables
            "R20K20": center-of-mass diagnostic observables (TODO)

        # wavefunction storage
        "save_wavefunctions" (bool): whether or not to save smwf files in (separate) archive

        # version parameters
        "h2_format" (int): h2 file format to use (values include: 0, 15099)
        "mfdn_executable" (string): mfdn executable name
        "mfdn_driver" (module): mfdn driver module

        # natural orbital parameters
        "natural_orbitals" (bool): enable/disable natural orbitals
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
    - 3/17/17 (pjf): Split mfdn/__init__.py into submodules.
    - 7/31/17 (pjf): Add mfdn_driver field.
    - 8/11/17 (pjf):
      + Replace TruncationMode with SingleParticleTruncationMode and ManyBodyTruncationMode.
      + Replace truncation_mode key with sp_truncation_mode and mb_truncation_mode.
      + Fix FCI truncation.
    - 09/22/17 (pjf): Take "observables" as list of tuples instead of dict.
    - 09/24/17 (pjf): Add option to save wavefunctions for postprocessing.
"""

__ALL__ = ['descriptors', 'handlers', 'radial', 'operators', 'tbme', 'utils', 'config']
from . import descriptors, handlers, radial, operators, tbme, utils, modes

if (__name__ == "__MAIN__"):
    pass
