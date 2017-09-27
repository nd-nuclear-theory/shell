"""modes.py -- define configuration options for MFDn scripting

Patrick Fasano
University of Notre Dame

- 3/22/17 (pjf): Created, split from __init__.py.
- 4/7/17 (pjf): Rename Configuration -> Environment.
- 6/3/17 (pjf): Remove dependence of filenames on natural orbital iteration.
- 6/5/17 (pjf): Clean up formatting.
- 8/11/17 (pjf): Split TruncationMode into SingleParticleTruncationMode and
    ManyBodyTruncationMode.
- 08/26/17 (pjf): Add parity flag for WeightMax many-body truncation mode.
- 9/12/17 (mac): Put mfdn executable filename under common mcscript install directory.
- 09/12/17 (pjf): Split config.py -> mode.py + environ.py.
- 09/27/17 (pjf): Add MFDnRunMode options for counting-only modes.
"""

import enum


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
# MFDn run modes
################################################################

@enum.unique
class MFDnRunMode(enum.IntEnum):
    """MFDn run mode (IFLAG_mode)

    kNormal: 0, normal MFDn diagonalization run

    kDimension: 1, count basis dimension

    kNonzeros: 3, count nonzero matrix elements

    """

    kNormal = 0
    kDimension = 1
    kNonzeros = 3
