/****************************************************************
  tbme_separable.h

  Evaluates TBMEs of identity and certain separable operators.

  Normalization convention: All matrix elements are stored as NAS
  RMEs.  These RMEs are stored under the group theory Wigner-Eckart
  normalization convention (i.e., "no dimension factor out front, just
  the Clebsch"), but, for scalar operators, note that this RME is
  equivalently, and more simply, the branched ME (with M'=M).

  Symmetrization convention: The full square matrix is populated on
  diagonal sectors.

  Mark A. Caprio
  University of Notre Dame

  + 10/26/16 (mac): Created.  Supplants code from shell_2body.cpp
    (created 4/25/11) for identity operator, h2gen.cpp (created
    3/13/12) for kinematic operators, and shell_separable.cpp (created
    5/24/15) for angular momentum operators.
  + 11/1/16 (mac):
    - Make identity operator the A-dependent many-body identity.
    - Convert from AS to NAS storage.
  + 11/4/16 (mac):
    - Complete implementation of kinematic operators.
    - Attempt OpenMP parallelization of kinematic operators.
  + 11/6/16 (mac):
    - Remove use of size() in omp for loop limit.
    - Refactor implementation of KinematicScalarTBME.
  + 09/20/17 (pjf): Add initial support for isospin operators (T^2,
    Tz, T+, T-).
  + 09/20/17 (mac): Fix normalization of NAS matrix elements for
    new isospin operators.
  + 09/21/17 (pjf): Fix missing Hat() factors in Racah's reduction
    formula for isospin operators.
  + 02/23/18 (pjf): Implement kinematic operators in terms of generic
    Racah's reduction formula code.
  + 12/10/18 (pjf): Deprecate AngularMomentumOperatorFamily and
    AngularMomentumOperatorSpecies in favor of am::AngularMomentumOperatorType
    and basis::OperatorTypePN, respectively.
  + 01/24/19 (pjf):
    - Rewrite UpgradeOneBodyOperatorJJJPN and RacahReduceTensorProductJJJPN
      for generic nonscalar operators.
    - Remove RacahReduceDotProductJJJPN (redundant to RacahReduceTensorProductJJJPN).

****************************************************************/

#ifndef TBME_SEPARABLE_H_
#define TBME_SEPARABLE_H_

#include "eigen3/Eigen/Dense"

#include "am/rme.h"
#include "basis/jjjpn_scheme.h"
#include "basis/proton_neutron.h"
#include "mcutils/deprecated.h"
#include "obme/obme_io.h"

namespace shell {
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  // upgrade one-body operator
  ////////////////////////////////////////////////////////////////

  basis::OperatorBlock<double>
  UpgradeOneBodyOperatorJJJPN(
      const basis::OrbitalSpaceLJPN& ob_orbital_space,
      const basis::OrbitalSectorsLJPN& ob_sectors,
      const basis::OperatorBlocks<double>& ob_matrices,
      const basis::TwoBodySectorsJJJPN::SectorType& sector,
      int A
    );
  // Populate matrix (for a given JJJPN sector) with two-body matrix
  // elements for the A-body "upgraded" one-body operator.
  //
  // See intrinsic equation (TBD).
  //
  // Arguments:
  //   ob_orbital_space (basis::OrbitalSpaceLJPN) : orbital space for
  //     one-body operator
  //   ob_sectors (basis::OrbitalSectorsLJPN) : sectors for one-body
  //     operator
  //   ob_matrices (basis::OperatorBlocks) : matrices for one-body
  //     operator
  //   sector (basis::TwoBodySectorsJJJPN::SectorType) : the sector to
  //     populate
  //   A (int): atomic mass number
  //
  // Returns:
  //   (basis::OperatorBlock<double>) : the matrix for this sector

  ////////////////////////////////////////////////////////////////
  // two-system tensor product
  ////////////////////////////////////////////////////////////////

  basis::OperatorBlock<double>
  RacahReduceTensorProductJJJPN(
      const basis::OrbitalSpaceLJPN& ob_orbital_space,
      const basis::OrbitalSectorsLJPN& ob_sectors1,
      const basis::OperatorBlocks<double>& ob_matrices1,
      const basis::OrbitalSectorsLJPN& ob_sectors2,
      const basis::OperatorBlocks<double>& ob_matrices2,
      const basis::TwoBodySectorsJJJPN::SectorType& sector,
      int J0
    );
  // Populate matrix (for a given JJJPN sector) with two-body matrix
  // elements for the symmetrized tensor-product of two one-body operators.
  //
  // See intrinsic equation (TBD).
  //
  // Arguments:
  //   ob_orbital_space (basis::OrbitalSpaceLJPN) : orbital space for
  //     one-body operator
  //   ob_sectors1 (basis::OrbitalSectorsLJPN) : sectors for one-body
  //     operator 1
  //   ob_matrices1 (basis::OperatorBlocks) : matrices for one-body
  //     operator 1
  //   ob_sectors2 (basis::OrbitalSectorsLJPN) : sectors for one-body
  //     operator 2
  //   ob_matrices2 (basis::OperatorBlocks) : matrices for one-body
  //     operator 2
  //   sector (basis::TwoBodySectorsJJJPN::SectorType) : the sector to
  //     populate
  //   J0 (int): angular momentum of coupled operator
  //
  // Returns:
  //   (basis::OperatorBlock<double>) : the matrix for this sector

  ////////////////////////////////////////////////////////////////
  // two-body identity operator
  ////////////////////////////////////////////////////////////////

  Eigen::MatrixXd
  IdentityOperatorMatrixJJJPN(
      const basis::TwoBodySectorsJJJPN::SectorType& sector,
      int A
    );
  // Populate matrix (for a given JJJPN sector) with two-body matrix
  // elements for the A-body identity matrix.
  //
  // See csbasis [PRC 86, 034312 (2012)] equation (D4).
  //
  // Precondition: It is assumed that the sector is valid for a scalar
  // (and isoscalar and positive parity), i.e., is a diagonal sector.
  //
  // Arguments:
  //   sector (basis::TwoBodySectorsJJJPN::SectorType) : the sector to
  //     populate
  //   A (int): atomic mass number
  //
  // Returns:
  //   (Eigen::MatrixXd) : the matrix for this sector

  ////////////////////////////////////////////////////////////////
  // loop timing test
  ////////////////////////////////////////////////////////////////

  Eigen::MatrixXd
  TimingTestMatrixJJJPN(
      const basis::TwoBodySectorsJJJPN::SectorType& sector,
      bool loop,
      bool store
    );
  // Provides function call with same control structure as the matrix
  // generating functions below but without any compute load, for
  // profiling purposes.
  //
  // Arguments:
  //   sector (basis::TwoBodySectorsJJJPN::SectorType) : the sector to
  //     populate
  //   loop (bool): do main loop over entries
  //   store (bool): do memory access to write entry

  ////////////////////////////////////////////////////////////////
  // kinematic operators
  ////////////////////////////////////////////////////////////////

  enum class CoordinateSqrOperatorType {kUTSqr,kVT1T2};

  Eigen::MatrixXd
  CoordinateSqrMatrixJJJPN(
      const basis::OrbitalSpaceLJPN& radial_orbital_space,
      const basis::OrbitalSectorsLJPN& radial_sectors,
      const basis::OperatorBlocks<double>& radial_matrices,
      const RadialOperatorType& radial_operator_type,
      CoordinateSqrOperatorType coordinate_sqr_operator_type,
      const basis::TwoBodySectorsJJJPN::SectorType& sector,
      int A, int J0
    );
  // Populate matrix with two-body matrix elements of coordinate-squared
  // operator.
  //
  // Case kUTSqr:
  //
  //   Populate matrix with two-body matrix elements for the
  //   *one-body* operator U_(T^2) obtained from a scalar one-body
  //   operator T^2, i.e., r^2 or k^2, on a given JJJPN sector.
  //
  //   This one-body operator is related to the two-body operator
  //   V_(T^2), also defined in csbasis, by U_(T^2) = 1/(A-1)*V_(T^2).
  //   See csbasis (51).
  //
  //   Obtained by csbasis (52)-(54).
  //
  // Case kUTSqr:
  //
  //   Populate matrix with two-body matrix elements for the
  //   *two-body* operator V_(T1.T2) obtained from a a dot product of
  //   one-body vector operators (i.e., T = r or k), on a given JJJPN
  //   sector.
  //
  //   Obtained by csbasis (55)-(60).
  //
  // Precondition: It is assumed that the sector is valid for an isoscalar.
  //
  // Arguments:
  //   radial_orbital_space, radial_sectors, radial_matrices (...):
  //      definition of radial matrix elements (T^2 for UTSqr or T for VT1T2)
  //   radial_operator_type (shell::RadialOperatorType): whether
  //     coordinate (kR) or momentum (kK) space (only affects phase factor
  //     calculation)
  //   coordinate_sqr_operator_type (CoordinateSqrOperatorType): whether UTSqr or VT1T2 matrix
  //   sector (basis::TwoBodySectorsJJJPN::SectorType) : the sector to
  //     populate
  //   A (int): atomic mass number
  //   J0 (int): two-body operator J0
  //
  // Returns:
  //   (Eigen::MatrixXd) : the matrix for this sector

  ////////////////////////////////////////////////////////////////
  // angular momentum operators
  ////////////////////////////////////////////////////////////////

  // Note: The implementation of AngularMomentumMatrixJJJPN is
  // redundant to CoordinateSqrMatrixJJJPN and could be reimplemented more
  // simply in terms of CoordinateSqrMatrixJJJPN by definition of
  // appropriate "radial" matrix realizations for the various angular
  // momentum operators.

  // enum class AngularMomentumOperatorFamily {kOrbital,kSpin,kTotal};
  DEPRECATED("use am::AngularMomentumOperatorType instead")
  typedef am::AngularMomentumOperatorType AngularMomentumOperatorFamily;
  // enum class AngularMomentumOperatorSpecies {kP,kN,kTotal};
  DEPRECATED("use basis::OperatorTypePN instead")
  typedef basis::OperatorTypePN AngularMomentumOperatorSpecies ;

  Eigen::MatrixXd
  AngularMomentumMatrixJJJPN(
      am::AngularMomentumOperatorType operator_family,
      basis::OperatorTypePN operator_species,
      const basis::TwoBodySectorsJJJPN::SectorType& sector,
      int A
    );
  // Populate matrix with two-body matrix elements for a squared
  // angular momentum operator, on a given JJJPN sector.
  //
  // Precondition: It is assumed that the sector is valid for a scalar
  // (and isoscalar and positive parity), i.e., is a diagonal sector.
  //
  // Only the upper triangle of the returned matrix is populated (this
  // operator only has diagonal sectors).
  //
  // Arguments:
  //   operator_family (shell::AngularMomentumOperatorFamily):
  //     identifies momentum operator type (kOrbital, kSpin, kTotal)
  //   operator_species (shell::AngularMomentumOperatorSpecies):
  //     whether operator is total or restricted (kP, kN, kTotal)
  //   sector (basis::TwoBodySectorsJJJPN::SectorType) : the sector to
  //     populate
  //   A (int): atomic mass number
  //
  // Returns:
  //   (Eigen::MatrixXd) : the matrix for this sector

  ////////////////////////////////////////////////////////////////
  // isospin operators
  ////////////////////////////////////////////////////////////////

  enum class IsospinOperatorType {kSquared, kTz, kRaising, kLowering};

  Eigen::MatrixXd
  IsospinMatrixJJJPN(
      // overlap matrix element data
      const basis::OrbitalSpaceLJPN& overlap_orbital_space,
      const basis::OrbitalSectorsLJPN& overlap_sectors,
      const basis::OperatorBlocks<double>& overlap_matrices,
      IsospinOperatorType operator_type,
      const basis::TwoBodySectorsJJJPN::SectorType& sector,
      int A
    );
  // Populate matrix with two-body matrix elements for an isospin
  // operator, on a given JJJPN sector.
  //
  // Precondition: It is assumed that the sector is valid for a scalar
  // (and positive parity).
  //
  // Only the upper triangle of the returned matrix is populated if
  // the sector is diagonal.
  //
  // Arguments:
  //   overlap_orbital_space, overlap_sectors, overlap_matrices (...):
  //      definition of pn-overlap matrix elements
  //   operator_type (shell::IsospinOperatorType): isospin operator to
  //     calculate
  //   sector (basis::TwoBodySectorsJJJPN::SectorType) : the sector to
  //     populate
  //   A (int): atomic mass number
  //
  // Returns:
  //   (Eigen::MatrixXd) : the matrix for this sector

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
