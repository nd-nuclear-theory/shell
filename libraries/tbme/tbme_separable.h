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

****************************************************************/

#ifndef TBME_SEPARABLE_H_
#define TBME_SEPARABLE_H_

#include "eigen3/Eigen/Dense"

#include "basis/jjjpn_scheme.h"
#include "basis/nlj_operator.h"
#include "radial/radial_io.h"

namespace shell {
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

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

  enum class KinematicOperatorType {kUTSqr,kVT1T2};

  Eigen::MatrixXd 
  KinematicMatrixJJJPN(
      const basis::OrbitalSpaceLJPN& radial_orbital_space,
      const basis::OrbitalSectorsLJPN& radial_sectors,
      const basis::MatrixVector& radial_matrices,
      KinematicOperatorType kinematic_operator_type,
      shell::RadialOperatorType radial_operator_type,
      const basis::TwoBodySectorsJJJPN::SectorType& sector,
      int A
    );
  // Populate matrix with two-body matrix elements of kinematic
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
  // Case kUTSwr:
  //
  //   Populate matrix with two-body matrix elements for the
  //   *two-body* operator V_(T1.T2) obtained from a a dot product of
  //   one-body vector operators (i.e., T = r or k), on a given JJJPN
  //   sector.
  //
  //   Obtained by csbasis (55)-(60).
  //
  // Precondition: It is assumed that the sector is valid for a scalar
  // (and isoscalar and positive parity), i.e., is a diagonal sector.
  //
  // Only the upper triangle of the returned matrix is populated (this
  // operator only has diagonal sectors).
  //
  // Arguments:
  //   radial_orbital_space, radial_sectors, radial_matrices (...):
  //      definition of radial matrix elements (T^2 for UTSqr or T for VT1T2)
  //   kinematic_operator_type (KinematicOperatorType): whether UTSqr or VT1T2 matrix
  //   radial_operator_type (shell::RadialOperatorType): whether
  //     coordinate (kR) or momentum (kK) space (only affects phase factor
  //     calculation)
  //   sector (basis::TwoBodySectorsJJJPN::SectorType) : the sector to
  //     populate
  //   A (int): atomic mass number
  //
  // Returns:
  //   (Eigen::MatrixXd) : the matrix for this sector

  ////////////////////////////////////////////////////////////////
  // angular momentum operators
  ////////////////////////////////////////////////////////////////

  // Note: The implementation of AngularMomentumMatrixJJJPN is
  // redundant to KinematicMatrixJJJPN and could be reimplemented more
  // simply in terms of KinematicMatrixJJJPN by definition of
  // appropriate "radial" matrix realizations for the various angular
  // momentum operators.

  enum class AngularMomentumOperatorFamily {kOrbital,kSpin,kTotal};
  enum class AngularMomentumOperatorSpecies {kP,kN,kTotal};

  Eigen::MatrixXd 
  AngularMomentumMatrixJJJPN(
      shell::AngularMomentumOperatorFamily operator_family, 
      shell::AngularMomentumOperatorSpecies operator_species, 
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
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
