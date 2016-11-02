/****************************************************************
  tbme_separable.h

  Evaluates TBMEs of certain separable operators.

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

****************************************************************/

#ifndef TBME_SEPARABLE_H_
#define TBME_SEPARABLE_H_

#include "eigen3/Eigen/Dense"

#include "basis/jjjpn_scheme.h"
#include "basis/nlj_operator.h"

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
  //   sector (basis::TwoBodySectorsJJJPN::SectorType) : The sector to
  //     populate.
  //   A (int): atomic mass number
  //
  // Returns:
  //   (Eigen::MatrixXd) : The matrix for this sector.

  ////////////////////////////////////////////////////////////////
  // kinematic operators
  ////////////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////
  // angular momentum operators
  ////////////////////////////////////////////////////////////////

  enum class AngularMomentumOperatorFamily {kOrbital,kSpin,kTotal};
  enum class AngularMomentumOperatorSpecies {kP,kN,kTotal};

  Eigen::MatrixXd 
  AngularMomentumMatrixJJJPN(
      shell::AngularMomentumOperatorFamily operator_family, 
      shell::AngularMomentumOperatorSpecies operator_species, 
      const basis::TwoBodySectorsJJJPN::SectorType& sector,
      int A
    );
  // Populate matrix (for a given JJJPN sector) with two-body matrix
  // elements for a squared angular momentum operator.
  //
  // Precondition: It is assumed that the sector is valid for a scalar
  // (and isoscalar and positive parity), i.e., is a diagonal sector.
  //
  // Arguments:
  //   operator_family (shell::AngularMomentumOperatorFamily):
  //     identifies momentum operator type (kOrbital, kSpin, kTotal)
  //   operator_species (shell::AngularMomentumOperatorSpecies):
  //     whether operator is total or restricted (kP, kN, kTotal)
  //   sector (basis::TwoBodySectorsJJJPN::SectorType) : The sector to
  //     populate.
  //   A (int): atomic mass number
  //
  // Returns:
  //   (Eigen::MatrixXd) : The matrix for this sector.


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
