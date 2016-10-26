/****************************************************************
  separable.h

  Evaluates TBMEs of certain separable operators.

  Normalization convention: All matrix elements are stored as AS
  (rather than NAS) RMEs.  And RMEs are stored under the group theory
  Wigner-Eckart normalization convention (i.e., "no dimension factor
  out front, just the Clebsch").  For scalar operators, note that this
  RME is more simply the branched (M'=M) ME.

  Symmetrization convention: The full square matrix is populated on
  diagonal sectors. 

  Mark A. Caprio
  University of Notre Dame

  10/26/16 (mac): Created.  Supplants code from shell_2body.cpp
    (created 4/25/11) for identity operator, h2gen.cpp (created
    3/13/12) for kinematic operators, and shell_separable.cpp (created
    5/24/15) for angular momentum operators.
     
****************************************************************/

#ifndef SEPARABLE_H_
#define SEPARABLE_H_

#include "eigen3/Eigen/Core"

#include "basis/jjjpn_scheme.h"

namespace shell {
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  Eigen::MatrixXd 
  IdentityOperatorSectorMatrixJJJPN(
      const basis::TwoBodySectorsJJJPN::SectorType& sector
    );
  // Populate sector of two-body identity matrix.
  //
  // This function should only be invoked on diagonal sectors, since
  // the identity operator is a scalar operator.  The returned matrix
  // is *twice* the identity matrix, on these diagonal sectors, due to
  // the AS normalization convention.
  //
  // Arguments:
  //   sector (basis::TwoBodySectorsJJJPN::SectorType) : The sector to
  //     populate.
  //
  // Returns:
  //   (Eigen::MatrixXd) : The matrix for this sector.



  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
