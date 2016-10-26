/****************************************************************
  separable.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "tbme/separable.h"

namespace shell {
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  Eigen::MatrixXd 
  IdentityOperatorSectorMatrixJJJPN(
      const basis::TwoBodySectorsJJJPN::SectorType& sector
    )
  {

    assert(sector.IsDiagonal());

    // set up aliases
    const basis::TwoBodySubspaceJJJPN& bra_subspace = sector.bra_subspace();
    const basis::TwoBodySubspaceJJJPN& ket_subspace = sector.ket_subspace();

    // generate matrix for sector
    Eigen::MatrixXd matrix;
    const double normalization_factor = 2.;
    matrix = normalization_factor*Eigen::MatrixXd::Identity(bra_subspace.size(),ket_subspace.size());

    return matrix;
  }
    

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
