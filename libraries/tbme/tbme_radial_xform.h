/****************************************************************
  tbme_radial_xform.h

  Carries out radial basis transformation on two-body matrix elements.

  Normalization convention: All matrix elements are stored as NAS
  RMEs.  These RMEs are stored under the group theory Wigner-Eckart
  normalization convention (i.e., "no dimension factor out front, just
  the Clebsch"), but, for scalar operators, note that this RME is
  equivalently, and more simply, the branched ME (with M'=M).

  (However, actually, the Wigner-Eckart normalization convention for
  these RMEs is actually irrelevant as far as similarity
  transformation is concerned.)

  Symmetrization convention: The full square matrix is populated on
  diagonal sectors.

  Mark A. Caprio
  University of Notre Dame

  + 10/30/16 (mac): Created.  Supplants code from h2xform.cpp.
  + 11/6/16 (mac): Finish implementing TwoBodyTransformedMatrix.
  + 10/20/18 (pjf): Generalize TwoBodyTransformedMatrix for
      nonscalar operators.
  + 05/09/19 (pjf): Use std::size_t for basis indices and sizes.
  + 12/22/23 (zz/mac): Rename from tbme_xform to tbme_radial_xform.

****************************************************************/

#ifndef TBME_XFORM_H_
#define TBME_XFORM_H_

#include "eigen3/Eigen/Dense"

#include "basis/jjjpn_scheme.h"
#include "basis/operator.h"
#include "tbme/tbme_mapping.h"


namespace shell {
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  Eigen::MatrixXd TwoBodyTransformedMatrix(
      // radial overlap data
      const basis::OrbitalSpaceLJPN& radial_source_orbital_space,
      const basis::OrbitalSpaceLJPN& radial_target_orbital_space,
      const basis::OrbitalSectorsLJPN& radial_sectors,
      const basis::OperatorBlocks<double>& radial_matrices,
      // two-body indexing
      const typename basis::TwoBodySectorsJJJPN::SectorType& source_sector,
      const typename basis::TwoBodySectorsJJJPN::SectorType& target_sector,
      // matrix data
      const Eigen::MatrixXd& source_matrix
    );
  // Carry out radial basis transformation on two-body matrix.
  //
  // Precondition: The given source matrix must obey the
  // symmetrization condition, i.e., full square matrices must have
  // been populated for diagonal sectors.

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
