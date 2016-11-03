/****************************************************************
  tbme_xform.h

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

  10/30/16 (mac): Created.  Supplants code from h2xform.cpp.
     
****************************************************************/

#ifndef TBME_XFORM_H_
#define TBME_XFORM_H_

#include "eigen3/Eigen/Dense"

#include "basis/jjjpn_scheme.h"
#include "basis/operator.h"
#include "tbme/two_body_mapping.h"


namespace shell {
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  Eigen::MatrixXd TwoBodyTransformedMatrix(
      const typename basis::TwoBodySectorsJJJPN::SectorType& source_sector,
      const Eigen::MatrixXd& source_matrix,
      const typename basis::TwoBodySectorsJJJPN::SectorType& target_sector,
      const basis::OrbitalSpaceLJPN& radial_source_orbital_space,
      const basis::OrbitalSpaceLJPN& radial_target_orbital_space,
      const basis::OrbitalSectorsLJPN& radial_sectors,
      const basis::MatrixVector& radial_matrices,
      const shell::TwoBodyMapping& two_body_mapping
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
