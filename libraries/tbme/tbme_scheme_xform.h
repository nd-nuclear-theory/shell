/****************************************************************
  tbme_scheme_xform.h

  Carries out scheme transformation on two-body matrix elements.

  Normalization convention: All matrix elements are stored as AS
  RMEs.  These RMEs are stored under the group theory Wigner-Eckart
  normalization convention (i.e., "no dimension factor out front, just
  the Clebsch"), but, for scalar operators, note that this RME is
  equivalently, and more simply, the branched ME (with M'=M).

  Zhou Zhou
  University of Notre Dame

  + 12/22/23 (zz): Created using function from TTz_T_TFilter.

****************************************************************/

#ifndef TBME_SCHEME_XFORM_H_
#define TBME_SCHEME_XFORM_H_

#include "eigen3/Eigen/Dense"

#include "basis/jjjt_operator.h"
#include "basis/jjjttz_operator.h"
#include "basis/operator.h"

namespace shell {
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  void TransformOperatorTwoBodyJJJTToTwoBodyJJJTTz( // assumptions: J, T, g, Tz are the same between bra and ket for a given matrix element
      const basis::TwoBodySpaceJJJT& two_body_jjjt_space,
      const std::array<basis::TwoBodySectorsJJJT,3>& two_body_jjjt_component_sectors,
      const std::array<basis::OperatorBlocks<double>,3>& two_body_jjjt_component_matrices,
      const basis::TwoBodySpaceJJJTTz& two_body_jjjttz_space,
      basis::TwoBodySectorsJJJTTz& two_body_jjjttz_sectors,
      basis::OperatorBlocks<double>& two_body_jjjttz_matrices
    );
  // Carry out radial basis transformation on two-body matrix.
  //
  // Precondition: The given source matrix must obey the
  // symmetrization condition, i.e., full square matrices must have
  // been populated for diagonal sectors.

  void TransformOperatorTwoBodyJJJTTzToTwoBodyJJJT(
      const basis::TwoBodySpaceJJJTTz& two_body_jjjttz_space,
      const basis::TwoBodySectorsJJJTTz& two_body_jjjttz_sectors,
      const basis::OperatorBlocks<double>& two_body_jjjttz_matrices,
      const basis::TwoBodySpaceJJJT& two_body_jjjt_space,
      std::array<basis::TwoBodySectorsJJJT,3>& two_body_jjjt_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& two_body_jjjt_component_matrices
    );

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
