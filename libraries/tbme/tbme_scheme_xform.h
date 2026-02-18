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

  + 12/22/23 (zz): Created using functions from TTz_T_TFilter.
  + 02/04/25 (mac): Provide jjJpn to jjJTTz conversion.

****************************************************************/

#ifndef TBME_SCHEME_XFORM_H_
#define TBME_SCHEME_XFORM_H_

#include <Eigen/Dense>

#include "basis/jjjpn_operator.h"
#include "basis/jjjt_operator.h"
#include "basis/jjjttz_operator.h"
#include "basis/operator.h"

namespace shell {
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  // JJJTT <-> JJJTTz
  
  void TransformOperatorTwoBodyJJJTToTwoBodyJJJTTz(
      const basis::TwoBodySpaceJJJT& two_body_jjjt_space,
      const std::array<basis::TwoBodySectorsJJJT,3>& two_body_jjjt_component_sectors,
      const std::array<basis::OperatorBlocks<double>,3>& two_body_jjjt_component_matrices,
      const basis::TwoBodySpaceJJJTTz& two_body_jjjttz_space,
      basis::TwoBodySectorsJJJTTz& two_body_jjjttz_sectors,
      basis::OperatorBlocks<double>& two_body_jjjttz_matrices
    );
  // Branch operator to two-body jjJTTz scheme representation (in TwoBodyJJJTTz
  // basis), from two-body JJJT scheme representation (in TwoBodyJJJT basis).
  //
  // Preconditions: Only "Hamiltonian-like" operators are supported, that is,
  // operators for which the only matrix elements connecting bras and kets
  // sharing the same (J, T, g, Tz) can be nonvanishing.  This means that the
  // provided JJJT operator is restricted to (J0,g0)=(0,0), matrix elements
  // between T=0 and T=1 subspaces are assumed vanishing (and ignored in the
  // conversion), and the resultant operator is taken to have Tz=0.
  //
  // Arguments:
  //   ...
  
  void TransformOperatorTwoBodyJJJTTzToTwoBodyJJJT(
      const basis::TwoBodySpaceJJJTTz& two_body_jjjttz_space,
      const basis::TwoBodySectorsJJJTTz& two_body_jjjttz_sectors,
      const basis::OperatorBlocks<double>& two_body_jjjttz_matrices,
      const basis::TwoBodySpaceJJJT& two_body_jjjt_space,
      std::array<basis::TwoBodySectorsJJJT,3>& two_body_jjjt_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& two_body_jjjt_component_matrices
    );
  // Upcouple operator to two-body JJJT scheme representation (in TwoBodyJJJT
  // basis), from two-body jjJTTz scheme representation (in TwoBodyJJJTTz
  // basis).
  //
  // Preconditions: Only "Hamiltonian-like" operators are supported, that is,
  // operators for which the only matrix elements connecting bras and kets
  // sharing the same (J, T, g, Tz) can be nonvanishing.  This means that the
  // provided JJJTTz operator is restricted to (J0,g0,Tz0)=(0,0,0), and matrix
  // elements between T=0 and T=1 subspaces are assumed vanishing (and ignored
  // in the conversion).
  //
  // Arguments: ...


  // JJJTTz <-> JJJPN
  
  //void TransformOperatorTwoBodyJJJTTzToTwoBodyJJJPN(
  //    const basis::TwoBodySpaceJJJTTz& two_body_jjjttz_space,
  //    const basis::TwoBodySectorsJJJTTz& two_body_jjjttz_sectors,
  //    const basis::OperatorBlocks<double>& two_body_jjjttz_matrices,
  //    //const basis::TwoBodySpaceJJJT& two_body_jjjt_space,
  //    //std::array<basis::TwoBodySectorsJJJT,3>& two_body_jjjt_component_sectors,
  //    //std::array<basis::OperatorBlocks<double>,3>& two_body_jjjt_component_matrices
  //  );

  void TransformOperatorTwoBodyJJJPNToTwoBodyJJJTTz(
      const basis::TwoBodySpaceJJJPN& two_body_jjjpn_space,
      const basis::TwoBodySectorsJJJPN& two_body_jjjpn_sectors,
      const basis::OperatorBlocks<double>& two_body_jjjpn_matrices,
      const basis::TwoBodySpaceJJJTTz& two_body_jjjttz_space,
      basis::TwoBodySectorsJJJTTz& two_body_jjjttz_sectors,
      basis::OperatorBlocks<double>& two_body_jjjttz_matrices
    );
  // Upcouple operator to two-body JJJTTz representation (in TwoBodyJJJTz
  // basis), from two-body jjJpn scheme representation (in TwoBodyJJJPN basis).
  //
  // Preconditions: Transformation from proton-neutron scheme back to isospin
  // scheme is only meaningful if proton and neutron orbital sets are identical
  // (both in labeling and underlying physical meaning, i.e., radial wave
  // functions).  To the extent that the orbital labels in JJJTTz representation
  // are taken to be oscillator orbital labels (N,j), the proton-neutron labels
  // must specifically be oscillator orbitals in standard (N,j) lexicographic
  // indexing.
  //
  // Arguments:
  //   ...

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
