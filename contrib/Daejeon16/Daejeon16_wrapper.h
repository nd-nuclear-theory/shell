/****************************************************************
  Daejeon16_wrapper.h

  Transformation (dilation) of relative operators.

  Language: C++17

  Patrick J. Fasano
  University of Notre Dame
  SPDX-License-Identifier: MIT

  + 10/13/20 (pjf): Created.

****************************************************************/

#ifndef DAEJEON16_WRAPPER_H_
#define DAEJEON16_WRAPPER_H_

#include "basis/lsjt_operator.h"
#include "basis/proton_neutron.h"

namespace contrib {

  void ConstructDaejeon16Operator(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_blocks
    );
  // Construct Daejeon16 operator in relative LSJT basis.
  //
  // These matrix elements are calculated for oscillator length
  // parameter b_rel = 2^(1/2) in the relative coordinate.  (See "Note
  // on oscillator length" at start of this header file.)
  //
  // See notes on "internal representation of an operator in JT
  // scheme" in lsjt_operator.h for the general principles of how the
  // operators are represented.
  //
  // Arguments:
  //   operator_labels (input): tensorial properties of operator
  //   relative_space (input): target space
  //   relative_component_sectors (output): target sectors
  //   relative_component_matrices (output): target matrices

}  // namespace contrib

#endif  // DAEJEON16_WRAPPER_H_