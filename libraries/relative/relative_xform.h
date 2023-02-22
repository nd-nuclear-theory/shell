/****************************************************************
  relative_xform.h

  Transformation (dilation) of relative operators.

  Language: C++17

  Patrick J. Fasano
  University of Notre Dame

  + 10/13/20 (pjf): Created.

****************************************************************/

#ifndef RELATIVE_XFORM_H_
#define RELATIVE_XFORM_H_

#include "basis/lsjt_operator.h"
#include "basis/proton_neutron.h"

namespace relative {

  ////////////////////////////////////////////////////////////////
  // relative dilation
  ////////////////////////////////////////////////////////////////

  void TransformRelativeOperatorLSJT(
      const basis::RelativeOperatorParametersLSJT& source_operator_parameters,
      const basis::RelativeSpaceLSJT& source_relative_space,
      const std::array<basis::RelativeSectorsLSJT,3>& source_relative_component_sectors,
      const std::array<basis::OperatorBlocks<double>,3>& source_relative_component_blocks,
      const basis::RelativeOperatorParametersLSJT& target_operator_parameters,
      const basis::RelativeSpaceLSJT& target_relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& target_relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& target_relative_component_blocks,
      double b_ratio,
      int num_steps,
      bool verbose=false
    );
  // Transform relative operator to new basis (via numerical integration).
  //
  // Makes use of overlaps from numerical integration.
  //
  // See notes on "internal representation of an operator in JT
  // scheme" in lsjt_operator.h for the general principles of how the
  // operators are represented.
  //
  // Arguments:
  //   operator_labels (input) : tensorial properties of operator
  //   source_relative_space (input) : source space
  //   source_relative_component_sectors (input) : source sectors
  //   source_relative_component_blocks (input) : source blocks
  //   target_relative_space (input) : target space
  //   target_relative_component_sectors (output) : target sectors
  //   target_relative_component_blocks (output) : target blocks
  //   b_ratio (input): ratio of target oscillator length parameter to source
  //     oscillator length parameter
  //   num_steps (input): number of steps for integration
  //   verbose (optional): print verbose output

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif  // RELATIVE_XFORM_H_
