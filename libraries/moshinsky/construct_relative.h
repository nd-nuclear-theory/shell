/****************************************************************
  construct_relative.h

  Construction of relative operators.

  (The identity operator is provided in basis/lsjt_operator.h.)

  WARNING: This header has its own namespace "relative".  It may be
  moved out of the moshinsky library someday, into a relative
  interaction/operator library.

  Language: C++11
                                 
  Mark A. Caprio
  University of Notre Dame

  7/12/16 (mac): Created.

****************************************************************/

#ifndef CONSTRUCT_RELATIVE_H_
#define CONSTRUCT_RELATIVE_H_

#include "basis/lsjt_operator.h"

namespace relative {

  ////////////////////////////////////////////////////////////////
  // relative LSJT operator construction
  ////////////////////////////////////////////////////////////////

  void ConstructDiagonalConstantOperator(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::MatrixVector,3>& relative_component_matrices,
      double c
    );
  // Construct diagonal constant operator in relative LSJT basis.  The
  // zero and identity operators are special cases.
  //
  // Only the T0=0 component is nonzero.
  //
  // See notes on "internal representation of an operator in JT
  // scheme" in lsjt_operator.h for the general principles of how the
  // operators are represented.
  //
  // Arguments:
  //   operator_labels (basis::OperatorLabelsJT) : tensorial properties of operator
  //   relative_space (...) : target space
  //   relative_component_sectors (..., output) : target sectors
  //   relative_component_matrices (..., output) : target matrices
  //   c (double) : the constant diagonal value

  enum class KinematicOperator {kKSqr,kRSqr};
  void ConstructKinematicOperator(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::MatrixVector,3>& relative_component_matrices,
      relative::KinematicOperator kinematic_operator
    );
  // Construct simple oscillator-basis operators in relative LSJT
  // basis.
  //
  // See notes on "internal representation of an operator in JT
  // scheme" in lsjt_operator.h for the general principles of how the
  // operators are represented.
  //
  // Arguments:
  //   operator_labels (basis::OperatorLabelsJT) : tensorial properties of operator
  //   relative_space (...) : target space
  //   relative_component_sectors (..., output) : target sectors
  //   relative_component_matrices (..., output) : target matrices
  //   relative::KinematicOperator kinematic_operator : which operator
  //     to construct

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
