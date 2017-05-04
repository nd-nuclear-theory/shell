/****************************************************************
  relative_me.h

  Construction of constant and kinematic (r^2 and k^2) relative
  operators.

  (The identity operator is provided in basis/lsjt_operator.h.)

  Language: C++11
                                 
  Mark A. Caprio
  University of Notre Dame

  7/12/16 (mac): Created (as construct_relative).
  3/26/17 (mac): Add ConstructCoulombOperator.
  3/28/17 (mac): Remove ConstructDiagonalConstantOperator.  Rename to relative_me.

****************************************************************/

#ifndef RELATIVE_ME_H_
#define RELATIVE_ME_H_

#include "basis/lsjt_operator.h"

namespace relative {

  enum class KinematicOperator {kKSqr,kRSqr};
  void ConstructKinematicOperator(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::MatrixVector,3>& relative_component_matrices,
      relative::KinematicOperator kinematic_operator
    );
  // Construct simple "kinematic" (r^2 and k^2) operators in relative
  // LSJT basis.
  //
  // Makes use of analytic matrix elements for r^2 and k^2 in
  // oscillator basis.
  //
  // See notes on "internal representation of an operator in JT
  // scheme" in lsjt_operator.h for the general principles of how the
  // operators are represented.
  //
  // Arguments:
  //   operator_labels (input) : tensorial properties of operator
  //   relative_space (input) : target space
  //   relative_component_sectors (output) : target sectors
  //   relative_component_matrices (output) : target matrices
  //   kinematic_operator (input): which operator to construct

  void ConstructCoulombOperator(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::MatrixVector,3>& relative_component_matrices,
      int num_steps
    );
  // Construct coulomb operator in relative LSJT basis.
  //
  // Makes use of spline integration.  An offset inner integration
  // point x0>0 is required to prevent the integrand from exploding.
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
  //   num_steps (input): number of integration steps

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
