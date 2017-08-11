/****************************************************************
  relative_me.h

  Construction of constant and kinematic (r^2 and k^2) relative
  operators.

  (The identity operator is provided in basis/lsjt_operator.h.)

  See notes on "internal representation of an operator in JT
  scheme" in lsjt_operator.h for the general principles of how the
  operators are represented.

  Language: C++11
                                 
  Mark A. Caprio
  University of Notre Dame

  + 7/12/16 (mac): Created (as construct_relative).
  + 3/26/17 (mac): Add ConstructCoulombOperator.
  + 3/28/17 (mac):
    - Rename to relative_me.
    - Remove ConstructDiagonalConstantOperator.
  + 7/1/17 (mac): Add ConstructQuadrupoleOperator.
****************************************************************/

#ifndef RELATIVE_ME_H_
#define RELATIVE_ME_H_

#include "basis/lsjt_operator.h"
#include "basis/proton_neutron.h"

namespace relative {

  enum class KinematicOperator {kKSqr,kRSqr};
  void ConstructKinematicOperator(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
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
  //   kinematic_operator (input): whether to construct coordinate
  //     or momentum operator

  void ConstructQuadrupoleOperator(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
      relative::KinematicOperator kinematic_operator,
      basis::TwoBodySpeciesPN operator_species
    );
  // Construct quadrupole operator in relative LSJT basis.
  //
  // Arguments:
  //   operator_labels (input) : tensorial properties of operator
  //   relative_space (input) : target space
  //   relative_component_sectors (output) : target sectors
  //   relative_component_matrices (output) : target matrices
  //   kinematic_operator (input): whether to construct coordinate
  //     or momentum operator (coopted from ConstructKinematicOperator)
  //   operator_species (input): operator for protons, neutrons, or both

  void ConstructCoulombOperator(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
      basis::TwoBodySpeciesPN operator_species,
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
  //   operator_species (input): operator for protons, neutrons, or both
  //   num_steps (input): number of integration steps

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
