/****************************************************************
  lenpic_relative_me.h

  Construction of LENPIC consistent Gamow-Teller operators in relative basis.

  See relative_me.h for note on oscillator length.

  Language: C++17

  Patrick J. Fasano
  University of Notre Dame

  + 05/25/22 (pjf): Created.

****************************************************************/

#ifndef LENPIC_RELATIVE_ME_H_
#define LENPIC_RELATIVE_ME_H_

#include "basis/lsjt_operator.h"
#include "basis/proton_neutron.h"

namespace relative::lenpic
{

////////////////////////////////////////////////////////////////
// N2LO Gamow-Teller operator
////////////////////////////////////////////////////////////////

void ConstructN2LOGTOperator(
    const basis::OperatorLabelsJT& operator_labels,
    const basis::RelativeSpaceLSJT& relative_space,
    std::array<basis::RelativeSectorsLSJT, 3>& relative_component_sectors,
    std::array<basis::OperatorBlocks<double>, 3>& relative_component_matrices,
    double regulator_R,
    double single_particle_b,
    std::size_t num_points=3000
  );
// Construct N2LO Gamow-Teller operator in relative LSJT basis.
//
// The matrix elements are calculated for *single-particle* oscillator length
// parameter b, which is related to the relative and center-of-mass oscillator
// lengths by b_rel=2^(1/2)*b and b_cm=2^(-1/2)*b.
//
// See notes on "internal representation of an operator in JT
// scheme" in lsjt_operator.h for the general principles of how the
// operators are represented.
//
// Arguments:
//   operator_labels (input): tensorial properties of operator
//   relative_space (input): target space
//   relative_cm_component_sectors (output): target sectors
//   relative_cm_component_matrices (output): target matrices
//   regulator_R (input): regulator "R" in fm
//   single_particle_b (input): single-particle oscillator length parameter in fm

};  // namespace relative::lenpic

#endif  // LENPIC_RELATIVE_ME_H_
