/****************************************************************
  shell_separable.h                       

  Defines two-body matrix elements of separable tensor-product
  operators, given OBMEs.
                                  
  Created by Mark A. Caprio, University of Notre Dame.

  5/24/15 (mac): Initiated, providing squared angular momentum
  function, modeling on functions in h2gen.cpp.

****************************************************************/

#ifndef shell_separable_h
#define shell_separable_h

#include <shell/shell_2body.h>
#include <shell/shell_indexing_nlj.h>
#include <shell/pair_indexing.h>

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <halfint/halfint.h>
#include <am/angular_momentum.h>

namespace shell {

  ////////////////////////////////////////////////////////////////
  // (Nlj)^2 matrix sector operations
  ////////////////////////////////////////////////////////////////

  typedef int AngularMomentumType;
  const AngularMomentumType kOrbital = 0;
  const AngularMomentumType kSpin = 1;
  const AngularMomentumType kTotal = 2;

  // void TwoBodyMatrixSectorAddAngularMomentum (double scale, AngularMomentumType op, TwoBodyMatrixNljTzJP& destination_matrix, TwoSpeciesStateType state_type, int J, int g);
  void TwoBodyMatrixSectorAddAngularMomentum (double scale, AngularMomentumType op, TwoSpeciesStateType operator_state_type, int A, TwoBodyMatrixNljTzJP& destination_matrix, const SectorNljTzJP& sector);

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
