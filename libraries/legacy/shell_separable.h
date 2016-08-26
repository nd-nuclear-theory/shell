/****************************************************************
  shell_separable.h                       

  Defines two-body matrix elements of separable tensor-product
  operators, given OBMEs.
                                  
  Created by Mark A. Caprio, University of Notre Dame.

  5/24/15 (mac): Initiated, providing squared angular momentum
    function, modeling on functions in h2gen.cpp.
  8/26/16 (mac): Minimal patches to serve as deprecated legacy 
    library w/in shell project.

****************************************************************/

#ifndef shell_separable_h
#define shell_separable_h

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "am/halfint.h"

#include "legacy/shell_2body.h"
#include "legacy/shell_indexing_nlj.h"
#include "legacy/pair_indexing.h"

namespace legacy {

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
