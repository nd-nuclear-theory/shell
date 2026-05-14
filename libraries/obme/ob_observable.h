/****************************************************************
  @file ob_observable.h

  Defines contraction of OBMEs with densities for evaluation of one-body
  observable RME.

  Language: C++11

  Mark A. Caprio and Patrick J. Fasano
  University of Notre Dame

  + 05/14/26 (mac): Created, refactoring from obscalc-ob.
****************************************************************/

#ifndef OB_OBSERVABLE_H_
#define OB_OBSERVABLE_H_

#include "obme/obme_operator.h"

namespace shell {

  double CalculateOneBodyObservableMatrixElement(
      const basis::OrbitalSpaceLJPN& space,
      const basis::OrbitalSectorsLJPN& sectors,
      const basis::OperatorBlocks<double>& blocks,
      const std::unique_ptr<shell::InOBDMEStream>& density_stream
    );
  // Contract OBMEs with OBDMEs to evaluate RME.
  //
  // Arguments:
  //   space (input): one-body space
  //   sectors (output): one-body sectors
  //   matrices (input): operator matrices
  //   density_stream (input): densities
  //
  // Returns:
  //   rme

}
#endif  // OB_OBSERVABLE_H_
