/****************************************************************
  shell_xform.h                       

  Defines TBME unitary transformation based on nlj orbitals.
                                  
  Created by Mark A. Caprio, University of Notre Dame.

  9/1/11 (mac): Extracted from h2xform.c.
  3/15/12 (mac): Updated to use sector iterator.
  3/17/12 (mac): Updated to take dummy delta-l parameter in input file.
  ReadRadialTransformation is DEPRECATED in favor (soon) of using
  shell_radial generic radial ME code.

  4/25/15 (mac): Reformatted source file.
  8/26/16 (mac): Minimal patches to serve as deprecated legacy 
    library w/in shell project.

****************************************************************/

#ifndef shell_xform_h
#define shell_xform_h

#include <iostream>
#include <vector>

#include "legacy/shell_2body.h"

namespace legacy {


  typedef std::vector< PairLookupArray< double > > TransformationContainer;

  void ReadRadialTransformation(std::istream& xs, TransformationContainer& transformation);

  void TwoBodyMatrixSectorTransform (const TwoBodyMatrixNljTzJP& source_matrix, TwoBodyMatrixNljTzJP& destination_matrix, 
				     const TransformationContainer& transformation, int N1b_ho_max,
				     const SectorNljTzJP& sector
				     );

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
