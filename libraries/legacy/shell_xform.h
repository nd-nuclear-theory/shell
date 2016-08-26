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

****************************************************************/

#ifndef shell_xform_h
#define shell_xform_h

#include <shell/shell_2body.h>
#include <iostream>
#include <vector>


namespace shell {


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
