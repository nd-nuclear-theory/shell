/****************************************************************
  pn_io.h

  Read Petr Navratil format relative matrix element files.

  Language: C++11
                                 
  Anna E. McCoy
  TRIUMF

  3/26/17 (aem): Based on pn_io.h in SpNCCI (created 5/2/18).

****************************************************************/

// Petr Navratil relative interaction file vrel_ccm.int

// Notes on Petr Navratil _ccm.int relative file format:
//
// matrix elements and labels are given by 
//    n' l' n l s j tz rme
// where n is the radial quantum number (N=2n+l).
// t is determined by l+s+t odd antisymmetry constraint
// tz=+1 for neutron-neutron, 0 proton-neutron and -1 proton-proton

#ifndef PN_IO_H_
#define PN_IO_H_

#include "basis/lsjt_operator.h"

namespace relative
{
  
  void ReadPNOperatorPN(
        const std::string& source_filename,
        const basis::RelativeSpaceLSJT& relative_space,
        const basis::RelativeOperatorParametersLSJT& operator_parameters,
        const std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
        std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
        bool verbose
      );
  // Read a single isoscalar Hamiltonian-like relative operator in Petr Navartil
  // format.
  //
  // The input operator is represented in terms of non-reduced matrix
  // elements between good-JT (Tz=0) states, which, for an isoscalar
  // operator, are identical to group-theory convention JT-reduced
  // matrix elments.
  //
  //
  // Arguments:
  //   source_filename (input): input filename
  //   relative_space (input): target space
  //   operator_labels (input): operator labels
  //   relative_component_sectors (input): target sectors
  //   relative_component_matrices (output): target matrices to be populated
  //   verbose (input, optional): verbosity

  ////////////////////////////////////////////////////////////////


} // namespace

#endif
