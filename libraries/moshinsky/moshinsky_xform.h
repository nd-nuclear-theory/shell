/****************************************************************
  moshinsky_xform.h

  Defines Moshinsky transformation of general relative operator in
  LSJT scheme.

  See notes "Moshinsky xform for operators" (2015) for derivation.

  Language: C++11
                                 
  Mark A. Caprio
  University of Notre Dame

  11/13/15 (mac): Created (as program moshinsky.cpp).
  11/21/15 (mac): Extract interaction I/O and manipulation functions
    to interaction_lsjt.  Allow input of generic relative interaction.
  11/26/15 (mac): Initial running version.
  1/6/16 (mac): Documentation updates.
  7/4/16 (mac): Overhaul:
    - Revise for restructured shell directory structure.
    - Revise to use block-N structure, with updated lsjt_scheme indexing.
    - Update to group theory Wigner-Eckart convention.
    - Extract transformation routines from program file into header.

****************************************************************/

#ifndef MOSHINSKY_XFORM_H_
#define MOSHINSKY_XFORM_H_

#include "am/wigner_gsl.h"
#include "basis/lsjt_scheme.h"
#include "basis/lsjt_operator.h"

namespace moshinsky {

////////////////////////////////////////////////////////////////
// Racah reduction factor
////////////////////////////////////////////////////////////////

double RacahReductionFactorFirstSystemGT(
    const HalfInt& J1p, const HalfInt& Jp, 
    const HalfInt& J1, const HalfInt& J, 
    const HalfInt& J2, const HalfInt& J0
  );
// Prefactor for first-system operator in Racah two-system 
// reduction formula.
//
// Follows group theoretic (Rose) conventions for RME normalization:
//
//   <J1',J2;J||A_1^J0||J1,J2;J>
//     =(-)^(J1'+J2+J+J0)*Hat(J1')*Hat(J)*{J1',J',J2;J,J1,J0} 
//
// Note: Under Edmonds convention, the Hat(J1') would be a Hat(J').
//
// Arguments (all HalfInt):
//   J1p, Jp: bra first-system and total a.m.
//   J1, J: ket first-system and total a.m.
//   J2: bra/ket second-system a.m.
//   J0: operator a.m.
//   

////////////////////////////////////////////////////////////////
// Moshinsky transformation
////////////////////////////////////////////////////////////////

Eigen::MatrixXd
MoshinskyMatrixNLSJT(
    const basis::RelativeCMSubspaceNLSJT& relative_cm_subspace,
    const basis::TwoBodySubspaceNLSJT& two_body_subspace
  );
// Generate Moshinsky transformation matrix in given NLSJT sector.
//
// The result consists of Moshinsky brackets.
//
// The spectator SJT labels for the subspaces are assumed to be identical.
//
// Arguments:
//   relative_cm_subspace (basis::RelativeCMSubspaceNLSJT) : the relative-cm subspace
//   two_body_subspace (basis::TwoBodySubspaceNLSJT) : the two-body subspace
//
// Returns:
//   (matrix) : the transformation brackets


Eigen::MatrixXd
RelativeCMMatrixNLSJT(
    const basis::RelativeSpaceLSJT& relative_space,
    const basis::RelativeSectorsLSJT& relative_sectors,
    const basis::MatrixVector& relative_matrices,
    const typename basis::RelativeCMSectorsNLSJT::SectorType& relative_cm_sector,
    int J0
  );
// Generate JT-reduced matrix elements of an operator A^(J0,T0,g0) on
// a single fixed-N relative-cm sector (NLSJT), from the relative
// JT-reduced matrix elements.
//
// Since we are dealing with JT-reduced matrix elements, we work with
// each individual isospin component (T0) of the operator separately.
//
// This function implements the relation
//
//   <Nr',lr',Nc',lc';L',S',J',T',g'||A||Nr,lr,Nc,lc;L,S,J,T,g>
//     = Z(lr',lc,S';L',Jr';J') * Z(lr,lc,S;L,J';J)
//       * R1(Jr',lc,J';Jr,lc,J;J0)
//       *  <Nr';lr',S',Jr',T',gr'||A_r||Nr;lr,S,Jr,T,gr>
//
//  with a delta on (Nc',lc')=(Nc,lc), where the Z factors are
//  (12)3-(13)2 unitary recoupling coefficients, and R1 is the Racach
//  reduction formula factor for a spherical tensor operator which
//  only acts on the "first system".
//
// Arguments:
//   relative_space (basis::RelativeSpaceLSJT) : relative space
//     on which relative operator is defined
//   relative_sectors (basis::RelativeSectorsLSJT) : source sectors 
//     on which relative operator is defined
//   relative_matrices (basis::MatrixVector) : source matrix elements 
//     defining relative operator
//   relative_cm_sector (basis::RelativeCMSectorsNLSJT::SectorType) : target sector
//     for relative-cm operator
//   J0 (int) : the angular momentum carried by the operator (needed
//     in recoupling coefficient)
//
// Returns:
//   (Eigen::MatrixXd) : the matrix representation of this sector

Eigen::MatrixXd 
TransformedSector(
    const basis::TwoBodySubspaceLSJT& two_body_subspace2,
    const basis::TwoBodySubspaceLSJT& two_body_subspace1,
    const basis::RelativeSpaceLSJT& relative_space,
    const basis::RelativeSectorsLSJT& relative_sectors,
    const basis::MatrixVector& relative_matrices,
    int J0
  );
// Carries out Moshinsky transform to generate the reduced matrix
// elements on a LSTJ-coupled two-body sector.  
// 
// Important: The input relative matrices must contain *reduced*
// matrix elements, not plain matrix elements.  The present code is
// oriented towards the tranformation of general spherical tensor
// relative operators, for which it is appropriate to work with
// reduced metrix elements.  In contrast, for simple scalar
// (interaction-type) relative operators, common convention is to
// specify the relative interaction in terms of its plain, unreduced
// matrix elements (which are M-independent) in the relative subspace.
//
// The output matrix element are antisymmetrized (AS) matrix elements,
// rather than normalized antisymmetrized (NAS) matrix elements (on
// the one hand) or distinguishable-particle matrix elements (on the
// other).
//
// Note inclusion of factor of 2 relative two "naive" application of
// Moshinsky transform on distinguishable-particle states, arising
// from the bookkeeping since our two-body states are antisymmetrized,
// rather than distinguishable-particle states.


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
