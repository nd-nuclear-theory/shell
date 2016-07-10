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
  7/9/16 (mac): Overhaul:
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
#include "basis/jjjt_scheme.h"

namespace moshinsky {

  ////////////////////////////////////////////////////////////////
  // Obtaining relative-cm matrix elements from relative
  ////////////////////////////////////////////////////////////////

  double RacahReductionFactorFirstSystemGT(
      const HalfInt& J1p, const HalfInt& J2p, const HalfInt& Jp, 
      const HalfInt& J1, const HalfInt& J2, const HalfInt& J, 
      const HalfInt& J0
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
  // Assertion: Arguments J2p and J2 must be equal.
  //
  // Arguments (all HalfInt):
  //   J1p, J2p, Jp: bra angular momenta
  //   J1, J2, J: ket angular momenta
  //   J0: operator angular momentum

  Eigen::MatrixXd
    RelativeCMMatrixNLSJT(
        const basis::RelativeSpaceLSJT& relative_space,
        const basis::RelativeSectorsLSJT& relative_sectors,
        const basis::MatrixVector& relative_matrices,
        const typename basis::RelativeCMSectorsNLSJT::SectorType& relative_cm_sector,
        int J0, int T0, int g0,
        basis::SymmetryPhaseMode symmetry_phase_mode
      );
  // Generate JT-reduced matrix elements of an operator A^(J0,T0,g0) on
  // a single fixed-N relative-cm sector (NLSJT), from the relative
  // JT-reduced matrix elements.
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
  // Symmetry of source matrix elements: Only canonical (upper
  // triangular) relative matrix elements are referenced.  The
  // canonicalization phase provided by CanonicalizeIndicesRelativeLSJT
  // is used to generate any noncanonical matrix elements that are
  // needed.
  //
  // Symmetry of target matrix elements: Only canonical (upper
  // triangular) "sectors" are generate.  However, within diagonal
  // sectors, the full square matrix is populated, to facilitate the
  // subsequent similarity transform.
  //
  // Since we are dealing with JT-reduced matrix elements, we work with
  // each individual isospin component (T0) of the operator separately.
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
  //    J0, T0, g0 (int) : operator properties
  //    symmetry_phase_mode (basis::SymmetryPhaseMode) : specification of
  //      matrix element conjugation properties of the operator
  //
  // Returns:
  //   (Eigen::MatrixXd) : the matrix representation of this sector

  void TransformOperatorRelativeLSJTToRelativeCMNLSJT(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      const std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      const std::array<basis::MatrixVector,3>& relative_component_matrices,
      const basis::RelativeCMSpaceNLSJT& relative_cm_nlsjt_space,
      std::array<basis::RelativeCMSectorsNLSJT,3>& relative_cm_nlsjt_component_sectors,
      std::array<basis::MatrixVector,3>& relative_cm_nlsjt_component_matrices
    );
  // Construct relative-cm representation of operator in NLSJT basis,
  // from relative representation.
  //
  // See notes on "internal representation of an operator in JT
  // scheme" in lsjt_operator.h for the general principles of how the
  // operators are represented.
  //
  // The Nmax truncation on the source relative space should be at
  // least as high as that on the target relative-cm space.
  //
  // Arguments:
  //   operator_labels (basis::OperatorLabelsJT) : tensorial properties of operator
  //   relative_space (...) : source space
  //   relative_component_sectors (...) : source sectors
  //   relative_component_matrices (...) : source matrices
  //   relative_cm_nlsjt_space (...) : target space
  //   relative_cm_nlsjt_component_sectors (..., output) : target sectors
  //   relative_cm_nlsjt_component_matrices (..., output) : target matrices


  ////////////////////////////////////////////////////////////////
  // Moshinsky transformation
  ////////////////////////////////////////////////////////////////

  Eigen::MatrixXd
    TransformationMatrixRelativeCMTwoBodyNLSJT(
        const basis::RelativeCMSubspaceNLSJT& relative_cm_subspace,
        const basis::TwoBodySubspaceNLSJT& two_body_subspace
      );
  // Generate Moshinsky transformation matrix in given NLSJT sector.
  //
  // The result consists of Moshinsky brackets, but multiplied by
  // sqrt(2.) to account for antisymmetry.  The resulting matrix
  // transforms from antisymmetry-allowed relative-cm states to
  // *antisymmetrized* (AS) two-body states, which one must recall are
  // not the same as normalized antisymmetrized (NAS) two-body states.
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
    TwoBodyMatrixNLSJT(
        const basis::RelativeCMSectorsNLSJT::SectorType& relative_cm_sector,
        const basis::TwoBodySectorsNLSJT::SectorType& two_body_sector,
        const Eigen::MatrixXd& relative_cm_matrix
      );
  // Obtain two-body LSJT sector by Moshinsky transformation on
  // relative-cm LSJT sector.
  //
  // It is assumed that all matrix elements are JT-reduced matrix
  // elements under group theoretical conventions, although this
  // Moshinsky transformation step is actually independent of
  // Wigner-Eckart convention.
  //
  // The output matrix element are antisymmetrized (AS) matrix elements,
  // rather than normalized antisymmetrized (NAS) matrix elements.
  //
  // Arguments:
  //   relative_cm_sector (basis::RelativeCMSectorNLSJT::SectorType) :
  //     source sector information
  //   two_body_sector (basis::TwoBodySectorNLSJT::SectorType) :
  //     target sector information
  //    relative_cm_matrix (Eigen::MatrixXd) : source sector matrix

  void TransformOperatorRelativeCMNLSJTToTwoBodyNLSJT(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeCMSpaceNLSJT& relative_cm_nlsjt_space,
      const std::array<basis::RelativeCMSectorsNLSJT,3>& relative_cm_nlsjt_component_sectors,
      const std::array<basis::MatrixVector,3>& relative_cm_nlsjt_component_matrices,
      const basis::TwoBodySpaceNLSJT& two_body_nlsjt_space,
      std::array<basis::TwoBodySectorsNLSJT,3>& two_body_nlsjt_component_sectors,
      std::array<basis::MatrixVector,3>& two_body_nlsjt_component_matrices
    );
  // Construct two-body representation of operator in NLSJT basis,
  // from relative-cm representation.
  //
  // See notes on "internal representation of an operator in JT
  // scheme" in lsjt_operator.h for the general principles of how the
  // operators are represented.
  //
  // Arguments:
  //   operator_labels (basis::OperatorLabelsJT) : tensorial properties of operator
  //   relative_cm_nlsjt_space (...) : source space
  //   relative_cm_nlsjt_component_sectors (...) : source sectors
  //   relative_cm_nlsjt_component_matrices (...) : source matrices
  //   two_body_nlsjt_space (...) : target space
  //   two_body_nlsjt_component_sectors (..., output) : target sectors
  //   two_body_nlsjt_component_matrices (..., output) : target matrices

  ////////////////////////////////////////////////////////////////
  // recoupling to jjJT scheme
  ////////////////////////////////////////////////////////////////

  void TransformOperatorTwoBodyNLSJTToTwoBodyNJJJT(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::TwoBodySpaceNLSJT& two_body_nlsjt_space,
      const std::array<basis::TwoBodySectorsNLSJT,3>& two_body_nlsjt_component_sectors,
      const std::array<basis::MatrixVector,3>& two_body_nlsjt_component_matrices,
      const basis::TwoBodySpaceNJJJT& two_body_njjjt_space,
      std::array<basis::TwoBodySectorsNJJJT,3>& two_body_njjjt_component_sectors,
      std::array<basis::MatrixVector,3>& two_body_njjjt_component_matrices
    );
  // Recouple operator to two-body jjJT scheme representation (in
  // TwoBodyNJJJT basis), from two-body LSJT representation (in
  // TwoBodyNLSJT basis).
  //
  // See notes on "internal representation of an operator in JT
  // scheme" in lsjt_operator.h for the general principles of how the
  // operators are represented.
  //
  // Arguments:
  //   operator_labels (basis::OperatorLabelsJT) : tensorial properties of operator
  //   two_body_nlsjt_space (...) : source space
  //   two_body_nlsjt_component_sectors (...) : source sectors
  //   two_body_nlsjt_component_matrices (...) : source matrices
  //   two_body_njjjt_space (...) : target space
  //   two_body_njjjt_component_sectors (..., output) : target sectors
  //   two_body_njjjt_component_matrices (..., output) : target matrices


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
