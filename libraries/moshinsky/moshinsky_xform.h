/****************************************************************
  moshinsky_xform.h

  Defines Moshinsky transformation of general relative operator in
  LSJT scheme.

  See notes "Moshinsky xform for operators" (2015) for derivation.

  Language: C++11

  Mark A. Caprio
  University of Notre Dame

  + 11/13/15 (mac): Created (as program moshinsky.cpp).
  + 11/21/15 (mac): Extract interaction I/O and manipulation functions
    to interaction_lsjt.  Allow input of generic relative interaction.
  + 11/26/15 (mac): Initial running version.
  + 01/06/16 (mac): Documentation updates.
  + 07/09/16 (mac): Overhaul:
    - Revise for restructured shell directory structure.
    - Revise to use block-N structure, with updated lsjt_scheme indexing.
    - Update to group theory Wigner-Eckart convention.
    - Extract transformation routines from program file into header.
  + 07/22/16 (mac): Implement isospin branching.
  + 10/25/16 (mac): Update to use TwoBodySectorsJJJPN constructor with
    explicit Tz=0.
  + 09/27/17 (mac): Change normalization for TwoBodyMatrixJJJPN results from
    AS to NAS.
  + 04/28/18 (mac): Extract RacahReductionFactorFirstSystemGT to am package.
  + 07/31/18 (pjf):
    - Fix ket phase in branching to jjJpn scheme.
    - Sanitize naming convention for quantum numbers; rename Jp->J_bra and
      J->J_ket, etc.
  + 05/09/19 (pjf): Use std::size_t for basis indices and sizes.
  + 06/19/19 (pjf): Generalize TwoBodyMatrixJJJPN and
    TransformOperatorTwoBodyJJJTToTwoBodyJJJPN for branching to Tz0 > 0.

****************************************************************/

#ifndef MOSHINSKY_XFORM_H_
#define MOSHINSKY_XFORM_H_

#include "eigen3/Eigen/Sparse"

#include "basis/lsjt_scheme.h"
#include "basis/lsjt_operator.h"
#include "basis/jjjt_scheme.h"
#include "basis/jjjpn_scheme.h"

namespace moshinsky {

  ////////////////////////////////////////////////////////////////
  // Obtaining relative-cm matrix elements from relative
  ////////////////////////////////////////////////////////////////

  Eigen::MatrixXd
    RelativeCMMatrixLSJTN(
        const basis::RelativeSpaceLSJT& relative_space,
        const basis::RelativeSectorsLSJT& relative_sectors,
        const basis::OperatorBlocks<double>& relative_matrices,
        const typename basis::RelativeCMSectorsLSJTN::SectorType& relative_cm_sector,
        int J0, int T0, int g0,
        basis::SymmetryPhaseMode symmetry_phase_mode
      );
  // Generate JT-reduced matrix elements of an operator A^(J0,T0,g0) on
  // a single fixed-N relative-cm sector (LSJTN), from the relative
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
  // triangular) "sectors" are generated.  However, within diagonal
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
  //   relative_matrices (basis::OperatorBlocks<double>) : source matrix elements
  //     defining relative operator
  //   relative_cm_sector (basis::RelativeCMSectorsLSJTN::SectorType) : target sector
  //     for relative-cm operator
  //   J0, T0, g0 (int) : operator properties
  //   symmetry_phase_mode (basis::SymmetryPhaseMode) : specification of
  //     matrix element conjugation properties of the operator
  //
  // Returns:
  //   (Eigen::MatrixXd) : the matrix representation of this sector

  void TransformOperatorRelativeLSJTToRelativeCMLSJTN(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      const std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      const std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
      const basis::RelativeCMSpaceLSJTN& relative_cm_lsjtn_space,
      std::array<basis::RelativeCMSectorsLSJTN,3>& relative_cm_lsjtn_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_cm_lsjtn_component_matrices
    );
  // Construct relative-cm representation of operator in LSJTN basis,
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
  //   relative_cm_lsjtn_space (...) : target space
  //   relative_cm_lsjtn_component_sectors (..., output) : target sectors
  //   relative_cm_lsjtn_component_matrices (..., output) : target matrices


  ////////////////////////////////////////////////////////////////
  // Moshinsky transformation
  ////////////////////////////////////////////////////////////////

  Eigen::MatrixXd
    TransformationMatrixRelativeCMTwoBodyLSJTN(
        const basis::RelativeCMSubspaceLSJTN& relative_cm_subspace,
        const basis::TwoBodySubspaceLSJTN& two_body_subspace
      );
  // Generate Moshinsky transformation matrix in given LSJTN sector.
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
  //   relative_cm_subspace (basis::RelativeCMSubspaceLSJTN) : the relative-cm subspace
  //   two_body_subspace (basis::TwoBodySubspaceLSJTN) : the two-body subspace
  //
  // Returns:
  //   (matrix) : the transformation brackets

  Eigen::MatrixXd
    TwoBodyMatrixLSJTN(
        const basis::RelativeCMSectorsLSJTN::SectorType& relative_cm_sector,
        const basis::TwoBodySectorsLSJTN::SectorType& two_body_sector,
        const Eigen::MatrixXd& relative_cm_matrix
      );
  // Obtain two-body LSJTN sector by Moshinsky transformation on
  // relative-cm LSJTN sector.
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
  //   relative_cm_sector (basis::RelativeCMSectorsLSJTN::SectorType) :
  //     source sector information
  //   two_body_sector (basis::TwoBodySectorsLSJTN::SectorType) :
  //     target sector information
  //   relative_cm_matrix (Eigen::MatrixXd) : source sector matrix

  void TransformOperatorRelativeCMLSJTNToTwoBodyLSJTN(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeCMSpaceLSJTN& relative_cm_lsjtn_space,
      const std::array<basis::RelativeCMSectorsLSJTN,3>& relative_cm_lsjtn_component_sectors,
      const std::array<basis::OperatorBlocks<double>,3>& relative_cm_lsjtn_component_matrices,
      const basis::TwoBodySpaceLSJTN& two_body_lsjtn_space,
      std::array<basis::TwoBodySectorsLSJTN,3>& two_body_lsjtn_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& two_body_lsjtn_component_matrices
    );
  // Construct two-body representation of operator in LSJTN basis,
  // from relative-cm representation.
  //
  // See notes on "internal representation of an operator in JT
  // scheme" in lsjt_operator.h for the general principles of how the
  // operators are represented.
  //
  // Arguments:
  //   operator_labels (basis::OperatorLabelsJT) : tensorial properties of operator
  //   relative_cm_lsjtn_space (...) : source space
  //   relative_cm_lsjtn_component_sectors (...) : source sectors
  //   relative_cm_lsjtn_component_matrices (...) : source matrices
  //   two_body_lsjtn_space (...) : target space
  //   two_body_lsjtn_component_sectors (..., output) : target sectors
  //   two_body_lsjtn_component_matrices (..., output) : target matrices



  ////////////////////////////////////////////////////////////////
  // recoupling to jjJT scheme
  ////////////////////////////////////////////////////////////////

  Eigen::SparseMatrix<double>
    TransformationMatrixTwoBodyLSJTNToTwoBodyJJJTN(
        const basis::TwoBodySubspaceLSJTN& two_body_lsjtn_subspace,
        const basis::TwoBodySubspaceJJJTN& two_body_jjjtn_subspace
      );
  // Generate recoupling transformation matrix between given LSJTN and
  // JJJTN sectors.
  //
  // The state
  //
  //   |(N1,l1)(N2,l1);LSJTN>
  //
  // in a given LSJTN subspace connects only to the states
  //
  //   |(N1,l1,j1)(N2,l1,j2);NJT>
  //
  // of the same (N1,l1)(N2,l1) and various j1 and j2 in the target
  // subspace, so the resulting transformation matrix is relatively
  // sparse.  There are at most four nonzero entries per row, from the
  // coupling of l with spin to give j as Delta(l,1/2,j).  There is at
  // most one entry per column, since the target state
  // (N1,l1,j1)(N2,l2,j2) derives only from the source state
  // (N1,l1)(N2,l2).
  //
  // Therefore, we use a sparse matrix, and rely upon Eigen's
  // sparse-dense matrix multiplication.
  //
  // The overlap bracket is the *unitary* 9-J recoupling coefficient
  //
  //     [ l1 s1 j1 ]
  //     [ l2 s2 j2 ]
  //     [ L  S  J  ]
  //
  // Arguments:
  //   two_body_lsjtn_subspace (...) : the source subspace
  //   two_body_jjjtn_subspace (...) : the target subspace
  //
  // Returns:
  //   (sparse matrix) : the transformation brackets

  Eigen::MatrixXd
    TwoBodyMatrixJJJTN(
        const basis::TwoBodySpaceLSJTN& two_body_lsjtn_space,
        const basis::TwoBodySectorsLSJTN& two_body_lsjtn_sectors,
        const basis::OperatorBlocks<double>& two_body_lsjtn_matrices,
        const basis::TwoBodySectorsJJJTN::SectorType& two_body_jjjtn_sector,
        int J0, int T0, int g0,
        basis::SymmetryPhaseMode symmetry_phase_mode
      );
  // Obtain two-body JJJTN sector by Moshinsky transformation on
  // relative-cm LSJT sector.
  //
  // We need all source sectors since we are performing a sum over
  // source L and S labels, for both bra and ket.
  //
  // The output matrix element are antisymmetrized (AS) matrix elements,
  // rather than normalized antisymmetrized (NAS) matrix elements.
  //
  // Arguments:
  //   two_body_lsjtn_space (basis::TwoBodySpaceLSJTN) :
  //     source space
  //   two_body_lsjtn_sectors (basis::TwoBodySectorsLSJTN) :
  //     all source sectors in given isospin component
  //   two_body_lsjtn_matrices (basis::OperatorBlocks<double>) :
  //     all source matrices in given isospin component
  //   two_body_jjjtn_sector (basis::TwoBodySectorJJJTN::SectorType) :
  //     target sector information
  //   J0, T0, g0 (int) : operator properties
  //   symmetry_phase_mode (basis::SymmetryPhaseMode) : specification of
  //     matrix element conjugation properties of the operator


  void TransformOperatorTwoBodyLSJTNToTwoBodyJJJTN(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::TwoBodySpaceLSJTN& two_body_lsjtn_space,
      const std::array<basis::TwoBodySectorsLSJTN,3>& two_body_lsjtn_component_sectors,
      const std::array<basis::OperatorBlocks<double>,3>& two_body_lsjtn_component_matrices,
      const basis::TwoBodySpaceJJJTN& two_body_jjjtn_space,
      std::array<basis::TwoBodySectorsJJJTN,3>& two_body_jjjtn_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& two_body_jjjtn_component_matrices
    );
  // Recouple operator to two-body jjJT scheme representation (in
  // TwoBodyJJJTN basis), from two-body LSJT representation (in
  // TwoBodyLSJTN basis).
  //
  // See notes on "internal representation of an operator in JT
  // scheme" in lsjt_operator.h for the general principles of how the
  // operators are represented.
  //
  // Arguments:
  //   operator_labels (basis::OperatorLabelsJT) : tensorial properties of operator
  //   two_body_lsjtn_space (...) : source space
  //   two_body_lsjtn_component_sectors (...) : source sectors
  //   two_body_lsjtn_component_matrices (...) : source matrices
  //   two_body_jjjtn_space (...) : target space
  //   two_body_jjjtn_component_sectors (..., output) : target sectors
  //   two_body_jjjtn_component_matrices (..., output) : target matrices

  ////////////////////////////////////////////////////////////////
  // branching to jjJpn scheme
  ////////////////////////////////////////////////////////////////

  Eigen::MatrixXd
    TwoBodyMatrixJJJPN(
        const basis::OperatorLabelsJT& operator_labels,
        const basis::TwoBodySpaceJJJT& two_body_jjjt_space,
        const std::array<basis::TwoBodySectorsJJJT,3>& two_body_jjjt_component_sectors,
        const std::array<basis::OperatorBlocks<double>,3>& two_body_jjjt_component_matrices,
        const basis::TwoBodySectorsJJJPN::SectorType& two_body_jjjpn_sector
      );
  // Obtain two-body JJJPN sector by isospin branching, from isospin
  // component JJJT sectors.
  //
  // It is assumed that all input matrix elements are JT-reduced
  // matrix elements and output matrix elements are J-reduced matrix
  // elements, under group theoretical conventions.  Although th
  // Moshinsky transformation per se is independent of Wigner-Eckart
  // convention, the isospin branching does depend upon convention.
  //
  // The input and output matrix elements are antisymmetrized (AS)
  // matrix elements, rather than normalized antisymmetrized (NAS)
  // matrix elements.
  //
  // We need all source isospin components, since we are performing a sum over
  // source isospin components.
  //
  // Isospin branching relations (for AS matrix elements):
  //
  //   Conventions:
  //     - proton is Tz=+1/2
  //     - group theory convention for W-E theorem
  //
  //   Reference: See notes "relation of JT and pn TBME (AS)".
  //
  //   <ab;Jab||W||cd;Jcd>_pp,AS = (1,+1,T0,0|1,+1) <ab;Jab,1||W||cd;Jcd,1>_AS
  //
  //   <ab;Jab||W||cd;Jcd>_pn,AS = 1/2*[
  //         (1,0,T0,0|1,0) <ab;Jab,1||W||cd;Jcd,1>_AS
  //       + (0,0,T0,0|1,0) <ab;Jab,1||W||cd;Jcd,0>_AS
  //       + (1,0,T0,0|0,0) <ab;Jab,0||W||cd;Jcd,1>_AS
  //       + (0,0,T0,0|0,0) <ab;Jab,0||W||cd;Jcd,0>_AS
  //      ]
  //
  //   <ab;Jab||W||cd;Jcd>_nn,AS = (1,-1,T0,0|1,-1) <ab;Jab,1||W||cd;Jcd,1>_AS
  //
  // PRECONDITION: The proton and neutron orbital label sets must be
  // defined identically, and the indexing order must be
  // lexicographical in (N,j).  This precondition is relied upon when
  // (1) proton and neutrons orbital labels are exchanged to place
  // them in canonical order when looking up isospin-scheme matrix
  // elements (2) orbital indices are compared to test canonicality.
  //
  // Arguments:
  //   operator_labels (basis::OperatorLabelsJT) : tensorial properties of operator
  //   two_body_jjjt_space (...) : source space
  //   two_body_jjjt_component_sectors (...) :
  //     all source sectors, for all isospin components
  //   two_body_jjjt_component_matrices (...) :
  //     all source matrices, for all isospin components
  //   two_body_jjjpn_sector (basis::TwoBodySectorJJJPN::SectorType) :
  //     target sector information

  void TransformOperatorTwoBodyJJJTToTwoBodyJJJPN(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::TwoBodySpaceJJJT& two_body_jjjt_space,
      const std::array<basis::TwoBodySectorsJJJT,3>& two_body_jjjt_component_sectors,
      const std::array<basis::OperatorBlocks<double>,3>& two_body_jjjt_component_matrices,
      const basis::TwoBodySpaceJJJPN& two_body_jjjpn_space,
      basis::TwoBodySectorsJJJPN& two_body_jjjpn_sectors,
      basis::OperatorBlocks<double>& two_body_jjjpn_matrices,
      int Tz0 = 0
    );
  // Branch operator to two-body jjJpn scheme representation (in
  // TwoBodyJJJPN basis), from two-body JJJT representation (in
  // TwoBodyJJJT basis).
  //
  // This transformation could have been defined directly from the
  // JJJTN representation, but it is conceptually simpler to transform
  // between more closely matched source and target sectors.
  //
  // See notes on "internal representation of an operator in JT
  // scheme" in lsjt_operator.h for the general principles of how the
  // operators are represented.
  //
  // Arguments:
  //   operator_labels (basis::OperatorLabelsJT) : tensorial properties of operator
  //   two_body_jjjt_space (...) : source space
  //   two_body_jjjt_component_sectors (...) : source sectors
  //   two_body_jjjt_component_matrices (...) : source matrices
  //   two_body_jjjpn_space (...) : target space
  //   two_body_jjjpn_sectors (..., output) : target sectors
  //   two_body_jjjpn_matrices (..., output) : target matrices


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
