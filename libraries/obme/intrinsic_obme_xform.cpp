/****************************************************************
  intrinsic_obme_xform.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include <Eigen/Core>

#include "am/am.h"
#include "moshinsky/moshinsky_bracket.h"

#include "obme/obme_operator.h"

namespace shell
{

  ///////////////////////////////////////////////////////////////
  // Indexing for one-body intrinsic density transformation matrix
  //
  // For the M matrix (or its inverse) is defined in (31) of Navratil 2021, we
  // define an indexing which recognizes that the matrix has a block sparsity
  // structure when written in terms of blocks by Delta_N.  We adopt a symmetric
  // labeling scheme
  //
  //     M(n1, l1, j1, n2, l2, j2; n1', l1', j1', n2', l2', j2')
  //
  // for the matrix elements.  Compare original
  //
  //     M(n, l, j, n', l', j'; n1, l1, j1, n2, l2, j2)
  //
  // in (31) of Navratil 2021.
  //
  // For purposes of the indexing scheme, within the nomenclature of the basis
  // package, the row indices are "bra" state indices, and the column indices
  // are "ket" state indices.  However, physically, we are defining indexing,
  // not for a basis of states, but rather a basis of multipole operators on the
  // single-particle space (a.k.a. fundamental multipole operators).
  //
  // The overall "state" indexing is provided by (n1,l1,j1,n2,l2,j2), where:
  //
  //  - The tuples (n1,l1,j1) and (n2,l2,j2) separately index nlj orbitals.
  //
  //  - We may equivalently label these orbitals by (N1,l1,j1) and (N2,l2,j2),
  //    where N=2*n+l is the oscillator principal quantum number.
  //
  //  - The two orbitals must combine to give specified multipolarity J0 and
  //    parity grade g0 for the one-body operator.
  //
  //  - We must impose some truncation on the orbitals or, more generally,
  //    orbital pairs included.  To cover all indices appearing for a one-body
  //    operator defined on a set of orbital truncated according to N<=Nmax, we
  //    must at least include N1<=Nmax and N2<=Nmax.
  //
  // Then, the natural subspace structure for present purposes arises since the
  // Moshinsky brackets in (13) of Navratil 2021 enforce that the M matrix is
  // block diagonal in the oscillator "shell shift" Delta_N=N1-N2 characterizing
  // the fundamental one-body operators.  As for what truncation we can impose
  // on Delta_N, if the one-body operators on which we will be acting with M in
  // (31) of Navratil are defined on a given set of orbitals with some one-body
  // Nmax, we need to include at least the maximal Delta_N=Nmax arising in the
  // OBMEs for this operator, while additional Delta_N subspaces are
  // superfluous.
  //
  // Within each Delta_N subspace, we need to include orbitals at least through
  // N1<=Nmax and N2<=Nmax, as noted above.  However, we *might* want to include
  // additional orbitals, to potentially reduce truncation error in deducing the
  // inverse matrix.  From the viewpoint of evaluating the required Moshinsky
  // brackets, a truncation based on N_tot=N1+N2 would be convenient, in which
  // case N_tot=2*Nmax would be the minimum sufficient truncation.
  //
  // However, OBMEs also involve a species s (or, equivalently, tz) dependence.
  // If we include such a dependence in the orbital indexing for M, then the
  // full labeling scheme for M becomes
  //
  //     M(s1, n1, l1, j1, s2, n2, l2, j2; s1', n1', l1', j1', s2', n2', l2', j2')
  //
  // the overall "state" indexing becomes (s1,n1,l1,j1,s2,n2,l2,j2), and the
  // orbitals (s1,n1,l1,j1) and (s2,n2,l2,j2) must also now combine to give
  // specified Tz0 for the one-body operator.
  //
  ///////////////////////////////////////////////////////////////
  //
  // ## Labeling ##
  //
  // The space is naturally defined by the (J0, g0) of the one-body operator of
  // interest:
  //
  //   * J0 (int): Multipolarity.
  //
  //   * g0 (int): Grade (=0,1) for the parity P.
  //
  // subspace labels: (Delta_N)
  //
  //   * Delta_N (int): Difference in principal quantum number (N2-N1).
  //
  //   Note that Delta_N is constrained by its relation g0~Delta_N to the
  //   parity grade.
  //
  // state labels within subspace: (n1, l1, j1, n2, l2, j2)
  //
  //   * n1 (int): Radial quantum number (0,1,...).
  //   * l1 (int): Orbital angular momentum.
  //   * j1 (HalfInt): Total angular momentum.
  //   * n2 (int): Radial quantum number (0,1,...).
  //   * l2 (int): Orbital angular momentum.
  //   * j2 (HalfInt): Total angular momentum.
  //
  //   These labels exhibit some redundancy, as n2 may be deduced from the
  //   others for given Delta_N.
  //
  //   The oscillator quantum number is deduced from the
  //   n and l quantum numbers:
  //
  //   * N1 (int): Oscillator quanta (N1=2*n1+l1).
  //   * N2 (int): Oscillator quanta (N2=2*n2+l2).
  //
  ///////////////////////////////////////////////////////////////
  //
  // ## Subspaces ##
  //
  // Within a full space defined by fixed (J0, g0), and subject to
  // single-particle truncation Nmax, subspaces are ordered by:
  //
  //    * Increasing Delta_N (Delta_N=0,1,...,Nmax=Nmax_op), either all even or all odd,
  //      as given by constraint delta_N~g0.
  //
  // For purposes of truncating Delta_N, we may, more specifically, use the
  // one-body truncaiton Nmax=Nmax_op applied to the orbitals on which the OBMEs
  // and densities appearing in (31) are defined.  Compare Nmax_mat below, used
  // for the intermediate step of matrix inversion.
  //
  ///////////////////////////////////////////////////////////////
  //
  // ## States ##
  //
  // Within a subspace, the states are ordered by:
  //
  //   * Increasing N_tot=N1+N2 (Ntot=0,1,...,2*Nmax-Delta_N),
  //
  //   * Increasing N1 (N1=0,1,...,Ntot).
  //
  //   * Then N2 is determined by N2=N1+Delta_N.
  //
  //   * Increasing j1 (implies increasing l1).
  //
  //   * Increasing j2 (implies increasing l2).
  //
  // Constraints:
  //
  //   * Orbital triangularity: triangle(l1,1/2,j1), triangle(l2,1/2,j2).
  //
  //   * Orbital parity: l1+l2~g0, already enforced by Delta_N.
  //
  //   * Multipole triangularity: triangle(j1,j2,J0).
  //
  // For purposes of truncating states, we may truncate by some combination of
  // N_tot<=N_tot_max ("triangular truncation"), as motivated above for
  // convenience in truncation of the required set of Moshinsky coefficients,
  // and N<=Nmax=Nmax_mat for the orbitals ("square truncation").  The minimal
  // truncation providing the matrix elements appearing (31) is given by
  // Nmax_mat=Nmax_op and N_tot_max=2*Nmax_op.  However, we may increase these
  // cutoffs for improved accuracy in matrix inversion.
  //
  ///////////////////////////////////////////////////////////////
  

}  // namespace shell
