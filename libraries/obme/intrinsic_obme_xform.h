/****************************************************************
  @file intrinsic_obme_xform.h

  Defines transformation of one-body matrix elements to obtain intrinsic matrix
  element.

  This transformation is defined by (14) of P. Navratil, Phys. Rev. C 104,
  064322 (2021) [https://doi.org/10.1103/PhysRevC.104.064322].

  Language: C++11

  Mark A. Caprio
  University of Notre Dame

  + 05/14/26 (mac): Created.
****************************************************************/

#ifndef OBME_INTRINSIC_OBME_XFORM_H_
#define OBME_INTRINSIC_OBME_XFORM_H_

#include "am/halfint.h"

#include "basis/basis.h"

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
  ///////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////
  //
  // OneBodyOperatorDeltaN: fundamental one-body operators organized by Delta N
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
  //    * Increasing Delta_N
  //      (Delta_N=-Delta_N_max,-Delta_N_max+2,...,Delta_N_max), either all even
  //      or all odd, as given by constraint delta_N~g0.
  //
  // The |Delta_N| of an operator can be no larger than the single-particle
  // N1max.  This maximal value is obtained by destroying a N=0 particle and
  // creating an N=N1max particle, or vice versa.  All operations we wish to
  // carry out on the matrix M (in particular, inversion) will preserve Delta_N
  // subspaces, and in the end we are only interested in operators and densities
  // appearing in (31).  Therefore, for purposes of truncating Delta_N, we may,
  // more specifically, use the one-body truncation Nmax=Nmax_op applied to the
  // orbitals on which the OBMEs and densities appearing in (31) are defined.
  // Compare Nmax_mat below, used for the intermediate step of matrix inversion.
  //
  // The truncation parameter for subspaces within the space, regardless of how
  // its value may be chosen, is specified as:
  //
  //   * Delta_N_max (int): Maximum Delta_N for orbital pair (N1, N2).
  //     Typically to be chosen as the maxinum single-particle N for the set of
  //     orbitals on which the one-body operators are defined.
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
  // The truncation parameters for states within a subspace, regardless of how
  // their values may be chosen, are specified as:
  //
  //   * N1max (int): Maximum single-particle N for orbitals individually (N1, N2).
  //   * N2max (int): Maximum two-particle N for both orbitals combined (N1+N2).
  //
  ///////////////////////////////////////////////////////////////
  
  // declarations

  class OneBodyOperatorDeltaNSubspace;
  class OneBodyOperatorDeltaNState;  // not really a "state" in the physical sense!
  class OneBodyOperatorDeltaNSpace;

  // labels

  typedef std::tuple<int> OneBodyOperatorDeltaNSubspaceLabels;
  typedef std::tuple<int,int,HalfInt,int,int,HalfInt> OneBodyOperatorDeltaNStateLabels;

  // subspace

  class OneBodyOperatorDeltaNSubspace
    : public basis::BaseSubspace<
    OneBodyOperatorDeltaNSubspace,OneBodyOperatorDeltaNSubspaceLabels,
    OneBodyOperatorDeltaNState,OneBodyOperatorDeltaNStateLabels
    >
    {

      public:

      // constructor

      OneBodyOperatorDeltaNSubspace() = default;
      // default constructor -- provided since required for certain
      // purposes by STL container classes (e.g., std::vector::resize)

      OneBodyOperatorDeltaNSubspace(int J0, int g0, int Delta_N, int N1max, int N2max);
      // Set up indexing with truncation by oscillator quanta.

      // accessors -- subspace labels
      int Delta_N() const {return std::get<0>(labels());}

      // accessors -- space-based constraint or truncation parameters
      int J0() const {return J0_;}
      int g0() const {return g0_;}
      int N1max() const {return N1max_;}
      int N2max() const {return N2max_;}

      // diagnostic strings

      std::string LabelStr() const;
      // Provide string representation of subspace labels.
      std::string DebugStr() const;
      // Dump subspace contents.

      private:

      //validation
      bool ValidLabels() const;

      // truncation
      int J0_;
      int g0_;
      int N1max_;
      int N2max_;

    };

  // state

  class OneBodyOperatorDeltaNState
    : public basis::BaseState<OneBodyOperatorDeltaNSubspace>
  {

    public:

    // pass-through constructors

    OneBodyOperatorDeltaNState(const SubspaceType& subspace, std::size_t index)
      // Construct state by index.
      : basis::BaseState<OneBodyOperatorDeltaNSubspace> (subspace, index) {}

    OneBodyOperatorDeltaNState(const SubspaceType& subspace, const StateLabelsType& state_labels)
      // Construct state by reverse lookup on labels.
      : basis::BaseState<OneBodyOperatorDeltaNSubspace> (subspace, state_labels) {}

    // pass-through accessors (from subspace)
    int J0() const {return subspace().J0();}
    int g0() const {return subspace().g0();}
    int Delta_N() const {return subspace().Delta_N();}

    // state label accessors
    int n1() const {return std::get<0>(labels());}
    int l1() const {return std::get<1>(labels());}
    HalfInt j1() const {return std::get<2>(labels());}
    int g1() const {return l1()%2;}
    int N1() const {return 2*n1()+l1();}

    int n2() const {return std::get<3>(labels());}
    int l2() const {return std::get<4>(labels());}
    HalfInt j2() const {return std::get<5>(labels());}
    int g2() const {return l2()%2;}
    int N2() const {return 2*n2()+l2();}
      

    // derived quantum numbers
    int Ntot() const {return N1()+N2();}
    // Calculate oscillator quantum number.

    // diagnostic strings
    std::string LabelStr() const;
    // Provide string representation of state labels.

    // comparison
    friend bool operator == (const OneBodyOperatorDeltaNState& a1, const OneBodyOperatorDeltaNState& a2)
    // Equality test based on labels (so permits comparison across different subspace indexings).
      {
        return (a1.labels() == a2.labels()) && (a1.subspace().labels() == a2.subspace().labels());
      }

  };

  // space

  class OneBodyOperatorDeltaNSpace
    : public basis::BaseSpace<OneBodyOperatorDeltaNSpace,OneBodyOperatorDeltaNSubspace>
  {

    public:

    // constructor

    OneBodyOperatorDeltaNSpace() = default;
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    explicit OneBodyOperatorDeltaNSpace(int J0, int g0, int Delta_N_max, int N1max, int N2max);
    // Set up indexing and weights in traditional oscillator Nmax
    // truncation.

    // accessors
    int J0() const {return J0_;}
    int g0() const {return g0_;}
    int Delta_N_max() const {return Delta_N_max_;}
    int N1max() const {return N1max_;}
    int N2max() const {return N2max_;}

    // diagnostic string
    std::string DebugStr() const;

    private:

    // truncation
    int J0_;
    int g0_;
    int Delta_N_max_;
    int N1max_;
    int N2max_;

  };

  // sectors

  class OneBodyOperatorDeltaNSectors
    : public basis::BaseSectors<OneBodyOperatorDeltaNSpace>
  {

    public:

    // constructor

    OneBodyOperatorDeltaNSectors() = default;
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    OneBodyOperatorDeltaNSectors(
        const OneBodyOperatorDeltaNSpace& space
      );
    // Enumerate diagonal sectors, as needed for M matrix.
    //
    // Arguments:
    //
    //   space (const OneBodyOperatorDeltaNSpace& space): Underlying space.

   private:

  };

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  
};      // namespace shell
#endif  // OBME_INTRINSIC_OBME_XFORM_H
