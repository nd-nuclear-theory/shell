/****************************************************************
  shell_indexing_nlj.h

  Defines single-particle (Nlj) indexing scheme and related indexing
  formulas.
                                  
  Mark A. Caprio, University of Notre Dame.

  2/22/11 (mac): Created.
  3/11/11 (mac): Enclose in namespace.
  4/14/11 (mac): 
    -- Rename *Nj to *Nlj.
    -- Add relational operators.
  4/17/11 (mac):
    -- Extract (Nl) scheme code to shell_indexing_nl.
    -- Two-body space indexing and matrix element storage. 
  4/25/11 (mac): Extract TBME definitions.
  10/20/13 (mac): Add simple function for jmax(N).
  4/25/15 (mac): 
    -- Reformat source file.
    -- Upgrade documentation.  
    -- Adde String() member function.
    -- Rename from shell_indexing.h to shell_indexing_nlj.h.

TODO: rename accessors to Google conventions, e.g., Getl() -> l(), and
in general update indexing to scheme in shell_indexing_lstj

****************************************************************/

#ifndef shell_indexing_nlj_h
#define shell_indexing_nlj_h

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <halfint/halfint.h>

namespace shell {

  ////////////////////////////////////////////////////////////////
  // Dimension and quantum number formulas
  ////////////////////////////////////////////////////////////////

  // LevelCountNlj() calculates the number of Nlj orbitals through a
  // given one-body shell cutoff.
  //
  // N: maximum oscillator quantum number (0-based)
  //
  // Returns number orbitals.

  inline
    int LevelCountNlj(int N)
  {
    return (N+1)*(N+2)/2;
  }

  // Shelljmax() calculates the maximum j arising under a given one-body
  // shell cutoff.
  //
  // N: maximum oscillator quantum number (0-based)
  //
  // Returns maximum j as a HalfInt.

  inline
    HalfInt Shelljmax(int N)
  {
    return N + HalfInt(1,2);
  }

  ////////////////////////////////////////////////////////////////
  // SPOrbitalNlj -- (Nlj) orbital class
  ////////////////////////////////////////////////////////////////

  // SPOrbitalNlj -- class for (N,l,j) orbital indexing arithmetic
  //
  // The (N,l,j) compound label for an orbital is mapped to a 0-based
  // index:
  //
  // 0 <-> (0,0,1/2) = "0s1/2"    (in spectroscopic nl_j notation)
  // 1 <-> (1,1,1/2) = "0p1/2"
  // 2 <-> (1,1,3/2) = "0p3/2"
  // ...
  //
  // Here N=2*n+l (=0,1,...) is the oscillator principal quantum number,
  // where n is the radial node number (=0,1,...), l is the orbital
  // angular momentum (l=N,N-2,...,N mod 2), and j is the total angular
  // momentum (half integer) obtained by coupling l with a spin of 1/2.
  //
  // Example construction:
  //
  //     SPOrbitalNlj();    // initialized to 0s1/2
  //     SPOrbitalNlj(2);   // initialized to 0p3/2
  //     SPOrbitalNlj(1,HalfInt(3,2));   // initialized to 0p3/2

  class SPOrbitalNlj {
  public:
    // constructors
    SPOrbitalNlj();
    explicit SPOrbitalNlj(int);  // 0-based
    SPOrbitalNlj(int, const HalfInt&);

    // label extraction
    int GetN() const;
    int Getl() const;
    HalfInt Getj() const;
    int Getn() const;

    // 0-based access
    int GetIndex() const;
    void SetIndex(int);

    // 1-based access
    int GetIndex1() const;
    void SetIndex1(int);

    // formatting
    std::string String() const;

    // incrementor
    SPOrbitalNlj& operator ++ (); // prefix ++, increments then returns reference
    SPOrbitalNlj operator ++ (int); // postfix ++, returns value then increments

  private:

    // SPLabelsNlj -- primitive label storage structure -- for internal use
    struct SPLabelsNlj {
      SPLabelsNlj(int N0, HalfInt j0) {N = N0; j = j0;};
      int N;
      HalfInt j;
    };

    // static label cache
    static std::vector<SPLabelsNlj> label_cache_;
    static int max_cached_N_;
    void ExtendCache() const;

    // index for this instance
    //   internal index is always 0-based
    int k_;
  };

  // SPOrbitalNlj implementation (inline)

  inline
    void SPOrbitalNlj::ExtendCache() const
  // ExtendCache() updates internal label cache to cover current
  //   instance's k value
  {
    // extend cache vector if necessary so labels are available
    for (; k_ >= LevelCountNlj(max_cached_N_); ++max_cached_N_)
      {
	int N = max_cached_N_ + 1;
	for (HalfInt j = HalfInt(1,2); j <= N + HalfInt(1,2); ++j)
	  {
	    label_cache_.push_back(SPLabelsNlj(N,j));
	  }
      }
  }

  // construct default 0s_1/2 orbital

  inline
    SPOrbitalNlj::SPOrbitalNlj()
  {
    // store index
    k_ = 0;

    // extend cache
    ExtendCache();
  }

  // construct orbital by index (0-based)

  inline
    SPOrbitalNlj::SPOrbitalNlj(int k)
  {
    // validate argument
    // ASSERTION: assert (0<k);
    if (k < 0)
      {
	std::cerr << "SPOrbitalNlj::SPOrbitalNlj: invalid index " << k << std::endl;
	std::exit(EXIT_FAILURE);
      }

    // store index
    k_ = k;

    // extend cache
    ExtendCache();
  }

  // construct orbital by (N,j) labels

  inline
    SPOrbitalNlj::SPOrbitalNlj(int N, const HalfInt& j)
  {
    // establish corresponding index (0-based)
    k_ = LevelCountNlj(N) - (2*N+1-TwiceValue(j))/2 -1;

    // extend cache
    ExtendCache();
  }

  // label retrieval

  inline
    int SPOrbitalNlj::GetN() const
  {
    return label_cache_[k_].N;
  }

  inline
    HalfInt SPOrbitalNlj::Getj() const
  {
    return label_cache_[k_].j;
  }

  inline
    int SPOrbitalNlj::Getl() const
  {
    int N = GetN();
    HalfInt j = Getj();
    return (TwiceValue(j)-1)/2 + (N+(TwiceValue(j)-1)/2)%2;
  }

  inline
    int SPOrbitalNlj::Getn() const
  {
    return (GetN() - Getl())/2;
  }

  inline
    int SPOrbitalNlj::GetIndex() const
  {
    return k_;
  }

  inline
    void SPOrbitalNlj::SetIndex(int i)
  {
    k_ = i;
  }

  inline
    int SPOrbitalNlj::GetIndex1() const
  {
    return k_ + 1;
  }

  inline
    void SPOrbitalNlj::SetIndex1(int i)
  {
    k_ = i - 1;
  }


  inline SPOrbitalNlj& SPOrbitalNlj::operator ++ ()
  {
    ++k_;
    ExtendCache();

    return *this;
  }

  inline SPOrbitalNlj SPOrbitalNlj::operator ++ (int)
  {
    SPOrbitalNlj x = *this;
    ++k_;
    ExtendCache();
    return x;
  }


  ////////////////////////////////////////////////////////////////
  // SPOrbitalNlj relational operators (inline)
  ////////////////////////////////////////////////////////////////
 
  inline bool operator == (const SPOrbitalNlj& a1, const SPOrbitalNlj& a2)
  {
    return a1.GetIndex() == a2.GetIndex();
  }

  inline bool operator < (const SPOrbitalNlj& a1, const SPOrbitalNlj& a2)
  {
    return a1.GetIndex() < a2.GetIndex();
  }

  inline bool operator > (const SPOrbitalNlj& a1, const SPOrbitalNlj& a2)
  {
    return !((a1==a2)||(a1<a2));
  }

  inline bool operator >= (const SPOrbitalNlj& a1, const SPOrbitalNlj& a2)
  {
    return !(a1 < a2);
  }

  inline bool operator <= (const SPOrbitalNlj& a1, const SPOrbitalNlj& a2)
  {
    return ((a1 < a2)||(a1 == a2));
  }

  inline bool operator != (const SPOrbitalNlj& a1, const SPOrbitalNlj& a2)
  {
    return !(a1 == a2);
  }


  ////////////////////////////////////////////////////////////////
  // SPOrbitalNlj stream operator
  ////////////////////////////////////////////////////////////////

  std::ostream& operator<< (std::ostream&, const SPOrbitalNlj&);

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
