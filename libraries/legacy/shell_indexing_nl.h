/****************************************************************
  shell_indexing_nl.h                       

  Defines single-particle index in (Nl) scheme, i.e., without spin.
                                  
  Mark A. Caprio, University of Notre Dame.

  4/17/11 (mac): Extracted from shell_indexing, originally coded 2/22/11.
  4/25/15 (mac): Reformat source file.

TODO: rename accessors to Google conventions, e.g., Getl() -> l(), and
in general update indexing to scheme in shell_indexing_lstj

****************************************************************/

#ifndef shell_indexing_nl_h
#define shell_indexing_nl_h

#include <cstdlib>
#include <iostream>
#include <vector>

#include <halfint/halfint.h>

namespace shell {

  ////////////////////////////////////////////////////////////////
  // SPOrbitalNl declaration
  ////////////////////////////////////////////////////////////////

  // TODO: revise from 1-based indexing to 0-based indexing with
  // 1-based alternate accessors


  // LevelCountNl(int N) returns number of l orbitals through shell N
  inline
    int LevelCountNl(int N)
  {
    return ( (N+4)*(N+1)/2 - (N+1)/2)/2;
  }

  // SPLabelsNl -- primative label storage structure

  //   meant primarily for internal use
  struct SPLabelsNl {
    SPLabelsNl(int N0, int l0) {N = N0; l = l0;};
    int N;
    int l;
  };

  // SPOrbitalNl -- class for (N,l) orbital indexing arithmetic

  class SPOrbitalNl {
  public:
    // constructors
    SPOrbitalNl();
    explicit SPOrbitalNl(int);
    SPOrbitalNl(int, int);

    // label extraction
    int GetN() const;
    int Getl() const;
    int Getn() const;

    // accessor for index
    int GetIndex() const;

    // incrementor
    SPOrbitalNl& operator ++ (); // prefix ++, increments then returns reference
    SPOrbitalNl operator ++ (int); // postfix ++, returns value then increments

  private:
    // static indexing convention
    static const int base_ = 1;

    // static label cache
    static std::vector<SPLabelsNl> label_cache_;
    static int max_cached_N_;
    void ExtendCache() const;

    // index for this instance
    //   internal index is always 0-based
    int k_;
  };

  // stream output

  std::ostream& operator<< (std::ostream&, const SPOrbitalNl&);

  ////////////////////////////////////////////////////////////////
  // SPOrbitalNl implementation
  ////////////////////////////////////////////////////////////////

  // ExtendCache() updates internal label cache to cover current
  //   instance's k value

  inline
    void SPOrbitalNl::ExtendCache() const
  {
    // extend cache vector if necessary so labels are available
    for (; k_ >= LevelCountNl(max_cached_N_); ++max_cached_N_)
      {
	int N = max_cached_N_ + 1;
	for (int l = N%2; l <= N; l+=2)
	  {
	    label_cache_.push_back(SPLabelsNl(N,l));
	  }
      }
  }

  // construct default 0s orbital

  inline
    SPOrbitalNl::SPOrbitalNl()
  {
    // store index
    k_ = 0;

    // extend cache
    ExtendCache();
  }

  // construct orbital by index

  inline
    SPOrbitalNl::SPOrbitalNl(int k)
  {
    // validate argument
    if (k < base_)
      {
	std::cerr << "SPOrbitalNl given index below base" << std::endl;
	std::exit(EXIT_FAILURE);
      }

    // store index
    k_ = k - base_;

    // extend cache
    ExtendCache();
  }

  // construct orbital by (N,j) labels

  inline
    SPOrbitalNl::SPOrbitalNl(int N, int l)
  {
    // establish corresponding index (0-based)
    k_ = LevelCountNl(N) - (N-l)/2 -1;

    // extend cache
    ExtendCache();
  }

  // label retrieval

  inline
    int SPOrbitalNl::GetN() const
  {
    return label_cache_[k_].N;
  }

  inline
    int SPOrbitalNl::Getl() const
  {
    return label_cache_[k_].l;
  }

  inline
    int SPOrbitalNl::Getn() const
  {
    return (GetN() - Getl())/2;
  }

  inline
    int SPOrbitalNl::GetIndex() const
  {
    return k_ + base_;
  }


  inline SPOrbitalNl& SPOrbitalNl::operator ++ ()
  {
    ++k_;
    ExtendCache();

    return *this;
  }

  inline SPOrbitalNl SPOrbitalNl::operator ++ (int)
  {
    SPOrbitalNl x = *this;
    ++k_;
    ExtendCache();
    return x;
  }

  ////////////////////////////////////////////////////////////////
  // two-body labeling in Nl scheme
  ////////////////////////////////////////////////////////////////

  struct TwoBodyStateNl {
    SPOrbitalNl a1, a2;
    int L;
  };

  typedef std::vector<TwoBodyStateNl> TwoBodyStateSetNl;

  TwoBodyStateSetNl TwoBodySpaceNL(int, int);


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
