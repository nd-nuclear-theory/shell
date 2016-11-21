/****************************************************************
  shell_indexing_lj.h                       

  Defines single-particle (lj) space indexing scheme.

  Mark A. Caprio, University of Notre Dame.

  4/25/15 (mac): Created, modeled on shell_indexing_nlj module.

****************************************************************/

#ifndef shell_indexing_lj_h
#define shell_indexing_lj_h

#include <algorithm>
#include <cstdlib>
#include <iostream>

#include "halfint.h"

namespace shell {

  ////////////////////////////////////////////////////////////////
  // SPSpacelj -- (lj) orbital space class
  ////////////////////////////////////////////////////////////////

  // SPSpacelj -- class for (l,j) single-particle space indexing
  // arithmetic
  //
  // The (l,j) compound label for an orbital is mapped to a 0-based
  // index:
  //
  // 0 <-> (0,1/2) = "s1/2"    (in spectroscopic l_j notation)
  // 1 <-> (1,1/2) = "p1/2"
  // 2 <-> (1,3/2) = "p3/2"
  // 3 <-> (2,3/2) = "d3/2"
  // 4 <-> (2,5/2) = "d5/2"
  // ...
  //
  // Indexing formulas:
  //   k = 2*l + [(j-l) - 1/2]   // if supposed normal mathematical division
  //     = l + (2*j-1)/2         // works with integer arithmetic 
  //   l = (k+1)/2               // integer division
  //   j = Fraction(2*(k-l)+1,2)
  //
  // Example construction:
  //
  //     SPSpacelj();    // initialized to 0s1/2
  //     SPSpacelj(2);   // initialized to p3/2
  //     SPSpacelj(1,HalfInt(3,2));   // initialized to p3/2

  class SPSpacelj {
  public:
    // constructors
    SPSpacelj();
    explicit SPSpacelj(int);  // 0-based
    SPSpacelj(int, const HalfInt&);

    // label extraction
    int Getl() const;
    HalfInt Getj() const;

    // index access
    int GetIndex() const;
    void SetIndex(int);

    // formatting
    std::string String() const;

    // incrementor
    SPSpacelj& operator ++ (); // prefix ++, increments then returns reference
    SPSpacelj operator ++ (int); // postfix ++, returns value then increments

  private:
    // index for this instance
    //   internal index is always 0-based
    int k_;
  };

  // construct default s_1/2 space

  inline
    SPSpacelj::SPSpacelj()
  {
    // store index
    k_ = 0;
  }

  // construct space by index (0-based)

  inline
    SPSpacelj::SPSpacelj(int k)
  {
    // validate argument
    // ASSERTION: assert (0<k);
    if (k < 0)
      {
	std::cerr << "SPSpacelj::SPSpacelj: invalid index " << k << std::endl;
	std::exit(EXIT_FAILURE);
      }

    // store index
    k_ = k;
  }

  // construct space by (l,j) labels

  inline
    SPSpacelj::SPSpacelj(int l, const HalfInt& j)
  {
    // establish corresponding index (0-based)
    k_ = l + (TwiceValue(j)-1)/2;
  }

  // label retrieval

  inline
    HalfInt SPSpacelj::Getj() const
  {
    int jj;
    jj=2*(k_-Getl())+1;
    return HalfInt(jj,2);
  }

  inline
    int SPSpacelj::Getl() const
  {
    return (k_+1)/2;
  }

  inline
    int SPSpacelj::GetIndex() const
  {
    return k_;
  }

  inline
    void SPSpacelj::SetIndex(int i)
  {
    k_ = i;
  }

  inline SPSpacelj& SPSpacelj::operator ++ ()
  {
    ++k_;
    return *this;
  }

  inline SPSpacelj SPSpacelj::operator ++ (int)
  {
    SPSpacelj x = *this;
    ++k_;
    return x;
  }

  ////////////////////////////////////////////////////////////////
  // SPSpacelj relational operators (inline)
  ////////////////////////////////////////////////////////////////
 
  inline bool operator == (const SPSpacelj& a1, const SPSpacelj& a2)
  {
    return a1.GetIndex() == a2.GetIndex();
  }

  inline bool operator < (const SPSpacelj& a1, const SPSpacelj& a2)
  {
    return a1.GetIndex() < a2.GetIndex();
  }

  inline bool operator > (const SPSpacelj& a1, const SPSpacelj& a2)
  {
    return !((a1==a2)||(a1<a2));
  }

  inline bool operator >= (const SPSpacelj& a1, const SPSpacelj& a2)
  {
    return !(a1 < a2);
  }

  inline bool operator <= (const SPSpacelj& a1, const SPSpacelj& a2)
  {
    return ((a1 < a2)||(a1 == a2));
  }

  inline bool operator != (const SPSpacelj& a1, const SPSpacelj& a2)
  {
    return !(a1 == a2);
  }

  ////////////////////////////////////////////////////////////////
  // SPSpacelj stream operator
  ////////////////////////////////////////////////////////////////

  std::ostream& operator<< (std::ostream&, const SPSpacelj&);


  ////////////////////////////////////////////////////////////////
  // Quantum number formulas
  ////////////////////////////////////////////////////////////////

  // SPSpaceljMax() evaluates the maximum lj space arising under a given
  // one-body shell cutoff.
  //
  // N: maximum oscillator quantum number (0-based)
  //
  // Returns maximum (l,j) as an SPSpacelj.

  inline
    SPSpacelj SPSpaceljMax(int N)
  {
    return SPSpacelj(N,N+HalfInt(1,2));
  }


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
