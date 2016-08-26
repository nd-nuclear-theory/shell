/****************************************************************
  shell_indexing_nl.cpp                       

  Defines single-particle index in (Nl) scheme.
                                  
  Mark A. Caprio, University of Notre Dame.
  Last modified 4/25/15.

****************************************************************/

#include <shell/shell_indexing_nl.h>

namespace shell {

  ////////////////////////////////////////////////////////////////
  // SPOrbitalNl implementation
  ////////////////////////////////////////////////////////////////

  using namespace std;

  // declare cache
  //   definition of vector outside class apparently necessary (under G++ 3.3)
  std::vector<SPLabelsNl> SPOrbitalNl::label_cache_;
  int SPOrbitalNl::max_cached_N_ = -1;

  // output operator

  ostream& operator<< (ostream& os, const SPOrbitalNl& orbital)
  {
    os << "(" << orbital.GetN() << "," << orbital.Getl() << ")";
	
    return os;
  }

  ////////////////////////////////////////////////////////////////
  // (Nl)^2 space -- fixed (N12, L)
  ////////////////////////////////////////////////////////////////

  TwoBodyStateSetNl TwoBodySpaceNL(int N, int L)
  {
    TwoBodyStateNl state;
    TwoBodyStateSetNl state_set;

    for (int N1 = 0; N1 <= N; ++N1)
      for (int l1 = (N1 % 2); l1 <= N1; l1 += 2)
	{
	  int N2 = N - N1;
	  int p = (N2 + l1 + L) % 2;
	  int l2_min = abs(l1 - L) + p;
	  int l2_max = min(N2, l1 + L - p);
	  for (int l2 = l2_min; l2 <= l2_max; l2 += 2)
	    {
	      state.a1 = SPOrbitalNl(N1,l1);
	      state.a2 = SPOrbitalNl(N2,l2);
	      state.L = L;
	      state_set.push_back(state);
	    }
	}

    return state_set;
  };


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
