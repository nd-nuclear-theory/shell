/******************************************************************************

  moshinsky_bracket_test.cpp

  M. A. Caprio
  University of Notre Dame

******************************************************************************/


#include <cmath>
#include <iostream>
#include <ostream>
#include <iomanip>
#include <algorithm>

#include "mcpp/arithmetic.h"
#include "mcpp/profiling.h"

#include "moshinsky/moshinsky_bracket.h"


int main(int argc, char **argv)
{

  ////////////////////////////////
  // spot checks
  ////////////////////////////////

  std::cout << "Spot checks" << std::endl;
  std::cout << std::setprecision(8) << std::fixed;

  // seed test -- from CM test
  // < 1 3 0 2 ; 4 | 0 2 0 5 ; 4 > 
  // expect -0.39086801 (TTB p. 3)
  std::cout << shell::MoshinskyBracket(1,3,0,2,0,2,0,5,4) << std::endl;


  // case generic -- from CM test
  // < 5 2 0 0 ; 2 | 3 0 2 2 ; 2 > 
  // expect 0.22332586 (TTB p. 119)
  std::cout << shell::MoshinskyBracket(5,2,0,0,3,0,2,2,2) << std::endl;

  // trace_moshinsky = true;

  // < 2 1 1 1 ; 1 | 2 1 1 1 ; 1> 
  // expect 0.27500000 (TTB p. 94)
  std::cout << shell::MoshinskyBracket(2,1,1,1,2,1,1,1,1) << std::endl;
  // < 1 2 0 4 ; 2 | 2 1 1 1 ; 2> 
  // expect -0.16431679 (TTB p. 94)
  std::cout << shell::MoshinskyBracket(1,2,0,4,2,1,1,1,2) << std::endl;

  std::cout << "****" << std::endl;

  // <0 1 1 0 ; 1 | 0 0 1 1 ; 1>
  // expect -0.45643548 (TTB p. 47)
  trace_moshinsky = true;
  std::cout << shell::MoshinskyBracket(0,1,1,0,0,0,1,1,1) << std::endl;
  trace_moshinsky = false;

  std::cout << "****" << std::endl;

  ////////////////////////////////
  // loop by N_max
  ////////////////////////////////

  std::cout << "Normalization checks" << std::endl;

  int N_max = 50;
  for (int N = 0; N <= N_max; N++)
    for (int L = (N % 2); L <= N; ++L)
      {
	std::cout << "N " << std::setw(2) << N << ", "
		  << "L " << std::setw(3) << L << ": ";
			
	TwoBodyStateSetNl states = TwoBodySpaceNL(N,L);
	int two_body_dim = states.size();
	std::cout << "dim " << std::setw(5) << two_body_dim << std::endl;

	int diag_count = 0, off_diag_count = 0;
	double max_diag_error = 0., max_off_diag_error = 0.;

	Timer t;
	t.Start();
	for (int i_ket = 0; i_ket < two_body_dim; ++i_ket)
	  // for each ket
	  {
	    for (int i_ketp = i_ket; i_ketp < two_body_dim; ++i_ketp)
	      // for each ket primed
	      {
		// DEBUG: std::cerr << "[" << i_ket << " " << i_ketp <<"]" <<std::endl;

		// evaluate norm sum
		double norm_sum = 0;
		for (int i_bra = 0; i_bra < two_body_dim; ++i_bra)
		  {
		    norm_sum += shell::MoshinskyBracket(states[i_bra],states[i_ket])
		      * shell::MoshinskyBracket(states[i_bra],states[i_ketp]);
		  }
					
		// process norm sum
		if (i_ket == i_ketp)
		  // diagonal entry
		  {
		    ++diag_count;
		    max_diag_error = std::max(max_diag_error,fabs(norm_sum-1));
		  }
		else
		  // off-diagonal entry
		  {
		    ++off_diag_count;
		    max_off_diag_error = std::max(max_diag_error,fabs(norm_sum));
		  }
					
	      }
				
	  }
	t.Stop();
			
	std::cout << std::setprecision(4) << std::scientific;
	std::cout << "  " << "diagonal: " << "entries " << diag_count << ", max error " << max_diag_error << std::endl;
	std::cout << "  " << "off-diag: " << "entries " << off_diag_count << ", max error " << max_off_diag_error << std::endl;
	std::cout << std::setprecision(4) << std::fixed;
	std::cout << "     " << "time " << t.ElapsedTime() << std::endl;
	std::cout << std::endl;
			
      }

  // termination
  return 0;
}
