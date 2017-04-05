/******************************************************************************
  
  Created by M. A. Caprio, University of Notre Dame, 2/15/11.

  Last modified 12/24/10.

******************************************************************************/

#include "moshinsky.h"

#include <iostream>
#include <ostream>
#include <iomanip>
#include <algorithm>

#include <halfint/halfint.h>
#include <mcpp/arithmetic.h>

using std::cout;
using std::endl;
using std::setw;
using std::setprecision;
using std::fixed;
//using std::resetiosflags;
//using std::ostream;
using std::max;

int main(int argc, char **argv)
{

        ////////////////////////////////
	// spot checks
	////////////////////////////////

	cout << "Spot checks" << endl;
	// seed test
	// < 1 3 0 2 ; 4 | 0 2 0 5 ; 4 > 
	// expect -0.39086801 (TTB p. 3)
	cout << setprecision(8) << fixed;
	cout << SeedMoshinskyBracket(1,3,0,2,2,5,4) << endl;

	// CM zero case tests

	// case ket seed:
	// < 1 5 0 0 ; 5 | 0 2 0 5 ; 5 > 
	// expect -0.13176157 (TTB p. 3)
	cout << CMMoshinskyBracket(0,2,0,5,5) << endl;

	// case bra seed: ??

	// case n1=0:
	// < 3 1 0 0 ; 1 | 0 1 2 2 ; 1 > 
	// expect 0.17677670 (TTB p. 87)
	cout << CMMoshinskyBracket(0,1,2,2,1) << endl;

	// case generic:
	// < 5 2 0 0 ; 2 | 3 0 2 2 ; 2 > 
	// expect 0.22332586 (TTB p. 119)
	cout << CMMoshinskyBracket(3,0,2,2,2) << endl;

	// case invalid:
	cout << CMMoshinskyBracket(3,1,2,2,2) << endl;

	cout << "****" << endl;

        ////////////////////////////////
	// loop by N_max
	////////////////////////////////

	cout << "Normalization checks" << endl;
	
	int N_max = 50;
	for (int N = 0; N <= N_max; N++)
		for (int L = (N % 2); L <= N; L += 2)
			// for allowed L *even*
			//   for N_CM=0, only even L are allowed
		{
			cout << "N " << setw(2) << N << ", "
			     << "L " << setw(3) << L << ": ";
			
			TwoBodyStateSetNl states = TwoBodySpaceNL(N,L);
			int two_body_dim = states.size();
			double norm_sum = 0;
			double max_abs = 0;

			for (int i = 0; i < two_body_dim; i++)
			{
				double bracket = CMMoshinskyBracket(states[i]);
				norm_sum += sqr(bracket);
				max_abs = max(max_abs,fabs(bracket));
				
			}

			cout << "dim " << setw(5) << two_body_dim << ", " 
			     << "norm sum " << setprecision(14) << setw(16)  << norm_sum << ", "
			     << "max abs " << setprecision(8) << setw(10)  << max_abs 
				// << resetiosflags(ostream::floatfield)
			     << endl;

		}
					
	// termination
	return 0;
}
