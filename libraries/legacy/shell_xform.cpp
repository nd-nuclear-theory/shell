/****************************************************************
  shell_xform.cpp

  Mark A. Caprio, University of Notre Dame.

****************************************************************/


#include <sstream>

#include "legacy/shell_xform.h"

using namespace std;

namespace legacy {

  void ReadRadialTransformation(istream& xs, TransformationContainer& transformation)
  {
    // file parameters
    // Note: currently coding only supports Mp = M (square coefficient array)
    // for convenience of reusing PairLookupArray instead of, say, an external
    // library for general array arithmetic
    int l_max;  // max l stored
    int dl;
    int Mp; // dimension of transformed radial space (np = 0..Mp-1)
    int M; // dimension of oscillator radial space (n = 0..M-1)

    // read file parameters
    string line;
    getline(xs, line);
    istringstream(line) >> l_max >> dl >> Mp >> M;
    if (Mp != M)
      {
	std::cerr << "ReadRadialTransformation: transformation coefficient arrays" 
		  << " are presently required to be square" << std::endl;
	std::exit(EXIT_FAILURE);
      }

    // read in l-spaces
    transformation.resize(l_max+1);
    for (int l = 0; l <= l_max; ++l)
      {
	transformation[l].Initialize(kSquare,M,0.);
	for (int np = 0; np < Mp; ++np)
	  for (int n = 0; n < M; ++n)
	    xs >> transformation[l](np,n);
      }
  }


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
