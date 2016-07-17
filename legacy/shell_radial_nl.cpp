/****************************************************************
  shell_radial_nl.cpp

  Mark A. Caprio, University of Notre Dame.
  Last modified 4/25/15 (mac).

****************************************************************/

#include <sstream>

#include <shell/shell_radial_nl.h>

namespace shell {

  void RadialMatrices::Initialize(std::istream& is)
  {
    // file parameters
    // Note: currently coding only supports Mp = M (square coefficient array)
    // for convenience of reusing PairLookupArray instead of, say, an external
    // library for general array arithmetic
    int l_max;  // max l stored
    int dl;  // max delta in l stored
    int Mp; // dimension of transformed radial space (np = 0..Mp-1)
    int M; // dimension of oscillator radial space (n = 0..M-1)

    // read file parameters
    std::string line;
    std::getline(is, line);
    std::istringstream(line) >> l_max >> dl >> Mp >> M;
    if (Mp != M)
      {
	std::cerr << "RadialMatrices: header " << l_max << " " << dl << " " << Mp << " " << M 
		  << " --- coefficient arrays" 
		  << " are presently required to be square" << std::endl;
	std::exit(EXIT_FAILURE);
      }

    // save file parameters
    l_max_ = l_max;
    dl_ = dl;
    n_max_ = M - 1;

    // read in (lp,l) sectors
    radial_matrices_.resize(l_max+1);
    for (int lp = 0; lp <= l_max; ++lp)
      for (int l = lp + (dl % 2); l <= std::min<int>(l_max,lp+dl); l += 2)
	{
	  // convert l to array index
	  int i = (l - lp) /2;

	  // initialize matrix
	  radial_matrices_[lp].resize(i+1);
	  //std::cout << "Reading radial matrix " << lp << " " << l << " " << i << std::endl;
	  radial_matrices_[lp][i].Initialize(kSquare,n_max_+1,0.);

	  // read matrix elements
	  for (int np = 0; np <= n_max_; ++np)
	    for (int n = 0; n <= n_max_; ++n)
	      is >> radial_matrices_[lp][i](np,n);
	}
				
  }


  void RadialMatrices::Initialize(const std::string& is_name)
  {
    // diagnostic output
    std::cout << "Reading radial matrices " << is_name << "..." << std::endl;

    // open input stream
    std::ifstream is;
    is.open(is_name.c_str());
    if (!is.is_open()) 
      {
	std::cerr << "RadialMatrices: open failed on " << is_name << std::endl;
	std::exit(EXIT_FAILURE);
      }

    // call stream version of initializer
    Initialize(is);

    // close input stream
    is.close();
  }

  void RadialMatrices::Free (){
    radial_matrices_.clear();
  }


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
