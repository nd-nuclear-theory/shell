/****************************************************************
  shell_radial_nl.h                       

  Defines radial matrix element storage and retrieval for nl orbitals.

  Uses same number of radial functions for all l sectors.  

  Input file format

    The header was realized in Mathematica to support a one-body
    truncation N1bMax, but with vast overkill on the number of radial
    matrix elements at high l, via:

      lMax = N1bMax;
      Mp = N1bMax/2 + 1;
      M = N1bMax/2 + 1;
      TableHeader = {{lMax, dl, Mp, M}};

    The header is: 

      l_max dl Mp M

      l_max: highest l space for bra
      dl: maximum difference of ket space l from bra
      Mp (=M): total number of radial states in each sector

    Confusingly, for the kets this means l runs up to l_max+dl, which
    is in general higher than l_max!

    The iteration scheme is defined in Mathematica via:

      {lp, 0, lMax}, {l, lp + Mod[dl, 2], Min[lp + dl, lMax], 2}

    This provides the matrix elements needed for all multipolarities
    up to dl, in steps of 2.

  Internal storage scheme

    vector lp = 0..l_max
      vector l = lp+(dl%2)..min<int>(lp+dl,l_max) step 2
        PairLookupArray<double>  
                                  
  Mark A. Caprio, University of Notre Dame.

  3/12/12 (mac): Created based on shell_xform.h.
  2/14/13 (mac): Added initializer for input given filename.
  2/25/13 (mac): Fig bugs in input for case of nonzero dl.
  4/25/15 (mac): 
    -- Source file reformatted.
    -- Renamed from shell_radial to shell_radial_nl.
  DEPRECATED 4/25/13 in favor of shell_radial_nlj module.

****************************************************************/

#ifndef shell_radial_nl_h
#define shell_radial_nl_h

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <shell/pair_indexing.h>
// #include <shell/shell_indexing_lj.h>

namespace shell {


  class RadialMatrices{

  public:

    typedef double MatrixElementType;

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    // None

    ////////////////////////////////////////////////////////////////
    // initialization
    ////////////////////////////////////////////////////////////////

    void Initialize (std::istream& is);
    void Initialize (const std::string& is_name);

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////
	
    // GetMatrixElement -- returns NAS ME, arguments must be canonical
    MatrixElementType GetMatrixElement (int lp, int l, int np, int n) const;


    ////////////////////////////////////////////////////////////////
    // deallocation
    ////////////////////////////////////////////////////////////////

    void Free ();

  private:

    ////////////////////////////////////////////////////////////////
    // internal storage
    ////////////////////////////////////////////////////////////////

    // dimension parameters
    int l_max_;
    int dl_;
    int n_max_;

    // matrix element storage
    typedef std::vector< std::vector< PairLookupArray< MatrixElementType > > > RadialMatrixContainer;
    RadialMatrixContainer radial_matrices_;
  };

  ////////////////////////////////////////////////////////////////
  // inline implementation
  ////////////////////////////////////////////////////////////////

  inline RadialMatrices::MatrixElementType RadialMatrices::GetMatrixElement (int lp, int l, int np, int n) const
  {
    double matrix_element;

    if (lp > l)
      {
	// canonicalize indices (lp <= l)
	matrix_element= GetMatrixElement(l,lp,n,np);
      }
    else
      {
	// lp <= l
	// validate indices
	if ( !( (0 <= lp) && ( (lp + l + dl_) % 2 == 0) && ((l - lp) <= dl_) && (l <= l_max_) 
		&& (0 <= np) && (np <= n_max_) && (0 <= n) && (n <= n_max_) ) )
	  {
	    std::cerr << __FILE__ << ":" << __LINE__ << ": "
		      << "RadialMatrices::GetMatrixElement: requested radial matrix element indices invalid"
		      << " (" << lp << "," << l << "," << np << "," << n <<")" 
		      << std::endl;
	    std::exit(EXIT_FAILURE);

	  }

	// extract matrix element
	int i = ((l - lp) - (dl_ % 2) ) / 2;
	// std::cout << " (" << lp << "," << l << "," << np << "," << n << ")" << i 
	// 	  << std::endl;
	// std::cout << " base size " << radial_matrices_.size()
	// 	  << std::endl;
	// std::cout << " lp size " << radial_matrices_[lp].size()
	// 	  << std::endl;

	matrix_element = radial_matrices_[lp][i](np,n);
      }

    return matrix_element;
  }



  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
