/****************************************************************
  shell_radial_nlj.h                       

  Defines radial matrix element storage and retrieval for nlj orbitals.
                                  
  Mark A. Caprio, University of Notre Dame.

  4/25/15 (mac): Created, based loosely on shell_radial_nl.
  TODO IMPLEMENTATION

****************************************************************/

#ifndef shell_radial_nlj_h
#define shell_radial_nlj_h

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <eigen3/Eigen/Core>

#include <shell/pair_indexing.h>
#include <shell/shell_indexing_lj.h>

namespace shell {


  class RadialMatrices{

  public:

    ////////////////////////////////////////////////////////////////
    // type definitions for internal storage
    ////////////////////////////////////////////////////////////////

    typedef double MatrixElementType;
    typedef std::pair< ShellSpacelj, ShellSpacelj >  RadialSectorType;
    typedef Eigen::MatrixXd  RadialMatrixType;

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
    RadialMatrixType& RadialMatrix(const ShellSpacelj& ljb, const ShellSpacelj& lja);

    ////////////////////////////////////////////////////////////////
    // text-mode I/O
    ////////////////////////////////////////////////////////////////

    void Read(const std::string& is_name);
    void Write(const std::string& is_name);

    ////////////////////////////////////////////////////////////////
    // deallocation
    ////////////////////////////////////////////////////////////////

    void Free ();

  private:

    ////////////////////////////////////////////////////////////////
    // internal storage
    ////////////////////////////////////////////////////////////////

    // indexing
    int N_max_;
    std::map< RadialSectorType, int > sector_indices_;

    // matrix storage 
    // Stored in vector indexed by sector index.
    std::vector< RadialMatrixType > radial_matrices_;
  };

#if 0
  class RadialMatriceslj{

  public:

        
    ////////////////////////////////////////////////////////////////
    // type definitions for internal storage
    ////////////////////////////////////////////////////////////////

    typedef  PairLookupArray< double >  RadialMatrixType;
    typedef RadialMatrixContainer std::map< std::pair< ljPair, ljPair>, RadialMatrixType >;

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    // construct with given truncation parameters
    //   N1b_max: maximum HO N quantum number for single-particle truncation
    //   nu: is power of radial operator (determines parity selection)
    //   lambda_max: determines multipolarity selection on l and j
    //
    // Constraints on (lb,jb)x(la,ja) sectors are:
    //    ...
    // Note: nu=0 supports overlaps
    RadialMatriceslj(int N1b_max, int nu);

    //TODO: will construct sector_list_ based on selection rules,
    // will allocate a matrix for each sector within map

    ////////////////////////////////////////////////////////////////
    // input/output
    ////////////////////////////////////////////////////////////////

    void Read(const std::string& is_name);
    void Write(const std::string& is_name);

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////
  
    RadialMatrixType& RadialMatrix(int lb, HalfInt jb, int la, HalfInt ja);

    ////////////////////////////////////////////////////////////////
    // deallocation
    ////////////////////////////////////////////////////////////////
  
    void Free ();

  private:

    ////////////////////////////////////////////////////////////////
    // internal storage
    ////////////////////////////////////////////////////////////////

    std::list,ljPair> sector_list_;
    RadialMatrixContainer radial_matrices_;
  };

#endif

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
