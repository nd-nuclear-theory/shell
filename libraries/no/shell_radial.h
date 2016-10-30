/****************************************************************
  shell_radial.h                       

  Defines radial matrix element storage and retrieval for nl orbitals.
                                  
  Created by Mark A. Caprio, University of Notre Dame.

  3/12/12 (mac): Created based on shell_xform.h/cpp.
  2/14/13 (mac): Added initializer for input given filename.
  2/25/13 (mac): Fig bugs in input for case of nonzero dl.
  Last modified 2/25/13 (mac).

  [Date stamp shows vc modified 1/16/16.  See vc comment in cpp file.]

****************************************************************/

#ifndef shell_radial_h
#define shell_radial_h

#include <halfint/halfint.h>
#include <am/angular_momentum.h>
#include <shell/shell_indexing_lj.h>

// map contains the map class and utility the pair class..
// Eigen is used to store the matrices as Eigen Matrices..

#include <map>
#include <utility>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include </global/homes/c/cconsta1/eigen/Eigen/Dense>
#include </global/homes/c/cconsta1/eigen/Eigen/Eigenvalues>

namespace shell {


// Mathematica radial matrix file specification
//
// lMax = N1bMax;
// Mp = N1bMax/2 + 1;
// M = N1bMax/2 + 1;
// 
// TableHeader = {{lMax, dl, Mp, M}};
//    {lp, 0, lMax}, {l, lp + Mod[dl, 2], Min[lp + dl, lMax], 2}

// storage scheme
//   vector lp = 0..l_max
//     vector l = lp+(dl%2)..min<int>(lp+dl,l_max) step 2
//       PairLookupArray<double>  


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
	MatrixElementType GetMatrixElement (int lp, int l, int np, int n);

	int Getlmax() const {return l_max_;}

	int Getdl() const {return dl_;}

	int Getnmax() const {return n_max_;}

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
	typedef std::map <std::pair <int,int>, Eigen::MatrixXd> RadialMatrixContainer;
	RadialMatrixContainer radial_matrices_;
};

////////////////////////////////////////////////////////////////
// lj-indexed radial matrices
//
////////////////////////////////////////////////////////////////



// Note: Using class SpSpacelj and map to store the matrices..

class RadialMatriceslj{

public:

        
  ////////////////////////////////////////////////////////////////
  // type definitions for internal storage
  ////////////////////////////////////////////////////////////////

  typedef double MatrixElement;
  typedef Eigen::MatrixXd  RadialMatrixType;
  typedef std::map< std::pair< SPSpacelj, SPSpacelj>, Eigen::MatrixXd > RadialMatrixContainer;

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

  //RadialMatriceslj(int N1b_max, int nu);

  //TODO: will construct sector_list_ based on selection rules,
  // will allocate a matrix for each sector within map

  ////////////////////////////////////////////////////////////////
  // input/output
  ////////////////////////////////////////////////////////////////

  void Read(const std::string& is_name);
  void Write(const std::string& is_name);

  //////////////////////////////////////////////////////////////
  // accessors
  ////////////////////////////////////////////////////////////////
	
  // GetMatrixElement -- returns NAS ME, arguments must be canonical
  
  int Getlmax() const {return l_max_;}
 
  int Getdl() const {return dl_;}

  int Getnmax() const {return n_max_;}

  // Returns an Eigen Matrix  

  //RadialMatrixType GetRadialMatrixlj(int la, HalfInt ja, int lb, HalfInt jb);

  MatrixElement GetMatrixElementlj(int na, int la, HalfInt ja, int nb, int lb, HalfInt jb);

  ////////////////////////////////////////////////////////////////
  // deallocation
  ////////////////////////////////////////////////////////////////
  
  void Free ();

private:

///////////////////////////////////////////////////////////////
// internal storage
////////////////////////////////////////////////////////////////

  int l_max_;
  int dl_;
  int n_max_;
  //int nu_;

  RadialMatrixContainer radial_matrices_;
  
};


////////////////////////////////////////////////////////////////
// inline implementation
////////////////////////////////////////////////////////////////

inline RadialMatrices::MatrixElementType RadialMatrices::GetMatrixElement (int lp, int l, int np, int n) 
{
	double matrix_element;
	std::pair <int, int> lpair;
  lpair = std::make_pair(lp,l);
     
 	if (lp > l)
		{
			// canonicalize indices (lp <= l)
			matrix_element = GetMatrixElement(l,lp,n,np);
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
                
      matrix_element = radial_matrices_[lpair](np,n);
		}

	return matrix_element;
}


inline RadialMatriceslj::MatrixElement RadialMatriceslj::GetMatrixElementlj (int na, int la, HalfInt ja, int nb, int lb, HalfInt jb) 
{
	double matrix_element;
  std::pair <SPSpacelj, SPSpacelj> ljpair;
  ljpair = std::make_pair( SPSpacelj(la,ja), SPSpacelj(lb,jb) );     
         
	if (la > lb)
		{
			// canonicalize indices (lp <= l)
			matrix_element = GetMatrixElementlj(nb,lb,jb,na,la,ja);
		}
	else
		{
			// lp <= l
			// validate indices
			if ( !( (0 <= la) && ( (la + lb + dl_) % 2 == 0) && ((lb - la) <= dl_) && (lb <= l_max_+ dl_) 
				&& (0 <= na) && (na <= n_max_) && (0 <= nb) && (nb <= n_max_) /*&& AllowedTriangle(ja,jb,dl_)*/ ) )
				{
					std::cerr << __FILE__ << ":" << __LINE__ << ": "
				  					<< "RadialMatriceslj::GetMatrixElementlj: requested radial matrix element indices invalid"
				  					<< " (" << na << "," << la << "," << ja << "," << nb << "," << lb << "," << jb << ")" 
				  					<< std::endl;
					std::exit(EXIT_FAILURE);
				}
                
			matrix_element = radial_matrices_[ljpair](na,nb);
		}

	return matrix_element;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
} // namespace

#endif
