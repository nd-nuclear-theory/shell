/****************************************************************
  shell_radial.cpp

  Created by Mark A. Caprio, University of Notre Dame.
  Last modified 2/25/13 (mac).

  V. Constantinou. Modifying the implementation to read the new radial-files
  in the lj format.

  The standard library classes map and pair will be used to help with indexing 
  instead of the old implementation which used 

****************************************************************/

#include "shell_radial.h"
#include <sstream>

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
  MatrixElementType integral;
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

  // old implementation
	// radial_matrices_.resize(l_max+1);

	std::pair <int,int> lpair;
   
  // map <std::pair <int,int>, Eigen::MatrixXd> something;

	for (int lp = 0; lp <= l_max; ++lp)
		for (int l = lp + (dl % 2); l <= std::min<int>(l_max,lp+dl); l += 2)
			{
      	lpair = std::make_pair (lp,l);
        Eigen::MatrixXd input(M,M);  
			
				//std::cout << "Reading radial matrix " << lp << " " << l << " " << i << std::endl;
				// read matrix elements
				for (int np = 0; np <= n_max_; ++np)
					{
          	std::getline(is,line);
            std::istringstream buffer(line);

						for (int n = 0; n <= n_max_; ++n)
            	{
			        	buffer >> integral;
                input(np,n) = integral;                   
              }
					}

				radial_matrices_[lpair] = input;
			}
}

void RadialMatriceslj::Read(const std::string& is_name)
{
	// file parameters
	// Note: currently coding only supports Mp = M (square coefficient array)
	// for convenience of reusing PairLookupArray instead of, say, an external
	// library for general array arithmetic

	std::ifstream is(is_name.c_str(),std::ios::in);

	int l_max;  // max l stored
	int dl;  // max delta in l stored
	int Mp; // dimension of transformed radial space (np = 0..Mp-1)
	int M; // dimension of oscillator radial space (n = 0..M-1)

	// read file parameters

	MatrixElement matrix_element;
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

	// read in ((lp,jp),(l,j)) sectors

  std::pair <SPSpacelj, SPSpacelj> ljpair;
   
	for (int lp = 0; lp <= l_max; ++lp)
		{       
			HalfInt jp_min = abs(lp - HalfInt(1,2));  
      HalfInt jp_max = (lp + HalfInt(1,2)); 
 
		  for (int l = lp + (dl % 2); l <= std::min<int>(l_max+dl,lp+dl); l += 2)
		  	{
        	HalfInt j_min =abs(l-HalfInt(1,2));
          HalfInt j_max = l + HalfInt(1,2); 
			
					//std::cout << "Reading radial matrix " << lp << " " << l << " " << i << std::endl;
			
					for(HalfInt jp = jp_min; jp <= jp_max; jp++)
						for(HalfInt j = j_min; j <= j_max; j++)
							{
							//if((AllowedTriangle(jp,j,dl)))
							//	{
              	ljpair = std::make_pair( SPSpacelj(lp,jp), SPSpacelj(l,j) );
								Eigen::MatrixXd input(M,M); 

								// read matrix elements
								for (int np = 0; np <= n_max_; ++np)
                	{
                  	std::getline(is,line);
                    std::istringstream buffer(line);

										for (int n = 0; n <= n_max_; ++n)
                    	{
			            			buffer >> matrix_element;
                        input(np,n) = matrix_element;                   
											}
                  }

                radial_matrices_[ljpair] = input;
						//	}
							}
				}	
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


void RadialMatrices::Free ()
{
	radial_matrices_.clear();
}


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
} // namespace
