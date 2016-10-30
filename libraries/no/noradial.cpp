// noradial.cpp produces the radial files which document the matrix elements 
// for the operators r^{1}, r^{2} and k^{1} and k^{2}. The matrix elements 
// are taken between natural orbitals, which are linear combinations of harmonic 
// oscillator or CS functions. The natural orbitals radial functions have the form
// NO_{na, la, ja}. 

/*
  8/26/16 (mac): Minimal patches to compile under shell project.
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Eigenvalues"

#include "no/no_filename_retrieval.h"

#include "am/halfint.h"

#include "legacy/shell_radial_nl.h"


// Input for the noradial code
// path_and_radial_file is the path to a radial file (including the radial file)
// ex. /home/valentino/radial_files/radial-me-HO-r1.dat

struct NORadialParameters
{
	int nmax, nv;
	std::string species;
	std::string path_and_radial_file;
};

void ReadNORadialParameters(NORadialParameters& input)
{
	std::string input_line;

	{
		std::getline(std::cin,input_line);
		std::istringstream line_stream(input_line);
		line_stream >> input.nmax >> input.nv >> input.species;
	}

	std::cout << input.nmax << input.nv << input.species << std::endl;

	{
		std::getline(std::cin,input_line);
		std::istringstream line_stream(input_line);
		line_stream >> input.path_and_radial_file;
	}
}

// NOBasisRadialIntegralFilel produces radial files that depend only on the highest 
// value of j (recall that l is coupled to spin to produce j, and j = l + 1/2). So if say l = 1
// then j = 3/2. This function produces radial files that have the same form as the radial files
// that were used with h2utils for conventional runs with the HO basis and the CSFN basis.
// 
// The kinematic radial matrix elements in the conventional files had the form <na la| Op | nb lb>, 
// because the matrix elements depended on radial functions that had the form R_{na,la}(r;b=1)(Op is 
// an operator such as r^{1} or k^{2}). The natural orbitals radial functions have a j-dependence
// so they have the form NO_{na, la, ja}, thus the matrix elements have the form <na la ja| Op | nb lb jb>.
// To calculate the new radial files we perform a similarity transformation on the old radial matrices. If
// X = <na la| Op | nb lb> is the matrix representation in the R_{na la} basis, and if U is the matrix whose
// rows represent the natural orbitals then X' = U^{\dagger} X U, where X' is the matrix representation in the 
// NO_{na la ja} basis.
//
// For the matrix multiplication we use eigen..
// 
// In NOBasisRadialIntegralFilel, only the highest value of j is considered so we effectively change basis from
// R_{na la} to NO_{na la ja+1/2}. This allows for the new files to be used directly with h2utils without any modifications.
// This was done for testing purposes..

void NOBasisRadialIntegralFilel ( NORadialParameters& input )
{
	int max_occupied_l = (input.nmax) + (input.nv);
	std::string species_type = input.species;  

	legacy::RadialMatrices radial_matrix;
	std::ifstream open_integrals_file( (input.path_and_radial_file).c_str(),std::ios::in);   
	radial_matrix.Initialize(open_integrals_file);

	int filetype_position = (input.path_and_radial_file).find_last_of(".");
	std::string filetype = (input.path_and_radial_file).substr( filetype_position - 2 );

	std::string output_filename = "radial-me-NO-"+ input.species + "-" + filetype;

	std::ofstream output_file(output_filename);

	int lmax = radial_matrix.l_max();
	int dl = radial_matrix.dl();
	int nmax = radial_matrix.n_max();   

	Eigen::VectorXi header(4);    
   
	header[0] = max_occupied_l;
	header[1] = dl;
	header[2] = nmax + 1;
	header[3] = nmax + 1;

	output_file << header.transpose() << std::endl;
        
        std::cout << "max occupied l is: " << max_occupied_l << std::endl;
	std::cout << "The header is: " << header.transpose() << std::endl;

	for(int la=0; la <= max_occupied_l; la++)
		{
			HalfInt ja_max = la+HalfInt(1,2);

			for(int lb=(la + dl%2); lb <= std::min(la+dl,max_occupied_l+dl) ; lb+=2)
				{  
					Eigen::MatrixXd old_basis_matrix(nmax+1,nmax+1); 
					HalfInt jb_max = lb + HalfInt(1,2);

					for(int i=0; i <= nmax; i++)
						{
							for(int j=0; j <= nmax; j++)
								{                               
									old_basis_matrix(i,j) = radial_matrix.GetMatrixElement(la,lb,i,j);
								}
						}	

					if( (lb <= max_occupied_l) )
						{
							Eigen::MatrixXd left_side_laja_transformation = NaturalOrbital(la, ja_max, species_type).transpose();
							Eigen::MatrixXd right_side_lbjb_transformation = NaturalOrbital(lb, jb_max, species_type);

              std::cout << "(la,ja) = (" << la << ", " << ja_max << ")" << std::endl;
              std::cout << "(lb,jb) = (" << lb << ", " << jb_max << ")" << std::endl;
              std::cout << std::endl << left_side_laja_transformation << std::endl;
              std::cout << std::endl << right_side_lbjb_transformation << std::endl;
              std::cout << std::endl << old_basis_matrix << std::endl;

		          output_file << left_side_laja_transformation * old_basis_matrix * right_side_lbjb_transformation << std::endl;
						}
					else if(lb > max_occupied_l)
						{
            	// Note to self: The matrix that transforms the radial integrals from the right is an identity matrix because  
              // the natural orbitals involve spaces with l <= max_occupied_l only!

              Eigen::MatrixXd left_side_laja_transformation = NaturalOrbital(la, ja_max, species_type).transpose();
              Eigen::MatrixXd right_side_lbjb_transformation = Eigen::MatrixXd::Identity(nmax+1,nmax+1);
              output_file << left_side_laja_transformation * old_basis_matrix * right_side_lbjb_transformation << std::endl;
						}                                                          
				}
		}
}

// NOBasisRadialIntegralFilelj produces radial files with both l and j dependence..

void NOBasisRadialIntegralFilelj ( NORadialParameters& input )
{
	int max_occupied_l = (input.nmax) + (input.nv);
  std::cout << "max occupied l is: " << max_occupied_l << std::endl;
	std::string species_type = input.species;  

	legacy::RadialMatrices radial_matrix;
	std::ifstream open_integrals_file( (input.path_and_radial_file).c_str(),std::ios::in);   
	radial_matrix.Initialize(open_integrals_file);

	int filetype_position = (input.path_and_radial_file).find_last_of(".");
	std::string filetype = (input.path_and_radial_file).substr( filetype_position - 2 );

	std::string output_filename = "radial-me-NO-"+ input.species + "-" + filetype;

	std::ofstream output_file(output_filename);

	int lmax = radial_matrix.l_max();
	int dl = radial_matrix.dl();
	int nmax = radial_matrix.n_max();   

	Eigen::VectorXi header(4);       
	header[0] = max_occupied_l;
	header[1] = dl;
	header[2] = nmax + 1;
	header[3] = nmax + 1;

	output_file << header.transpose() << std::endl;

	for(int la=0; la <= max_occupied_l; la++)
		{
			HalfInt ja_min = abs(la-HalfInt(1,2));
      HalfInt ja_max = la+HalfInt(1,2);
         
      //std::cout << "la is: " << la << std::endl;

			for(int lb=(la + dl%2); lb <= std::min(la+dl,max_occupied_l+dl) ; lb+=2)
				{  
					Eigen::MatrixXd old_basis_matrix(nmax+1,nmax+1); 
          HalfInt jb_min = abs(lb-HalfInt(1,2));
					HalfInt jb_max = lb + HalfInt(1,2);
                
	        //std::cout << "lb is: " << lb << std::endl;

					for(int i=0; i <= nmax; i++)
						for(int j=0; j <= nmax; j++)
							{                               
								old_basis_matrix(i,j) = radial_matrix.GetMatrixElement(la,lb,i,j);
							}
		                        	

		      //std::cout << "The old matrix is: " << std::endl << old_basis_matrix << std::endl;

					for(HalfInt ja = ja_min; ja <= ja_max; ja++)
						for(HalfInt jb = jb_min; jb <= jb_max; jb++)
							{  
							//if((AllowedTriangle(ja,jb,dl)))
							//  { 
							if( (lb <= max_occupied_l) )
              	{ 
									Eigen::MatrixXd left_side_laja_transformation = NaturalOrbital(la, ja, species_type).transpose();	
    							Eigen::MatrixXd right_side_lbjb_transformation = NaturalOrbital(lb, jb, species_type);
									Eigen::MatrixXd transformed_matrix = left_side_laja_transformation * old_basis_matrix * right_side_lbjb_transformation;

    							std::cout << "(la,ja) = (" << la << ", " << ja << ")" << std::endl;
									std::cout << std::endl << left_side_laja_transformation << std::endl;
                  std::cout << "(lb,jb) = (" << lb << ", " << jb << ")" << std::endl;                  
                  std::cout << std::endl << right_side_lbjb_transformation << std::endl;
									std::cout << "Matrix in old basis: " << std::endl;
                  std::cout << std::endl << old_basis_matrix << std::endl;
                  std::cout << "Transformed matrix: " << std::endl;
                  std::cout << std::endl << transformed_matrix << std::endl;
	
		              output_file << transformed_matrix << std::endl;
                }
              else if(lb > max_occupied_l)
							  {
                  // Note to self: The matrix that transforms the radial integrals from the right is an identity matrix because  
                  // the natural orbitals involve spaces with l <= max_occupied_l only!
												
  								Eigen::MatrixXd left_side_laja_transformation = NaturalOrbital(la, ja, species_type).transpose();
        					Eigen::MatrixXd right_side_lbjb_transformation = Eigen::MatrixXd::Identity(nmax+1,nmax+1);
									Eigen::MatrixXd transformed_matrix = left_side_laja_transformation * old_basis_matrix * right_side_lbjb_transformation;

									std::cout << "(la,ja) = (" << la << ", " << ja << ")" << std::endl;
									std::cout << std::endl << left_side_laja_transformation << std::endl;
                  std::cout << "(lb,jb) = (" << lb << ", " << jb << ")" << std::endl;                  
                  std::cout << std::endl << right_side_lbjb_transformation << std::endl;
									std::cout << "Matrix in old basis: " << std::endl;
                  std::cout << std::endl << old_basis_matrix << std::endl;
                  std::cout << "Transformed matrix: " << std::endl;
                  std::cout << std::endl << transformed_matrix << std::endl;
	
		              output_file << transformed_matrix << std::endl;
                }
							//}                                                          
              }                       
				}			
		}
}

// NOBasisRadialIntegralIdentityFilelj produces radial files with lj dependence only that these
// files are the same for each different ja, jb combination. These files are to be used to test that
// the new version of h2utils works properly..

void NOBasisRadialIntegralIdentityFilelj ( NORadialParameters& input )
{
	int max_occupied_l = (input.nmax) + (input.nv);
  std::cout << "max occupied l is: " << max_occupied_l << std::endl;
	std::string species_type = input.species;  

	legacy::RadialMatrices radial_matrix;
	std::ifstream open_integrals_file( (input.path_and_radial_file).c_str(),std::ios::in);   
	radial_matrix.Initialize(open_integrals_file);

	int filetype_position = (input.path_and_radial_file).find_last_of(".");
	std::string filetype = (input.path_and_radial_file).substr( filetype_position - 2 );

	std::string output_filename = "radial-me-NO-"+ input.species + "-" + filetype;

	std::ofstream output_file(output_filename);

	int lmax = radial_matrix.l_max();
	int dl = radial_matrix.dl();
	int nmax = radial_matrix.n_max();   

	Eigen::VectorXi header(4);       
	header[0] = max_occupied_l;
	header[1] = dl;
	header[2] = nmax + 1;
	header[3] = nmax + 1;

	output_file << header.transpose() << std::endl;

	for(int la=0; la <= max_occupied_l; la++)
		{
			HalfInt ja_min = abs(la-HalfInt(1,2));
      HalfInt ja_max = la+HalfInt(1,2);
         
      //std::cout << "la is: " << la << std::endl;

			for(int lb=(la + dl%2); lb <= std::min(la+dl,max_occupied_l+dl) ; lb+=2)
				{  
					Eigen::MatrixXd old_basis_matrix(nmax+1,nmax+1); 
          HalfInt jb_min = abs(lb-HalfInt(1,2));
					HalfInt jb_max = lb + HalfInt(1,2);
                
	        //std::cout << "lb is: " << lb << std::endl;

					for(int i=0; i <= nmax; i++)
						for(int j=0; j <= nmax; j++)
							{                               
								old_basis_matrix(i,j) = radial_matrix.GetMatrixElement(la,lb,i,j);
							}
		                        	

		      //std::cout << "The old matrix is: " << std::endl << old_basis_matrix << std::endl;

					for(HalfInt ja = ja_min; ja <= ja_max; ja++)
						for(HalfInt jb = jb_min; jb <= jb_max; jb++)
							{
								Eigen::MatrixXd left_side_laja_transformation = Eigen::MatrixXd::Identity(nmax+1,nmax+1).transpose();
        				Eigen::MatrixXd right_side_lbjb_transformation = Eigen::MatrixXd::Identity(nmax+1,nmax+1);
								Eigen::MatrixXd transformed_matrix = left_side_laja_transformation * old_basis_matrix * right_side_lbjb_transformation;
								
								std::cout << "(la,ja) = (" << la << ", " << ja << ")" << std::endl;
								std::cout << std::endl << left_side_laja_transformation << std::endl;
                std::cout << "(lb,jb) = (" << lb << ", " << jb << ")" << std::endl;                  
                std::cout << std::endl << right_side_lbjb_transformation << std::endl;
								std::cout << "Matrix in old basis: " << std::endl;
                std::cout << std::endl << old_basis_matrix << std::endl;
                std::cout << "Transformed matrix: " << std::endl;
                std::cout << std::endl << transformed_matrix << std::endl;
	
		            output_file << transformed_matrix << std::endl;                                                                                  
              }                       
				}			
		}
}


// main...

int main()
{
	NORadialParameters input;
  ReadNORadialParameters(input);

  NOBasisRadialIntegralFilelj(input);
	//NOBasisRadialIntegralIdentityFilelj(input);

	return 0;
}








