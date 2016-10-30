// noradial.cpp produces the radial files which document the matrix elements 
// for the operators r^{1}, r^{2} and k^{1} and k^{2}. The matrix elements 
// are taken between natural orbitals, which are linear combinations of harmonic 
// oscillator or CS functions. The natural orbitals radial functions have the form
// NO_{na, la, ja}. 

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "halfint.h"
#include "wigner_gsl.h"
#include "angular_momentum.h"
#include "no_filename_retrieval.h"
#include "shell_radial.h"

#include </global/homes/c/cconsta1/eigen/Eigen/Dense>
#include </global/homes/c/cconsta1/eigen/Eigen/Eigenvalues>

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

// 08/15/15
// NOBasisRadialIntegralFileljiter produces radial files with both l and j dependence..
// The difference with the function NOBasisRadialIntegralFilelj is that this time 
// the radial file already has an lj dependence from an initial NO change of basis
// but now we are finding the change of basis for the second time hence the word iter.
// We are iterating the NOs..

void NOBasisRadialIntegralFileljiter ( NORadialParameters& input )
{
	int max_occupied_l = (input.nmax) + (input.nv);
  std::cout << "max occupied l is: " << max_occupied_l << std::endl;
	std::string species_type = input.species;  

	shell::RadialMatriceslj radial_matrix;  
	radial_matrix.Read(input.path_and_radial_file);

	int filetype_position = (input.path_and_radial_file).find_last_of(".");
	std::string filetype = (input.path_and_radial_file).substr( filetype_position - 2 );

	std::string output_filename = "radial-me-NO-"+ input.species + "-" + filetype;

	std::ofstream output_file(output_filename);

	int lmax = radial_matrix.Getlmax();
	int dl = radial_matrix.Getdl();
	int nmax = radial_matrix.Getnmax();   

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
          HalfInt jb_min = abs(lb-HalfInt(1,2));
					HalfInt jb_max = lb + HalfInt(1,2);
                
	        //std::cout << "lb is: " << lb << std::endl;                      

					for(HalfInt ja = ja_min; ja <= ja_max; ja++)
						for(HalfInt jb = jb_min; jb <= jb_max; jb++)
							{ 

							Eigen::MatrixXd old_basis_matrix(nmax+1,nmax+1);
				
									for(int i=0; i <= nmax; i++)
										for(int j=0; j <= nmax; j++)
											{                               
												old_basis_matrix(i,j) = radial_matrix.GetMatrixElementlj(i,la,ja,j,lb,jb);
											}
 
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

// main...

int main()
{
	NORadialParameters input;
  ReadNORadialParameters(input);

  NOBasisRadialIntegralFileljiter(input);
	//NOBasisRadialIntegralIdentityFilelj(input);

	return 0;
}








