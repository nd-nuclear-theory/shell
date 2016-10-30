#include <iostream>
#include <vector>
#include <string>

#include "halfint.h"
#include "wigner_gsl.h"
#include "angular_momentum.h"
#include "no_filename_retrieval.h"
#include "shell_radial.h"

#include </global/homes/c/cconsta1/eigen/Eigen/Dense>
#include </global/homes/c/cconsta1/eigen/Eigen/Eigenvalues>


struct NOxformParameters
{
	int nmax, nv;
	std::string species;
	std::string path_and_xform_file;
};

void ReadNOxformParameters(NOxformParameters& input)
{
	std::string input_line;

	{
		std::getline(std::cin,input_line);
		std::istringstream line_stream(input_line);
		line_stream >> input.nmax >> input.nv >> input.species;
	}

	std::cout << input.nmax << " " << input.nv << " " << input.species << std::endl;

	{
		std::getline(std::cin,input_line);
		std::istringstream line_stream(input_line);
		line_stream >> input.path_and_xform_file;
	}
	std::cout << input.path_and_xform_file << std::endl;
}

// Agust 15th 2015
// Iterating the NOs..

void NOxformFilesljiter(NOxformParameters& input)
{
	int max_occupied_l = (input.nmax) + (input.nv);
  std::string species_type = input.species;  

  shell::RadialMatriceslj radial_matrix;   
  radial_matrix.Read(input.path_and_xform_file);

  int filetype_position = (input.path_and_xform_file).find_last_of(".");
  std::string filetype = (input.path_and_xform_file).substr( filetype_position - 5 );

	std::string output_filename = "radial-xform-NO-"+ input.species + "-" + filetype;

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

      for(int lb=(la + dl%2); lb <= std::min(la+dl,max_occupied_l+dl) ; lb+=2)
      	{  
          HalfInt jb_min = abs(lb-HalfInt(1,2));
					HalfInt jb_max = lb+HalfInt(1,2);                       					
                     				
					for(HalfInt ja = ja_min; ja <= ja_max; ja++)
			    	for(HalfInt jb = jb_min; jb <= jb_max; jb++)
							{ 
              //if((AllowedTriangle(ja,jb,dl)))
              //  {
								Eigen::MatrixXd old_basis_matrix(nmax+1,nmax+1);
				
									for(int i=0; i <= nmax; i++)
										for(int j=0; j <= nmax; j++)
											{                               
												old_basis_matrix(i,j) = radial_matrix.GetMatrixElementlj(i,la,ja,j,lb,jb);
											}


              	Eigen::MatrixXd right_side_lbjb_transformation = NaturalOrbital(lb, jb, species_type);

              	std::cout << "(la,ja) = (" << la << ", " << ja << ")" << std::endl;
                std::cout << "(lb,jb) = (" << lb << ", " << jb << ")" << std::endl;
                std::cout << "The old matrix is: " << std::endl << old_basis_matrix << std::endl << std::endl;
                std::cout << "The NO is: " << std::endl << right_side_lbjb_transformation << std::endl << std::endl;
								std::cout << "The transformed matrix is: " << std::endl << old_basis_matrix * right_side_lbjb_transformation << std::endl << std::endl;
								
                output_file << old_basis_matrix * right_side_lbjb_transformation << std::endl;
								
							//  }
							}       
				}
		}
}

// main...
int main()
{
	NOxformParameters input;
	ReadNOxformParameters(input); 

	NOxformFilesljiter(input);

	return 0;
}












