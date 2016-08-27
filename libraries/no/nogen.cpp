//
// nogen.cpp 05/03/15
//
// nogen.cpp diagonalizes the density matrix which couples pairs of orbitals like c_{a}^{\dagger} c_{b},
// where a = (nalaja) and b =(nblbjb). We don't diagonalize the whole matrix but only submatrices which are
// created using matrix elements for which la = lb and ja = jb (read below)..

// To generate the natural orbitals up to a given N_{max}, we create 
// square matrices with la = lb and ja = jb. Therefore the density sub-matrix for a family of orbitals
// that share the same l and the same j is symmetric. Furthermore recall that N <= N_{max} + N_{0} is the 
// total number of oscillator quanta in the wave function. Therefore the highest possible (lj) pair has l = N
// and j = N + 1/2. This is because N = 2n + l, thus if n = 0, then l = N. Also j = l + 1/2 or l - 1/2 therefore
// j_{max} = N + 1/2. For example if the nucleus is a p-shell nucleus then N_{0} = 1. If we perform an N_{max} = 2
// run, then N = 3, and the highest (lj) pair has (l,j) = (3, 7/2). 

// If we take the combination l = 0 and j = 1/2, then since N <= 3 and because N = 2n + l, we conclude that the only
// l = 0 and j = 1/2 orbitals in our basis have (n,l,j) = (0, 0, 1/2) and (1, 0, 1/2). The matrix that "mixes" these two orbitals
// will therefore be a 2 x 2 matrix, from the diagonalization of which we get the natural orbitals. The natural orbitals will be a linear
// combination of the two orbitals 0s1/2 and 1s1/2 (in spectroscopic notation). If the two natural orbitals are |0>, and |1>, then 
// |0> = a |0s1/2> + b |1s1/2>, where the coefficients are determined from the diagonalization. The eigenvalues that correspond to the two natural orbitals 
// |0> and |1> are the occupancies of the two natural orbitals.

// The diagonalization is done using methods from the library eigen. The caveat is that when we create the mixing matrix we zero-pad it, so that it will have the
// same dimensions as the current radial-files.. This is because we use the matrix of eigenvectors for the change of basis from say the harmonic 
// oscillator radial integrals, to the natural orbitals radial integrals.. The old radial integral files contained matrices with dimensions 11 x 11, with n running
// from (0, 1, 2, .., 10) for each l. Therefore our mixing matrices (which are submatrices of the density files) will have size 11 x 11. In the above example of the 
// s1/2 space the only non-zero matrix elements are the elements a_{0,0}, a_{0,1}, a_{1,0} and a_{1,1}.

// Diagonalization: The eigenvectors are outputted as columns of a square matrix.. 

// July 4th 2015: Change input to number operator function from PostProcessingParameters that is 
// mostly used when running emcalc to NaturalOrbitalIngredients.
// Rename number_operator to NumberOperator

/*
  8/26/16 (mac): Minimal patches to compile under shell project.
*/


#include <iostream>
#include <vector>
#include <string>

#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Eigenvalues"

#include "am/halfint.h"

#include "no/mfdn_file_processing.h"
#include "no/no_filename_retrieval.h"



// Input parameters:
// nmax and nv were descriped above, species determines whether we will get the natural orbitals for protons or neutrons
// Then we need the full path to the info file, and the full path to the statrobdme file..

struct NaturalOrbitalsIngredients
{
   int nmax, nv;
   std::string species;
   std::string reduced_one_body_density_path; // path to reduced one body density
   std::string orbital_info_path; // path to orbital file
};

void ReadNaturalOrbitalsIngredients(NaturalOrbitalsIngredients& input_paths)
{
  std::string input_line;

  {
    std::getline(std::cin,input_line);
    std::istringstream line_stream(input_line);
    line_stream >> input_paths.nmax >> input_paths.nv >> input_paths.species;
  }

  // std::cout << input_paths.nmax << input_paths.nv << input_paths.species << std::endl;

  {
    std::getline(std::cin,input_line);
    std::istringstream line_stream(input_line);
    line_stream >> input_paths.orbital_info_path;
  }

  // std::cout << input_paths.orbital_info_path << std::endl;

  {
    std::getline(std::cin,input_line);
    std::istringstream line_stream(input_line);
    line_stream >> input_paths.reduced_one_body_density_path;
  }
  
   // std::cout << input_paths.reduced_one_body_density_path << std::endl;
}  

// Included here as a quick diagnostic test in case we need to check the files integrity..

double NumberOperator(NaturalOrbitalsIngredients& input_file_parameters)
{ 
	double occupation_number = 0;
	string type = input_file_parameters.species;
  std::vector <OrbitalInformation> orbitals;
  std::vector <double> densities;

  ReadOrbitalInformationFile(orbitals, input_file_parameters.orbital_info_path);  
	ReadOneBodyDensitiesFile(densities,input_file_parameters.reduced_one_body_density_path);

  int half_size = (densities.size())/2;
    

  //proton sector

  if(type == "p")
  	{
    	for(int i=0; i < half_size; i++)
      	{
        	if( (orbitals[i].na == orbitals[i].nb) && (orbitals[i].la == orbitals[i].lb) && 
            (orbitals[i].ja == orbitals[i].jb) && (orbitals[i].k == 0)  )
          	{
            	occupation_number = occupation_number + Hat(orbitals[i].ja) * densities[i];
            }
				}
		}
	else if(type == "n") //neutron sector
  	{
    	for(int i=half_size; i < 2*half_size; i++)
      	{
        	if( (orbitals[i-half_size].na == orbitals[i-half_size].nb) && (orbitals[i-half_size].la == orbitals[i-half_size].lb) && 
            (orbitals[i-half_size].ja == orbitals[i-half_size].jb) && (orbitals[i-half_size].k == 0)  )
          	{
            	occupation_number = occupation_number + Hat(orbitals[i-half_size].ja) * densities[i];
            }
				}
		}  
     
	//cout << "The sum is: " << protons + neutrons << endl;

	return (occupation_number);
}                       

Eigen::MatrixXd MixingMatrix(/*const int& n1, const int& n2, */const int& l, const HalfInt& j, const std::string& type,
                      std::vector <OrbitalInformation>& orbital_data_vector, std::vector <double>& one_body_densities)
{
	Eigen::MatrixXd mixing_matrix;

	// mixing_matrix = Eigen::MatrixXd::Zero(n1,n2);

  vector <double> mixing_elements;

  //cout << "m is: " << endl << mixing_matrix << endl;

  std::string species = type;   
   
  int file_size = orbital_data_vector.size();

  if(species == "p")
		{
    	for (int i=0; i < file_size; i++)
      	{
        	if( (orbital_data_vector[i].k == 0) && (orbital_data_vector[i].la == l) && (orbital_data_vector[i].ja == j) &&
                 (orbital_data_vector[i].lb == l) && (orbital_data_vector[i].jb == j) )
          	{
            	mixing_elements.push_back(one_body_densities[i]);
            }
        } 
     }
	else if(species == "n") 
		{
    	for (int i=0; i < file_size; i++)
      	{
        	if( (orbital_data_vector[i].k == 0) && (orbital_data_vector[i].la == l) && (orbital_data_vector[i].ja == j) &&
                 (orbital_data_vector[i].lb == l) && (orbital_data_vector[i].jb == j) )
        	  {
            	mixing_elements.push_back(one_body_densities[i+file_size]);
            }
         } 
    }

	int mixing_matrix_size = sqrt(mixing_elements.size());
	mixing_matrix = Eigen::MatrixXd(mixing_matrix_size,mixing_matrix_size);

  for(int i=0; i < mixing_matrix_size; i++)
  	{
    	for(int j=0;j<mixing_matrix_size; j++)
      	{
        	mixing_matrix(i,j) =  mixing_elements[ i * mixing_matrix_size + j ];
        }
    }

	return mixing_matrix;
}


// Takes the NaturalOrbitalsIngredients struct and returns all the possible natural orbitals (for given N1b = Nmax + Nv)
// as text files, as well as the occupation probabilities files.. 


void NaturalOrbitalsFiles(NaturalOrbitalsIngredients& input)
{
	std::vector <OrbitalInformation> orbitals;
  ReadOrbitalInformationFile(orbitals, input.orbital_info_path);

  std::vector <double> densities;
  ReadOneBodyDensitiesFile(densities,input.reduced_one_body_density_path);

  int max_occupied_l = (input.nmax) + (input.nv);

  double number_operator = NumberOperator(input);

  std::string type = input.species;
   
  // Report occupation sum..
  double occupation_sum = 0;


  // This must be equal to the size of the matrices stored in the radial files..

  // int matrix_dimension = 11;

  // Note the size of the MixingMatrix must coincide with the size of the radial integral files.. Maybe this should be done internally here by opening the 
  // radial matrix file and extracting all the header information..


  // By the way in eigen the eigenvectors of a matrix are stored as columns of a matrix... So eigenvector one is in column one etc..  

  // In nogen.out we report the mixing matrix for each orbital, we output the NOs, the occupations, we test the orthonormality
  // of the eigenvectors, and we output the number operator from the sum over densities.
  std::string diagnostic_output_filename = "nogen-" + type + ".out";
  std::ofstream diagnostic_output(diagnostic_output_filename); 

  for (int l=0; l <= max_occupied_l ; l++)
		{
	    HalfInt j_min = abs(l - HalfInt(1,2));
      HalfInt j_max = l + HalfInt(1,2);
          
      for(HalfInt j=j_min; j<=j_max; j++)
      	{
        	std::string no_output_filename = NaturalOrbitalsFilename(l,j,type);
          std::string occupations_output_filename = OccupationsFilename(l,j,type);

          std::ofstream no_output(no_output_filename);
          std::ofstream occupations_output(occupations_output_filename);

          Eigen::MatrixXd orbital_mixing = MixingMatrix(/*matrix_dimension, matrix_dimension,*/ l, j, type, orbitals, densities);
          Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensystem(orbital_mixing);

					Eigen::VectorXd eigenvalues = eigensystem.eigenvalues().real().reverse();

					Eigen::MatrixXd natural_orbitals = eigensystem.eigenvectors().real().rowwise().reverse();

					Eigen::VectorXd::Index max_index[natural_orbitals.cols()];
					Eigen::VectorXd max_value(natural_orbitals.cols());

					for(int i=0; i< natural_orbitals.cols(); i++)
						{
							max_value(i) = natural_orbitals.col(i).cwiseAbs().maxCoeff( &max_index[i] );

							int sign_of_max_coeff = ((natural_orbitals.col(i))(max_index[i]) > 0) ? 1 : (((natural_orbitals.col(i))(max_index[i]) < 0) ? -1 : 0);

							natural_orbitals.col(i) = sign_of_max_coeff*natural_orbitals.col(i);

							//std::cout << "Max coeff is " <<  max_value(i) << "@ " << max_index[i] << "and the sign is: " << sign_of_max_coeff << std::endl;
						}

					Eigen::MatrixXd eigenvectors = Eigen::MatrixXd::Identity(11,11);

				  eigenvectors.topLeftCorner(natural_orbitals.rows(),natural_orbitals.cols()) = natural_orbitals;

          // Eigen::MatrixXd eigenvectors = eigensystem.eigenvectors().real();
          // Eigen::MatrixXd eigenvalues  = eigensystem.eigenvalues().real();

          occupation_sum = occupation_sum + Hat(j)*eigenvalues.sum();

          no_output << eigenvectors << std::endl;
          occupations_output << eigenvalues << std::endl;  

          // Report everything in a nice diagnostic output file that contains the mixing matrix, the natural
          // orbitals with their occupations, a check for orthonormality etc..
          // Added on July 5th 2015..

          diagnostic_output << ">>> Mixing matrix:" << std::endl << "(lj): (" << l << ", " << j << "), species: " << type << std::endl;
          diagnostic_output << std::endl << orbital_mixing << std::endl << std::endl;

					//debugging

					diagnostic_output << "Debugging.. Natural Orbitals Before reversing.. " << std::endl;
  				diagnostic_output << std::endl << eigensystem.eigenvectors().real() << std::endl << std::endl;

          diagnostic_output << "Natural orbitals matrix:" << std::endl;
          diagnostic_output << std::endl << eigenvectors << std::endl << std::endl;
          diagnostic_output << "Orthonormality of the natural orbital matrix check:" << std::endl;
          diagnostic_output << std::endl << eigenvectors*(eigenvectors.transpose()) << std::endl << std::endl;
          diagnostic_output << "Occupation numbers for each natural orbital: " << std::endl << std::endl;
          diagnostic_output << std::endl << eigenvalues << std::endl << std::endl;         
        }         
    }
       
    diagnostic_output << "Natural orbitals occupation sum:" << std::endl;   
    diagnostic_output << occupation_sum << std::endl;
    diagnostic_output << "Number operator expectation value (calculated by summing the initial robdme): " << std::endl << number_operator << std::endl;
}

// Created on July 23rd 2015 to make sure that the density matrix does not mix states with different angular momenta. 
// This was not necessary since we are dealing with a matrix coupled to zero therefore the momenta better be eequal..

/*void MixingMatrixOffDiagonal(const int& l1, const HalfInt& j1, const int& l2, const HalfInt& j2, const std::string& type,
                      std::vector <OrbitalInformation>& orbital_data_vector, std::vector <double>& one_body_densities)
{
  std::string species = type;   
   
  int file_size = orbital_data_vector.size();

  if(species == "p")
		{
    	for (int i=0; i < file_size; i++)
      	{
        	if( (orbital_data_vector[i].k == 0) && (orbital_data_vector[i].la == l1) && (orbital_data_vector[i].ja == j1) &&
                 (orbital_data_vector[i].lb == l2) && (orbital_data_vector[i].jb == j2) )
          	{
            	std::cout << "The mixing element for (l1,j1) = (" << l1 << "," << j1 << ") and (l2,j2) = (" 
												<< l2 << "," << j2 << ") is: " << one_body_densities[i] << std::endl;
            }
        } 
     }
	else if(species == "n") 
		{
    	for (int i=0; i < file_size; i++)
      	{
        	if( (orbital_data_vector[i].k == 0) && (orbital_data_vector[i].la == l1) && (orbital_data_vector[i].ja == j1) &&
                 (orbital_data_vector[i].lb == l2) && (orbital_data_vector[i].jb == j2) )
        	  {
            	std::cout << "The mixing element for (l1,j1) = (" << l1 << "," << j1 << ") and (l2,j2) = (" 
												<< l2 << "," << j2 << ") is: " << one_body_densities[i] << std::endl;
            }
         } 
    }
}*/

// main. Creates the no files..

int main()
{
	NaturalOrbitalsIngredients input;
	ReadNaturalOrbitalsIngredients(input);

  //double number_operator = NumberOperator(input);

  //cout << "The number operator expectation value is equal to: " << number_operator << endl;

  NaturalOrbitalsFiles(input);

	// The following lines (commented out with an asterisk) were created on July 23rd 2015 
	// to make sure that the density matrix does not mix states with different l1, j1 and
	// l2, j2. In retrospect there is no way that two spherical tensors coupled to a total angular momentum
	// zero can have different l or j so it was a waste of time but oh well..	

	/*std::vector <OrbitalInformation> orbitals;
  ReadOrbitalInformationFile(orbitals, input.orbital_info_path);

  std::vector <double> densities;
  ReadOneBodyDensitiesFile(densities,input.reduced_one_body_density_path);
	std::string type = input.species;

	MixingMatrixOffDiagonal(1, HalfInt(1,2), 1 , HalfInt(3,2), type, orbitals, densities);*/

	return 0;
}


