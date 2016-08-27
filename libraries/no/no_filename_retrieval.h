/*
  8/26/16 (mac): Minimal patches to compile under shell project.
*/

// Simple functions that return the filename for the natural orbitals
// and their occupancies. The functions are used both for the creation of the 
// natural orbitals and their retrieval.. 
//

#ifndef no_filename_retrieval_h
#define no_filename_retrieval_h

#include <iostream>
#include <fstream>
#include <string>

#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Eigenvalues"

#include "am/halfint.h"


// Returns the string with the NO filename

inline
std::string NaturalOrbitalsFilename(const int& l, const HalfInt& j, std::string& species)
{
   return "natural-orbital-2l" + std::to_string(2*l) + "-2j" + std::to_string(j.TwiceValue()) + "-" + species + ".dat";
}

// Returns the string with the occupations file filename

inline
std::string OccupationsFilename(const int& l, const HalfInt& j, std::string& species)
{
   return "occupations-orbital-2l" + std::to_string(2*l) + "-2j" + std::to_string(j.TwiceValue()) + "-" + species + ".dat";
}

// Opens up the already created NO file and returns it as a matrix
// returns the matrix (in eigen form) of the natural orbitals

inline
Eigen::MatrixXd NaturalOrbital(int& l, HalfInt& j, std::string& type)
{
    std::string filename = NaturalOrbitalsFilename(l,j,type);
    std::ifstream input_file(filename);

    // The matrix size is 11 because the radial files that we use as our basis to construct
    // radial matrices in the NO basis have dimensio 11 x 11

    int matrix_size = 11;

    Eigen::MatrixXd m(matrix_size, matrix_size);
    double element;

    for(int i = 0; i < matrix_size; i++)
       {
          for(int j =0; j < matrix_size; j++)
             {
                input_file >> element;
                m(i,j) = element;
             }
       }

    return m;
}

#endif
