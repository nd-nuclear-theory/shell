/******************************************************************************/
/**
  @file natorb-gen.cpp

  generate natural orbitals

  Syntax:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    natorb-gen input_orbitals robdme_info robdme_stat output_xform output_orbitals
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Patrick J. Fasano
  University of Notre Dame

  + 11/13/16 (pjf): Created, based on orbital-gen.
  + 10/12/17 (pjf): Update for changes to radial_io.
  + 11/28/17 (pjf): Print header with version.
  + 12/29/17 (pjf): Use input orbital indexing for output orbitals.
  + 12/30/17 (pjf): Ensure orbital file is properly sorted.
  + 07/27/18 (pjf): Update for new OBDME input routines.
  + 05/09/19 (pjf): Use std::size_t for basis indices and sizes.
  + 08/16/19 (pjf): Remove radial operator type and power from OutOBMEStream.

******************************************************************************/

#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>

#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Eigenvalues"

#include "basis/nlj_orbital.h"
#include "density/obdme_io.h"
#include "obme/obme_io.h"
#include "basis/operator.h"
#include "mcutils/eigen.h"

////////////////////////////////////////////////////////////////
// define a function mapping eigenvalue to weight
/////////////////////////////////////////////////////////////////

// traditional oscillator weight
void ComputeOrbitalWeights(std::vector<basis::OrbitalPNInfo>& orbitals) {
  for (auto& orbital: orbitals) {
    orbital.weight = 2 * orbital.n + orbital.l;
  }
}

// // weight = 0.25 * |log(eigenvalue/max_eigenvalue)|
// //   -- max_eigenvalue specific to species
// void ComputeOrbitalWeights(std::vector<basis::OrbitalPNInfo>& orbitals) {
//   double max_eigenvalue[2];
//   for (auto& orbital : orbitals) {
//     if (orbital.weight > max_eigenvalue[static_cast<int>(orbital.orbital_species)]) {
//       max_eigenvalue[static_cast<int>(orbital.orbital_species)] = orbital.weight;
//     }
//   }
//
//   for (auto& orbital : orbitals) {
//     if (orbital.weight > 0 ) {
//       orbital.weight = 0.25 *
//         std::abs(std::log(orbital.weight / max_eigenvalue[static_cast<int>(orbital.orbital_species)]));
//     } else {
//       orbital.weight = 999;
//     }
//   }
// }

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

// Stores simple parameters for run
struct RunParameters {
  // filenames
  std::string input_orbital_file;
  std::string robdme_info;
  std::string robdme_stat;
  std::string output_xform_file;
  std::string output_orbital_file;
};

void PrintUsage(const char **argv) {
  std::cout << "Usage: " << argv[0]
            << " input_orbitals robdme_info robdme_stat output_xform output_orbitals"
            << std::endl;
}

void ProcessArguments(int argc, const char *argv[], RunParameters& run_parameters) {
  // argument counter
  int arg = 0;

  // usage message
  if (argc-1 != 5) {
    PrintUsage(argv);
    std::cerr << "ERROR: Incorrect number of arguments." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // input orbital file
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.input_orbital_file;
    struct stat st;
    if (stat(run_parameters.input_orbital_file.c_str(), &st) != 0) {
      std::cerr << "ERROR: file " << run_parameters.input_orbital_file << " does not exist!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // input ROBDME info file
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.robdme_info;
    struct stat st;
    if (stat(run_parameters.robdme_info.c_str(), &st) != 0) {
      std::cerr << "ERROR: file " << run_parameters.robdme_info << " does not exist!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // input static ROBDME file
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.robdme_stat;
    struct stat st;
    if (stat(run_parameters.robdme_stat.c_str(), &st) != 0) {
      std::cerr << "ERROR: file " << run_parameters.robdme_stat << " does not exist!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // output xform file
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.output_xform_file;
    struct stat st;
    if (stat(run_parameters.output_xform_file.c_str(), &st) == 0) {
      std::cerr << "WARN: overwriting file " << run_parameters.output_xform_file << std::endl;
    }
  }

  // output orbital file
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.output_orbital_file;
    struct stat st;
    if (stat(run_parameters.output_orbital_file.c_str(), &st) == 0) {
      std::cerr << "WARN: overwriting file " << run_parameters.output_orbital_file << std::endl;
    }
  }
}

///////////////////////////////////////////////////////////////////////
// Take an EigenSolver and extract sorted eigenvalues and eigenvectors.
///////////////////////////////////////////////////////////////////////
void SortEigensystem(Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>& eigensolver,
                     Eigen::VectorXd& eigenvalues,
                     Eigen::MatrixXd& eigenvectors)
{
  // initialize Eigen objects
  std::size_t dimension = eigensolver.eigenvalues().size();
  eigenvalues.resize(dimension);
  Eigen::MatrixXd permutation_matrix = Eigen::MatrixXd::Zero(dimension, dimension);
  // generate a vector of pairs (sortable by std::sort) which keeps track of
  // the eigenvalue and the eigenvalue's original index
  std::vector<std::pair<double,std::size_t>> eigenvalue_keys;
  for (std::size_t i=0; i < dimension; ++i) {
    eigenvalue_keys.emplace_back(eigensolver.eigenvalues()(i), i);
  }

  #if __cplusplus > 201103L
  // C++14 has a nice way to sort in descending lexicographical order
  std::sort(eigenvalue_keys.begin(), eigenvalue_keys.end(), std::greater<>());
  #else
  // C++11's way of sorting is not as pretty, using reverse iterators
  std::sort(eigenvalue_keys.rbegin(), eigenvalue_keys.rend());
  #endif

  for (std::size_t i=0; i<dimension; ++i) {
    eigenvalues(i) = eigenvalue_keys[i].first;
    permutation_matrix(eigenvalue_keys[i].second, i) = 1;
  }

  eigenvectors = eigensolver.eigenvectors() * permutation_matrix;
}

////////////////////////////////////////////////////////////////
// main program
////////////////////////////////////////////////////////////////
int main(int argc, const char *argv[]) {
  // header
  std::cout << std::endl;
  std::cout << "natorb-gen -- natural orbital generation" << std::endl;
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  RunParameters run_parameters;
  ProcessArguments(argc, argv, run_parameters);

  // Read input orbitals
  std::ifstream input_orbital_s(run_parameters.input_orbital_file);
  basis::OrbitalSpaceLJPN input_space(basis::ParseOrbitalPNStream(input_orbital_s, true));
  input_orbital_s.close();

  // Get OBDMEs
  basis::OperatorBlocks<double> density_matrices;
  shell::InOBDMEStreamMulti obdme_reader(run_parameters.robdme_info, run_parameters.robdme_stat, input_space);

  // Get sectors and blocks
  basis::OrbitalSectorsLJPN sectors;
  obdme_reader.GetMultipole(0, sectors, density_matrices);

  // Here we loop through the density matrices and diagonalize each sector.
  basis::OperatorBlocks<double> xform_matrices;
  // TODO(pjf) generalize to use eigenvalues for orbital weights
  // std::vector<basis::OrbitalPNInfo> output_orbitals;
  for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index) {
    // get next sector
    const basis::OrbitalSectorsLJPN::SectorType sector = sectors.GetSector(sector_index);
    basis::OrbitalSpeciesPN species = sector.bra_subspace().orbital_species();
    int l = sector.bra_subspace().l();
    HalfInt j = sector.bra_subspace().j();

    // the one-body static density matrix should be real and symmetric, so
    // we can use the SelfAdjointEigenSolver from Eigen
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(density_matrices[sector_index]);

    // use sort code to ge eigen vectors and eigenvalues
    Eigen::VectorXd eigenvalues;
    Eigen::MatrixXd eigenvectors;
    SortEigensystem(eigensolver, eigenvalues, eigenvectors);

    // add to output
    xform_matrices.push_back(eigenvectors);
    // TODO(pjf) generalize to use eigenvalues for orbital weights
    // for (int i=0; i < eigenvalues.size(); ++i) {
    //   output_orbitals.emplace_back(species, i, l, j, eigenvalues(i));
    // }
  }

  // TODO(pjf) generalize to use eigenvalues for orbital weights
  // ComputeOrbitalWeights(output_orbitals);

  // write xform out to file
  // TODO(pjf) construct new indexing for output
  basis::OrbitalSpaceLJPN& output_space = input_space;
  basis::OrbitalSectorsLJPN output_sectors(input_space, output_space, 0, 0, 0);
      // note hard-coded l0max=0, Tz0=0
  shell::OutOBMEStream xs(run_parameters.output_xform_file,
                          input_space, output_space, output_sectors,
                          basis::OneBodyOperatorType::kRadial);
  xs.Write(xform_matrices);

  // sort orbitals and write out
  std::vector<basis::OrbitalPNInfo> output_orbitals = output_space.OrbitalInfo();
  std::sort(output_orbitals.begin(), output_orbitals.end(), basis::OrbitalSortCmpWeight);
  std::ofstream output_orbital_s(run_parameters.output_orbital_file);
  output_orbital_s << basis::OrbitalDefinitionStr(output_orbitals, true);
  output_orbital_s.close();

  /* return code */
  return EXIT_SUCCESS;
}
