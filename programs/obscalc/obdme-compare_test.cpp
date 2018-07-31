/******************************************************************************/
/**
  @file natorb-gen.cpp

  generate natural orbitals

  Syntax:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obdme-compare_test input_orbitals robdme_info robdme robdme_info robdme
    obdme-compare_test input_orbitals robdme_info robdme robdme
    obdme-compare_test input_orbitals robdme robdme
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Patrick J. Fasano
  University of Notre Dame

  + 07/25/18 (pjf): Created, based on natorb-gen.

******************************************************************************/

#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>

#include "basis/nlj_orbital.h"
#include "density/obdme_io.h"
#include "basis/operator.h"
#include "mcutils/eigen.h"

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

// Stores simple parameters for run
struct RunParameters {
  // filenames
  std::string input_orbital_file;
  std::string robdme_info1;
  std::string robdme1;
  std::string robdme_info2;
  std::string robdme2;
};

void PrintUsage(const char **argv) {
  std::cout << "Usage: " << argv[0]
                         << " input_orbitals [robdme_info] robdme [robdme_info] robdme"
                         << std::endl;
}

void ProcessArguments(int argc, const char *argv[], RunParameters& run_parameters) {
  // argument counter
  int arg = 0;

  // usage message
  if ((argc-1 < 3) || (argc-1 > 5)) {
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
  if (argc-1 >= 4)
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.robdme_info1;
    struct stat st;
    if (stat(run_parameters.robdme_info1.c_str(), &st) != 0) {
      std::cerr << "ERROR: file " << run_parameters.robdme_info1 << " does not exist!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  else
  {
    run_parameters.robdme_info1 = "";
  }

  // input ROBDME file
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.robdme1;
    struct stat st;
    if (stat(run_parameters.robdme1.c_str(), &st) != 0) {
      std::cerr << "ERROR: file " << run_parameters.robdme1 << " does not exist!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // input ROBDME info file
  if (argc-1 == 5)
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.robdme_info2;
    struct stat st;
    if (stat(run_parameters.robdme_info1.c_str(), &st) != 0) {
      std::cerr << "ERROR: file " << run_parameters.robdme_info2 << " does not exist!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  else
  {
    run_parameters.robdme_info2 = "";
  }

  // input static ROBDME file
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.robdme2;
    struct stat st;
    if (stat(run_parameters.robdme2.c_str(), &st) != 0) {
      std::cerr << "ERROR: file " << run_parameters.robdme2 << " does not exist!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
}


////////////////////////////////////////////////////////////////
// main program
////////////////////////////////////////////////////////////////
int main(int argc, const char *argv[]) {
  // header
  std::cout << std::endl;
  std::cout << "obdme-compare_test -- natural orbital generation" << std::endl;
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  RunParameters run_parameters;
  ProcessArguments(argc, argv, run_parameters);

  // Read input orbitals
  std::ifstream input_orbital_s(run_parameters.input_orbital_file);
  basis::OrbitalSpaceLJPN input_space(basis::ParseOrbitalPNStream(input_orbital_s, true));
  input_orbital_s.close();

  // Get OBDMEs
  shell::InOBDMEStream in_obdme1, in_obdme2;
  if (run_parameters.robdme_info1 == "") {
    std::cout << "Hello ";
    in_obdme1 = shell::InOBDMEStreamSingle(run_parameters.robdme1, input_space);
  } else {
    in_obdme1 = shell::InOBDMEStreamMulti(run_parameters.robdme_info1, run_parameters.robdme1, input_space);
  }
  if (run_parameters.robdme_info2 == "") {
    std::cout << "world!" << std::endl;
    in_obdme2 = shell::InOBDMEStreamSingle(run_parameters.robdme2, input_space);
  } else {
    in_obdme2 = shell::InOBDMEStreamMulti(run_parameters.robdme_info2, run_parameters.robdme2, input_space);
  }

  assert(in_obdme1.g0() == in_obdme2.g0());
  assert(in_obdme1.Tz0() == in_obdme2.Tz0());

  int min_K = std::max(in_obdme1.min_K(), in_obdme2.min_K());
  int max_K = std::min(in_obdme1.max_K(), in_obdme2.max_K());
  for (int K = min_K; K <= max_K; K++) {
    std::cout << "Comparing K=" << K << std::endl;
    basis::OrbitalSectorsLJPN sectors;
    basis::OperatorBlocks<double> matrices1, matrices2, ratio_matrices;
    in_obdme1.GetMultipole(K, sectors, matrices1);
    in_obdme2.GetMultipole(K, sectors, matrices2);
    basis::SetOperatorToZero(sectors, ratio_matrices);

    for (int sector_index = 0; sector_index < sectors.size(); ++sector_index) {
      // get next sector
      const basis::OrbitalSectorsLJPN::SectorType sector = sectors.GetSector(sector_index);
      // get sizes
      const int bra_subspace_size = sector.bra_subspace().size();
      const int ket_subspace_size = sector.ket_subspace().size();

      // main loop
      #pragma omp parallel for collapse(2)
      for (int j=0; j < bra_subspace_size; ++j) {
        for (int k=0; k < ket_subspace_size; ++k) {
          if ((matrices2[sector_index](j,k) < 1e-12) || (matrices1[sector_index](j,k) < 1e-12))
            continue;
          ratio_matrices[sector_index](j,k) = matrices2[sector_index](j,k)/matrices1[sector_index](j,k);
        }
      }
      std::cout << mcutils::FormatMatrix(ratio_matrices[sector_index], "12.4e") << std::endl << std::endl;
    }
  }

  /* return code */
  return EXIT_SUCCESS;
}
