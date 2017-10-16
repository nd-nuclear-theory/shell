/******************************************************************************
  obscalc-ob.cpp

  Calculate one-body observable from one-body operator matrix elements and
one-body
  density matrix elements.

  Syntax:
    + obscalc-ob operator_filename robdme_info robdme_file

  Patrick J. Fasano
  University of Notre Dame

  + 10/08/17 (pjf): Created.

******************************************************************************/

#include <sys/stat.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "basis/nlj_orbital.h"
#include "basis/proton_neutron.h"
#include "cppformat/format.h"
#include "density/obdme_io.h"
#include "radial/radial_io.h"

// Stores simple parameters for run
struct RunParameters {
  // filenames
  std::string operator_filename;
  std::string robdme_info;
  std::string robdme_file;
};

void PrintUsage(char** argv) {
  std::cout << "Usage: " << argv[0]
            << " operator_filename robdme_info robdme_file" << std::endl;
}

void ProcessArguments(int argc, char** argv, RunParameters& run_parameters) {
  // argument counter
  int arg = 1;

  // usage message (must at least provide mode)
  if (argc - 1 < 3) {
    PrintUsage(argv);
    std::exit(EXIT_FAILURE);
  }

  // operator file
  {
    std::istringstream parameter_stream(argv[arg++]);
    parameter_stream >> run_parameters.operator_filename;
    if (!parameter_stream) {
      PrintUsage(argv);
      std::cerr << "Invalid operator filename." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    struct stat st;
    if (stat(run_parameters.operator_filename.c_str(), &st) != 0) {
      std::cerr << "ERROR: file " << run_parameters.operator_filename
                << " does not exist!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // input ROBDME info file
  {
    std::istringstream parameter_stream(argv[arg++]);
    parameter_stream >> run_parameters.robdme_info;
    std::cout << run_parameters.robdme_info << std::endl;
    struct stat st;
    if (stat(run_parameters.robdme_info.c_str(), &st) != 0) {
      std::cerr << "ERROR: file " << run_parameters.robdme_info
                << " does not exist!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // input static ROBDME file
  {
    std::istringstream parameter_stream(argv[arg++]);
    parameter_stream >> run_parameters.robdme_file;
    std::cout << run_parameters.robdme_file << std::endl;
    struct stat st;
    if (stat(run_parameters.robdme_file.c_str(), &st) != 0) {
      std::cerr << "ERROR: file " << run_parameters.robdme_file
                << " does not exist!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
}

int main(int argc, char** argv) {
  // header
  std::cout << std::endl;
  std::cout << "obscalc-ob -- one-body observable evaluation" << std::endl;
  std::cout << std::endl;

  // process arguments
  RunParameters run_parameters;
  ProcessArguments(argc, argv, run_parameters);

  // parallel performance diagnostic
  std::cout << fmt::format("INFO: OMP max_threads {}, num_procs {}",
                           omp_get_max_threads(), omp_get_num_procs())
            << std::endl
            << std::endl;

  // read operator
  shell::InRadialStream is(run_parameters.operator_filename);
  basis::OperatorBlocks<double> operator_matrices;
  is.Read(operator_matrices);

  // get indexing
  basis::OrbitalSpaceLJPN bra_space, ket_space;
  basis::OrbitalSectorsLJPN sectors;
  is.SetToIndexing(bra_space, ket_space, sectors);
  assert(bra_space.OrbitalInfo() == ket_space.OrbitalInfo());
  basis::OrbitalSpaceLJPN& space = bra_space;
  is.Close();

  // get OBDMEs
  basis::OperatorBlocks<double> density_matrices;
  shell::InOBDMEReader obdme_reader(run_parameters.robdme_info, space,
                                    sectors.g0(), sectors.Tz0());
  obdme_reader.ReadMultipole(run_parameters.robdme_file, sectors.j0(),
                             density_matrices);

  // loop and sum over \sum_{a,b} rho_{ab} T_{ba}
  double value = 0;
  for (int subspace_index_a = 0; subspace_index_a < space.size();
       ++subspace_index_a) {
    for (int subspace_index_b = 0; subspace_index_b < space.size();
         ++subspace_index_b) {
      auto subspace_a = space.GetSubspace(subspace_index_a);
      auto subspace_b = space.GetSubspace(subspace_index_b);
      auto sector_index =
          sectors.LookUpSectorIndex(subspace_index_a, subspace_index_b);
      if (sector_index == basis::kNone) continue;

      for (int state_index_a = 0; state_index_a < subspace_a.size();
           ++state_index_a) {
        for (int state_index_b = 0; state_index_b < subspace_b.size();
             ++state_index_b) {
          value += operator_matrices[sector_index](state_index_a, state_index_b)
                   * density_matrices[sector_index](state_index_a, state_index_b);
        }
      }
    }
  }

  std::cout << "value: " << value << std::endl;
}
