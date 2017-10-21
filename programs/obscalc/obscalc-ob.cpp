/******************************************************************************
  obscalc-ob.cpp

  Calculate one-body observable from one-body operator matrix elements and
one-body
  density matrix elements.

  Syntax:
    + obscalc-ob

  Input format:

    set-output-file output_filename
    set-indexing orbital_filename
    define-operator name operator_filename
    define-static-densities J g n robdme_info_filename robdme_filename
    define-transition-densities Jf gf nf Ji gi fi robdme_info_filename
robdme_filename

  Patrick J. Fasano
  University of Notre Dame

  + 10/08/17 (pjf): Created.

******************************************************************************/

#include <sys/stat.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "am/halfint.h"
#include "basis/nlj_orbital.h"
#include "basis/proton_neutron.h"
#include "cppformat/format.h"
#include "density/obdme_io.h"
#include "radial/radial_io.h"

// Store one-body operators
struct OneBodyOperator {
  std::string name;
  basis::OrbitalSpaceLJPN space;
  basis::OrbitalSectorsLJPN sectors;
  basis::OperatorBlocks<double> operator_matrices;

  explicit OneBodyOperator(const std::string& name__, const std::string& filename)
      : name(name__) {
    // read operator
    shell::InRadialStream is(filename);
    is.Read(operator_matrices);

    // get indexing
    basis::OrbitalSpaceLJPN ket_space;
    is.SetToIndexing(space, ket_space, sectors);
    assert(space.OrbitalInfo() == ket_space.OrbitalInfo());
    is.Close();
  }
}

// store one-body static densities
struct StaticOneBodyDensities {
  HalfInt J;
  int g, n;
  std::string robdme_filename;
  shell::InOBDMEReader obdme_reader;

  StaticOneBodyDensities(HalfInt J__, int g__, int n__,
                         const basis::OrbitalSpaceLJPN& space,
                         const std::string& robdme_info,
                         const std::string& robdme_file)
      : J(J__),
        g(g__),
        n(n__),
        obdme_reader(robdme_info, space, 0, 0),
        robdme_filename(robdme_file) {}
  // TODO(pjf): generalize to parity- and isospin-changing operators

  void ReadMultipole(int j0, basis::OperatorBlocks<double> density_matrices) {
    obdme_reader.ReadMultipole(robdme_file, j0, density_matrices);
  }
}

// store one-body static densities
struct TransitionOneBodyDensities {
  HalfInt Ji, Jf;
  int gi, gf, ni, nf;
  std::string robdme_filename;
  shell::InOBDMEReader obdme_reader;

  TransitionOneBodyDensities(HalfInt Jf__, int gf__, int nf__, HalfInt Ji__,
                             int gi__, int ni__,
                             const basis::OrbitalSpaceLJPN& space,
                             const std::string robdme_info,
                             const std::string robdme_file)
      : Jf(Jf__),
        gf(gf__),
        nf(nf__),
        Ji(Ji__),
        gi(gi__),
        ni(ni__),
        obdme_reader(robdme_info, space, 0, 0),
        robdme_filename(robdme_file) {
    assert(gf == gi);
    // TODO(pjf): generalize to parity- and isospin-changing operators
  }

  void ReadMultipole(int j0, basis::OperatorBlocks<double> density_matrices) {
    // assert(am::AllowedTriangle(Jf, Ji, j0));
    obdme_reader.ReadMultipole(robdme_file, j0, density_matrices);
  }
}

// Stores parameters for run
struct RunParameters {
  // filenames
  std::string output_filename;
  basis::OrbitalSpaceLJPN space;
  std::vector<OneBodyOperator> operators;
  std::vector<StaticOneBodyDensities> static_densities;
  std::vector<TransitionOneBodyDensities> transition_densities;
};

void PrintUsage(char** argv) { std::cout << "Usage: " << argv[0] << std::endl; }

void ReadParameters(RunParameters& run_parameters) {
  std::string line;
  int line_count = 0;

  while (std::getline(std::cin, line)) {
    ++line_count;

    // set up for line parsing
    std::istringstream line_stream(line);
    std::string keyword;
    line_stream >> keyword;

    // skip blank line or hash comment line
    if ((keyword == "") || (keyword == "#")) continue;

    // select action based on keyword
    if (keyword == "set-output-file") {
      line_stream >> run_parameters.output_filename;
      ParsingCheck(line_stream, line_count, line);
      struct stat st;
      if (stat(run_parameters.output_filename.c_str(), &st) == 0) {
        std::cerr << "WARN: overwriting file " << run_parameters.output_filename
                  << std::endl;
      }
    } else if (keyword == "set-indexing") {
      std::string orbital_filename;
      line_stream >> transition_type;
      ParsingCheck(line_stream, line_count, line);

      std::ifstream orbital_stream(orbital_filename);
      std::vector<basis::OrbitalPNInfo> input_orbitals =
          basis::ParseOrbitalPNStream(orbital_stream, true);
      run_parameters.space = basis::OrbitalSpaceLJPN(input_orbitals);
    } else if (keyword == "define-operator") {
      std::string name, filename;
      line_stream >> name >> filename;
      ParsingCheck(line_stream, line_count, line);
      struct stat st;
      if (stat(filename, &st) != 0) {
        std::cerr << "ERROR: file " << filename << " does not exist!" << std::endl;
        std::exit(EXIT_FAILURE);
      }

      run_parameters.operators.emplace_back(name, filename);
    } else if (keyword == "set-basis-scale-factor") {
      line_stream >> run_parameters.scale_factor;
      ParsingCheck(line_stream, line_count, line);
    } else if (keyword == "define-radial-source") {
      line_stream >> run_parameters.radial_me_filename;
      ParsingCheck(line_stream, line_count, line);
      struct stat st;
      if (stat(run_parameters.radial_me_filename.c_str(), &st) != 0) {
        std::cerr << "ERROR: file " << run_parameters.radial_me_filename
                  << " does not exist!" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    } else if (keyword == "set-charge") {
      std::string charge, species_code;
      double value;
      line_stream >> charge >> species_code >> value;
      ParsingCheck(line_stream, line_count, line);

      int species_index;
      if (species_code == "p" || species_code == "P") {
        species_index = 0;
      } else if (species_code == "n" || species_code == "N") {
        species_index = 1;
      }

      if (charge == "q") {
        run_parameters.electric_charge[species_index] = value;
      } else if (charge == "gl") {
        run_parameters.orbital_g[species_index] = value;
      } else if (charge == "gs") {
        run_parameters.g[species_index] = value;
      }
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
