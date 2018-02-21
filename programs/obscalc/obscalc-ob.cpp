/******************************************************************************
  obscalc-ob.cpp

  Calculate one-body observable from one-body operator matrix elements and
  one-body density matrix elements.

  Syntax:
    + obscalc-ob

  Input format:

    set-output-file output_filename
    set-indexing orbital_filename
    set-robdme-info robdme_info_filename
    define-operator name operator_filename
    define-static-densities 2J g n robdme_filename
    define-transition-densities 2Jf gf nf 2Ji gi fi robdme_filename

  Patrick J. Fasano
  University of Notre Dame

  + 10/08/17 (pjf): Created.
  + 10/23/17 (pjf): Rewrite for reading/writing many observables at one time.
  + 10/25/17 (pjf): Make robdme info filename global to run.
  + 11/28/17 (pjf): Include version in header.

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
#include "obme/obme_operator.h"
#include "obme/obme_io.h"

struct OneBodyOperator : public shell::OneBodyOperator {
  std::string filename;

  OneBodyOperator(const std::string& name_, const std::string& filename_) : filename(filename_)
  {
    name = name_;
  }
};

struct OneBodyDensities {
  HalfInt Ji, Jf;
  int gi, gf, ni, nf;
  std::string robdme_filename;
};

// Stores parameters for run
struct RunParameters {
  // filenames
  std::string output_filename;
  basis::OrbitalSpaceLJPN space;
  std::string robdme_info_filename;
  std::vector<OneBodyOperator> operators;
  std::vector<OneBodyDensities> static_densities;
  std::vector<OneBodyDensities> transition_densities;
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
      line_stream >> orbital_filename;
      ParsingCheck(line_stream, line_count, line);
      struct stat st;
      if (stat(orbital_filename.c_str(), &st) != 0) {
        std::cerr << "ERROR: file " << orbital_filename << " does not exist!"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }

      std::ifstream orbital_stream(orbital_filename);
      std::vector<basis::OrbitalPNInfo> input_orbitals =
          basis::ParseOrbitalPNStream(orbital_stream, true);
      run_parameters.space = basis::OrbitalSpaceLJPN(input_orbitals);
    } else if (keyword == "set-robdme-info") {
        std::string robdme_info_filename;
        line_stream >> robdme_info_filename;
        ParsingCheck(line_stream, line_count, line);
        struct stat st;
        if (stat(robdme_info_filename.c_str(), &st) != 0) {
          std::cerr << "ERROR: file " << robdme_info_filename
                    << " does not exist!" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        run_parameters.robdme_info_filename = robdme_info_filename;
    } else if (keyword == "define-operator") {
      std::string name, filename;
      line_stream >> name >> filename;
      ParsingCheck(line_stream, line_count, line);
      struct stat st;
      if (stat(filename.c_str(), &st) != 0) {
        std::cerr << "ERROR: file " << filename << " does not exist!" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      run_parameters.operators.emplace_back(name, filename);
    } else if (keyword == "define-static-densities") {
      OneBodyDensities densities;
      int twiceJ, g, n;
      std::string robdme_filename;
      line_stream >> twiceJ >> g >> n >> robdme_filename;
      ParsingCheck(line_stream, line_count, line);
      struct stat st;
      if (stat(robdme_filename.c_str(), &st) != 0) {
        std::cerr << "ERROR: file " << robdme_filename << " does not exist!"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      densities.Jf = densities.Ji = HalfInt(twiceJ, 2);
      densities.gf = densities.gi = g;
      densities.nf = densities.ni = n;
      densities.robdme_filename = robdme_filename;
      run_parameters.static_densities.push_back(densities);
    } else if (keyword == "define-transition-densities") {
      OneBodyDensities densities;
      int twiceJf, gf, nf, twiceJi, gi, ni;
      std::string robdme_filename;
      line_stream >> twiceJf >> gf >> nf >> twiceJi >> gi >> ni >> robdme_filename;
      ParsingCheck(line_stream, line_count, line);
      struct stat st;
      if (stat(robdme_filename.c_str(), &st) != 0) {
        std::cerr << "ERROR: file " << robdme_filename << " does not exist!"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      densities.Jf = HalfInt(twiceJf, 2);
      densities.gf = gf;
      densities.nf = nf;
      densities.Ji = HalfInt(twiceJi, 2);
      densities.gi = gi;
      densities.ni = ni;
      densities.robdme_filename = robdme_filename;
      run_parameters.transition_densities.push_back(densities);
    }
  }
}

double CalculateMatrixElement(const RunParameters& run_parameters,
                              const shell::OneBodyOperator& op,
                              const OneBodyDensities& densities) {
  basis::OperatorBlocks<double> density_matrices;
  const basis::OrbitalSpaceLJPN& space = run_parameters.space;
  const basis::OrbitalSectorsLJPN& sectors = op.sectors;
  const basis::OperatorBlocks<double>& operator_matrices = op.matrices;
  shell::InOBDMEReader obdme_reader(run_parameters.robdme_info_filename, space, 0, 0);
  obdme_reader.ReadMultipole(densities.robdme_filename, op.sectors.j0(),
                             density_matrices);

  // loop and sum over \sum_{a,b} rho_{ab} T_{ba}
  double value = 0;
  for (int subspace_index_a = 0; subspace_index_a < space.size();
       ++subspace_index_a) {
    for (int subspace_index_b = 0; subspace_index_b < space.size();
         ++subspace_index_b) {
      const auto subspace_a = space.GetSubspace(subspace_index_a);
      const auto subspace_b = space.GetSubspace(subspace_index_b);
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
  return (1./Hat(op.sectors.j0())) * value;
}

int main(int argc, char** argv) {
  // header
  std::cout << std::endl;
  std::cout << "obscalc-ob -- one-body observable evaluation" << std::endl;
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  // process input
  RunParameters run_parameters;
  ReadParameters(run_parameters);

  // parallel performance diagnostic
  std::cout << fmt::format("INFO: OMP max_threads {}, num_procs {}",
                           omp_get_max_threads(), omp_get_num_procs())
            << std::endl
            << std::endl;

  // read operators
  for (auto op : run_parameters.operators)
  {
    shell::InOBMEStream stream(op.filename);
    stream.GetOneBodyOperator(op);
    stream.Close();

    assert(op.bra_orbital_space.OrbitalInfo() == op.ket_orbital_space.OrbitalInfo());
    assert(op.operator_type == basis::OneBodyOperatorType::kSpherical);
  }

  // open output
  std::ios_base::openmode mode_argument = std::ios_base::trunc;
  std::ofstream out_stream(run_parameters.output_filename, mode_argument);
  StreamCheck(bool(out_stream), run_parameters.output_filename,
              "Failure opening file for output");

  out_stream << "[Observables]" << std::endl;
  out_stream << "names =";
  for (const auto& op : run_parameters.operators) {
    out_stream << " " << op.name;
  }
  out_stream << std::endl << std::endl;

  // static observables
  out_stream << "[Static one-body observables]" << std::endl;
  out_stream << fmt::format("# {:>2} {:>2} {:>2} ", "J", "g", "n");
  for (const auto& op : run_parameters.operators) {
    out_stream << fmt::format(" {:>15}", op.name);
  }
  out_stream << std::endl;

  for (const auto& densities : run_parameters.static_densities) {
    out_stream << fmt::format("  {:2d} {:2d} {:2d} ", float(densities.Jf),
                              densities.gf, densities.nf);

    for (const auto& op : run_parameters.operators) {
      out_stream << fmt::format(
          " {:15.8E}", CalculateMatrixElement(run_parameters, op, densities));
    }
    out_stream << std::endl;
  }
  out_stream << std::endl;

  // transition observables
  out_stream << "[Transition one-body observables]" << std::endl;
  out_stream << fmt::format("# {:>3} {:>3} {:>3}  {:>3} {:>3} {:>3} ", "Jf",
                            "gf", "nf", "Ji", "gi", "ni");
  for (const auto& op : run_parameters.operators) {
    out_stream << fmt::format(" {:>15}", op.name);
  }
  out_stream << std::endl;

  for (const auto& densities : run_parameters.transition_densities) {
    out_stream << fmt::format("  {:3d} {:3d} {:3d}  {:3d} {:3d} {:3d} ",
                              float(densities.Jf),
                              densities.gf,
                              densities.nf,
                              float(densities.Ji),
                              densities.gi,
                              densities.ni);

    for (const auto& op : run_parameters.operators) {
      out_stream << fmt::format(
          " {:15.8E}", CalculateMatrixElement(run_parameters, op, densities));
    }
    out_stream << std::endl;
  }
}
