/******************************************************************************
  obscalc-ob.cpp

  Calculate one-body observable from one-body operator matrix elements and
  one-body density matrix elements.

  Note: Output RMEs are in Edmonds convention, in order to match MFDn's output.

  Syntax:
    + obscalc-ob

  Input format:

    set-output-file output_filename
    set-indexing orbital_filename
    define-operator name operator_filename
    define-static-densities 2J g n robdme_filename [robdme_info_filename]
    define-transition-densities 2Jf gf nf 2Ji gi ni robdme_filename [robdme_info_filename]

  Patrick J. Fasano
  University of Notre Dame

  + 10/08/17 (pjf): Created.
  + 10/23/17 (pjf): Rewrite for reading/writing many observables at one time.
  + 10/25/17 (pjf): Make robdme info filename global to run.
  + 11/28/17 (pjf): Include version in header.
  + 07/27/18 (pjf): Update for new OBDME input routines.
  + 03/30/19 (pjf): Add support for single-file ROBDME formats.
  + 04/03/19 (pjf): Do all calculations in Rose convention, *except for final
    output*.
  + 05/09/19 (pjf): Use std::size_t for basis indices and sizes.
  + 05/28/19 (pjf): Output zero matrix element if densities missing from file.
  + 05/30/19 (pjf): Reduce file I/O requirements by refactoring to pass  compute
    all operators for a given set of densities.
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
#include "fmt/format.h"
#include "density/obdme_io.h"
#include "obme/obme_operator.h"
#include "obme/obme_io.h"

// Store one-body operators
struct OneBodyOperator {
  std::string name;
  basis::OrbitalSpaceLJPN space;
  basis::OrbitalSectorsLJPN sectors;
  basis::OperatorBlocks<double> blocks;

  explicit OneBodyOperator(const std::string& name__, const std::string& filename)
      : name(name__)
  {
    // read operator
    shell::InOBMEStream is(filename);
    is.Read(blocks);

    // get indexing
    basis::OrbitalSpaceLJPN ket_space;
    is.SetToIndexing(space, ket_space, sectors);
    assert(space.OrbitalInfo() == ket_space.OrbitalInfo());
    is.Close();
  }
};

struct OneBodyDensities {
  HalfInt Ji, Jf;
  int gi, gf, ni, nf;
  std::string robdme_filename;
  std::string robdme_info_filename;
};

// Stores parameters for run
struct RunParameters {
  // filenames
  std::string output_filename;
  basis::OrbitalSpaceLJPN space;
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
      mcutils::ParsingCheck(line_stream, line_count, line);
      struct stat st;
      if (stat(run_parameters.output_filename.c_str(), &st) == 0) {
        std::cerr << "WARN: overwriting file " << run_parameters.output_filename
                  << std::endl;
      }
    } else if (keyword == "set-indexing") {
      std::string orbital_filename;
      line_stream >> orbital_filename;
      mcutils::ParsingCheck(line_stream, line_count, line);
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
    } else if (keyword == "define-operator") {
      std::string name, filename;
      line_stream >> name >> filename;
      mcutils::ParsingCheck(line_stream, line_count, line);
      mcutils::FileExistCheck(filename, true, false);
      run_parameters.operators.emplace_back(name, filename);
    } else if (keyword == "define-static-densities") {
      OneBodyDensities densities;
      int twiceJ, g, n;
      std::string robdme_filename, robdme_info_filename="";
      line_stream >> twiceJ >> g >> n >> robdme_filename;
      mcutils::ParsingCheck(line_stream, line_count, line);
      mcutils::FileExistCheck(robdme_filename, true, false);

      if (line_stream.good()) {
        line_stream >> robdme_info_filename;
        mcutils::ParsingCheck(line_stream, line_count, line);
        mcutils::FileExistCheck(robdme_info_filename, true, false);
      }
      densities.Jf = densities.Ji = HalfInt(twiceJ, 2);
      densities.gf = densities.gi = g;
      densities.nf = densities.ni = n;
      densities.robdme_filename = robdme_filename;
      run_parameters.static_densities.push_back(densities);
    } else if (keyword == "define-transition-densities") {
      OneBodyDensities densities;
      int twiceJf, gf, nf, twiceJi, gi, ni;
      std::string robdme_filename, robdme_info_filename="";
      line_stream >> twiceJf >> gf >> nf >> twiceJi >> gi >> ni >> robdme_filename;
      mcutils::ParsingCheck(line_stream, line_count, line);
      mcutils::FileExistCheck(robdme_filename, true, false);
      if (line_stream.good()) {
        line_stream >> robdme_info_filename;
        mcutils::ParsingCheck(line_stream, line_count, line);
        mcutils::FileExistCheck(robdme_info_filename, true, false);
      }

      densities.Jf = HalfInt(twiceJf, 2);
      densities.gf = gf;
      densities.nf = nf;
      densities.Ji = HalfInt(twiceJi, 2);
      densities.gi = gi;
      densities.ni = ni;
      densities.robdme_filename = robdme_filename;
      densities.robdme_info_filename = robdme_info_filename;
      run_parameters.transition_densities.push_back(densities);
    }
  }
}

std::vector<double> CalculateMatrixElements(
    const RunParameters& run_parameters,
    const std::vector<OneBodyOperator>& operators,
    const OneBodyDensities& densities
  )
{
  basis::OrbitalSectorsLJPN density_sectors;
  basis::OperatorBlocks<double> density_blocks;
  const basis::OrbitalSpaceLJPN& space = run_parameters.space;
  shell::InOBDMEStream obdme_s;
  if (densities.robdme_info_filename == "")
  {
    // TODO(pjf): check that density quantum numbers match input quantum numbers
    obdme_s = shell::InOBDMEStreamSingle(densities.robdme_filename, space);
  }
  else
  {
    obdme_s = shell::InOBDMEStreamMulti(
        densities.robdme_info_filename,
        densities.robdme_filename,
        space, /*g0=*/0, /*Tz0=*/0
      );
  }

  std::vector<double> return_values;
  for (const auto& op : operators)
  {
    // convenience variables
    const basis::OrbitalSectorsLJPN sectors = op.sectors;
    const basis::OperatorBlocks<double>& operator_blocks = op.blocks;
    if ((op.sectors.j0() < obdme_s.j0_min()) || (op.sectors.j0() > obdme_s.j0_max()))
    {
      // output zero if obdmes missing
      return_values.push_back(0.);
      continue;
    }
    obdme_s.GetMultipole(op.sectors.j0(), density_sectors, density_blocks);

    // loop and sum over \sum_{a,b} rho_{ab} T_{ba}
    double value = 0.;
    for (std::size_t subspace_index_a = 0; subspace_index_a < space.size();
        ++subspace_index_a) {
      for (std::size_t subspace_index_b = 0; subspace_index_b < space.size();
          ++subspace_index_b) {
        const auto subspace_a = space.GetSubspace(subspace_index_a);
        const auto subspace_b = space.GetSubspace(subspace_index_b);
        auto sector_index =
            sectors.LookUpSectorIndex(subspace_index_a, subspace_index_b);
        if (sector_index == basis::kNone) continue;

        // dimension factor only present in Rose convention
        double dimension_factor = double(2*subspace_b.j()+1);

        for (std::size_t state_index_a = 0; state_index_a < subspace_a.size(); ++state_index_a) {
          for (std::size_t state_index_b = 0; state_index_b < subspace_b.size(); ++state_index_b) {
            value += dimension_factor * operator_blocks[sector_index](state_index_a, state_index_b)
                    * density_blocks[sector_index](state_index_a, state_index_b);
          }
        }
      }
    }
    // convert to Edmonds convention
    value /= Hat(op.sectors.j0());
    // store value for return
    return_values.push_back(value);
  }
  return return_values;
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
  out_stream << fmt::format("# {:>4} {:>2} {:>2} ", "J", "g", "n");
  for (const auto& op : run_parameters.operators) {
    out_stream << fmt::format(" {:>15}", op.name);
  }
  out_stream << std::endl << std::flush;

  for (const auto& densities : run_parameters.static_densities) {
    out_stream << fmt::format(
        "  {:4.1f} {:2d} {:2d} ",
        float(densities.Jf), densities.gf, densities.nf
      );

    auto matrix_elements = CalculateMatrixElements(
        run_parameters,
        run_parameters.operators,
        densities
      );
    for (const auto& matrix_element : matrix_elements)
    {
      out_stream << fmt::format(" {:15.8E}", matrix_element);
    }
    out_stream << std::endl;
  }
  out_stream << std::endl;

  // transition observables
  out_stream << "[Transition one-body observables]" << std::endl;
  out_stream << fmt::format(
      "# {:>4} {:>3} {:>3}  {:>4} {:>3} {:>3} ",
      "Jf", "gf", "nf",
      "Ji", "gi", "ni"
    );
  for (const auto& op : run_parameters.operators) {
    out_stream << fmt::format(" {:>15}", op.name);
  }
  out_stream << std::endl;

  for (const auto& densities : run_parameters.transition_densities) {
    out_stream << fmt::format(
        "  {:4.1f} {:3d} {:3d}  {:4.1f} {:3d} {:3d} ",
        float(densities.Jf), densities.gf, densities.nf,
        float(densities.Ji), densities.gi, densities.ni
      );

    auto matrix_elements = CalculateMatrixElements(
        run_parameters,
        run_parameters.operators,
        densities
      );
    for (const auto& matrix_element : matrix_elements)
    {
      out_stream << fmt::format(" {:15.8E}", matrix_element);
    }
    out_stream << std::endl;
  }
}
