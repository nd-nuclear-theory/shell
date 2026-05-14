/******************************************************************************
  obscalc-ob.cpp

  Calculate one-body observable from one-body operator matrix elements and
  one-body density matrix elements.

  Note: Output RMEs are in Edmonds convention, in order to match MFDn's output.

  Syntax:
    + obscalc-ob

  Input format:

    set-output-file output_filename [append]
    set-indexing orbital_filename
    define-operator name operator_filename
    define-densities Jf gf nf Ji gi ni robdme_filename [robdme_info_filename]

  Example output:

    [One-body observable]
    #  J0  g0 Tz0  name
        2   0   0  E2p
    #   Jf  gf  nf    Ji  gi  ni              rme
       0.5   1   1   1.5   1   1   2.10359485e+00
       0.5   1   1   1.5   1   2   9.99786880e+00
    ...

    [One-body observable]
    #  J0  g0 Tz0  name
        1   0   0  Dlp
    #   Jf  gf  nf    Ji  gi  ni              rme
       0.5   1   1   0.5   1   1   2.33091440e-01
       0.5   1   1   1.5   1   1   2.52027931e-01
    ...

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
  + 08/17/19 (pjf):
    - Fix Rose convention.
    - Fix input parsing.
  + 04/08/20 (pjf): Rewrite.
    - Remove distinction between "static" and "transition" densities.
    - Suppress output if transition cannot be calculated due to Clebsch zero.
    - New output format to use simple listings rather than tables.
    - Reorder loops over operators and densities.
  + 08/14/20 (pjf): Add append mode to allow repeated calls to obscalc-ob.
  + 10/12/23 (pjf): Fix for corrected normalization of OBDMEs.
  + 05/14/26 (mac): Refactor matrix element calculation to obme/ob_observable.

******************************************************************************/

#include <sys/stat.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "am/halfint.h"
#include "am/wigner_gsl.h"
#include "basis/nlj_orbital.h"
#include "basis/proton_neutron.h"
#include "fmt/format.h"
#include "density/obdme_io.h"
#include "mcutils/parsing.h"
#include "obme/obme_operator.h"
#include "obme/obme_io.h"
#include "obme/ob_observable.h"

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

// Stores parameters for run
struct RunParameters {
  // filenames
  std::string output_filename;
  std::ios_base::openmode file_mode;
  basis::OrbitalSpaceLJPN space;
  std::vector<OneBodyOperator> operators;
  std::vector<std::unique_ptr<shell::InOBDMEStream>> density_streams;
};

void PrintUsage(char** argv) { std::cout << "Usage: " << argv[0] << std::endl; }

void ReadParameters(RunParameters& run_parameters) {
  std::string line;
  int line_count = 0;

  while (mcutils::GetLine(std::cin, line, line_count)) {
    // set up for line parsing
    std::istringstream line_stream(line);
    std::string keyword;
    line_stream >> keyword;

    // select action based on keyword
    if (keyword == "set-output-file") {
      line_stream >> run_parameters.output_filename;
      run_parameters.file_mode = std::ios_base::trunc;
      if (!line_stream.eof()) {
        std::string mode;
        line_stream >> mode;
        if (mode == "append")
          run_parameters.file_mode = std::ios_base::app;
      }
      mcutils::ParsingCheck(line_stream, line_count, line);
      mcutils::FileExistCheck(
        run_parameters.output_filename, false,
        run_parameters.file_mode==std::ios_base::trunc
        );
    } else if (keyword == "set-indexing") {
      std::string orbital_filename;
      line_stream >> orbital_filename;
      mcutils::ParsingCheck(line_stream, line_count, line);
      mcutils::FileExistCheck(orbital_filename, true, false);

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
    } else if (keyword == "define-densities") {
      std::unique_ptr<shell::InOBDMEStream> density_stream;
      float Jf, Ji;
      int gf, nf, gi, ni;
      std::string robdme_filename, robdme_info_filename="";
      line_stream >> Jf >> gf >> nf >> Ji >> gi >> ni >> robdme_filename;
      mcutils::ParsingCheck(line_stream, line_count, line);
      mcutils::FileExistCheck(robdme_filename, true, false);
      if (!line_stream.eof()) {
        line_stream >> robdme_info_filename;
        mcutils::ParsingCheck(line_stream, line_count, line);
        mcutils::FileExistCheck(robdme_info_filename, true, false);

        // construct multi-file stream
        run_parameters.density_streams.emplace_back(new shell::InOBDMEStreamMulti(
            robdme_info_filename, robdme_filename, run_parameters.space,
            HalfInt(2*Jf,2), gf, nf, HalfInt(2*Ji, 2), gi, ni
          ));
      } else {
        // construct single-file stream
        run_parameters.density_streams.emplace_back(new shell::InOBDMEStreamSingle(
            robdme_filename, run_parameters.space,
            HalfInt(2*Jf,2), gf, nf, HalfInt(2*Ji, 2), gi, ni
          ));

      }
    } else {
        mcutils::ParsingError(line_count,line,"Unrecognized keyword");
    }


  }
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
  std::ofstream out_stream(run_parameters.output_filename, run_parameters.file_mode);
  mcutils::StreamCheck(
      bool(out_stream), run_parameters.output_filename,
      "Failure opening file for output"
    );

  for (const auto& op : run_parameters.operators)
  {
    std::ostringstream section_stream;
    section_stream << "[One-body observable]" << std::endl;
    section_stream << fmt::format(
        "# {:>3} {:>3} {:>3}  {:s}",
        "J0", "g0", "Tz0", "name"
      ) << std::endl;
    section_stream << fmt::format(
        "  {:>3d} {:>3d} {:>3d}  {:s}",
        op.sectors.J0(), op.sectors.g0(), op.sectors.Tz0(), op.name
      ) << std::endl;
    section_stream << fmt::format(
        "# {:>4} {:>3} {:>3}  {:>4} {:>3} {:>3}  {:>15s}",
        "Jf", "gf", "nf",
        "Ji", "gi", "ni",
        "rme"
      ) << std::endl;

    std::size_t count = 0;
    for (const auto& density_stream : run_parameters.density_streams)
    {
      auto matrix_element = shell::CalculateOneBodyObservableMatrixElement(
          // run_parameters.space, op.sectors, op.blocks, density_stream
        );
      if (std::isnan(matrix_element)) continue;
      ++count;
      section_stream << fmt::format(
          "  {:>4.1f} {:>3d} {:>3d}  {:>4.1f} {:>3d} {:>3d}  {:15.8e}",
          float(density_stream->J_bra()), density_stream->g_bra(), density_stream->n_bra(),
          float(density_stream->J_ket()), density_stream->g_ket(), density_stream->n_ket(),
          matrix_element
        ) << std::endl;
    }

    // extra newline separating sections
    section_stream << std::endl;

    // if count is greater than 0, write section to file
    if (count > 0)
      out_stream << section_stream.str() << std::flush;
    else
      std::cerr << fmt::format(
        "WARN: no observables calculated for operator \"{:s}\"", op.name
       ) << std::endl;

  }
}
