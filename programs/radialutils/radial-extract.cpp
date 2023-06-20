/******************************************************************************/
/**
  @file radial-extract.cpp

  extract orbital transformation for consumption by Mathematica

  Syntax:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    radial-extract xform_file output_file
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Patrick J. Fasano
  University of Notre Dame

  + 03/04/19 (pjf): Created, based on radial-plot.
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

#include "am/halfint.h"
#include "am/halfint_fmt.h"
#include "fmt/format.h"
#include "basis/nlj_orbital.h"
#include "basis/operator.h"
#include "mcutils/parsing.h"
#include "obme/obme_io.h"

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

// Stores simple parameters for run
struct RunParameters {
  // filenames
  std::string xform_file;
  std::string output_file;
};

void PrintUsage(const char **argv) {
  std::cout << "Usage: " << argv[0]
            << " xform_file output_file"
            << std::endl;
}

void ProcessArguments(int argc, const char *argv[], RunParameters& run_parameters) {
  // argument counter
  int arg = 0;

  // usage message
  if (argc-1 != 2) {
    PrintUsage(argv);
    std::cerr << "ERROR: Incorrect number of arguments." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // xform file
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.xform_file;
    mcutils::FileExistCheck(run_parameters.xform_file, true, false);
  }

  // output file
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.output_file;
    mcutils::FileExistCheck(run_parameters.output_file, false, true);
  }
}

////////////////////////////////////////////////////////////////
// main program
////////////////////////////////////////////////////////////////
int main(int argc, const char *argv[]) {
  // header
  std::cout << std::endl;
  std::cout << "radial-extract -- radial function extraction" << std::endl;
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  RunParameters run_parameters;
  ProcessArguments(argc, argv, run_parameters);

  // get indexing and blocks
  basis::OrbitalSpaceLJPN bra_space, ket_space;
  basis::OrbitalSectorsLJPN sectors;
  basis::OperatorBlocks<double> xform_matrices;
  shell::InOBMEStream xform_s(run_parameters.xform_file);
  xform_s.SetToIndexing(bra_space, ket_space, sectors);
  xform_s.Read(xform_matrices);
  xform_s.Close();

  // check that xform is valid
  assert(xform_s.operator_type() == basis::OneBodyOperatorType::kRadial);
  assert(sectors.J0() == 0);
  assert(sectors.g0() == 0);
  assert(sectors.Tz0() == 0);

  // open output file for writing
  std::ofstream os(run_parameters.output_file);

  // loop over space
  for (std::size_t subspace_index=0; subspace_index<ket_space.size(); ++subspace_index)
  {
    // get indexing
    const auto& subspace = ket_space.GetSubspace(subspace_index);
    const std::size_t bra_subspace_index = bra_space.LookUpSubspaceIndex(subspace.labels());
    const auto& bra_subspace = bra_space.GetSubspace(bra_subspace_index);
    const std::size_t sector_index = sectors.LookUpSectorIndex(bra_subspace_index, subspace_index);
    const auto& xform_matrix = xform_matrices[sector_index];

    // output state information
    for (std::size_t state_index=0; state_index<subspace.size(); ++state_index)
    {
      const basis::OrbitalStateLJPN state(subspace, state_index);

      os << fmt::format(
          " {:4d} {:4d} {:4.1f} {:4.1f} {:6d}",
          state.n(), state.l(), double(state.j()), double(state.Tz()), bra_subspace.size()
        );
      for (std::size_t bra_state_index=0; bra_state_index<bra_subspace.size(); ++bra_state_index)
        os << fmt::format(" {:16.8e}", xform_matrix(bra_state_index, state_index));
      os << std::endl;
    }
  }

  os.close();

  /* return code */
  return EXIT_SUCCESS;
}
