/**************************************************************************//**
  @file radial-xform.cpp

  transform single-particle matrix elements

  Syntax:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    radial-xform input_olap_file input_me_file output_me_file
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  @author Patrick J. Fasano
  University of Notre Dame

  + 11/7/16 (pjf): Created, based on radial-scale.cpp.

******************************************************************************/

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "mcutils/profiling.h"
#include "basis/nlj_orbital.h"
#include "radial/radial_io.h"

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

// Stores simple parameters for run
struct RunParameters {
  // filenames
  std::string olap_filename;
  std::string input_filename;
  std::string output_filename;
};

void PrintUsage(const char *argv[]) {
  std::cout << "Usage: "
            << argv[0] << " input_olap_file input_me_file output_me_file"
            << std::endl;
}

void ProcessArguments(const int argc, const char *argv[], RunParameters& run_parameters) {
  // usage message
  if (argc-1 != 3) {
    PrintUsage(argv);
    std::cerr << "Wrong number of arguments." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // olap file
  {
    std::istringstream parameter_stream(argv[1]);
    parameter_stream >> run_parameters.olap_filename;
    if (!parameter_stream) {
      PrintUsage(argv);
      std::cerr << "Invalid input filename." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    /** @todo check for file existence */
  }

  // input file
  {
    std::istringstream parameter_stream(argv[2]);
    parameter_stream >> run_parameters.input_filename;
    if (!parameter_stream) {
      PrintUsage(argv);
      std::cerr << "Invalid input filename." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    /** @todo check for file existence */
  }

  // output file
  {
    std::istringstream parameter_stream(argv[3]);
    parameter_stream >> run_parameters.output_filename;
    if (!parameter_stream) {
      PrintUsage(argv);
      std::cerr << "Invalid output filename." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    /** @todo check and warn if overwriting existing file */
  }
}

int main(int argc, char **argv) {
  RunParameters run_parameters;
  ProcessArguments(argc, argv, run_parameters);

  // Read input
  shell::InRadialStream is(run_parameters.input_filename);
  shell::InRadialStream olaps(run_parameters.olap_filename);

  // get indexing
  basis::OrbitalSpaceLJPN in_bra_space, in_ket_space;
  basis::OrbitalSectorsLJPN in_sectors;
  shell::RadialOperatorType in_operator_type = is.radial_operator_type();
  is.SetToIndexing(in_bra_space, in_ket_space, in_sectors);

  basis::OrbitalSpaceLJPN olap_bra_space, olap_ket_space;
  basis::OrbitalSectorsLJPN olap_sectors;
  shell::RadialOperatorType olap_type = olaps.radial_operator_type();
  xs.SetToIndexing(olap_bra_space, olap_ket_space, olap_sectors);

  // Eigen initialization
  basis::MatrixVector olap_matrices, input_matrices, output_matrices;
  is.Read(input_matrices);
  olaps.Read(olap_matrices);

  if ((olap_type != shell::RadialOperatorType::kOverlaps)
      (in_bra_space.OrbitalInfo() != in_ket_space.OrbitalInfo()) ||
      (in_ket_space.OrbitalInfo() != olap_bra_space.OrbitalInfo()))
  {
    std::cerr << "Invalid olap for this input." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // construct new sectors
  basis::OrbitalSectorsLJPN out_sectors(olap_ket_space, olap_ket_space, in_sectors.l0max(), in_sectors.Tz0());

  // main loop
  for (int sector_index=0; sector_index < out_sectors.size(); ++sector_index) {
    output_matrices.push_back(
      olap_matrices[sector_index].transpose() * input_matrices[sector_index] * olap_matrices[sector_index]
    );
  }

  // write out to file
  shell::OutRadialStream os(run_parameters.output_filename,
                            olap_ket_space, olap_ket_space, out_sectors,
                            in_operator_type);
  os.Write(output_matrices);

  is.Close();
  os.Close();

  return EXIT_SUCCESS;
}
