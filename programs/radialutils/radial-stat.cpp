/**************************************************************************/
/**
  @file radial-stat.cpp

  produce annotated radial files

  Syntax:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    radial-stat input_file output_file
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  @author Patrick J. Fasano
  University of Notre Dame

  + 08/11/17 (pjf): Created, based on radial-compose.cpp.
  + 10/12/17 (pjf): Update for changes to radial_io.
  + 11/28/17 (pjf): Print header with version.

******************************************************************************/

#include <sys/stat.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "mcutils/profiling.h"
#include "basis/nlj_orbital.h"
#include "basis/nlj_operator.h"
#include "obme/obme_io.h"

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

// Stores simple parameters for run
struct RunParameters {
  // filenames
  std::string input_filename;
  std::string output_filename;
};

void PrintUsage(const char *argv[]) {
  std::cout << "Usage: "
            << argv[0] << " input_file output_file"
            << std::endl;
}

void ProcessArguments(const int argc, const char *argv[], RunParameters& run_parameters) {
  // argument counter
  int arg = 0;

  // usage message
  if (argc-1 != 2) {
    PrintUsage(argv);
    std::cerr << "ERROR: Wrong number of arguments." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // input file
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.input_filename;
    if (!parameter_stream) {
      PrintUsage(argv);
      std::cerr << "ERROR: Invalid input filename." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    struct stat st;
    if (stat(run_parameters.input_filename.c_str(), &st) != 0) {
      std::cerr << "ERROR: file " << run_parameters.input_filename << " does not exist!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // output file
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.output_filename;
    if (!parameter_stream) {
      PrintUsage(argv);
      std::cerr << "ERROR: Invalid output filename." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    struct stat st;
    if (stat(run_parameters.output_filename.c_str(), &st) == 0) {
      std::cerr << "WARN: overwriting file " << run_parameters.output_filename << std::endl;
    }
  }
}

int main(int argc, const char *argv[]) {
  // header
  std::cout << std::endl;
  std::cout << "radial-stat -- radial matrix element file statistics" << std::endl;
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  RunParameters run_parameters;
  ProcessArguments(argc, argv, run_parameters);

  // Read input
  shell::InOBMEStream inputs(run_parameters.input_filename);

  // get indexing
  basis::OrbitalSpaceLJPN bra_space, ket_space;
  basis::OrbitalSectorsLJPN sectors;
  basis::OneBodyOperatorType operator_type = inputs.operator_type();
  shell::RadialOperatorType radial_operator_type = inputs.radial_operator_type();
  int operator_power = inputs.radial_operator_power();
  inputs.SetToIndexing(bra_space, ket_space, sectors);

  // Eigen initialization
  basis::OperatorBlocks<double> matrices;
  inputs.Read(matrices);
  inputs.Close();

  // write out to file
  std::cout << "INFO: Writing to file " << run_parameters.output_filename << std::endl;
  shell::OutOBMEStream os(run_parameters.output_filename,
                          bra_space, ket_space, sectors,
                          operator_type, radial_operator_type, operator_power, true);
                            // verbose_mode = true
  os.Write(matrices);
  os.Close();


  return EXIT_SUCCESS;
}
