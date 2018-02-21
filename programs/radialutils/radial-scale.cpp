/**************************************************************************//**
  @file radial-scale.cpp

  rescale radial matrix elements with new oscillator length

  Syntax:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    radial-scale proton_length_ratio neutron_length_ratio input_file output_file
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  @author Patrick J. Fasano
  University of Notre Dame

  + 11/4/16 (pjf): Created, based on radial-gen.cpp.
  + 10/12/17 (pjf): Update for changes to radial_io.
  + 11/28/17 (pjf): Print header with version.

******************************************************************************/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "mcutils/profiling.h"
#include "basis/nlj_orbital.h"
#include "obme/obme_io.h"

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

// Stores simple parameters for run
struct RunParameters {
  // mode
  float proton_scale;
  float neutron_scale;
  // filenames
  std::string input_filename;
  std::string output_filename;
};

void PrintUsage(char **argv) {
  std::cout << "Usage: " << argv[0]
            << " proton_length_ratio neutron_length_ratio input_file output_file"
            << std::endl;
}

void ProcessArguments(int argc, char **argv, RunParameters& run_parameters) {
  // usage message
  if (argc-1 != 4) {
    PrintUsage(argv);
    std::exit(EXIT_FAILURE);
  }

  // proton scale
  {
    std::istringstream parameter_stream(argv[1]);
    parameter_stream >> run_parameters.proton_scale;
    if (!parameter_stream) {
      PrintUsage(argv);
      std::cerr << "Invalid proton scaling factor" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // neutron scale
  {
    std::istringstream parameter_stream(argv[2]);
    parameter_stream >> run_parameters.neutron_scale;
    if (!parameter_stream) {
      PrintUsage(argv);
      std::cerr << "Invalid neutron scaling factor" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // input file
  {
    std::istringstream parameter_stream(argv[3]);
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
    std::istringstream parameter_stream(argv[4]);
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
  // header
  std::cout << std::endl;
  std::cout << "radial-scale -- analytic radial matrix element scaling" << std::endl;
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  RunParameters run_parameters;
  ProcessArguments(argc, argv, run_parameters);

  // Read input
  shell::InOBMEStream is(run_parameters.input_filename);

  // get indexing
  basis::OrbitalSpaceLJPN bra_orbital_space, ket_orbital_space;
  basis::OrbitalSectorsLJPN sectors;
  basis::OneBodyOperatorType operator_type = is.operator_type();
  shell::RadialOperatorType radial_operator_type = is.radial_operator_type();
  int radial_operator_power = is.radial_operator_power();
  is.SetToIndexing(bra_orbital_space, ket_orbital_space, sectors);

  // check that this is an operator we can scale
  if (radial_operator_type == shell::RadialOperatorType::kO) {
    std::cerr << "ERROR: Overlaps cannot be scaled. Exiting." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (radial_operator_type == shell::RadialOperatorType::kGeneric) {
    std::cerr << "ERROR: Generic operators cannot be scaled. Exiting." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (sectors.Tz0() != 0) {
    std::cerr << "ERROR: Cannot scale pn sectors. Exiting." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  float proton_scale_factor = 1.;
  float neutron_scale_factor = 1.;
  if (radial_operator_type == shell::RadialOperatorType::kR) {
    proton_scale_factor = std::pow(run_parameters.proton_scale, radial_operator_power);
    neutron_scale_factor = std::pow(run_parameters.neutron_scale, radial_operator_power);
  } else if (radial_operator_type == shell::RadialOperatorType::kK) {
    proton_scale_factor = std::pow(run_parameters.proton_scale, -1*radial_operator_power);
    neutron_scale_factor = std::pow(run_parameters.neutron_scale, -1*radial_operator_power);
  }

  // Eigen initialization
  basis::OperatorBlocks<double> matrices;
  is.Read(matrices);

  // main loop
  const int sectors_size = sectors.size();
  #pragma omp parallel for
  for (int sector_index=0; sector_index < sectors_size; ++sector_index) {
    if (sectors.GetSector(sector_index).bra_subspace().orbital_species() == basis::OrbitalSpeciesPN::kP) {
      matrices[sector_index] *= proton_scale_factor;
    } else if (sectors.GetSector(sector_index).bra_subspace().orbital_species() == basis::OrbitalSpeciesPN::kN) {
      matrices[sector_index] *= neutron_scale_factor;
    }
  }

  // write out to file
  shell::OutOBMEStream os(run_parameters.output_filename,
                          bra_orbital_space, ket_orbital_space, sectors,
                          operator_type, radial_operator_type, radial_operator_power);
  os.Write(matrices);

  is.Close();
  os.Close();

  return EXIT_SUCCESS;
}
