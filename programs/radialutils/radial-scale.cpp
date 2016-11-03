/// @file
/******************************************************************************

  radial-scale.cpp -- rescale radial matrix elements with new oscillator length

  Syntax:
    radial-scale proton_length_ratio neutron_length_ratio input_file output_file

  Patrick J. Fasano
  University of Notre Dame

  11/4/16 (pjf): Created, based on radial-gen.cpp.

******************************************************************************/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "mcutils/profiling.h"
#include "basis/nlj_orbital.h"
#include "radial/radial_io.h"
#include "spline/wavefunction_class.h"

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
  RunParameters run_parameters;
  ProcessArguments(argc, argv, run_parameters);

  // Read input
  shell::InRadialStream is(run_parameters.input_filename);

  // check that this is an operator we can scale
  if (is.radial_operator_type() == shell::RadialOperatorType::kO) {
    std::cerr << "ERROR: Overlaps cannot be scaled. Exiting." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (is.sectors().Tz0() != 0) {
    std::cerr << "ERROR: Cannot scale pn sectors. Exiting." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  float proton_scale_factor = 1.;
  float neutron_scale_factor = 1.;
  if (is.radial_operator_type() == shell::RadialOperatorType::kR) {
    proton_scale_factor = std::pow(run_parameters.proton_scale, is.sectors().l0max());
    neutron_scale_factor = std::pow(run_parameters.neutron_scale, is.sectors().l0max());
  } else if (is.radial_operator_type() == shell::RadialOperatorType::kK) {
    proton_scale_factor = std::pow(run_parameters.proton_scale, -1*is.sectors().l0max());
    neutron_scale_factor = std::pow(run_parameters.neutron_scale, -1*is.sectors().l0max());
  }

  // get indexing
  const basis::OrbitalSectorsLJPN& sectors = is.sectors();

  // Eigen initialization
  basis::MatrixVector matrices;
  is.Read(matrices);

  // main loop
  #pragma omp parallel for
  for (int sector_index=0; sector_index < sectors.size(); ++sector_index) {
    if (sectors.GetSector(sector_index).bra_subspace().orbital_species() == basis::OrbitalSpeciesPN::kP) {
      matrices[sector_index] *= proton_scale_factor;
    } else if (sectors.GetSector(sector_index).bra_subspace().orbital_species() == basis::OrbitalSpeciesPN::kN) {
      matrices[sector_index] *= neutron_scale_factor;
    }
  }

  // write out to file
  shell::OutRadialStream os(run_parameters.output_filename,
                            is.bra_orbital_space(), is.ket_orbital_space(), sectors,
                            is.radial_operator_type());
  os.Write(matrices);

  is.Close();
  os.Close();

  return EXIT_SUCCESS;
}
