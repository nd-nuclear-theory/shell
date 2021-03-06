/******************************************************************************

  h2conv.cpp -- copy and convert MFDn H2 TBME file

  Syntax:
    h2conv format input_filename output_filename

  Mark A. Caprio
  University of Notre Dame

  + 10/19/16 (mac): Created, based on h2stat.
  + 10/22/16 (mac): Update syntax.
  + 11/28/17 (pjf): Include version in header.
  + 02/21/19 (pjf): Add H2 Version15200 support.
  + 05/09/19 (pjf): Use std::size_t for basis indices and sizes.

******************************************************************************/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>

#include "mcutils/profiling.h"
#include "tbme/h2_io.h"

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

struct RunParameters
// Stores simple parameters for run
{
  // filenames
  std::string input_filename;
  std::string output_filename;
  // mode
  shell::H2Format output_h2_format;
};

void ProcessArguments(int argc, char **argv, RunParameters& run_parameters)
{
  // usage message
  if (argc-1 < 3)
    {
      std::cout << "Syntax: h2conv format input_filename output_filename" << std::endl;
      std::exit(EXIT_SUCCESS);
    }

  // format
  std::istringstream parameter_stream(argv[1]);
  parameter_stream >> run_parameters.output_h2_format;
  if (!parameter_stream)
    {
      std::cerr << "Expecting numeric value for H2 format argument" << std::endl;
      std::exit(EXIT_FAILURE);
    }

  // input filename
  run_parameters.input_filename = argv[2];

  // output filename
  run_parameters.output_filename = argv[3];
}

////////////////////////////////////////////////////////////////
// main program
/////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  ////////////////////////////////////////////////////////////////
  // initialization
  ////////////////////////////////////////////////////////////////

  // header
  std::cout << std::endl;
  std::cout << "h2conv -- MFDn H2 file conversion" << std::endl;
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  // read parameters
  RunParameters run_parameters;
  ProcessArguments(argc,argv,run_parameters);

  // start timing
  mcutils::SteadyTimer total_time;
  total_time.Start();

  // stream initialization

  std::cout << "Input stream" << std::endl;
  shell::InH2Stream input_stream(run_parameters.input_filename);
  std::cout << input_stream.DiagnosticStr();
  std::cout << std::endl;

  const basis::OrbitalSpacePN& orbital_space = input_stream.orbital_space();
  const basis::TwoBodySpaceJJJPN& input_space = input_stream.space();
  const basis::TwoBodySectorsJJJPN& input_sectors = input_stream.sectors();
  basis::TwoBodySpaceJJJPNOrdering space_ordering =
    shell::kH2SpaceOrdering.at(run_parameters.output_h2_format);
  const basis::TwoBodySpaceJJJPN space = basis::TwoBodySpaceJJJPN(
      orbital_space,
      input_space.weight_max(),
      space_ordering
    );
  const basis::TwoBodySectorsJJJPN sectors = basis::TwoBodySectorsJJJPN(
      space, input_sectors.J0(), input_sectors.g0(), input_sectors.Tz0()
    );

  std::cout << "Output stream" << std::endl;
  shell::OutH2Stream output_stream(
      run_parameters.output_filename,
      orbital_space,space,sectors,
      run_parameters.output_h2_format
    );
  std::cout << output_stream.DiagnosticStr();
  std::cout << std::endl;

  ////////////////////////////////////////////////////////////////
  // copy sectors
  ////////////////////////////////////////////////////////////////

  // iterate over sectors
  for (std::size_t sector_index = 0; sector_index < input_stream.num_sectors(); ++sector_index)
    {
      Eigen::MatrixXd matrix;
      // construct target sector
      const auto& sector = sectors.GetSector(sector_index);
      // locate corresponding input sector
      std::size_t input_bra_subspace_index
        = input_space.LookUpSubspaceIndex(sector.bra_subspace().labels());
      std::size_t input_ket_subspace_index
        = input_space.LookUpSubspaceIndex(sector.ket_subspace().labels());
      // Note: We cannot simply look up by target_sector's Key, since that uses target subspace indices.
      std::size_t input_sector_index
        = input_sectors.LookUpSectorIndex(input_bra_subspace_index,input_ket_subspace_index);
      input_stream.ReadSector(input_sector_index, matrix);
      output_stream.WriteSector(sector_index, matrix);
      std::cout << "." << std::flush;
    }
  std::cout << std::endl;


  ////////////////////////////////////////////////////////////////
  // termination
  ////////////////////////////////////////////////////////////////

  // explicitly force stream closure
  //   for neatness
  input_stream.Close();
  output_stream.Close();

  // end timing
  total_time.Stop();
  std::cout << "(Total time: " << total_time.ElapsedTime() << ")" << std::endl;
  std::cout << std::endl;

  // exit
  return EXIT_SUCCESS;
}
