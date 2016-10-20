/******************************************************************************

  h2conv.cpp -- copy and convert MFDn H2 TBME file

  Syntax:
    h2stat input_filename output_filename format

  Mark A. Caprio
  University of Notre Dame

  10/19/16 (mac): Created, based on h2stat.

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
      std::cout << "Syntax: h2stat input_filename" << std::endl;
      std::exit(EXIT_SUCCESS);
    }

  // input filename
  run_parameters.input_filename = argv[1];
  run_parameters.output_filename = argv[2];
  std::istringstream parameter_stream(argv[3]);
  parameter_stream >> run_parameters.output_h2_format;
  if (!parameter_stream)
    {
      std::cerr << "Expecting numeric value for H2 format argument" << std::endl;
      std::exit(EXIT_FAILURE);
    }
}

int main(int argc, char **argv)
{

  ////////////////////////////////////////////////////////////////
  // initialization
  ////////////////////////////////////////////////////////////////

  // header
  std::cout << std::endl;
  std::cout << "h2conv -- MFDn H2 file conversion" << std::endl;
  std::cout << std::endl;

  // read parameters
  RunParameters run_parameters;
  ProcessArguments(argc,argv,run_parameters);

  // start timing
  Timer total_time;
  total_time.Start();

  // stream initialization

  std::cout << "Input stream" << std::endl;
  shell::InH2Stream input_stream(run_parameters.input_filename);
  std::cout << input_stream.DiagnosticStr() << std::endl;
  std::cout << std::endl;

  const basis::OrbitalSpacePN& orbital_space = input_stream.orbital_space();
  const basis::TwoBodySpaceJJJPN& space = input_stream.space();
  const basis::TwoBodySectorsJJJPN& sectors = input_stream.sectors();

  std::cout << "Output stream" << std::endl;
  shell::OutH2Stream output_stream(
      run_parameters.output_filename,
      orbital_space,space,sectors,
      run_parameters.output_h2_format
    );
  std::cout << output_stream.DiagnosticStr() << std::endl;
  std::cout << std::endl;

  ////////////////////////////////////////////////////////////////
  // copy sectors
  ////////////////////////////////////////////////////////////////

  // iterate over sectors
  for (int sector_index = 0; sector_index < input_stream.num_sectors(); ++sector_index)
    {
      Eigen::MatrixXd matrix;
      input_stream.ReadSector(matrix);
      output_stream.WriteSector(matrix);
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
