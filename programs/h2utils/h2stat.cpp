/******************************************************************************

  h2stat.cpp -- generate statistics for MFDn H2 TBME file

  Syntax:
    h2stat input_filename

  Mark A. Caprio
  University of Notre Dame

  10/19/16 (mac): Created.

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
};

void ProcessArguments(int argc, char **argv, RunParameters& run_parameters)
{
  // usage message
  if (argc-1 < 1)
    {
      std::cout << "Syntax: h2stat input_filename" << std::endl;
      std::exit(EXIT_SUCCESS);
    }

  // input filename
  //std::string input_filename(argv[1]);
  run_parameters.input_filename = argv[1];
}

int main(int argc, char **argv)
{

  ////////////////////////////////////////////////////////////////
  // initialization
  ////////////////////////////////////////////////////////////////

  // header
  std::cout << std::endl;
  std::cout << "h2stat  -- MFDn H2 file statistics" << std::endl;
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

  ////////////////////////////////////////////////////////////////
  // statistics
  ////////////////////////////////////////////////////////////////

  std::cout << "Orbitals" << std::endl;
  std::cout << basis::OrbitalDefinitionStr(input_stream.orbital_space().OrbitalInfo());
  std::cout << std::endl;

  // iterate over sectors
  for (int sector_index = 0; sector_index < input_stream.num_sectors(); ++sector_index)
    {
      // progress indicator
      input_stream.SkipSector();
      std::cout << "." << std::flush;
    }
  std::cout << std::endl;


  ////////////////////////////////////////////////////////////////
  // termination
  ////////////////////////////////////////////////////////////////

  // explicitly force stream closure
  //   for neatness
  input_stream.Close();

  // end timing
  total_time.Stop();
  std::cout << "(Total time: " << total_time.ElapsedTime() << ")" << std::endl;
  std::cout << std::endl;

  // exit
  return EXIT_SUCCESS;
}
