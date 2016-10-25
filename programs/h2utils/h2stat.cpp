/******************************************************************************

  h2stat.cpp -- generate statistics for MFDn H2 TBME file

  Syntax:
    h2stat mode input_filename

  Mark A. Caprio
  University of Notre Dame

  10/19/16 (mac): Created.
  10/22/16 (mac): Add mode argument.

******************************************************************************/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>

#include "cppformat/format.h"
#include "mcutils/profiling.h"
#include "tbme/h2_io.h"

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

struct RunParameters
// Stores simple parameters for run
{
  // statistics mode
  std::string mode;
  // filenames
  std::string input_filename;
};

void ProcessArguments(int argc, char **argv, RunParameters& run_parameters)
{
  // usage message
  if (argc-1 < 2)
    {
      std::cout << "Syntax: h2stat mode input_filename" << std::endl;
      std::exit(EXIT_SUCCESS);
    }

  // statistics mode
  run_parameters.mode = argv[1];

  // input filename
  run_parameters.input_filename = argv[2];
}

////////////////////////////////////////////////////////////////
// statistics modes
/////////////////////////////////////////////////////////////////

void DoVerify(shell::InH2Stream& input_stream)
{
  std::cout << "Verification scan" << std::endl;

  // iterate over sectors
  for (int sector_index = 0; sector_index < input_stream.num_sectors(); ++sector_index)
    {
      // read sector
      input_stream.SkipSector();

      // progress indicator
      std::cout << "." << std::flush;
    }
  std::cout << std::endl;
}

void DoOrbitals(shell::InH2Stream& input_stream)
{
  std::cout << "Orbitals" << std::endl;
  std::cout << basis::OrbitalDefinitionStr(input_stream.orbital_space().OrbitalInfo());
  std::cout << std::endl;
}

void DoIndexing(shell::InH2Stream& input_stream)
{
  std::cout << "Orbital space" << std::endl
            << std::endl;
  const basis::OrbitalSpacePN& orbital_space = input_stream.orbital_space();

  std::cout << " Subspaces" << std::endl;
  std::cout << orbital_space.DebugStr();
  std::cout << std::endl;

  for (int subspace_index=0; subspace_index<orbital_space.size(); ++subspace_index)
    {
      std::cout << fmt::format(" Subspace {} {}",subspace_index,orbital_space.GetSubspace(subspace_index).LabelStr()) << std::endl;
      std::cout << orbital_space.GetSubspace(subspace_index).DebugStr();
      std::cout << std::endl;
    }

  std::cout << "Two-body space" << std::endl
            << std::endl;
  const basis::TwoBodySpaceJJJPN& space = input_stream.space();

  std::cout << " Subspaces" << std::endl;
  std::cout << space.DebugStr();
  std::cout << std::endl;

  for (int subspace_index=0; subspace_index<space.size(); ++subspace_index)
    {
      std::cout << fmt::format(" Subspace {} {}",subspace_index,space.GetSubspace(subspace_index).LabelStr()) << std::endl;
      std::cout << space.GetSubspace(subspace_index).DebugStr();
      std::cout << std::endl;
    }

  // std::cout << "Two-body sectors" << std::endl;
  // std::cout << input_stream.sectors().DebugStr();
  // std::cout << std::endl;
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
  std::cout << input_stream.DiagnosticStr();
  std::cout << std::endl;

  ////////////////////////////////////////////////////////////////
  // statistics
  ////////////////////////////////////////////////////////////////

  std::cout << fmt::format("Statistics mode: {}",run_parameters.mode)
            << std::endl
            << std::endl;

  std::cout << "----------------------------------------------------------------"
            << std::endl
            << std::endl;

  if (run_parameters.mode=="verify")
    DoVerify(input_stream);
  else if (run_parameters.mode=="orbitals")
    DoOrbitals(input_stream);
  else if (run_parameters.mode=="indexing")
    DoIndexing(input_stream);
  else
    {
      std::cout << "Unrecognized statistics mode specified" << std::endl;
      std::exit(EXIT_FAILURE);
    }


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
