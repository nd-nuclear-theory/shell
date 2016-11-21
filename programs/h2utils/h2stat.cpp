/******************************************************************************

  h2stat.cpp -- generate statistics for MFDn H2 TBME file

  Syntax:
    h2stat mode [mode] input_filename

  Mark A. Caprio
  University of Notre Dame

  10/19/16 (mac): Created.
  10/22/16 (mac): Add mode argument.
  10/30/16 (mac): Update command line syntax.
  11/13/16 (mac): Add sector matrix output mode.

******************************************************************************/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

#include "cppformat/format.h"
#include "mcutils/eigen.h"
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
  if ((argc-1 < 1) || (std::string(argv[1])=="--help"))
    {
      std::cout
        << "Syntax: h2stat [mode] input_filename" << std::endl
        << std::endl
        // << "  --help: display this message" << std::endl
        << "  --header: display file header (default)" << std::endl
        << "  --orbitals: tabulate orbitals" << std::endl
        << "  --indexing: provide full orbital and two-body indexing information" << std::endl
        << "  --verify: read file to verify integrity" << std::endl
        << "  --matrices: display sectors as matrices" << std::endl
        << "      Sectors are labeled by [Tz J g] for bra and ket subspaces." << std::endl
        << std::endl;

      std::exit(EXIT_SUCCESS);
    }

  if (argc-1 == 1)
    // just filename
    {
      run_parameters.mode = "header";
      run_parameters.input_filename = argv[1];
    }
  else
    // mode + filename
    {
      // statistics mode
      if (std::string(argv[1]).substr(0,2)=="--")
        {
          run_parameters.mode = std::string(argv[1]).substr(2);
        }
      else
        {
          std::cerr << "Expected mode option beginning with '--'." << std::endl;
          std::exit(EXIT_FAILURE);
        }

      // input filename
      run_parameters.input_filename = argv[2];
    }
}

////////////////////////////////////////////////////////////////
// statistics modes
/////////////////////////////////////////////////////////////////

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

void DoMatrices(shell::InH2Stream& input_stream)
{
  std::cout << "Sector matrices" << std::endl;
  std::cout << std::endl;

  // iterate over sectors
  for (int sector_index = 0; sector_index < input_stream.num_sectors(); ++sector_index)
    {
      // read sector
      Eigen::MatrixXd matrix;
      input_stream.ReadSector(matrix);

      // head sector
      const typename basis::TwoBodySectorsJJJPN::SectorType& sector
        = input_stream.sectors().GetSector(sector_index);
      std::cout
        << fmt::format(
            "Sector {}: {} {}",
            sector_index,
            sector.bra_subspace().LabelStr(),
            sector.ket_subspace().LabelStr()
          )
        << std::endl
        << fmt::format("  Diagonal: {}",sector.IsDiagonal())
        << std::endl;
      std::cout << std::endl;

      // fill in lower triangle of matrix
      //
      // The full square matrix must be populated before remapping.
      //
      // Caution: We naively assume a symmetric matrix.  This is
      // appropriate for scalar operators, but this may need to change
      // to a more general phase relation for two-body nonscalar
      // operators.
      //
      // Limitation: For now we assume a scalar operator, until we
      // work out symmetry phase issues for nonscalar operators.
      assert(input_stream.sectors().J0()==0);
      assert(input_stream.sectors().g0()==0);
      assert(input_stream.sectors().Tz0()==0);
      if (sector.IsDiagonal())
        mcutils::CompleteLowerTriangle(matrix);

      // write matrix
      std::cout << mcutils::FormatMatrix(matrix,"+.8e") <<std::endl;
      std::cout << std::endl;
    }
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

  if (run_parameters.mode!="")
    {
      std::cout << fmt::format("Mode: {}",run_parameters.mode)
                << std::endl
                << std::endl;
    }

  if (run_parameters.mode=="header")
    {
      //pass
    }
  else if (run_parameters.mode=="orbitals")
    DoOrbitals(input_stream);
  else if (run_parameters.mode=="indexing")
    DoIndexing(input_stream);
  else if (run_parameters.mode=="verify")
    DoVerify(input_stream);
  else if (run_parameters.mode=="matrices")
    DoMatrices(input_stream);
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
