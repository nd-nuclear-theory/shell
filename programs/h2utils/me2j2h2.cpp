/******************************************************************************

  me2j2h2.cpp -- ME2J to H2  TBME file conversion

  Syntax:

    me2j2h2 rank cutoff input_filename output_filename

    The rank may be specified as "ob" or "tb", in which case cutoff represents
    N1max or N2max, respectively.

  Zhou Zhou
  University of Notre Dame

  + Created by zz ~01/16/24.
  + 02/04/25 (mac): Generalize to handle either one-body or two-body truncation.

******************************************************************************/

#include <fstream>

#include "basis/jjjt_operator.h"
#include "basis/jjjpn_scheme.h"
#include "basis/jjjpn_operator.h"
#include "mcutils/eigen.h"
#include "mcutils/parsing.h"
#include "mcutils/profiling.h"
#include "tbme/h2_io.h"
#include "tbme/me2j_io.h"
#include "tbme/tbme_scheme_xform.h"

#include "moshinsky/moshinsky_xform.h"

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

struct RunParameters
// Stores simple parameters for run
{
  // filenames
  std::string input_filename;
  std::string output_filename;

  // file format
  std::size_t float_size;

  // truncation
  basis::Rank truncation_rank;
  int truncation_cutoff;

  // default constructor
  RunParameters()
    : input_filename(""), output_filename(""), float_size(4)
  {}

};

void PrintUsage(const char **argv) {
  std::cout << "Usage: " << argv[0]
            << " [--float-size 4|8]"
            << " rank cutoff input_filename output_filename"
            << std::endl;
}

void ProcessArguments(int argc, const char *argv[], RunParameters& run_parameters)
{
  
  int arg = 1;

  // process options
  while (arg < argc && argv[arg][0] == '-')
    {
      std::istringstream parameter_stream(argv[arg++]);

      if (parameter_stream.str() == "--help" || parameter_stream.str() == "-h")
        {
          PrintUsage(argv);
          std::exit(EXIT_SUCCESS);
        }
      else if (parameter_stream.str() == "--float-size")
        {
          if (argc-arg < 1)
            {
              PrintUsage(argv);
              std::cerr << "Insufficient arguments for --float-size" << std::endl;
              std::exit(EXIT_FAILURE);
            }

          std::size_t float_size;
          std::istringstream float_size_stream(argv[arg++]);
          float_size_stream >> float_size;
          if (!float_size_stream || ! ((float_size == 4) || (float_size == 8))) {
            PrintUsage(argv);
            std::cerr << "Invalid float_size" << std::endl;
            std::exit(EXIT_FAILURE);
          }
          run_parameters.float_size = float_size;
        }
      else
        {
          PrintUsage(argv);
          std::cerr << "Unrecognized option '" << parameter_stream.str() << "'" << std::endl;
          std::exit(EXIT_FAILURE);
        }
    }
  
  // process fixed arguments
  if (argc-arg < 4)
    {
      PrintUsage(argv);
      std::cerr << "Insufficient arguments" << std::endl;
      std::exit(EXIT_FAILURE);
    }

  // rank
  std::string parameter(argv[arg++]);
  if (parameter=="ob")
    run_parameters.truncation_rank = basis::Rank::kOneBody;
  else if (parameter=="tb")
    run_parameters.truncation_rank = basis::Rank::kTwoBody;
  else
    {
      std::cerr << "ERROR: Expecting ob or tb for truncation rank" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    
  // cutoff
  std::istringstream parameter_stream(argv[arg++]);
  parameter_stream >> run_parameters.truncation_cutoff;
  if (!parameter_stream)
    {
      std::cerr << "ERROR: Expecting numeric value for truncation cutoff" << std::endl;
      std::exit(EXIT_FAILURE);
    }

  // input filename
  run_parameters.input_filename = argv[arg++];
  mcutils::FileExistCheck(run_parameters.input_filename, true, false);

  // output filename
  run_parameters.output_filename = argv[arg++];
  
}


int main(int argc, const char **argv)
{

  ////////////////////////////////////////////////////////////////
  // initialization
  ////////////////////////////////////////////////////////////////

  // header
  std::cout << std::endl;
  std::cout << "me2j2h2  -- ME2J to H2 TBME file conversion" << std::endl;
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  // read parameters
  RunParameters run_parameters;
  ProcessArguments(argc, argv, run_parameters);

  // start timing
  mcutils::SteadyTimer total_time;
  total_time.Start();

  ////////////////////////////////////////////////////////////////
  // input
  ////////////////////////////////////////////////////////////////
  
  // define operator labels
  int J0 = 0;
  int g0 = 0;
  int Tz0 = 0;

  // initialize
  basis::TwoBodySpaceJJJTTz two_body_jjjttz_space(run_parameters.truncation_rank, run_parameters.truncation_cutoff);
  basis::TwoBodySectorsJJJTTz two_body_jjjttz_sectors(two_body_jjjttz_space, J0, g0, Tz0);
  basis::OperatorBlocks<double> two_body_jjjttz_matrices;
  SetOperatorToZero(two_body_jjjttz_sectors, two_body_jjjttz_matrices);

  // read me2j file
  std::cout << "Input stream" << std::endl;
  shell::ReadMe2jFile(
      two_body_jjjttz_space, two_body_jjjttz_sectors, two_body_jjjttz_matrices,
      run_parameters.input_filename, run_parameters.float_size
    );
  std::cout << std::endl;

  ////////////////////////////////////////////////////////////////
  // branching to two-body JJJPN (and output)
  ////////////////////////////////////////////////////////////////

  std::cout << "Branch to two-body jjJpn..." << std::endl
            << std::endl;

  // transform
  basis::TwoBodySpaceJJJT two_body_jjjt_space(run_parameters.truncation_rank, run_parameters.truncation_cutoff);
  std::array<basis::TwoBodySectorsJJJT,3> two_body_jjjt_component_sectors;
  std::array<basis::OperatorBlocks<double>,3> two_body_jjjt_component_matrices;
  shell::TransformOperatorTwoBodyJJJTTzToTwoBodyJJJT(two_body_jjjttz_space,two_body_jjjttz_sectors,two_body_jjjttz_matrices,
    two_body_jjjt_space,two_body_jjjt_component_sectors,two_body_jjjt_component_matrices
    );

  // define space and operator containers
  int N1max = run_parameters.truncation_cutoff;
  basis::OrbitalSpacePN orbital_space(N1max);
  basis::TwoBodySpaceJJJPNOrdering two_body_jjjpn_space_ordering =
    shell::kH2SpaceOrdering.at(0);
  basis::TwoBodySpaceJJJPN two_body_jjjpn_space(
      orbital_space,
      basis::WeightMax(run_parameters.truncation_rank, run_parameters.truncation_cutoff),
      two_body_jjjpn_space_ordering
    );
  basis::TwoBodySectorsJJJPN two_body_jjjpn_sectors(
      two_body_jjjpn_space,
      J0, g0, Tz0
    );

  // stream initialization
  std::cout << "Output stream" << std::endl;
  shell::OutH2Stream output_stream(
      run_parameters.output_filename,
      orbital_space, two_body_jjjpn_space, two_body_jjjpn_sectors,
      shell::kVersion0
    );
  std::cout << output_stream.DiagnosticStr();

  // do branching
  for (std::size_t sector_index=0; sector_index<two_body_jjjpn_sectors.size(); ++sector_index)
    {
      // make reference to target sector
      const basis::TwoBodySectorsJJJPN::SectorType& two_body_jjjpn_sector
        = two_body_jjjpn_sectors.GetSector(sector_index);

      // transform
      basis::OperatorLabelsJT operator_labels(J0,g0,0,2,basis::SymmetryPhaseMode::kHermitian);
      auto matrix = moshinsky::TwoBodyMatrixJJJPN(
          operator_labels,
          two_body_jjjt_space,
          two_body_jjjt_component_sectors,
          two_body_jjjt_component_matrices,
          two_body_jjjpn_sector
        );
      // mcutils::ChopMatrix(matrix);
      output_stream.WriteSector(
          sector_index,
          matrix,
          basis::NormalizationConversion::kASToNAS
        );
      std::cout << "." << std::flush;
    }
  std::cout << std::endl;
  output_stream.Close();


  ////////////////////////////////////////////////////////////////
  // termination
  ////////////////////////////////////////////////////////////////

  // end timing
  total_time.Stop();
  std::cout << "(Total time: " << total_time.ElapsedTime() << ")" << std::endl;
  std::cout << std::endl;

  // exit
  return EXIT_SUCCESS;
}
