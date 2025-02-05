/******************************************************************************

  h22me2j.cpp -- H2 to ME2J TBME file conversion

  Restrictions:

    Matrix elements connecting T=0 and T=1 states are discarded in me2j format.    

  Syntax:

    h22me2j input_filename output_filename

  Mark A. Caprio
  University of Notre Dame

  + 10/25/24 (mac): Created, based on h2stat/xpn2h2/h22me2j.

******************************************************************************/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>

#include "basis/jjjpn_scheme.h"
#include "basis/jjjpn_operator.h"
#include "basis/jjjttz_scheme.h"
#include "basis/jjjttz_operator.h"
#include "fmt/format.h"
#include "mcutils/parsing.h"
#include "mcutils/profiling.h"

#include "tbme/h2_io.h"
#include "tbme/me2j_io.h"
#include "tbme/tbme_scheme_xform.h"

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

struct RunParameters
// Stores simple parameters for run
{
  // filenames
  std::string input_filename;
  std::string output_filename;

  // default constructor
  RunParameters()
    : input_filename(""), output_filename("")
  {}
  
};

void PrintUsage(const char **argv) {
  std::cout << "Usage: " << argv[0]
            << " input_filename output_filename"
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
      else
        {
          PrintUsage(argv);
          std::cerr << "Unrecognized option '" << parameter_stream.str() << "'" << std::endl;
          std::exit(EXIT_FAILURE);
        }
    }
  
  // process fixed arguments
  if (argc-arg < 2)
    {
      PrintUsage(argv);
      std::cerr << "Insufficient arguments" << std::endl;
      std::exit(EXIT_FAILURE);
    }

  // input filename
  run_parameters.input_filename = argv[arg++];
  mcutils::FileExistCheck(run_parameters.input_filename, true, false);

  // output filename
  run_parameters.output_filename = argv[arg++];
  
}


////////////////////////////////////////////////////////////////
// h2 input
/////////////////////////////////////////////////////////////////

void ReadH2File(
    const std::string& filename,
    basis::OrbitalSpacePN& orbital_space, basis::TwoBodySpaceJJJPN& two_body_space,
    basis::TwoBodySectorsJJJPN& two_body_sectors, basis::OperatorBlocks<double>& two_body_matrices,
    basis::NormalizationConversion conversion_mode
  )
// Read all data from h2 file.
//
// Arguments:
//   filename (std::string): filename
//   orbital_space (basis::OrbitalSpacePN, output): orbitals
//   two_body_space (basis::TwoBodySpaceJJJPN, output): two-body space
//   two_body_sectors (basis::TwoBodySectorsJJJPN, output): two-body sectors
//   two_body_matrices (basis::OperatorBlocks<double>, output): TBME matrices
//   conversion_mode (basis::NormalizationConversion, optional): selects AS/NAS conversion mode
{

  // initialize stream

  std::cout << "Input stream" << std::endl;
  shell::InH2Stream input_stream(filename);
  std::cout << input_stream.DiagnosticStr();
  std::cout << std::endl;

  // initialize data structures
  orbital_space = input_stream.orbital_space();
  two_body_space = input_stream.space();
  two_body_sectors = input_stream.sectors();
  const basis::TwoBodySpaceJJJPNOrdering space_ordering =
    shell::kH2SpaceOrdering.at(input_stream.h2_format());
  two_body_space = basis::TwoBodySpaceJJJPN(
      orbital_space,
      two_body_space.weight_max(),
      space_ordering
    );
  two_body_sectors = basis::TwoBodySectorsJJJPN(
      two_body_space, two_body_sectors.J0(), two_body_sectors.g0(), two_body_sectors.Tz0()
    );
  two_body_matrices.resize(two_body_sectors.size());
  
  // read sectors
  for (std::size_t sector_index = 0; sector_index < input_stream.num_sectors(); ++sector_index)
    {
      auto& matrix = two_body_matrices[sector_index];
      input_stream.ReadSector(sector_index, matrix, conversion_mode);
    }

  // close stream
  input_stream.Close();

}

////////////////////////////////////////////////////////////////
// main program
/////////////////////////////////////////////////////////////////

int main(int argc, const char **argv)
{

  ////////////////////////////////////////////////////////////////
  // initialization
  ////////////////////////////////////////////////////////////////

  // header
  std::cout << std::endl;
  std::cout << "h22me2j -- H2 to ME2J TBME file conversion" << std::endl;
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  // read parameters
  RunParameters run_parameters;
  ProcessArguments(argc, argv, run_parameters);

  // start timing
  mcutils::SteadyTimer total_time;
  total_time.Start();

  ////////////////////////////////////////////////////////////////
  // conversion
  ////////////////////////////////////////////////////////////////
  
  // read h2

  basis::OrbitalSpacePN orbital_space;
  basis::TwoBodySpaceJJJPN two_body_jjjpn_space;
  basis::TwoBodySectorsJJJPN two_body_jjjpn_sectors;
  basis::OperatorBlocks<double> two_body_jjjpn_matrices;
  ReadH2File(
      run_parameters.input_filename,
      orbital_space, two_body_jjjpn_space,
      two_body_jjjpn_sectors, two_body_jjjpn_matrices,
      basis::NormalizationConversion::kNASToAS
    );

  // extract source operator information
  
  // validate and extract orbital truncation
  if (!orbital_space.is_oscillator_like())
    {
      std::cerr << "ERROR: Input h2 file defines orbital set which is not oscillator-like" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  int Nmax_orb = int(orbital_space.weight_max());
  
  // validate and extract two-body space truncation
  const basis::WeightMax& weight_max = two_body_jjjpn_space.weight_max();
  if (!(
      (weight_max.one_body[0] == int(weight_max.one_body[0]))
      && (weight_max.one_body[0] == weight_max.one_body[1])
        ))
    {
      std::cerr << "ERROR: Input h2 file one-body truncation not pn-symmetric and oscillator-like" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  int N1max = int(weight_max.one_body[0]);
  if (N1max > Nmax_orb)
    {
      std::cerr << "ERROR: Input h2 file one-body truncation not consistent with orbital set" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  if (!(
          (weight_max.two_body[0] == int(weight_max.two_body[0]))
          && (weight_max.two_body[0] == weight_max.two_body[1])
          && (weight_max.two_body[0] == weight_max.two_body[2])
        ))
    {
      std::cerr << "ERROR: Input h2 file two-body truncation not pn-symmetric and oscillator-like" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  int N2max = int(weight_max.two_body[0]);
  basis::Rank truncation_rank;
  if (N2max==2*N1max)
    truncation_rank = basis::Rank::kOneBody;
  else if (N2max==N1max)
    truncation_rank = basis::Rank::kTwoBody;
  else
    {
      std::cerr << "ERROR: Unsupported combination of N1max and N2max" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  int truncation_cutoff = N1max;
  // std::cout << fmt::format("Input truncation: rank {:d} cutoff {:d}", int(truncation_rank), truncation_cutoff) << std::endl;
  // std::cout << std::endl;
  
  // extract and validate operator quantum numbers
  int J0 = two_body_jjjpn_sectors.J0();
  int g0 = two_body_jjjpn_sectors.g0();
  int Tz0 = two_body_jjjpn_sectors.Tz0();
  if (J0!=0 || g0!=0 || Tz0!=0)
    {
      std::cerr << "ERROR: Input h2 file must specify a Hamiltonian-like (J0,g0,Tz0)=(0,0,0) operator" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  
  // convert scheme

  std::cout << "Upcouple to two-body jjJTTz..." << std::endl
            << std::endl;
  
  basis::TwoBodySpaceJJJTTz two_body_jjjttz_space(truncation_rank, truncation_cutoff);
  basis::TwoBodySectorsJJJTTz two_body_jjjttz_sectors;
  basis::OperatorBlocks<double> two_body_jjjttz_matrices;
  shell::TransformOperatorTwoBodyJJJPNToTwoBodyJJJTTz(
      two_body_jjjpn_space,
      two_body_jjjpn_sectors,
      two_body_jjjpn_matrices,
      two_body_jjjttz_space,
      two_body_jjjttz_sectors,
      two_body_jjjttz_matrices
    );
  
  // write output
  std::cout << "Output stream" << std::endl;
  shell::WriteMe2jFile(
      two_body_jjjttz_space,
      two_body_jjjttz_sectors,
      two_body_jjjttz_matrices,
      run_parameters.output_filename
    );
  std::cout << std::endl;
  
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
