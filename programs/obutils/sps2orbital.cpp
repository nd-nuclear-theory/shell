/******************************************************************************

  sps2orbital.cpp -- convert BIGSTICK sps to shell orbital file

  Syntax:
    + sps2orbital input_filename output_filename

  Limitations:
    + Input sps file support is limited to "iso"-type sps files.
    + Output orbital file format is fixed as Version 15099.

  See Sec 4.2 "Defining the model space" in C. W. Johnson et al., "BIGSTICK: A
  flexible configuration-interaction shell-model code", arXiv:1801.08432.

  Mark A. Caprio
  University of Notre Dame

  + 03/01/23 (mac): Created, based on orbital-gen.
  + 08/20/24 (mac): Add option --clear-weights.

******************************************************************************/

#include <fstream>

#include "basis/nlj_orbital.h"
#include "fmt/format.h"
#include "mcutils/parsing.h"

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

struct RunParameters
// Stores simple parameters for run
{
  // filenames
  std::string input_filename;
  std::string output_filename;
  // filtering
  bool clear_weights;
  // output format
  basis::MFDnOrbitalFormat output_format;

  // default constructor
  RunParameters()
    : output_filename(""), input_filename(""),
      clear_weights(false),
      output_format(basis::MFDnOrbitalFormat::kVersion15099)
  {}
 
};

void PrintUsage(const char **argv) {
  std::cout << "Usage: " << argv[0]
            << " [--clear-weights] input_file output_file"
            << std::endl;
}

void ProcessArguments(int argc, const char *argv[], RunParameters& run_parameters)
{
  // usage message
  if (argc-1 == 0)
    {
      PrintUsage(argv);
      std::exit(EXIT_SUCCESS);
    }

  int arg = 1;
  while (arg < argc && argv[arg][0] == '-')
    {
      std::istringstream parameter_stream(argv[arg++]);

      if (parameter_stream.str() == "--help" || parameter_stream.str() == "-h")
        {
          PrintUsage(argv);
          std::exit(EXIT_SUCCESS);
        }
      else if (parameter_stream.str() == "--clear-weights")
        {
          run_parameters.clear_weights = true;
        }
      else
        {
          PrintUsage(argv);
          std::cerr << "Unrecognized option '" << parameter_stream.str() << "'" << std::endl;
          std::exit(EXIT_FAILURE);
        }
    }

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

basis::OrbitalPNList ReadSPS(const std::string& filename, bool clear_weights)
// Read BIGSTICK sps file for single-particle orbitals.
{

  // open sps file
  std::ifstream is(filename);
  if (!is)
    {
      std::cout << "ERROR: Failure opening sps file" << std::endl;
      std::exit(EXIT_SUCCESS);
    }

  // set up line counter for use in error messages
  std::string line;
  int line_count = 0;

  // parse header
  // line 1: mode
  // Only "iso" mode supported.
  {
    mcutils::GetLine(is,line,line_count);
    std::istringstream line_stream(line);
    std::string mode;
    line_stream >> mode;
    // std::cout << fmt::format("mode {}", mode) << std::endl;
    mcutils::ParsingCheck(line_stream,line_count,line);
    assert(mode=="iso");
  }

  // line 2: number of orbitals (same for p and n)
  int num_orbitals;
  {
    mcutils::GetLine(is,line,line_count);
    std::istringstream line_stream(line);
    std::string mode;
    line_stream >> num_orbitals;
    // std::cout << fmt::format("num_orbitals {}", num_orbitals) << std::endl;
    mcutils::ParsingCheck(line_stream,line_count,line);
    assert(num_orbitals>0);
  }

  // read in orbital info
  basis::OrbitalPNList orbitals;
  for (int orbital_index=0; orbital_index<num_orbitals; ++orbital_index)
    {
      mcutils::GetLine(is,line,line_count);
      std::istringstream line_stream(line);
      std::string mode;
      float n_raw, l_raw, j_raw, w_raw;
      line_stream >> n_raw >> l_raw >> j_raw >> w_raw;
      mcutils::ParsingCheck(line_stream,line_count,line);
      // std::cout << fmt::format("Raw input data: {} {} {} {}", n_raw, l_raw, j_raw, w_raw) << std::endl;

      // convert orbital parameters
      int n = int(n_raw);
      int l = int(l_raw);
      HalfInt j = HalfInt(2*j_raw, 2);
      float weight = w_raw;
      if (clear_weights)
        weight = 0.;

      // store as proton orbital
      basis::OrbitalPNInfo proton_orbital(basis::OrbitalSpeciesPN::kP, n, l, j, weight);
      orbitals.push_back(proton_orbital);
    }
  
  // append identical list of neutron orbitals
  for (int orbital_index=0; orbital_index<num_orbitals; ++orbital_index)
    {
      const basis::OrbitalPNInfo& proton_orbital = orbitals[orbital_index];
      basis::OrbitalPNInfo neutron_orbital(basis::OrbitalSpeciesPN::kN, proton_orbital.n, proton_orbital.l, proton_orbital.j, proton_orbital.weight);
      orbitals.push_back(neutron_orbital);
    }

  return orbitals;
}


int main(int argc, const char *argv[])
{
  // header
  std::cout << std::endl;
  std::cout << "sps2orbital  -- sps to orbital conversion" << std::endl;
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  // read parameters
  RunParameters run_parameters;
  ProcessArguments(argc,argv,run_parameters);

  // read orbitals from sps file
  std::cout << fmt::format("Reading orbital file {}...", run_parameters.input_filename) << std::endl;
  basis::OrbitalPNList orbitals = ReadSPS(run_parameters.input_filename, run_parameters.clear_weights);
  
  // write orbitals
  std::cout << fmt::format("Writing orbital file {}...", run_parameters.output_filename) << std::endl;
  std::ofstream os(run_parameters.output_filename);
  os << basis::OrbitalDefinitionStr(
      orbitals, true, run_parameters.output_format
    );
  os.close();

  return EXIT_SUCCESS;
}
