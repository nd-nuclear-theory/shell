/**************************************************************************//**
  @file orbital-gen.cpp

  generate MFDn v15 orbital files

  Syntax:
    + orbital-gen --Nmax Nmax output_file
    + orbital-gen --triangular weight_max a b output_file
      - weight = a*n + b*l

  Patrick J. Fasano
  University of Notre Dame

  + 11/06/16 (pjf): Created, based on radial-gen.
  + 11/28/17 (pjf): Print header with version.
  + 07/27/18 (pjf): Rename "oscillator" truncation to "Nmax"
  + 10/10/19 (pjf): Add "convert" mode.

******************************************************************************/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>

#include "basis/nlj_orbital.h"
#include "mcutils/io.h"

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

enum class TruncationMode {kConvert, kNmax, kTriangular};

// Stores simple parameters for run
struct RunParameters {
  // filename
  std::string output_filename;
  std::string input_filename;
  // format
  basis::MFDnOrbitalFormat input_format;
  basis::MFDnOrbitalFormat output_format;
  // mode
  TruncationMode mode;
  double weight_max;
  double a;
  double b;

  /** default constructor */
  RunParameters()
    : output_filename(""), input_filename(""),
      input_format(basis::MFDnOrbitalFormat::kVersion15099),
      output_format(basis::MFDnOrbitalFormat::kVersion15099),
      mode(TruncationMode::kNmax), weight_max(0), a(0), b(0)
  {}
};

void PrintUsage(const char **argv) {
  std::cout << "Usage: " << argv[0]
            << " --convert input_format input_file output_format output_file"
            << std::endl;
  std::cout << "Usage: " << argv[0]
            << " --Nmax Nmax output_file"
            << std::endl;
  std::cout << "       " << argv[0]
            << " --triangular weight_max a b output_file"
            << std::endl;
}

void ProcessArguments(int argc, const char *argv[], RunParameters& run_parameters) {
  // argument counter
  int arg = 1;

  // usage message
  if (argc-1 < 3) {
    PrintUsage(argv);
    std::cerr << "Incorrect number of arguments." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // operation mode
  {
    std::istringstream parameter_stream(argv[arg++]);
    if (parameter_stream.str() == "--convert") {
      run_parameters.mode = TruncationMode::kConvert;
    } else if (parameter_stream.str() == "--Nmax") {
      run_parameters.mode = TruncationMode::kNmax;
    } else if (parameter_stream.str() == "--triangular") {
      run_parameters.mode = TruncationMode::kTriangular;
    } else {
      PrintUsage(argv);
      std::cerr << "Invalid operation mode." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // mode-specific options
  if (run_parameters.mode == TruncationMode::kConvert) {
    if (argc-1 != 5)
    {
      PrintUsage(argv);
      std::cerr << "Incorrect number of arguments." << std::endl;
      std::exit(EXIT_FAILURE);
    }

    // input format
    {
      std::istringstream parameter_stream(argv[arg++]);
      int format_code;
      parameter_stream >> format_code;
      if (!parameter_stream) {
        PrintUsage(argv);
        std::cerr << "Invalid input format" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      run_parameters.input_format = static_cast<basis::MFDnOrbitalFormat>(format_code);
    }

    // input filename
    {
      std::istringstream parameter_stream(argv[arg++]);
      parameter_stream >> run_parameters.input_filename;
      if (!parameter_stream) {
        PrintUsage(argv);
        std::cerr << "Invalid input filename" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      mcutils::FileExistCheck(run_parameters.input_filename, true, false);
    }

    // output format
    {
      std::istringstream parameter_stream(argv[arg++]);
      int format_code;
      parameter_stream >> format_code;
      if (!parameter_stream) {
        PrintUsage(argv);
        std::cerr << "Invalid output format" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      run_parameters.output_format = static_cast<basis::MFDnOrbitalFormat>(format_code);
    }

  }
  else if (run_parameters.mode == TruncationMode::kNmax)
  {
    // check number of arguments
    if (argc-1 != 3) {
      PrintUsage(argv);
      std::cerr << "Incorrect number of arguments." << std::endl;
      std::exit(EXIT_FAILURE);
    }

    // Nmax
    {
      std::istringstream parameter_stream(argv[arg++]);
      int Nmax;
      parameter_stream >> Nmax;
      if (!parameter_stream) {
        PrintUsage(argv);
        std::cerr << "Invalid Nmax." << std::endl;
        std::exit(EXIT_FAILURE);
      }
      run_parameters.weight_max = Nmax;
    }
  } else if (run_parameters.mode == TruncationMode::kTriangular) {
    if (argc-1 != 5) {
      PrintUsage(argv);
      std::cerr << "Incorrect number of arguments" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    // weight max
    {
      std::istringstream parameter_stream(argv[arg++]);
      parameter_stream >> run_parameters.weight_max;
      if (!parameter_stream) {
        PrintUsage(argv);
        std::cerr << "Invalid weight_max" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }

    // a
    {
      std::istringstream parameter_stream(argv[arg++]);
      parameter_stream >> run_parameters.a;
      if (!parameter_stream) {
        PrintUsage(argv);
        std::cerr << "Invalid a." << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }

    // b
    {
      std::istringstream parameter_stream(argv[arg++]);
      parameter_stream >> run_parameters.b;
      if (!parameter_stream) {
        PrintUsage(argv);
        std::cerr << "Invalid b." << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
  }

  // output file
  {
    std::istringstream parameter_stream(argv[arg++]);
    parameter_stream >> run_parameters.output_filename;
    if (!parameter_stream) {
      PrintUsage(argv);
      std::cerr << "Invalid output filename." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    mcutils::FileExistCheck(run_parameters.output_filename, false, true);
  }
}

int main(int argc, const char *argv[]) {
  // header
  std::cout << std::endl;
  std::cout << "orbital-gen -- orbital file generation" << std::endl;
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  RunParameters run_parameters;
  ProcessArguments(argc, argv, run_parameters);
  std::vector<basis::OrbitalPNInfo> orbitals;

  if (run_parameters.mode == TruncationMode::kConvert)
  {
    std::ifstream is(run_parameters.input_filename);
    orbitals = basis::ParseOrbitalPNStream(is, true, run_parameters.input_format);
    is.close();
  }
  else if (run_parameters.mode == TruncationMode::kNmax) {
    orbitals = basis::OrbitalSpacePN(static_cast<int>(run_parameters.weight_max)).OrbitalInfo();
  } else if (run_parameters.mode == TruncationMode::kTriangular) {
    double& a = run_parameters.a;
    double& b = run_parameters.b;
    double& weight_max = run_parameters.weight_max;
    for (int species=0; species <= 1; ++species) {
      basis::OrbitalSpeciesPN orbital_species = static_cast<basis::OrbitalSpeciesPN>(species);
      for (int n=0; a*n <= weight_max; ++n) {
        for (int l=0; a*n + b*l <= weight_max; ++l) {
          if (l > 0) {
            basis::OrbitalPNInfo orbital(orbital_species, n, l, l-HalfInt(1, 2), a*n + b*l);
            orbitals.push_back(orbital);
          }
          basis::OrbitalPNInfo orbital(orbital_species, n, l, l+HalfInt(1, 2), a*n + b*l);
          orbitals.push_back(orbital);
        }
      }
    }
  }
  if (orbitals.size() == 0) {
    PrintUsage(argv);
    std::cerr << "ERROR: No orbitals generated." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // sort orbitals before output
  std::sort(orbitals.begin(), orbitals.end(), basis::OrbitalSortCmpWeight);

  std::ofstream os(run_parameters.output_filename);
  os << basis::OrbitalDefinitionStr(orbitals, true, run_parameters.output_format);

  /* return code */
  return EXIT_SUCCESS;
}
