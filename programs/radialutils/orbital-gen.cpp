/**************************************************************************//**
  @file orbital-gen.cpp

  generate MFDn v15 orbital files

  Syntax:
    + orbital-gen --oscillator Nmax output_file
    + orbital-gen --triangular weight_max a b output_file
      - weight = a*n + b*l

  Patrick J. Fasano
  University of Notre Dame

  + 11/6/16 (pjf): Created, based on radial-gen.

******************************************************************************/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>

#include "basis/nlj_orbital.h"

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

enum class AnalyticBasisType : int {
  kOscillator = 0, kLaguerre = 1
};

enum class TruncationMode {kOscillator, kTriangular};

// Stores simple parameters for run
struct RunParameters {
  // filename
  std::string output_filename;
  // mode
  TruncationMode mode;
  double weight_max;
  double a;
  double b;

  /** default constructor */
  RunParameters()
    : output_filename(""), mode(TruncationMode::kOscillator), weight_max(0), a(0), b(0) {}
};

void PrintUsage(const char **argv) {
  std::cout << "Usage: " << argv[0]
            << " --oscillator Nmax output_file"
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
    if (parameter_stream.str() == "--oscillator") {
      run_parameters.mode = TruncationMode::kOscillator;
    } else if (parameter_stream.str() == "--triangular") {
      run_parameters.mode = TruncationMode::kTriangular;
    } else {
      PrintUsage(argv);
      std::cerr << "Invalid operation mode." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // mode-specific options
  if (run_parameters.mode == TruncationMode::kOscillator) {
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
    /** @todo check and warn if overwriting existing file */
  }
}

int main(int argc, const char *argv[]) {
  RunParameters run_parameters;
  ProcessArguments(argc, argv, run_parameters);
  std::vector<basis::OrbitalPNInfo> orbitals;

  if (run_parameters.mode == TruncationMode::kOscillator) {
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
  os << basis::OrbitalDefinitionStr(orbitals, true);

  /* return code */
  return EXIT_SUCCESS;
}
