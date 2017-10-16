/******************************************************************************
  @file radial-gen.cpp

  compute radial matrix elements

  Syntax:
    + radial-gen --kinematic operator_type order analytic_basis_type orbital_file output_filename
      - operator_type={r,k}
    + radial-gen --radial order j0 g0 analytic_basis_type orbital_file output_filename
    + radial-gen --xform scale_ratio analytic_basis_type bra_orbital_file [ket_orbital_file] output_filename
    + radial-gen --pn-overlaps orbital_file output_filename
    + radial-gen --identity bra_orbital_file [ket_orbital_file] output_filename

  @note Currently only computes radial matrix elements between harmonic oscillator
  or Laguerre basis functions with identical bra and ket spaces.

  Patrick J. Fasano
  University of Notre Dame

  + 11/2/16 (pjf): Created, based on h2conv.
  + 11/4/16 (pjf): Implement different modes of operation:
    - Kinematic mode calculates kinematic matrix elements between states of a
     single space.
    - Overlaps mode calculates radial overlaps between states with different
     length parameter.
  + 11/6/16 (mac):
    - Fix radial_operator_type in overlap mode.
    - Overhaul control logic for matrix.
    - Implement transformation to between different basis function types.
    - Flag possible fencepost error in spline integration point/steps.
  + 12/29/16 (mac): Add OMP diagnostic.
  + 1/23/16 (pjf): Add identity mode.
  + 1/24/16 (pjf):
    - Add non-square overlap mode.
    - Add file existence checks.
  + 1/27/16 (pjf): Add identity for non-square matrices
  + 09/20/17 (pjf): Add support for generating pn overlaps.
  + 10/12/17 (pjf): Update for changes to radial_io:
    - Add radial mode. (ugly hack, clean up later)

******************************************************************************/

#include <sys/stat.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "basis/nlj_orbital.h"
#include "cppformat/format.h"
#include "mcutils/profiling.h"
#include "radial/radial_io.h"
#include "spline/wavefunction_class.h"

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

enum class AnalyticBasisType : int {
  kOscillator = 0, kLaguerre = 1
};

enum class OperationMode {kKinematic, kRadial, kXform, kPNOverlaps, kIdentity};

// Stores simple parameters for run
struct RunParameters {
  // filenames
  std::string bra_orbital_filename;
  std::string ket_orbital_filename;
  std::string output_filename;
  // mode
  OperationMode mode;
  shell::RadialOperatorType radial_operator;
  int order;
  int j0;
  int g0;
  int Tz0;
  float scale_ratio;
  AnalyticBasisType basis_type;
};

void PrintUsage(char **argv) {
  std::cout << "Usage: " << argv[0]
            << " --kinematic {r|k} order analytic_basis_type orbital_file output_filename"
            << std::endl;
  std::cout << "Usage: " << argv[0]
            << " --radial order j0 g0 Tz0 analytic_basis_type orbital_file output_filename"
            << std::endl;
  std::cout << "       " << argv[0]
            << " --xform scale_ratio analytic_basis_type bra_orbital_file [ket_orbital_file] output_filename"
            << std::endl;
  std::cout << "       " << argv[0]
            << " --pn-overlaps orbital_file output_filename"
            << std::endl;
  std::cout << "       " << argv[0]
            << " --identity bra_orbital_file [ket_orbital_file] output_filename"
            << std::endl;
}

void ProcessArguments(int argc, char **argv, RunParameters& run_parameters) {
  // argument counter
  int arg = 1;

  // usage message (must at least provide mode)
  if (argc-1 < 3) {
    PrintUsage(argv);
    std::exit(EXIT_FAILURE);
  }

  // operation mode
  {
    std::istringstream parameter_stream(argv[arg++]);
    if (parameter_stream.str() == "--kinematic") {
      run_parameters.mode = OperationMode::kKinematic;
      run_parameters.Tz0 = 0;
    } else if (parameter_stream.str() == "--radial") {
      run_parameters.mode = OperationMode::kRadial;
      run_parameters.radial_operator = shell::RadialOperatorType::kR;
    } else if (parameter_stream.str() == "--xform") {
      run_parameters.mode = OperationMode::kXform;
      run_parameters.order = 0;
      run_parameters.Tz0 = 0;
      run_parameters.radial_operator = shell::RadialOperatorType::kO;
    } else if (parameter_stream.str() == "--pn-overlaps") {
      run_parameters.mode = OperationMode::kPNOverlaps;
      run_parameters.order = 0;
      run_parameters.Tz0 = 1;
      run_parameters.radial_operator = shell::RadialOperatorType::kO;
    } else if (parameter_stream.str() == "--identity") {
      run_parameters.mode = OperationMode::kIdentity;
      run_parameters.order = 0;
      run_parameters.Tz0 = 0;
      run_parameters.radial_operator = shell::RadialOperatorType::kO;
    } else {
      PrintUsage(argv);
      std::cerr << "Invalid operation mode." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // mode-specific options
  if (run_parameters.mode == OperationMode::kKinematic) {
    // operator code
    {
      std::istringstream parameter_stream(argv[arg++]);
      char operator_code;
      parameter_stream >> operator_code;
      if (!parameter_stream) {
        PrintUsage(argv);
        std::cerr << "Valid operators: r, k" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      run_parameters.radial_operator = static_cast<shell::RadialOperatorType>(operator_code);
      if (run_parameters.radial_operator != shell::RadialOperatorType::kR &&
          run_parameters.radial_operator != shell::RadialOperatorType::kK)
      {
        PrintUsage(argv);
        std::cerr << "Valid operators: r, k" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }

    // operator order
    {
      std::istringstream parameter_stream(argv[arg++]);
      parameter_stream >> run_parameters.order;
      if (!parameter_stream) {
        PrintUsage(argv);
        std::cerr << "Order must be an integer." << std::endl;
        std::exit(EXIT_FAILURE);
      }

      if (run_parameters.order == 0) {
        std::cout << "WARN: Order is zero. Do you mean xform?" << std::endl;
        run_parameters.order = 0;
      }
    }

    // reference scale is 1.0 for kinematic operators
    run_parameters.scale_ratio = 1.0;
  } else if (run_parameters.mode == OperationMode::kRadial) {
    // operator order
    {
      std::istringstream parameter_stream(argv[arg++]);
      parameter_stream >> run_parameters.order;
      if (!parameter_stream) {
        PrintUsage(argv);
        std::cerr << "Order must be an integer." << std::endl;
        std::exit(EXIT_FAILURE);
      }

      if (run_parameters.order == 0) {
        std::cout << "WARN: Order is zero. Do you mean xform?" << std::endl;
        run_parameters.order = 0;
      }
    }

    // j0
    {
      std::istringstream parameter_stream(argv[arg++]);
      parameter_stream >> run_parameters.j0;
      if (!parameter_stream) {
        PrintUsage(argv);
        std::cerr << "j0 must be an integer." << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }

    // g0
    {
      std::istringstream parameter_stream(argv[arg++]);
      parameter_stream >> run_parameters.g0;
      if (!parameter_stream) {
        PrintUsage(argv);
        std::cerr << "g0 must be an integer." << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }

    // Tz0
    {
      std::istringstream parameter_stream(argv[arg++]);
      parameter_stream >> run_parameters.Tz0;
      if (!parameter_stream) {
        PrintUsage(argv);
        std::cerr << "Tz0 must be an integer." << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }

    // reference scale is 1.0 for kinematic operators
    run_parameters.scale_ratio = 1.0;
  } else if (run_parameters.mode == OperationMode::kXform) {
    // order for xform is 0 (i.e. r^0)
    run_parameters.radial_operator = shell::RadialOperatorType::kO;
    run_parameters.order = 0;

    // scale ratio
    {
      std::istringstream parameter_stream(argv[arg++]);
      parameter_stream >> run_parameters.scale_ratio;
      if (!parameter_stream) {
        PrintUsage(argv);
        std::cerr << "Invalid scale ratio." << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
  }

  // basis type
  if (run_parameters.mode == OperationMode::kKinematic || run_parameters.mode == OperationMode::kRadial || run_parameters.mode == OperationMode::kXform)
  {
    std::istringstream parameter_stream(argv[arg++]);
    if (parameter_stream.str() == "oscillator") {
      run_parameters.basis_type = AnalyticBasisType::kOscillator;
      std::cout << "INFO: Using oscillator basis functions" << std::endl;
    } else if (parameter_stream.str() == "laguerre") {
      run_parameters.basis_type = AnalyticBasisType::kLaguerre;
      std::cout << "INFO: Using Laguerre basis functions" << std::endl;
    } else {
      PrintUsage(argv);
      std::cerr << "Valid analytic basis types: oscillator|laguerre" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // bra orbital file
  {
    std::istringstream parameter_stream(argv[arg++]);
    parameter_stream >> run_parameters.bra_orbital_filename;
    if (!parameter_stream) {
      PrintUsage(argv);
      std::cerr << "Invalid orbital filename." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    struct stat st;
    if (stat(run_parameters.bra_orbital_filename.c_str(), &st) != 0) {
      std::cerr << "ERROR: file " << run_parameters.bra_orbital_filename << " does not exist!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // ket orbital file -- only parse an argument if there were 6 arguments
  if ((run_parameters.mode == OperationMode::kXform && (argc-1) == 6) ||
      (run_parameters.mode == OperationMode::kIdentity && (argc-1) == 4))
  {
    std::istringstream parameter_stream(argv[arg++]);
    parameter_stream >> run_parameters.ket_orbital_filename;
    if (!parameter_stream) {
      PrintUsage(argv);
      std::cerr << "Invalid ket orbital filename." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    struct stat st;
    if (stat(run_parameters.ket_orbital_filename.c_str(), &st) != 0) {
      std::cerr << "ERROR: file " << run_parameters.ket_orbital_filename << " does not exist!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  } else {
    run_parameters.ket_orbital_filename = run_parameters.bra_orbital_filename;
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
    struct stat st;
    if (stat(run_parameters.output_filename.c_str(), &st) == 0) {
      std::cerr << "WARN: overwriting file " << run_parameters.output_filename << std::endl;
    }
  }
}

void CalculateMatrixElements(
    spline::Basis bra_basis_type, spline::Basis ket_basis_type,
    double bra_scale_ratio, double ket_scale_ratio,
    int order,
    const basis::OrbitalSectorsLJPN& sectors,
    basis::OperatorBlocks<double>& matrices
  ) {
  for (int sector_index=0; sector_index < sectors.size(); ++sector_index) {
    // get sector
    const basis::OrbitalSectorsLJPN::SectorType sector = sectors.GetSector(sector_index);
    Eigen::MatrixXd sector_matrix(sector.bra_subspace().size(), sector.ket_subspace().size());

    // get sizes
    const int bra_subspace_size = sector.bra_subspace().size();
    const int ket_subspace_size = sector.ket_subspace().size();

    // main loop
    #pragma omp parallel for collapse(2)
    for (int j=0; j < bra_subspace_size; ++j) {
      for (int k=0; k < ket_subspace_size; ++k) {
        // get bra state
        basis::OrbitalStateLJPN bra_state(sector.bra_subspace(), j);
        spline::WaveFunction bra_wavefunction(bra_state.n(), bra_state.l(), bra_scale_ratio, bra_basis_type);
        // get ket state
        basis::OrbitalStateLJPN ket_state(sector.ket_subspace(), k);
        spline::WaveFunction ket_wavefunction(ket_state.n(), ket_state.l(), ket_scale_ratio, ket_basis_type);

        const int num_size = 3000;
        sector_matrix(j, k) = bra_wavefunction.MatrixElement(num_size, ket_wavefunction, order);
      }
    }
    matrices.push_back(sector_matrix);
  }
}

int main(int argc, char **argv) {

  // header
  std::cout << std::endl;
  std::cout << "radial-gen -- radial integral evaluation" << std::endl;
  std::cout << std::endl;

  // process arguments
  RunParameters run_parameters;
  ProcessArguments(argc, argv, run_parameters);

  // parallel performance diagnostic
  std::cout
    << fmt::format(
        "INFO: OMP max_threads {}, num_procs {}",
        omp_get_max_threads(), omp_get_num_procs()
      )
    << std::endl
    << std::endl;

  // Read orbitals
  std::ifstream bra_orbital_stream(run_parameters.bra_orbital_filename);
  std::vector<basis::OrbitalPNInfo> bra_input_orbitals =
    basis::ParseOrbitalPNStream(bra_orbital_stream, true);
  std::ifstream ket_orbital_stream(run_parameters.ket_orbital_filename);
  std::vector<basis::OrbitalPNInfo> ket_input_orbitals =
    basis::ParseOrbitalPNStream(ket_orbital_stream, true);

  // Eigen initialization
  Eigen::initParallel();
  basis::OperatorBlocks<double> matrices;

  // main control logic
  spline::Basis bra_basis_type, ket_basis_type;
  double bra_scale_ratio, ket_scale_ratio;
  if (run_parameters.mode == OperationMode::kKinematic) {
    bra_scale_ratio = ket_scale_ratio = 1.;
    if (run_parameters.radial_operator == shell::RadialOperatorType::kR) {
      if (run_parameters.basis_type == AnalyticBasisType::kOscillator)
        bra_basis_type = ket_basis_type = spline::Basis::HC;
      else if (run_parameters.basis_type == AnalyticBasisType::kLaguerre)
        bra_basis_type = ket_basis_type = spline::Basis::LC;
    } else if (run_parameters.radial_operator == shell::RadialOperatorType::kK) {
      if (run_parameters.basis_type == AnalyticBasisType::kOscillator)
        bra_basis_type = ket_basis_type = spline::Basis::HM;
      else if (run_parameters.basis_type == AnalyticBasisType::kLaguerre)
        bra_basis_type = ket_basis_type = spline::Basis::LM;
    }
  } else if (run_parameters.mode == OperationMode::kRadial) {
      bra_scale_ratio = ket_scale_ratio = 1.;
      if (run_parameters.basis_type == AnalyticBasisType::kOscillator)
        bra_basis_type = ket_basis_type = spline::Basis::HC;
      else if (run_parameters.basis_type == AnalyticBasisType::kLaguerre)
        bra_basis_type = ket_basis_type = spline::Basis::LC;
  } else if (run_parameters.mode == OperationMode::kXform) {
    bra_scale_ratio = 1.;
    ket_scale_ratio = run_parameters.scale_ratio;
    bra_basis_type = spline::Basis::HC;
    if (run_parameters.basis_type == AnalyticBasisType::kOscillator)
      ket_basis_type = spline::Basis::HC;
    else if (run_parameters.basis_type == AnalyticBasisType::kLaguerre)
      ket_basis_type = spline::Basis::LC;
  }

  // Construct indexing
  basis::OrbitalSpaceLJPN bra_space(bra_input_orbitals);
  basis::OrbitalSpaceLJPN ket_space(ket_input_orbitals);
  basis::OrbitalSectorsLJPN sectors;
  if (run_parameters.mode == OperationMode::kRadial) {
    sectors = basis::OrbitalSectorsLJPN(bra_space, ket_space, run_parameters.j0, run_parameters.g0, run_parameters.Tz0);
    std::cout << sectors.DebugStr() << std::endl;
  } else {
     sectors = basis::OrbitalSectorsLJPN(bra_space, ket_space, run_parameters.order, run_parameters.Tz0);
  }

  if (run_parameters.mode == OperationMode::kKinematic || run_parameters.mode == OperationMode::kRadial || run_parameters.mode == OperationMode::kXform)
  {
    CalculateMatrixElements(
        bra_basis_type, ket_basis_type,
        bra_scale_ratio, ket_scale_ratio,
        run_parameters.order,
        sectors, matrices
      );
  }
  else if (run_parameters.mode == OperationMode::kIdentity || run_parameters.mode == OperationMode::kPNOverlaps)
  {
    for (int sector_index = 0; sector_index < sectors.size(); ++sector_index) {
      const basis::OrbitalSectorsLJPN::SectorType sector = sectors.GetSector(sector_index);
      matrices.push_back(Eigen::MatrixXd::Identity(sector.bra_subspace().size(), sector.ket_subspace().size()));
    }
  }

  // write out to file
  shell::OutRadialStream os(run_parameters.output_filename,
                            bra_space, ket_space, sectors,
                            run_parameters.radial_operator, run_parameters.order);
  os.Write(matrices);

  return EXIT_SUCCESS;
}
