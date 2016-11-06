/**************************************************************************//**
  @file radial-gen.cpp

  compute radial matrix elements

  Syntax:
    + radial-gen --kinematic operator_type order analytic_basis_type orbital_file output_filename
      - operator_type={r,k}
    + radial-gen --overlaps scale_ratio analytic_basis_type orbital_file output_filename

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

******************************************************************************/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "mcutils/profiling.h"
#include "basis/nlj_orbital.h"
#include "radial/radial_io.h"
#include "spline/wavefunction_class.h"

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

enum class AnalyticBasisType : int {
  kOscillator = 0, kLaguerre = 1
};

enum class OperationMode {kKinematic, kOverlaps};

// Stores simple parameters for run
struct RunParameters {
  // filenames
  std::string orbital_filename;
  std::string output_filename;
  // mode
  OperationMode mode;
  shell::RadialOperatorType radial_operator;
  int order;
  float scale_ratio;
  AnalyticBasisType basis_type;
};

void PrintUsage(char **argv) {
  std::cout << "Usage: " << argv[0]
            << " --kinematic {r|k} order analytic_basis_type orbital_file output_filename"
            << std::endl;
  std::cout << "       " << argv[0]
            << " --overlaps scale_ratio analytic_basis_type orbital_file output_filename"
            << std::endl;
}

void ProcessArguments(int argc, char **argv, RunParameters& run_parameters) {
  // argument counter
  int arg = 1;

  // usage message
  if (argc-1 < 5) {
    PrintUsage(argv);
    std::exit(EXIT_FAILURE);
  }

  // operation mode
  {
    std::istringstream parameter_stream(argv[arg++]);
    if (parameter_stream.str() == "--kinematic") {
      run_parameters.mode = OperationMode::kKinematic;
    } else if (parameter_stream.str() == "--overlaps"){
      run_parameters.mode = OperationMode::kOverlaps;
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
        std::cout << "WARN: Order is zero. Do you mean overlaps?" << std::endl;
        run_parameters.order = 0;
      }
    }

    // reference scale is 1.0 for kinematic operators
    run_parameters.scale_ratio = 1.0;
  } else if (run_parameters.mode == OperationMode::kOverlaps) {
    // order for overlaps is 0 (i.e. r^0)
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

  // orbital file
  {
    std::istringstream parameter_stream(argv[arg++]);
    parameter_stream >> run_parameters.orbital_filename;
    if (!parameter_stream) {
      PrintUsage(argv);
      std::cerr << "Invalid orbital filename." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    /** @todo check for file existence */
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

void CalculateMatrixElements(
    spline::Basis bra_basis_type, spline::Basis ket_basis_type,
    double bra_scale_ratio, double ket_scale_ratio,
    int order,
    const basis::OrbitalSectorsLJPN& sectors,
    basis::MatrixVector& matrices
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

        const int num_integration_steps_or_maybe_points_or_whatever = 3000;
        sector_matrix(j, k) = bra_wavefunction.MatrixElement(num_integration_steps_or_maybe_points_or_whatever, ket_wavefunction, order);
      }
    }
    matrices.push_back(sector_matrix);
  }
}

int main(int argc, char **argv) {
  RunParameters run_parameters;
  ProcessArguments(argc, argv, run_parameters);

  // Read orbitals
  std::ifstream is(run_parameters.orbital_filename);
  std::vector<basis::OrbitalPNInfo> input_orbitals =
    basis::ParseOrbitalPNStream(is, true);

  // Construct indexing
  basis::OrbitalSpaceLJPN space(input_orbitals);
  basis::OrbitalSectorsLJPN sectors(space, space, run_parameters.order, 0);
  /** @note currently has hard-coded Tz0=0 */

  // Eigen initialization
  Eigen::initParallel();
  basis::MatrixVector matrices;

  // main control logic
  spline::Basis bra_basis_type, ket_basis_type;
  double bra_scale_ratio, ket_scale_ratio;
  int order;
  if (run_parameters.mode == OperationMode::kKinematic)
    {
      bra_scale_ratio = ket_scale_ratio = 1.;
      order = run_parameters.order;
      if (run_parameters.radial_operator == shell::RadialOperatorType::kR)
        {
          if (run_parameters.basis_type == AnalyticBasisType::kOscillator)
            bra_basis_type = ket_basis_type = spline::Basis::HC;
          else if (run_parameters.basis_type == AnalyticBasisType::kLaguerre)
            bra_basis_type = ket_basis_type = spline::Basis::LC;
        }
      else if (run_parameters.radial_operator == shell::RadialOperatorType::kK)
        {
          if (run_parameters.basis_type == AnalyticBasisType::kOscillator)
            bra_basis_type = ket_basis_type = spline::Basis::HM;
          else if (run_parameters.basis_type == AnalyticBasisType::kLaguerre)
            bra_basis_type = ket_basis_type = spline::Basis::LM;
        }
    }
  else if (run_parameters.mode == OperationMode::kOverlaps)
    {
      bra_scale_ratio = 1.;
      ket_scale_ratio = run_parameters.scale_ratio;
      bra_basis_type = spline::Basis::HC;
      order = 0;
      if (run_parameters.basis_type == AnalyticBasisType::kOscillator)
        ket_basis_type = spline::Basis::HC;
      else if (run_parameters.basis_type == AnalyticBasisType::kLaguerre)
        ket_basis_type = spline::Basis::LC;
    }

  CalculateMatrixElements(
      bra_basis_type, ket_basis_type,
      bra_scale_ratio, ket_scale_ratio,
      order,
      sectors, matrices
    );

  // write out to file
  shell::OutRadialStream os(run_parameters.output_filename,
                            space, space, sectors,
                            run_parameters.radial_operator);
  os.Write(matrices);

  return EXIT_SUCCESS;
}
