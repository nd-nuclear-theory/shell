/******************************************************************************
  @file radial-gen.cpp

  compute radial matrix elements

  Syntax:
    + radial-gen

  Input format:
    set-ket-indexing orbital_filename
    set-analytic-basis-type basis_type
      basis_type = oscillator|laguerre
    define-operator-target operator_type order j0 output_filename
      operator_type = r|k
      (order, j0) = (1,1)|(2,0)|(2,2)
    define-radial-target operator_type order j0 g0 Tz0 output_filename
      operator_type = r|k
    define-xform-target scale_factor Tz0 bra_orbital_file output_filename
  @note Currently only computes operator/radial matrix elements between harmonic oscillator
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
  + 11/28/17 (pjf): Include version in header.
  + 02/12/18 (pjf):
    - Rewrite to take input on stdin.

******************************************************************************/

#include <sys/stat.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <set>

#include "am/am.h"
#include "am/wigner_gsl.h"
#include "basis/nlj_orbital.h"
#include "basis/nlj_operator.h"
#include "cppformat/format.h"
#include "mcutils/profiling.h"
#include "obme/obme_io.h"
#include "spline/wavefunction_class.h"

const double kPi = 3.1415926535897;

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

enum class AnalyticBasisType : int {
  kOscillator = 0, kLaguerre = 1
};

enum class OperationMode {kOperator, kRadial, kXform};

std::set<std::pair<int, int>> AllowedOperators = {{1,1}, {2,0}, {2,2}};

// Stores simple parameters for run
struct OperatorParameters {
  // filenames
  std::string bra_orbital_filename;
  std::string ket_orbital_filename;
  std::string output_filename;

  // basis parameters
  float scale_factor;
  AnalyticBasisType basis_type;

  // mode
  OperationMode mode;

  // operator parameters
  shell::RadialOperatorType radial_operator;
  int order;
  int j0;
  int g0;
  int Tz0;
};

void PrintUsage(char **argv) {
  std::cout << "Usage: " << argv[0] << std::endl;
}

void ReadParameters(std::vector<OperatorParameters>& operator_parameters) {
  std::string line;
  int line_count = 0;
  /*****************************************************************/
  /* default parameters */
  std::string bra_orbital_filename = "orbitals.dat";
  std::string ket_orbital_filename = "orbitals.dat";
  AnalyticBasisType basis_type = AnalyticBasisType::kOscillator;
  /*****************************************************************/

  while (std::getline(std::cin, line)) {
    ++line_count;

    // set up for line parsing
    std::istringstream line_stream(line);
    std::string keyword;
    line_stream >> keyword;

    // skip blank line or hash comment line
    if ((keyword == "") || (keyword == "#")) continue;

    // select action based on keyword
    if (keyword == "set-ket-indexing") {
      line_stream >> ket_orbital_filename;
      ParsingCheck(line_stream, line_count, line);
      struct stat st;
      if (stat(ket_orbital_filename.c_str(), &st) != 0) {
        std::cerr << "ERROR: file " << ket_orbital_filename
                  << " does not exist!" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    } else if (keyword == "set-analytic-basis-type") {
      std::string basis_type_str;
      line_stream >> basis_type_str;
      ParsingCheck(line_stream, line_count, line);
      if (basis_type_str == "oscillator") {
        basis_type = AnalyticBasisType::kOscillator;
        std::cout << "INFO: Using oscillator basis functions" << std::endl;
      } else if (basis_type_str == "laguerre") {
        basis_type = AnalyticBasisType::kLaguerre;
        std::cout << "INFO: Using Laguerre basis functions" << std::endl;
      } else {
        std::cerr << "Valid analytic basis types: oscillator|laguerre" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    } else if (keyword == "define-operator-target") {
      // define-operator-target operator_type order j0 g0 Tz0 output_filename
      //   operator_type = r|k
      //   NOTE: order must equal j0 except for order=2, j0=0 (dot product)
      char operator_type;
      std::string output_filename;
      int order;
      int j0, g0, Tz0;

      line_stream >> operator_type >> order >> j0 >> output_filename;
      ParsingCheck(line_stream, line_count, line);

      if (!AllowedOperators.count({order,j0})) {
        std::cerr << "ERROR: Invalid operator." << std::endl;
        std::exit(EXIT_FAILURE);
      }

      operator_parameters.push_back({ket_orbital_filename, ket_orbital_filename, output_filename,
        1.0, basis_type, OperationMode::kOperator,
        static_cast<shell::RadialOperatorType>(operator_type), order, j0, j0%2, 0});
    } else if (keyword == "define-radial-target") {
      // define-radial-target operator_type order j0 g0 Tz0 output_filename
      //   operator_type = r|k
      char operator_type;
      std::string output_filename;
      int order;
      int j0, g0, Tz0;

      line_stream >> operator_type >> order >> j0 >> g0 >> Tz0 >> output_filename;
      ParsingCheck(line_stream, line_count, line);

      operator_parameters.push_back({ket_orbital_filename, ket_orbital_filename, output_filename,
        1.0, basis_type, OperationMode::kRadial,
        static_cast<shell::RadialOperatorType>(operator_type), order, j0, g0, Tz0});
    } else if (keyword == "define-xform-target") {
      // define-xform-target scale_factor Tz0 bra_orbital_file output_filename
      std::string bra_orbital_file, output_filename;
      float scale_factor;
      int Tz0;

      line_stream >> scale_factor >> Tz0 >> bra_orbital_file >> output_filename;
      ParsingCheck(line_stream, line_count, line);

      operator_parameters.push_back({bra_orbital_filename, ket_orbital_filename, output_filename,
        scale_factor, basis_type, OperationMode::kXform,
        shell::RadialOperatorType::kO, 0, 0, 0, Tz0});
    }
  }
}

void CalculateRadialMatrixElements(
    spline::Basis bra_basis_type, spline::Basis ket_basis_type,
    double bra_scale_factor, double ket_scale_factor,
    int order,
    const basis::OrbitalSectorsLJPN& sectors,
    basis::OperatorBlocks<double>& matrices
  )
{
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
        // get states
        basis::OrbitalStateLJPN bra_state(sector.bra_subspace(), j);
        basis::OrbitalStateLJPN ket_state(sector.ket_subspace(), k);
        if ((order == 0) && (bra_basis_type == ket_basis_type) && (bra_scale_factor == ket_scale_factor)) {
          // short circuit on trivial identity/permutation
          if ((bra_state.n() == ket_state.n()) && (bra_state.l() == ket_state.l())) {
            sector_matrix(j, k) = 1.;
          } else {
            sector_matrix(j, k) = 0.;
          }
          continue;
        }

        // get wave functions
        spline::WaveFunction bra_wavefunction(bra_state.n(), bra_state.l(), bra_scale_factor, bra_basis_type);
        spline::WaveFunction ket_wavefunction(ket_state.n(), ket_state.l(), ket_scale_factor, ket_basis_type);

        const int num_size = 3000;
        sector_matrix(j, k) = bra_wavefunction.MatrixElement(num_size, ket_wavefunction, order);
      }
    }
    matrices.push_back(sector_matrix);
  }
}

void BuildOperator(OperatorParameters operator_parameters) {
  // Read orbitals
  std::ifstream bra_orbital_stream(operator_parameters.bra_orbital_filename);
  auto bra_input_orbitals = basis::ParseOrbitalPNStream(bra_orbital_stream, true);
  std::ifstream ket_orbital_stream(operator_parameters.ket_orbital_filename);
  auto ket_input_orbitals = basis::ParseOrbitalPNStream(ket_orbital_stream, true);
  // Construct indexing
  basis::OrbitalSpaceLJPN bra_space(bra_input_orbitals);
  basis::OrbitalSpaceLJPN ket_space(ket_input_orbitals);
  basis::OrbitalSectorsLJPN sectors = basis::OrbitalSectorsLJPN(
    bra_space, ket_space,
    operator_parameters.j0, operator_parameters.g0, operator_parameters.Tz0
  );
  basis::OperatorBlocks<double> matrices;

  // main control logic
  spline::Basis bra_basis_type, ket_basis_type;
  double bra_scale_factor, ket_scale_factor;
  int order = operator_parameters.order;
  // map (operator, basis) to spline::Basis enum
  std::map<std::pair<shell::RadialOperatorType, AnalyticBasisType>, spline::Basis> operator_map = {
    {{shell::RadialOperatorType::kR, AnalyticBasisType::kOscillator}, spline::Basis::HC},
    {{shell::RadialOperatorType::kR, AnalyticBasisType::kLaguerre},   spline::Basis::LC},
    {{shell::RadialOperatorType::kK, AnalyticBasisType::kOscillator}, spline::Basis::HM},
    {{shell::RadialOperatorType::kK, AnalyticBasisType::kLaguerre},   spline::Basis::HM},
  };
  if (operator_parameters.mode == OperationMode::kOperator || operator_parameters.mode == OperationMode::kRadial) {
    auto it = operator_map.find({operator_parameters.radial_operator, operator_parameters.basis_type});
    if (it == operator_map.end()) {
      std::cerr << "Invalid operator type/basis combination!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    bra_basis_type = ket_basis_type = it->second;
    bra_scale_factor = ket_scale_factor = 1.;
  } else if (operator_parameters.mode == OperationMode::kXform) {
    bra_scale_factor = 1.;
    ket_scale_factor = operator_parameters.scale_factor;
    bra_basis_type = spline::Basis::HC;
    auto it = operator_map.find({shell::RadialOperatorType::kR, operator_parameters.basis_type});
    if (it == operator_map.end()) {
      std::cerr << "Invalid basis!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    ket_basis_type = it->second;
  }

  std::cout << bra_scale_factor << ket_scale_factor << order << std::endl;
  CalculateRadialMatrixElements(
      bra_basis_type, ket_basis_type,
      bra_scale_factor, ket_scale_factor,
      order,
      sectors,
      matrices
    );

  // angular factors
  basis::OneBodyOperatorType operator_type;
  if (operator_parameters.mode == OperationMode::kOperator)
  {
    operator_type = basis::OneBodyOperatorType::kSpherical;
    for (int sector_index=0; sector_index < sectors.size(); ++sector_index)
    {
      const basis::OrbitalSectorsLJPN::SectorType sector = sectors.GetSector(sector_index);
      HalfInt bra_j = sector.bra_subspace().j();
      HalfInt ket_j = sector.ket_subspace().j();
      if (operator_parameters.order == 1 && operator_parameters.j0 == 1)
      {
        // see csbasis (58) and (60)
        matrices[sector_index] *=
          ParitySign(ket_j - HalfInt(1,2) + 1) * Hat(bra_j) * Hat(ket_j)
          * am::Wigner3J(bra_j, ket_j, 1, HalfInt(1,2), HalfInt(-1,2), 0);
      }
      else if (operator_parameters.order == 2 && operator_parameters.j0 == 0)
      {
        // dimension factor -- Edmonds convention
        matrices[sector_index] *= Hat(ket_j);
      }
      else if (operator_parameters.order == 2 && operator_parameters.j0 == 2)
      {
        // see Varshalovich 3.2.1 eq. 23 and Suhonen (2.57)
        matrices[sector_index] *= std::sqrt(2/3)
          * ParitySign(ket_j - HalfInt(1, 2) + 2) * Hat(bra_j) * Hat(ket_j)
          * am::Wigner3J(bra_j, ket_j, 2, HalfInt(1, 2), HalfInt(-1, 2), 0);
      }
    }
  }
  else
  {
    operator_type = basis::OneBodyOperatorType::kRadial;
  }

  // write out to file
  shell::OutOBMEStream os(operator_parameters.output_filename,
                            bra_space, ket_space, sectors,
                            operator_type,
                            operator_parameters.radial_operator, operator_parameters.order);
  os.Write(matrices);
}

int main(int argc, char **argv) {
  // header
  std::cout << std::endl;
  std::cout << "radial-gen -- radial integral evaluation" << std::endl;
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  // process arguments
  std::vector<OperatorParameters> operators;
  ReadParameters(operators);

  // parallel performance diagnostic
  std::cout
    << fmt::format(
        "INFO: OMP max_threads {}, num_procs {}",
        omp_get_max_threads(), omp_get_num_procs()
      )
    << std::endl
    << std::endl;

  // Eigen initialization
  Eigen::initParallel();

  for (auto& operator_parameters : operators) {
    BuildOperator(operator_parameters);
  }

  return EXIT_SUCCESS;
}
