/******************************************************************************
  @file radial-gen.cpp

  compute radial matrix elements

  Syntax:
    + radial-gen

  Input format:
    set-ket-basis basis_type orbital_filename
      basis_type = oscillator|laguerre
    define-operator-target operator_type order j0 output_filename
      operator_type = r|k
      (order, j0) = (1,1)|(2,0)|(2,2)
    define-radial-target operator_type order j0 g0 Tz0 output_filename
      operator_type = r|k|o
    define-xform-target scale_factor bra_basis_type bra_orbital_file output_filename
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
  + 10/18/18 (pjf): Update to use new obme library.

******************************************************************************/

#include <sys/stat.h>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "am/am.h"
#include "am/wigner_gsl.h"
#include "basis/nlj_operator.h"
#include "basis/nlj_orbital.h"
#include "fmt/format.h"
#include "mcutils/parsing.h"
#include "mcutils/profiling.h"
#include "obme/obme_io.h"
#include "obme/radial.h"
#include "obme/obme.h"

const double kPi = 3.1415926535897;

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

enum class OperationMode { kOperator, kRadial, kXform };

std::map<std::string, shell::RadialBasisType> kBasisTypeDefinitions(
    {{"oscillator", shell::RadialBasisType::kOscillator},
     {"laguerre", shell::RadialBasisType::kLaguerre}});

// Stores simple parameters for run
struct OperatorParameters {
  // filenames
  std::string bra_orbital_filename;
  std::string ket_orbital_filename;
  std::string output_filename;

  // basis parameters
  shell::RadialBasisType bra_basis_type;
  shell::RadialBasisType ket_basis_type;
  double scale_factor;

  // mode
  OperationMode mode;

  // operator parameters
  shell::RadialOperatorType radial_operator;
  int order;
  int j0;
  int g0;
  int Tz0;
};

void PrintUsage(char** argv) { std::cout << "Usage: " << argv[0] << std::endl; }

void ReadParameters(std::vector<OperatorParameters>& operator_parameters)
{
  std::string line;
  int line_count = 0;
  /*****************************************************************/
  /* default parameters */
  std::string bra_orbital_filename = "orbitals.dat";
  std::string ket_orbital_filename = "orbitals.dat";
  shell::RadialBasisType ket_basis_type = shell::RadialBasisType::kOscillator;
  /*****************************************************************/

  while (mcutils::GetLine(std::cin, line, line_count)) {
    // set up for line parsing
    std::istringstream line_stream(line);
    std::string keyword;
    line_stream >> keyword;

    // select action based on keyword
    if (keyword == "set-ket-basis") {
      // set-ket-basis basis_type orbital_filename
      //   basis_type = oscillator|laguerre
      std::string basis_type_str;
      line_stream >> basis_type_str >> ket_orbital_filename;
      mcutils::ParsingCheck(line_stream, line_count, line);
      // convert basis
      if (kBasisTypeDefinitions.count(basis_type_str) == 0) {
        mcutils::ParsingError(line_count, line, "Valid analytic basis types: oscillator|laguerre");
      }
      ket_basis_type = kBasisTypeDefinitions.at(basis_type_str);
      mcutils::FileExistCheck(ket_orbital_filename, true, false);
    } else if (keyword == "define-operator-target") {
      // define-operator-target operator_type order j0 output_filename
      //   operator_type = r|k
      char operator_type;
      std::string output_filename;
      int order;
      int j0, g0, Tz0;

      line_stream >> operator_type >> order >> j0 >> output_filename;
      mcutils::ParsingCheck(line_stream, line_count, line);

      g0 = j0 % 2;
      Tz0 = 0;

      mcutils::FileExistCheck(output_filename, false, true);

      operator_parameters.push_back(
          {ket_orbital_filename, ket_orbital_filename, output_filename,
           ket_basis_type, ket_basis_type, 1.0, OperationMode::kOperator,
           static_cast<shell::RadialOperatorType>(operator_type), order, j0, g0,
           Tz0});
    } else if (keyword == "define-radial-target") {
      // define-radial-target operator_type order j0 g0 Tz0 output_filename
      //   operator_type = r|k|o
      char operator_type;
      std::string output_filename;
      int order;
      int j0, g0, Tz0;

      line_stream >> operator_type >> order >> j0 >> g0 >> Tz0 >> output_filename;
      mcutils::ParsingCheck(line_stream, line_count, line);

      // special case: treat order=0, j0=0, g0=0, Tz0=1 as overlaps
      if ((order == 0) && (j0 == 0) && (g0 == 0) && (Tz0 == 1)) {
        operator_type = 'o';
      }

      mcutils::FileExistCheck(output_filename, false, true);

      operator_parameters.push_back(
          {ket_orbital_filename, ket_orbital_filename, output_filename,
           ket_basis_type, ket_basis_type, 1.0, OperationMode::kRadial,
           static_cast<shell::RadialOperatorType>(operator_type), order, j0, g0,
           Tz0});
    } else if (keyword == "define-xform-target") {
      // define-xform-target scale_factor bra_basis_type bra_orbital_file output_filename
      std::string basis_type_str, output_filename;
      shell::RadialBasisType bra_basis_type;
      int order=0, j0 = 0, g0 = 0, Tz0 = 0;
      double scale_factor;

      line_stream >> scale_factor >> basis_type_str >> bra_orbital_filename
          >> output_filename;
      mcutils::ParsingCheck(line_stream, line_count, line);

      if (kBasisTypeDefinitions.count(basis_type_str) == 0) {
        mcutils::ParsingError(line_count, line, "Valid analytic basis types: oscillator|laguerre");
      }
      bra_basis_type = kBasisTypeDefinitions.at(basis_type_str);

      operator_parameters.push_back(
          {bra_orbital_filename, ket_orbital_filename, output_filename,
           bra_basis_type, ket_basis_type, scale_factor, OperationMode::kXform,
           shell::RadialOperatorType::kO, order, j0, g0, Tz0});
    } else {
      mcutils::ParsingError(line_count, line, "Unknown keyword");
    }
  }
}

void BuildOperator(OperatorParameters operator_parameters)
{
  // Read orbitals
  std::ifstream bra_orbital_stream(operator_parameters.bra_orbital_filename);
  auto bra_input_orbitals = basis::ParseOrbitalPNStream(bra_orbital_stream, true);
  std::ifstream ket_orbital_stream(operator_parameters.ket_orbital_filename);
  auto ket_input_orbitals = basis::ParseOrbitalPNStream(ket_orbital_stream, true);
  // Construct indexing
  basis::OrbitalSpaceLJPN bra_space(bra_input_orbitals);
  basis::OrbitalSpaceLJPN ket_space(ket_input_orbitals);
  basis::OrbitalSectorsLJPN sectors = basis::OrbitalSectorsLJPN(
      bra_space, ket_space, operator_parameters.j0, operator_parameters.g0,
      operator_parameters.Tz0);
  basis::OperatorBlocks<double> matrices;

  // main control logic
  if (operator_parameters.mode == OperationMode::kOperator) {
    shell::SolidHarmonicOneBodyOperator(
      operator_parameters.ket_basis_type,
      operator_parameters.radial_operator,
      operator_parameters.order,
      ket_space,
      sectors,
      matrices);
  } else if (operator_parameters.mode == OperationMode::kRadial) {
    shell::GenerateRadialOperator(
      operator_parameters.ket_basis_type,
      operator_parameters.radial_operator,
      operator_parameters.order,
      ket_space,
      sectors,
      matrices);
  } else if (operator_parameters.mode == OperationMode::kXform) {
    shell::GenerateRadialOverlaps(
      operator_parameters.bra_basis_type,
      operator_parameters.ket_basis_type,
      operator_parameters.scale_factor,
      bra_space,
      ket_space,
      sectors,
      matrices);
  }

  // file mode flag
  basis::OneBodyOperatorType operator_type;
  if (operator_parameters.mode == OperationMode::kOperator) {
    operator_type = basis::OneBodyOperatorType::kSpherical;
  } else {
    operator_type = basis::OneBodyOperatorType::kRadial;
  }

  // write out to file
  shell::OutOBMEStream os(
      operator_parameters.output_filename, bra_space, ket_space, sectors,
      operator_type, operator_parameters.radial_operator,
      operator_parameters.order);
  os.Write(matrices);
}

int main(int argc, char** argv)
{
  // header
  std::cout << std::endl;
  std::cout << "radial-gen -- radial integral evaluation" << std::endl;
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  // process arguments
  std::vector<OperatorParameters> operators;
  ReadParameters(operators);

  // parallel performance diagnostic
  std::cout << fmt::format(
                   "INFO: OMP max_threads {}, num_procs {}",
                   omp_get_max_threads(), omp_get_num_procs())
            << std::endl
            << std::endl;

  // Eigen initialization
  Eigen::initParallel();

  for (auto& operator_parameters : operators) {
    BuildOperator(operator_parameters);
  }

  return EXIT_SUCCESS;
}
