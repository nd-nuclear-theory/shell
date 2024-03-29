/******************************************************************************
  @file radial-gen.cpp

  compute radial and one-body operator matrix elements

  Syntax:
    + radial-gen

  Input format:
    set-ket-basis <basis_type> <orbital_filename>
      basis_type = oscillator|laguerre
    define-operator-target <mode> <id> ...
      define-operator-target kinematic <id> <output_filename>
        id = identity|r|ik|r.r|ik.ik
      define-operator-target solid-harmonic <id> <order> <output_filename>
        id = r|ik
      define-operator-source am <id> <output_filename>
        id = l|l2|s|s2|j|j2
      define-operator-target isospin <id> <output_filename>
        id = tz|t+|t-
    define-radial-target <operator_type> <order> <J0> <g0> <Tz0> <output_filename>
      operator_type = r|k
    define-xform-target <scale_factor> <bra_basis_type> <bra_orbital_file> <output_filename>

  @note `solid-harmonic` mode generates matrix elements of r^n Y_n (using the
  more common normalization of the spherical harmonic) rather than matrix
  elements of r^n C_n (a.k.a. Racah's normalization).

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
  + 08/14/29 (pjf): Add additional operator types.
  + 08/16/19 (pjf): Remove radial operator type and power from OutOBMEStream.
  + 08/27/19 (pjf): Separate out solid harmonic generation as a separate mode.
  + 08/28/19 (pjf): Remove operator species option from am mode.
  + 09/04/19 (pjf): Fix typo in kKinematicOneBodyOperatorDefinitions for "k".
  + 04/19/22 (pjf): Clarify solid harmonic normalization convention.

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
#include "basis/proton_neutron.h"
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

enum class OperationMode {
  kKinematic, kSolidHarmonic, kAngularMomentum, kIsospin, kRadial, kXform
};

std::unordered_map<std::string,std::tuple<shell::RadialOperatorType,int,int>>
kKinematicOneBodyOperatorDefinitions =
  {
    // {id, {operator type, order, J0}}
    {"identity", {shell::RadialOperatorType::kO, 0, 0}},
    {"r", {shell::RadialOperatorType::kR, 1, 1}},
    {"ik", {shell::RadialOperatorType::kK, 1, 1}},
    {"r.r", {shell::RadialOperatorType::kR, 2, 0}},
    {"ik.ik", {shell::RadialOperatorType::kK, 2, 0}},
  };

std::unordered_map<std::string,std::tuple<am::AngularMomentumOperatorType,int,int>>
kAngularMomentumOneBodyOperatorDefinitions =
  {
    {"l",  {am::AngularMomentumOperatorType::kOrbital, 1, 1}},
    {"l2", {am::AngularMomentumOperatorType::kOrbital, 2, 0}},
    {"s",  {am::AngularMomentumOperatorType::kSpin,    1, 1}},
    {"s2", {am::AngularMomentumOperatorType::kSpin,    2, 0}},
    {"j",  {am::AngularMomentumOperatorType::kTotal,   1, 1}},
    {"j2", {am::AngularMomentumOperatorType::kTotal,   2, 0}}
  };

std::unordered_map<std::string,int>
kIsospinOneBodyOperatorDefinitions =
  {
    {"tz",  0},
    {"t+", +1},
    {"t-", -1}
  };

std::unordered_map<std::string, shell::RadialBasisType>
kBasisTypeDefinitions =
  {
    {"oscillator", shell::RadialBasisType::kOscillator},
    {"laguerre", shell::RadialBasisType::kLaguerre}
  };

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
  union {
    shell::RadialOperatorType radial_operator_type;
    am::AngularMomentumOperatorType am_operator_type;
  };
  int order;
  int J0;
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
    if (keyword == "set-ket-basis")
    {
      // set-ket-basis basis_type orbital_filename
      //   basis_type = oscillator|laguerre
      std::string basis_type_str;
      line_stream >> basis_type_str >> ket_orbital_filename;
      mcutils::ParsingCheck(line_stream, line_count, line);
      // convert basis
      if (kBasisTypeDefinitions.count(basis_type_str) == 0)
      {
        mcutils::ParsingError(
            line_count, line,
            "Valid analytic basis types: oscillator|laguerre"
          );
      }
      ket_basis_type = kBasisTypeDefinitions.at(basis_type_str);
      mcutils::FileExistCheck(ket_orbital_filename, true, false);
    }
    else if (keyword == "define-operator-target")
    {
      // define-operator-target <mode> ...
      OperatorParameters parameters;
      std::string mode, id;

      // set basic parameters
      parameters.bra_orbital_filename = ket_orbital_filename;
      parameters.ket_orbital_filename = ket_orbital_filename;
      parameters.bra_basis_type = ket_basis_type;
      parameters.ket_basis_type = ket_basis_type;
      parameters.scale_factor = 1.0;

      line_stream >> mode >> id;
      mcutils::ParsingCheck(line_stream,line_count,line);

      if (mode == "kinematic")
      {
        if (!kKinematicOneBodyOperatorDefinitions.count(id))
          mcutils::ParsingError(line_count,line,"Unrecognized kinematic operator ID");
        parameters.mode = OperationMode::kKinematic;
        std::tie(parameters.radial_operator_type, parameters.order, parameters.J0)
          = kKinematicOneBodyOperatorDefinitions.at(id);

        parameters.g0 = parameters.J0%2;
        parameters.Tz0 = 0;
      }
      else if (mode == "solid-harmonic")
      {
        int order;
        line_stream >> order;
        mcutils::ParsingCheck(line_stream,line_count,line);

        parameters.mode = OperationMode::kSolidHarmonic;
        parameters.radial_operator_type = shell::kCharCodeRadialOperatorType.at(id);
        parameters.order = parameters.J0 = order;
        parameters.g0 = parameters.J0%2;
        parameters.Tz0 = 0;
      }
      else if (mode == "am")
      {
        std::string operator_species_s;
        if (!kAngularMomentumOneBodyOperatorDefinitions.count(id))
          mcutils::ParsingError(line_count,line,"Unrecognized angular momentum operator ID");
        parameters.mode = OperationMode::kAngularMomentum;
        std::tie(parameters.am_operator_type, parameters.order, parameters.J0)
          = kAngularMomentumOneBodyOperatorDefinitions.at(id);
        parameters.g0 = 0;
        parameters.Tz0 = 0;
      }
      else if (mode == "isospin")
      {
        if (!kIsospinOneBodyOperatorDefinitions.count(id))
          mcutils::ParsingError(line_count,line,"Unrecognized isospin operator ID");
        parameters.mode = OperationMode::kIsospin;
        parameters.radial_operator_type = shell::RadialOperatorType::kGeneric;
        parameters.order = 0;
        parameters.J0 = 0;
        parameters.g0 = 0;
        parameters.Tz0 = kIsospinOneBodyOperatorDefinitions.at(id);
      }

      line_stream >> parameters.output_filename;
      mcutils::ParsingCheck(line_stream,line_count,line);
      mcutils::FileExistCheck(parameters.output_filename, false, true);

      operator_parameters.push_back(parameters);
    }
    else if (keyword == "define-radial-target")
    {
      // define-radial-target <operator_type> <order> <J0> <g0> <Tz0> <output_filename>
      //   operator_type = r|k
      OperatorParameters parameters;
      std::string operator_type;

      // set basic parameters
      parameters.bra_orbital_filename = ket_orbital_filename;
      parameters.ket_orbital_filename = ket_orbital_filename;
      parameters.bra_basis_type = ket_basis_type;
      parameters.ket_basis_type = ket_basis_type;
      parameters.scale_factor = 1.0;
      parameters.mode = OperationMode::kRadial;

      // get parameters from input line
      line_stream >> operator_type >> parameters.order
                  >> parameters.J0 >> parameters.g0 >> parameters.Tz0
                  >> parameters.output_filename;
      mcutils::ParsingCheck(line_stream, line_count, line);

      // fill remaining parameters
      parameters.radial_operator_type = shell::kCharCodeRadialOperatorType.at(operator_type);

      mcutils::FileExistCheck(parameters.output_filename, false, true);

      operator_parameters.push_back(parameters);
    } else if (keyword == "define-xform-target") {
      // define-xform-target <scale_factor> <bra_basis_type> <bra_orbital_file> <output_filename>
      OperatorParameters parameters;
      std::string basis_type_str;

      // set basic parameters
      parameters.ket_orbital_filename = ket_orbital_filename;
      parameters.ket_basis_type = ket_basis_type;
      parameters.mode = OperationMode::kXform;
      parameters.order = 0;
      parameters.J0 = 0;
      parameters.g0 = 0;
      parameters.Tz0 = 0;

      line_stream >> parameters.scale_factor
                  >> basis_type_str
                  >> parameters.bra_orbital_filename
                  >> parameters.output_filename;
      mcutils::ParsingCheck(line_stream, line_count, line);

      if (kBasisTypeDefinitions.count(basis_type_str) == 0) {
        mcutils::ParsingError(line_count, line, "Valid analytic basis types: oscillator|laguerre");
      }
      parameters.bra_basis_type = kBasisTypeDefinitions.at(basis_type_str);

      mcutils::FileExistCheck(parameters.bra_orbital_filename, true, false);
      mcutils::FileExistCheck(parameters.output_filename, false, true);

      operator_parameters.push_back(parameters);
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
      bra_space, ket_space, operator_parameters.J0, operator_parameters.g0,
      operator_parameters.Tz0);
  basis::OperatorBlocks<double> matrices;

  // main control logic
  if (operator_parameters.mode == OperationMode::kKinematic)
  {
    shell::SolidHarmonicOneBodyOperator(
        operator_parameters.ket_basis_type,
        operator_parameters.radial_operator_type,
        operator_parameters.order,
        ket_space,
        sectors,
        matrices
      );
  }
  else if (operator_parameters.mode == OperationMode::kSolidHarmonic)
  {
    shell::SolidHarmonicOneBodyOperator(
        operator_parameters.ket_basis_type,
        operator_parameters.radial_operator_type,
        operator_parameters.order,
        ket_space,
        sectors,
        matrices
      );

    // difference between C and Y -- Brink & Satchler (1993), app. IV,  p.145
    basis::ScalarMultiplyOperator(
        sectors, matrices,
        Hat(operator_parameters.order) * am::kInvSqrt4Pi
      );
  }
  else if (operator_parameters.mode == OperationMode::kAngularMomentum)
  {
    if (operator_parameters.order == 1)
    {
      shell::AngularMomentumOneBodyOperator(
          operator_parameters.am_operator_type,
          ket_space, sectors, matrices
        );
    }
    else if (operator_parameters.order == 2)
    {
      shell::AngularMomentumSquaredOneBodyOperator(
          operator_parameters.am_operator_type,
          ket_space, sectors, matrices
        );
    }
  }
  else if (operator_parameters.mode == OperationMode::kIsospin)
  {
    shell::IsospinOneBodyOperator(ket_space, sectors, matrices);
  }
  else if (operator_parameters.mode == OperationMode::kRadial)
  {
    shell::GenerateRadialOperator(
        operator_parameters.ket_basis_type,
        operator_parameters.radial_operator_type,
        operator_parameters.order,
        ket_space,
        sectors,
        matrices
      );
  }
  else if (operator_parameters.mode == OperationMode::kXform)
  {
    shell::GenerateRadialOverlaps(
        operator_parameters.bra_basis_type,
        operator_parameters.ket_basis_type,
        operator_parameters.scale_factor,
        bra_space,
        ket_space,
        sectors,
        matrices
      );
  }

  // file mode flag
  basis::OneBodyOperatorType operator_type;
  if ((operator_parameters.mode == OperationMode::kKinematic)
      ||(operator_parameters.mode == OperationMode::kSolidHarmonic)
      ||(operator_parameters.mode == OperationMode::kAngularMomentum)
      ||(operator_parameters.mode == OperationMode::kIsospin)) {
    operator_type = basis::OneBodyOperatorType::kSpherical;
  } else {
    operator_type = basis::OneBodyOperatorType::kRadial;
  }

  // write out to file
  shell::OutOBMEStream os(
      operator_parameters.output_filename,
      bra_space, ket_space, sectors,
      operator_type
    );
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
