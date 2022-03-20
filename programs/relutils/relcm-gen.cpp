/****************************************************************
  relcm-gen.cpp

  Compute relative-cm operator matrix elements.

  See lsjt_operator.h for documentation of operator storage and the
  relative-cm operator file format.

  Standard input:
    J0 g0 T0_min T0_max
    Nmax
    operator_name [parameters]
    target_filename

  The parameters are

    J0 g0 T0_min T0_max : relative operator tensor properties

    Nmax : relative basis trunctation

  The operator_name line may be:

    LENPIC-NLOM1 regulator oscillator_length steps

      LENPIC NLO M1 operator
      Note: oscillator_length is the *single-particle* oscillator length

  Example (LENPICchi2bM11C_hw20.0_Nmax04_relcm.in):
    1 0 1 1
    4
    LENPIC-NLOM1 1.0 1.43998333278936 500
    LENPICchi2bM11C_hw20.0_Nmax04_relcm.dat

  Language: C++11

  Mark A. Caprio
  University of Notre Dame

  + 03/04/22 (pjf): Created (based on relative-gen.cpp).

****************************************************************/

#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#include "basis/lsjt_operator.h"
#include "fmt/format.h"
#include "mcutils/parsing.h"
#include "relative/relcm_lenpic_me.h"

#ifdef USE_DAEJEON16
#  include "Daejeon16/Daejeon16_wrapper.h"
#endif  // USE_DAEJEON16

////////////////////////////////////////////////////////////////
// parameter input
////////////////////////////////////////////////////////////////

struct Parameters
// Container for run input parameters.
{
  basis::RelativeCMOperatorParametersLSJT operator_parameters;
  std::string operator_name;
  std::string target_filename;

  // optional parameters for specific operators
  double regulator_parameter;
  double oscillator_length;
  int num_steps;
};

void ReadParameters(Parameters& parameters)
// Read run parameters from stdin.
//
// Arguments:
//   parameters (Parameters, output) :
//     container for input parameters
{
  // set up line counter for use in error messages
  std::string line;
  int line_count = 0;


  // line 1: operator tensor properties
  {
    ++line_count;
    std::getline(std::cin, line);
    std::istringstream line_stream(line);
    line_stream >> parameters.operator_parameters.J0
        >> parameters.operator_parameters.g0
        >> parameters.operator_parameters.T0_min
        >> parameters.operator_parameters.T0_max;
    mcutils::ParsingCheck(line_stream, line_count, line);
  }

  // line 2: operator basis parameters
  {
    ++line_count;
    std::getline(std::cin, line);
    std::istringstream line_stream(line);
    line_stream >> parameters.operator_parameters.Nmax;
    mcutils::ParsingCheck(line_stream, line_count, line);
  }

  // set miscellaneous operator_parameters fields
  parameters.operator_parameters.symmetry_phase_mode =
      basis::SymmetryPhaseMode::kHermitian;

  // impose current limitations on tensor character
  // assert(parameters.operator_parameters.g0==0);
  // assert(parameters.operator_parameters.T0_min==0);

  // line 3: operator choice and parameters
  {
    ++line_count;
    std::getline(std::cin, line);
    std::istringstream line_stream(line);
    line_stream >> parameters.operator_name;
    mcutils::ParsingCheck(line_stream, line_count, line);

    // parse additional operator parameters

    if (parameters.operator_name == "LENPIC-NLOM1")
    // LENPIC NLO M1
    {
      line_stream >> parameters.regulator_parameter
          >> parameters.oscillator_length >> parameters.num_steps;
      mcutils::ParsingCheck(line_stream, line_count, line);
    }
  }

  // line 4: output filename
  {
    ++line_count;
    std::getline(std::cin, line);
    std::istringstream line_stream(line);
    line_stream >> parameters.target_filename;
    mcutils::ParsingCheck(line_stream, line_count, line);
  }
}

////////////////////////////////////////////////////////////////
// operator work
////////////////////////////////////////////////////////////////

void PopulateOperator(
    const Parameters& parameters,
    basis::RelativeCMSpaceLSJT& relative_cm_space,
    std::array<basis::RelativeCMSectorsLSJT, 3>& relative_cm_component_sectors,
    std::array<basis::OperatorBlocks<double>, 3>& relative_cm_component_blocks
  )
// Define operator.
//
// Arguments:
//   parameters (Parameters) : includes tensorial properties of operator
//      choice of operator to use
//   relative_space (..., output) : target space
//   relative_component_sectors (..., output) : target sectors
//   relative_component_blocks (..., output) : target blocks
{
  // define shortcut reference to operator parameters
  const basis::RelativeCMOperatorParametersLSJT& operator_parameters =
      parameters.operator_parameters;

  // set up relative space
  relative_cm_space = basis::RelativeCMSpaceLSJT(operator_parameters.Nmax);

  // populate operator containers
  // if (parameters.operator_name == "zero")
  //   {
  //     basis::ConstructZeroOperatorRelativeLSJT(
  //         operator_parameters,
  //         relative_space,relative_component_sectors,relative_component_blocks
  //       );
  //   }
  /* else */ if (parameters.operator_name == "LENPIC-NLOM1")
  {
    relative::lenpic::ConstructNLOM1Operator(
        operator_parameters,
        relative_cm_space,
        relative_cm_component_sectors,
        relative_cm_component_blocks,
        parameters.regulator_parameter,
        parameters.oscillator_length
      );
  }
  // else if (parameters.operator_name == "identity")
  // {
  //   basis::ConstructIdentityOperatorRelativeLSJT(
  //       operator_parameters,
  //       relative_space,
  //       relative_component_sectors,
  //       relative_component_blocks
  //     );
  // }
  else
  {
    std::cerr << "ERROR: unknown operator name " << parameters.operator_name
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
  // header
  std::cout << std::endl;
  std::cout << "relativecm-gen -- relative operator matrix element generation"
            << std::endl;
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  // parameter input
  Parameters parameters;
  ReadParameters(parameters);

  // set up operator
  basis::RelativeCMSpaceLSJT relative_cm_space;
  std::array<basis::RelativeCMSectorsLSJT, 3> relative_cm_component_sectors;
  std::array<basis::OperatorBlocks<double>, 3> relative_cm_component_blocks;
  PopulateOperator(
      parameters, relative_cm_space, relative_cm_component_sectors, relative_cm_component_blocks
    );

  // write operator
  basis::WriteRelativeCMOperatorLSJT(
      parameters.target_filename,
      relative_cm_space,
      parameters.operator_parameters,  // only need operator labels
      relative_cm_component_sectors,
      relative_cm_component_blocks,
      true  // verbose
    );

  // termination
  return 0;
}
