/****************************************************************
  relative-gen.cpp

  Compute relative operator matrix elements.

  See lsjt_operator.h for documentation of operator storage and the
  relative operator file format.

  Standard input:
    J0 g0 T0_min T0_max
    Nmax
    operator_name [parameters]
    target_filename

  The parameters are

    J0 g0 T0_min T0_max : relative operator tensor properties

    Nmax : relative basis trunctation

  The operator_name line may be:

    zero

      Zero operator.

    identity

      Identity operator.

    ksqr

      Relative k^2 operator (~intrinsic kinetic energy).
      [DEPRECATED in favor of coordinate-sqr]

    rsqr

      Relative r^2 operator (~intrinsic r^2 operator).
      [DEPRECATED in favor of coordinate-sqr]

    coordinate-sqr r|k T0

      Relative r^2 or k^2 operator.
      [TODO: implement isovector variant]

    dipole r|k T0

      Relative electric dipole operator.
      [NOTE: only "r" version implemented; only T0=1 defined]

    quadrupole r|k T0

      Relative electric quadrupole operator.
      [NOTE: only "r" version implemented]

    orbital-am T0

      Relative orbital angular momentum (orbital part of magnetic dipole)
      operator.

    spin-am T0

      Spin operator (spin part of magnetic dipole), construed as relative two
      body operator.

    isospin

      Isospin operator, construed as relative two-body operator.

    su4casimir

      SU(4) quadratic Casimir operator, construed as relative two-body operator.

    coulomb p|n|total steps

      Coulomb potential.

    symmunit T0 Np Lp Sp Jp Tp N L S J T

      Symmetrized unit tensor.  The labels must specify a canonical
      (upper triangle) matrix element.

    interaction <interaction>

      Two-body interaction.
      interaction = Daejeon16

    LENPIC-LOGT

      LENPIC LO Gamow-Teller operator, construed as relative two-body operator

    LENPIC-N2LOGT regulator oscillator_length steps

      LENPIC N2LO Gamow-Teller operator
      Note: oscillator_length is the *single-particle* oscillator length

  Other operator parameter values are taken as:
    symmetry_phase_mode = kHermitian
    Jmax = Nmax+1

  Example (coulomb_Nmax20_p.in):
    0 0 0 2
    20
    coulomb p 500
    coulomb_Nmax6_rel.dat


  Language: C++11

  Mark A. Caprio
  University of Notre Dame

  + 07/16/16 (mac): Created (writerel.cpp).
  + 07/25/16 (mac): Update to use WriteRelativeOperatorLSJT.
  + 10/09/16 (pjf): Rename mcpp -> mcutils.
  + 03/06/17 (mac): Rough in Coulomb interaction code.
  + 03/26/17 (mac): Finish implementing Coulomb.  Rename to relative-gen.cpp.
  + 04/05/17 (mac): Update call to ConstructZeroOperatorRelativeLSJT.
  + 11/28/17 (pjf): Print header with version.
  + 05/04/17 (mac): Add support for isovector transition operators.
  + 05/09/19 (pjf): Use std::size_t for basis indices and sizes.
  + 05/15/19 (mac): Rename "matrices" to "blocks" in variable names.
  + 06/20/19 (pjf): Add isospin operator.
  + 10/30/20 (pjf): Add (optional) Daejeon16 interaction.
  + 10/31/20 (pjf): Remove Jmax option from interaction generation.
  + 06/02/22 (pjf): Add LENPIC N2LO Gamow-Teller operator.
  + 08/08/22 (pjf): Add LENPIC LO Gamow-Teller operator.

****************************************************************/

#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#include "basis/lsjt_operator.h"
#include "basis/proton_neutron.h"
#include "fmt/format.h"
#include "mcutils/parsing.h"
#include "relative/relative_me.h"
#include "relative/lenpic_relative_me.h"
#include "spline/wavefunction_class.h"

#ifdef USE_DAEJEON16
#include "Daejeon16/Daejeon16_wrapper.h"
#endif  // USE_DAEJEON16

////////////////////////////////////////////////////////////////
// parameter input
////////////////////////////////////////////////////////////////

struct UnitTensorLabels
// Labels for a single relative matrix element, as needed to define a
// unit tensor.
{
  int T0, Np, Lp, Sp, Jp, Tp, N, L, S, J, T;
};

struct Parameters
// Container for run input parameters.
{
  basis::RelativeOperatorParametersLSJT operator_parameters;
  std::string operator_name;
  std::string target_filename;

  // optional parameters for specific operators
  std::string interaction_name;
  UnitTensorLabels unit_tensor_labels;
  basis::OperatorTypePN operator_type_pn;
  relative::CoordinateType coordinate_type;
  double regulator_parameter;
  double oscillator_length;
  int num_steps;
  int T0;
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
    std::getline(std::cin,line);
    std::istringstream line_stream(line);
    line_stream >> parameters.operator_parameters.J0
                >> parameters.operator_parameters.g0
                >> parameters.operator_parameters.T0_min
                >> parameters.operator_parameters.T0_max;
    mcutils::ParsingCheck(line_stream,line_count,line);
  }

  // line 2: operator basis parameters
  {
    ++line_count;
    std::getline(std::cin,line);
    std::istringstream line_stream(line);
    line_stream >> parameters.operator_parameters.Nmax;
    mcutils::ParsingCheck(line_stream,line_count,line);
  }

  // set miscellaneous operator_parameters fields
  parameters.operator_parameters.Jmax = parameters.operator_parameters.Nmax+1;
  parameters.operator_parameters.symmetry_phase_mode = basis::SymmetryPhaseMode::kHermitian;

  // impose current limitations on tensor character
  // assert(parameters.operator_parameters.g0==0);
  // assert(parameters.operator_parameters.T0_min==0);

  // line 3: operator choice and parameters
  {
    ++line_count;
    std::getline(std::cin,line);
    std::istringstream line_stream(line);
    line_stream >> parameters.operator_name;
    mcutils::ParsingCheck(line_stream,line_count,line);

    // parse additional operator parameters

    if (
        (parameters.operator_name == "orbital-am")
        || (parameters.operator_name == "spin-am")
      )
      // arguments: T0
      {
        line_stream
          >> parameters.T0;
        mcutils::ParsingCheck(line_stream,line_count,line);
      }
    else if (
        (parameters.operator_name == "coordinate-sqr")
        || (parameters.operator_name == "dipole")
        || (parameters.operator_name == "quadrupole")
      )
      // arguments: r|k T0
      {
        std::string coordinate_type_string;
        line_stream
          >> coordinate_type_string
          >> parameters.T0;
        mcutils::ParsingCheck(line_stream,line_count,line);
        if (coordinate_type_string=="r")
          parameters.coordinate_type = relative::CoordinateType::kR;
        else if (coordinate_type_string=="k")
          parameters.coordinate_type = relative::CoordinateType::kK;
        else
          mcutils::ParsingError(line_count,line,"invalid operator type (r/k)");
      }
    else if (
        (parameters.operator_name == "coulomb")
      )
      // coulomb
      {
        std::string operator_type_pn_string;
        line_stream
          >> operator_type_pn_string
          >> parameters.num_steps;
        mcutils::ParsingCheck(line_stream,line_count,line);
        if (operator_type_pn_string=="p")
          parameters.operator_type_pn = basis::OperatorTypePN::kP;
        else if (operator_type_pn_string=="n")
          parameters.operator_type_pn = basis::OperatorTypePN::kN;
        else if (operator_type_pn_string=="total")
          parameters.operator_type_pn = basis::OperatorTypePN::kTotal;
        else
          mcutils::ParsingError(line_count,line,"invalid operator type (p/n/total)");
      }
    else if (parameters.operator_name == "symmunit")
      // special for "symmunit"
      {
        line_stream
          >> parameters.unit_tensor_labels.T0
          >> parameters.unit_tensor_labels.Np
          >> parameters.unit_tensor_labels.Lp
          >> parameters.unit_tensor_labels.Sp
          >> parameters.unit_tensor_labels.Jp
          >> parameters.unit_tensor_labels.Tp
          >> parameters.unit_tensor_labels.N
          >> parameters.unit_tensor_labels.L
          >> parameters.unit_tensor_labels.S
          >> parameters.unit_tensor_labels.J
          >> parameters.unit_tensor_labels.T;
        mcutils::ParsingCheck(line_stream,line_count,line);
      }
    else if (parameters.operator_name == "interaction")
      // arguments: interaction name
      {
        line_stream >> parameters.interaction_name;
        mcutils::ParsingCheck(line_stream,line_count,line);
#ifdef USE_DAEJEON16
        if (parameters.interaction_name == "Daejeon16")
          parameters.operator_parameters.Jmax = std::min(parameters.operator_parameters.Jmax, 6);
#endif  // USE_DAEJEON16
      }
    else if (parameters.operator_name == "LENPIC-N2LOGT")
      {
        line_stream >> parameters.regulator_parameter
                    >> parameters.oscillator_length
                    >> parameters.num_steps;
        mcutils::ParsingCheck(line_stream,line_count,line);
      }
  }

  // line 4: output filename
  {
    ++line_count;
    std::getline(std::cin,line);
    std::istringstream line_stream(line);
    line_stream >> parameters.target_filename;
    mcutils::ParsingCheck(line_stream,line_count,line);
  }

}

////////////////////////////////////////////////////////////////
// operator work
////////////////////////////////////////////////////////////////

void PopulateOperator(
    const Parameters& parameters,
    basis::RelativeSpaceLSJT& relative_space,
    std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
    std::array<basis::OperatorBlocks<double>,3>& relative_component_blocks
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
  const basis::RelativeOperatorParametersLSJT& operator_parameters
    = parameters.operator_parameters;

  // set up relative space
  relative_space = basis::RelativeSpaceLSJT(
      operator_parameters.Nmax,operator_parameters.Jmax
    );

  // populate operator containers
  if (parameters.operator_name == "zero")
    {
      basis::ConstructZeroOperatorRelativeLSJT(
          operator_parameters,
          relative_space,relative_component_sectors,relative_component_blocks
        );
    }
  else if (parameters.operator_name == "identity")
    {
      basis::ConstructIdentityOperatorRelativeLSJT(
          operator_parameters,
          relative_space,relative_component_sectors,relative_component_blocks
        );
    }
  else if (parameters.operator_name == "rsqr")  // DEPRECATED
    {
      relative::ConstructCoordinateSqr(
          operator_parameters,
          relative_space,relative_component_sectors,relative_component_blocks,
          relative::CoordinateType::kR,
          0  // T0
        );
    }
  else if (parameters.operator_name == "ksqr")  // DEPRECATED
    {
      relative::ConstructCoordinateSqr(
          operator_parameters,
          relative_space,relative_component_sectors,relative_component_blocks,
          relative::CoordinateType::kK,
          0  // T0
        );
    }
  else if (parameters.operator_name == "coordinate-sqr")
    {
      relative::ConstructCoordinateSqr(
          operator_parameters,
          relative_space,relative_component_sectors,relative_component_blocks,
          parameters.coordinate_type,
          parameters.T0
        );
    }
  else if (parameters.operator_name == "dipole")
    {
      relative::ConstructDipoleOperator(
          operator_parameters,
          relative_space,relative_component_sectors,relative_component_blocks,
          parameters.coordinate_type,
          parameters.T0
        );
    }
  else if (parameters.operator_name == "quadrupole")
    {
      relative::ConstructQuadrupoleOperator(
          operator_parameters,
          relative_space,relative_component_sectors,relative_component_blocks,
          parameters.coordinate_type,
          parameters.T0
        );
    }
  else if (parameters.operator_name == "orbital-am")
    {
      relative::ConstructOrbitalAMOperator(
          operator_parameters,
          relative_space,relative_component_sectors,relative_component_blocks,
          parameters.T0
        );
    }
  else if (parameters.operator_name == "spin-am")
    {
      relative::ConstructSpinAMOperator(
          operator_parameters,
          relative_space,relative_component_sectors,relative_component_blocks,
          parameters.T0
        );
    }
  else if (parameters.operator_name == "isospin")
    {
      relative::ConstructIsospinOperator(
          operator_parameters,
          relative_space,relative_component_sectors,relative_component_blocks
        );
    }
  else if (parameters.operator_name == "su4casimir")
    {
      relative::ConstructSU4CasimirOperator(
          operator_parameters,
          relative_space,relative_component_sectors,relative_component_blocks
        );
    }
  else if (parameters.operator_name == "coulomb")
    {
      // 500 steps seems to suffice for ~8 digits precision at Nmax20
      relative::ConstructCoulombOperator(
          operator_parameters,
          relative_space,relative_component_sectors,relative_component_blocks,
          parameters.operator_type_pn,
          parameters.num_steps
        );
    }
  else if (parameters.operator_name == "symmunit")
    {
      // start from zero operator
      basis::ConstructZeroOperatorRelativeLSJT(
          operator_parameters,
          relative_space,relative_component_sectors,relative_component_blocks
        );

      // define shortcut reference to unit tensor labels
      const UnitTensorLabels& unit_tensor_labels = parameters.unit_tensor_labels;

      // look up subspace indices
      int gp = unit_tensor_labels.Lp%2;
      std::size_t relative_subspace_index_bra = relative_space.LookUpSubspaceIndex(
          basis::RelativeSubspaceLSJTLabels(
              unit_tensor_labels.Lp,
              unit_tensor_labels.Sp,
              unit_tensor_labels.Jp,
              unit_tensor_labels.Tp,
              gp
            )
        );
      const basis::RelativeSubspaceLSJT& relative_subspace_bra
        = relative_space.GetSubspace(relative_subspace_index_bra);
      std::cout << "Bra subspace:"
                << " " << relative_subspace_bra.LabelStr()
                << " -> index " << relative_subspace_index_bra
                << std::endl;
      int g = unit_tensor_labels.L%2;
      std::size_t relative_subspace_index_ket = relative_space.LookUpSubspaceIndex(
          basis::RelativeSubspaceLSJTLabels(
              unit_tensor_labels.L,
              unit_tensor_labels.S,
              unit_tensor_labels.J,
              unit_tensor_labels.T,
              g
            )
        );
      const basis::RelativeSubspaceLSJT& relative_subspace_ket
        = relative_space.GetSubspace(relative_subspace_index_ket);
      std::cout << "Ket subspace:"
                << " " << relative_subspace_ket.LabelStr()
                << " -> index " << relative_subspace_index_ket
                << std::endl;

      // look up state indices
      std::size_t relative_state_index_bra = relative_subspace_bra.LookUpStateIndex(
          basis::RelativeStateLSJT::StateLabelsType(
              unit_tensor_labels.Np
            )
        );
      std::cout << "Bra state:"
                << " " << unit_tensor_labels.Np
                << " -> index " << relative_state_index_bra
                << std::endl;
      std::size_t relative_state_index_ket = relative_subspace_ket.LookUpStateIndex(
          basis::RelativeStateLSJT::StateLabelsType(
              unit_tensor_labels.N
            )
        );
      std::cout << "ket state:"
                << " " << unit_tensor_labels.N
                << " -> index " << relative_state_index_ket
                << std::endl;

      // select T0 component
      const basis::RelativeSectorsLSJT& relative_sectors = relative_component_sectors[unit_tensor_labels.T0];
      basis::OperatorBlocks<double>& relative_blocks = relative_component_blocks[unit_tensor_labels.T0];

      // look up sector
      std::size_t relative_sector_index
        = relative_sectors.LookUpSectorIndex(
            relative_subspace_index_bra,
            relative_subspace_index_ket
          );


      // check that matrix element is canonical
      if (relative_sectors.GetSector(relative_sector_index).IsDiagonal())
        assert(relative_state_index_bra<=relative_state_index_ket);

      // set matrix element
      relative_blocks[relative_sector_index](relative_state_index_bra,relative_state_index_ket) = 1.;

    }
  else if (parameters.operator_name == "interaction")
    {
      std::unordered_set<std::string> known_interactions;
      #ifdef USE_DAEJEON16
      known_interactions.insert("Daejeon16");
      if (parameters.interaction_name == "Daejeon16")
        {
          contrib::ConstructDaejeon16Operator(
              operator_parameters,
              relative_space,relative_component_sectors,relative_component_blocks
            );
        }
      #endif  // USE_DAEJEON16

      // exit with failure if interaction not recognized
      if (!known_interactions.count(parameters.interaction_name))
        {
          std::cerr << "ERROR: Unknown interaction " << parameters.interaction_name << std::endl;
          std::cerr << "  Known interactions:";
          for (const auto& interaction_name : known_interactions)
            std::cerr << " " << interaction_name;
          std::cerr << std::endl;
          std::exit(EXIT_FAILURE);
        }
    }
  else if (parameters.operator_name == "LENPIC-LOGT")
    {
      relative::lenpic::ConstructLOGTOperator(
          operator_parameters,
          relative_space, relative_component_sectors, relative_component_blocks
        );
    }
  else if (parameters.operator_name == "LENPIC-N2LOGT")
    {
      relative::lenpic::ConstructN2LOGTOperator(
          operator_parameters,
          relative_space, relative_component_sectors, relative_component_blocks,
          parameters.regulator_parameter, parameters.oscillator_length, parameters.num_steps
        );
    }
  else
    {
      std::cerr << "ERROR: unknown operator name " << parameters.operator_name << std::endl;
      std::exit(EXIT_FAILURE);
    }
}

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  // header
  std::cout << std::endl;
  std::cout << "relative-gen -- relative operator matrix element generation" << std::endl;
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  // parameter input
  Parameters parameters;
  ReadParameters(parameters);

  // set up operator
  basis::RelativeSpaceLSJT relative_space;
  std::array<basis::RelativeSectorsLSJT,3> relative_component_sectors;
  std::array<basis::OperatorBlocks<double>,3> relative_component_blocks;
  PopulateOperator(
      parameters,
      relative_space,
      relative_component_sectors,relative_component_blocks
    );

  // write operator
  basis::WriteRelativeOperatorLSJT(
      parameters.target_filename,
      relative_space,
      parameters.operator_parameters,  // only need operator labels
      relative_component_sectors,
      relative_component_blocks,
      true  // verbose
    );

  // termination
  return 0;
}
