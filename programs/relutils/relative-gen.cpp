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

    rsqr
      Relative r^2 operator (~intrinsic r^2 operator).
      [TODO: merge with ksqr?  provide pp|nn variants?]

    qrel r|k pp|nn|pn [TODO]

    coulomb pp|nn|pn [TODO]
      Coulomb potential.

    symmunit T0 Np Lp Sp Jp Tp N L S J T
      Symmetrized unit tensor.  The labels must specify a canonical
      (upper triangle) matrix element.

  Other operator parameter values are taken as:
    symmetry_phase_mode = kHermitian
    Jmax = Nmax+1

  Example (writerel_coulomb_Nmax6.in):
    0 0 0 2
    6
    coulomb
    coulomb_Nmax6_rel.dat


  Language: C++11
                                 
  Mark A. Caprio
  University of Notre Dame

  7/16/16 (mac): Created (writerel.cpp).
  7/25/16 (mac): Update to use WriteRelativeOperatorLSJT.
  10/9/16 (pjf): Rename mcpp -> mcutils.
  3/6/17 (mac): Rough in Coulomb interaction code.
  3/26/17 (mac): Finish implementing Coulomb.  Rename to relative-gen.cpp.
  4/5/17 (mac): Update call to ConstructZeroOperatorRelativeLSJT.

****************************************************************/

#include "basis/lsjt_operator.h"
#include "cppformat/format.h"
#include "mcutils/parsing.h"
#include "relative/relative_me.h"
#include "spline/wavefunction_class.h"

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
  UnitTensorLabels unit_tensor_labels;
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
    ParsingCheck(line_stream,line_count,line);
  }

  // line 2: operator basis parameters
  {
    ++line_count;
    std::getline(std::cin,line);
    std::istringstream line_stream(line);
    line_stream >> parameters.operator_parameters.Nmax;
    ParsingCheck(line_stream,line_count,line);
  }

  // set miscellaneous operator_parameters fields
  parameters.operator_parameters.Jmax = parameters.operator_parameters.Nmax+1;
  parameters.operator_parameters.symmetry_phase_mode = basis::SymmetryPhaseMode::kHermitian;

  // impose current limitations on tensor character
  assert(parameters.operator_parameters.J0==0);
  assert(parameters.operator_parameters.g0==0);
  assert(parameters.operator_parameters.T0_min==0);

  // line 3: operator choice and parameters
  {
    ++line_count;
    std::getline(std::cin,line);
    std::istringstream line_stream(line);
    line_stream >> parameters.operator_name;
    ParsingCheck(line_stream,line_count,line);

    // additional operator parameters
    if (parameters.operator_name == "symmunit")
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
        ParsingCheck(line_stream,line_count,line);
      }
  }

  // line 4: output filename
  {
    ++line_count;
    std::getline(std::cin,line);
    std::istringstream line_stream(line);
    line_stream >> parameters.target_filename;
    ParsingCheck(line_stream,line_count,line);
  }

}

////////////////////////////////////////////////////////////////
// operator work
////////////////////////////////////////////////////////////////

void PopulateOperator(
    const Parameters& parameters,
    basis::RelativeSpaceLSJT& relative_space,
    std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
    std::array<basis::MatrixVector,3>& relative_component_matrices
  )
// Define operator.
//
// Arguments:
//   parameters (Parameters) : includes tensorial properties of operator
//      choice of operator to use 
//   relative_space (..., output) : target space
//   relative_component_sectors (..., output) : target sectors
//   relative_component_matrices (..., output) : target matrices
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
          relative_space,relative_component_sectors,relative_component_matrices
        );
    }
  else if (parameters.operator_name == "identity")
    {
      basis::ConstructIdentityOperatorRelativeLSJT(
          operator_parameters,
          relative_space,relative_component_sectors,relative_component_matrices
        );
    }
  else if (parameters.operator_name == "rsqr")
    {
      relative::ConstructKinematicOperator(
          operator_parameters,
          relative_space,relative_component_sectors,relative_component_matrices,
          relative::KinematicOperator::kRSqr
        );
    }
  else if (parameters.operator_name == "ksqr")
    {
      relative::ConstructKinematicOperator(
          operator_parameters,
          relative_space,relative_component_sectors,relative_component_matrices,
          relative::KinematicOperator::kKSqr
        );
    }
  else if (parameters.operator_name == "coulomb")
    {
      int num_steps = 500;  // 500 steps seems to suffice for ~8 digits precision at Nmax20
      relative::ConstructCoulombOperator(
          operator_parameters,
          relative_space,relative_component_sectors,relative_component_matrices,
          basis::TwoBodySpeciesPN::kPP,
          num_steps
        );
    }
  else if (parameters.operator_name == "symmunit")
    {
      // start from zero operator
      basis::ConstructZeroOperatorRelativeLSJT(
          operator_parameters,
          relative_space,relative_component_sectors,relative_component_matrices
        );

      // define shortcut reference to unit tensor labels
      const UnitTensorLabels& unit_tensor_labels = parameters.unit_tensor_labels;

      // look up subspace indices
      int gp = unit_tensor_labels.Lp%2;
      int relative_subspace_index_bra = relative_space.LookUpSubspaceIndex(
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
      int relative_subspace_index_ket = relative_space.LookUpSubspaceIndex(
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
      int relative_state_index_bra = relative_subspace_bra.LookUpStateIndex(
          basis::RelativeStateLSJT::StateLabelsType(
              unit_tensor_labels.Np
            )
        );
      std::cout << "Bra state:"
                << " " << unit_tensor_labels.Np
                << " -> index " << relative_state_index_bra
                << std::endl;
      int relative_state_index_ket = relative_subspace_ket.LookUpStateIndex(
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
      basis::MatrixVector& relative_matrices = relative_component_matrices[unit_tensor_labels.T0];

      // look up sector
      int relative_sector_index
        = relative_sectors.LookUpSectorIndex(
            relative_subspace_index_bra,
            relative_subspace_index_ket
          );


      // check that matrix element is canonical
      if (relative_sectors.GetSector(relative_sector_index).IsDiagonal())
        assert(relative_state_index_bra<=relative_state_index_ket);

      // set matrix element
      relative_matrices[relative_sector_index](relative_state_index_bra,relative_state_index_ket) = 1.;

    }
}

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  // parameter input
  Parameters parameters;
  ReadParameters(parameters);

  // set up operator
  basis::RelativeSpaceLSJT relative_space;
  std::array<basis::RelativeSectorsLSJT,3> relative_component_sectors;
  std::array<basis::MatrixVector,3> relative_component_matrices;
  PopulateOperator(
      parameters,
      relative_space,
      relative_component_sectors,relative_component_matrices
    );

  // write operator
  basis::WriteRelativeOperatorLSJT(
      parameters.target_filename,
      relative_space,
      parameters.operator_parameters,  // only need operator labels
      relative_component_sectors,
      relative_component_matrices,
      true  // verbose
    );

  // termination
  return 0;
}
