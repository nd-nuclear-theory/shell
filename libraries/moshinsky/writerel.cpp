/****************************************************************
  writerel.cpp

  Write simple relative operator files.

  See lsjt_operator.h for documentation of operator storage and the
  relative operator file format.

  Standard input:
    J0 g0 T0_min T0_max
    Nmax
    operator_name [parameters]
    filename

  The parameters are

    J0 g0 T0_min T0_max : relative operator tensor properties

    Nmax : relative basis trunctation

  The operator_name line may be:

    zero
      Zero operator.

    identity
      Identity operator.

    ksqr
      Relative k^2 operator (~relative kinetic energy).

    rsqr
      Relative r^2 operator (~relative kinetic energy).

    symmunit T0 Np Lp Sp Jp Tp N L S J T
      Symmetrized unit tensor.  The labels must specify a canonical
      (upper triangle) matrix element.

  Other operator parameter values are taken as:
    symmetry_phase_mode = kHermitian
    Jmax = Nmax+1

  Language: C++11
                                 
  Mark A. Caprio
  University of Notre Dame

  7/16/16 (mac): Created.

****************************************************************/

#include <fstream>
#include <iomanip>

#include "basis/lsjt_operator.h"
#include "mcpp/parsing.h"
#include "mcpp/profiling.h"
#include "moshinsky/construct_relative.h"

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
  std::string filename;
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
  parameters.operator_parameters.version = 1;
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
    line_stream >> parameters.filename;
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
      relative::ConstructDiagonalConstantOperator(
          operator_parameters,
          relative_space,relative_component_sectors,relative_component_matrices,
          0.
        );
    }
  else if (parameters.operator_name == "identity")
    {
      relative::ConstructDiagonalConstantOperator(
          operator_parameters,
          relative_space,relative_component_sectors,relative_component_matrices,
          1.
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
  else if (parameters.operator_name == "symmunit")
    {
      // start from zero operator
      relative::ConstructDiagonalConstantOperator(
          operator_parameters,
          relative_space,relative_component_sectors,relative_component_matrices,
          0.
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

void WriteOperator(
    const Parameters& parameters,
    const basis::RelativeSpaceLSJT& relative_space,
    const std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
    const std::array<basis::MatrixVector,3>& relative_component_matrices
  )
// Write operator to output file.
{

  // define shortcut reference to operator parameters
  const basis::RelativeOperatorParametersLSJT& operator_parameters
    = parameters.operator_parameters;

  // set up stream for output
  std::ofstream os(parameters.filename.c_str());

  // write header parameters
  basis::WriteRelativeOperatorParametersLSJT(os,operator_parameters);

  // write matrices
  for (int T0=operator_parameters.T0_min; T0<=operator_parameters.T0_max; ++T0)
    {
      basis::WriteRelativeOperatorComponentLSJT(
          os,
          T0,
          relative_component_sectors[T0],relative_component_matrices[T0]
        );
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
  WriteOperator(
      parameters,
      relative_space,
      relative_component_sectors,relative_component_matrices
    );


  // termination
  return 0;
}
