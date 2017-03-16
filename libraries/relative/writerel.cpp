/****************************************************************
  writerel.cpp

  Write simple relative operator files.

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
      Relative k^2 operator (~relative kinetic energy).

    rsqr
      Relative r^2 operator (~relative kinetic energy).

    coulomb
      Coulomb potential (WIP!!!).

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

  7/16/16 (mac): Created.
  7/25/16 (mac): Update to use WriteRelativeOperatorLSJT.
  10/9/16 (pjf): Rename mcpp -> mcutils.
  3/6/17 (mac): Rough in coulomb interaction code.

****************************************************************/

#include "basis/lsjt_operator.h"
#include "mcutils/parsing.h"
#include "relative/construct_relative.h"

#include "cppformat/format.h"
#include "spline/wavefunction_class.h"


namespace relative {

  void ConstructCoulombOperator(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::MatrixVector,3>& relative_component_matrices
    )
  // TODO implement isospin stuff, move out to construct_relative.cpp
  {

    // zero initialize operator
    ConstructDiagonalConstantOperator(
        operator_labels,relative_space,relative_component_sectors,relative_component_matrices,0.
      );

    for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
      // for each isospin component
      {
        // select T0=0 component
        const basis::RelativeSectorsLSJT& sectors = relative_component_sectors[0];
        basis::MatrixVector& matrices = relative_component_matrices[0];

        // iterate over sectors
        for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
          {

            // extract sector
            const basis::RelativeSectorsLSJT::SectorType& sector = sectors.GetSector(sector_index);


            // short-circuit select diagonal sectors
            //
            // The a.m. scalar nature imposes diagonal on (L,S,J),
            // positive parity imposes diagonal on (g), and then there
            // is that apparent parity selection imposing diagonal in
            // T.
            if (!sector.IsDiagonal())
              continue;

            // short circuit a.m. scalar
            // book am_allowed = true;
            // am_allowed &= (bra_subspace.L()  == ket_subspace.L());
            // am_allowed &= (bra_subspace.S()  == ket_subspace.S());
            // am_allowed &= (bra_subspace.J()  == ket_subspace.J());


            // extract subspace and labels
            const basis::RelativeSubspaceLSJT& subspace = sector.ket_subspace();
            int Nmax = subspace.Nmax();
            int L = subspace.L();
            int S = subspace.S();
            int J = subspace.J();
            int T = subspace.T();
            int nmax = (Nmax-L)/2;   // max radial quantum number

            // Although actually we could get nmax just from the
            // subspace dimension...
            assert(nmax==subspace.size()-1);

            // populate nonzero entries
            //
            // We make use of the known indexing scheme for a
            // RelativeLSJT basis, that the radial quantum number n is
            // just the 0-based state index.

            Eigen::MatrixXd& matrix = matrices[sector_index];
            for (int bra_n=0; bra_n<=nmax; ++bra_n)
              for (int ket_n=0; ket_n<=nmax; ++ket_n)
                {

                  // get bra and ket states
                  spline::WaveFunction bra_wavefunction(bra_n,L,1,spline::Basis::HC);
                  spline::WaveFunction ket_wavefunction(ket_n,L,1,spline::Basis::HC);
                  std::cout << fmt::format("bra_n {} ket_n {} L {}") << std::endl;
                  
                  // evaluate radial integral
                  const int num_integration_steps_or_maybe_points_or_whatever = 1000;
                  double radial_integral = bra_wavefunction.MatrixElement(num_integration_steps_or_maybe_points_or_whatever,ket_wavefunction,-1);

                  // impose isospin factors -- TODO
                  matrix(bra_n,ket_n) = radial_integral;
                }
          }
      }

  }
}

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
  else if (parameters.operator_name == "coulomb")
    {
      relative::ConstructCoulombOperator(
          operator_parameters,
          relative_space,relative_component_sectors,relative_component_matrices
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
