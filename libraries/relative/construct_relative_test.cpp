/****************************************************************
  construct_relative_test.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "relative/construct_relative.h"

#include "cppformat/format.h"


void TestCoulombSpline()
{

  std::cout << "TestCoulombSpline" << std::endl;

  // set tensorial labels
  basis::OperatorLabelsJT operator_labels;
  operator_labels.J0=0;
  operator_labels.g0=0;
  operator_labels.symmetry_phase_mode=basis::SymmetryPhaseMode::kHermitian;
  operator_labels.T0_min=0;
  operator_labels.T0_max=2;

  // set basis parameters
  int Nmax=20;
  int Jmax = Nmax+1;

  // set up relative space
  basis::RelativeSpaceLSJT relative_space(Nmax,Jmax);

  // populate operator containers
  std::array<basis::RelativeSectorsLSJT,3> relative_component_sectors;
  std::array<basis::MatrixVector,3> relative_component_matrices;
  basis::ConstructIdentityOperatorRelativeLSJT(
      operator_labels,
      relative_space,
      relative_component_sectors,
      relative_component_matrices
    );

  std::vector<int> num_steps_list({50,100,200,500,1000});
  for (int num_steps :  num_steps_list)
    {
      relative::ConstructCoulombOperator(
          operator_labels,
          relative_space,
          relative_component_sectors,
          relative_component_matrices,
          num_steps
        );
  
      std::string filename = fmt::format("coulomb_test_Nmax{}_steps{}.dat",Nmax,num_steps);
      basis::WriteRelativeOperatorLSJT(
          filename,
          relative_space,
          operator_labels,relative_component_sectors,relative_component_matrices,
          true  // verbose
        );
    }
}

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  // TestCoulombSU11();
  TestCoulombSpline();

  // termination
  return 0;
}
