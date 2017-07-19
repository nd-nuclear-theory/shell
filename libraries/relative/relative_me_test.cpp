/****************************************************************
  relative_me_test.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "relative/relative_me.h"

#include "cppformat/format.h"

void TestQuadrupole()
{

  std::cout << "TestQuadrupole" << std::endl;

  // set tensorial labels
  basis::OperatorLabelsJT operator_labels;
  operator_labels.J0=2;
  operator_labels.g0=0;
  operator_labels.symmetry_phase_mode=basis::SymmetryPhaseMode::kHermitian;
  operator_labels.T0_min=0;
  operator_labels.T0_max=2;

  // set basis parameters
  int Nmax=6;
  int Jmax = Nmax+1;

  // set up relative space
  basis::RelativeSpaceLSJT relative_space(Nmax,Jmax);
 
  for (const basis::TwoBodySpeciesPN& operator_species : {basis::TwoBodySpeciesPN::kPP,basis::TwoBodySpeciesPN::kNN,basis::TwoBodySpeciesPN::kPN})
    {
      // populate operator containers
      std::array<basis::RelativeSectorsLSJT,3> relative_component_sectors;
      std::array<basis::MatrixVector,3> relative_component_matrices;
      basis::ConstructZeroOperatorRelativeLSJT(
          operator_labels,
          relative_space,
          relative_component_sectors,
          relative_component_matrices
        );

      relative::ConstructQuadrupoleOperator(
          operator_labels,
          relative_space,
          relative_component_sectors,
          relative_component_matrices,
          relative::KinematicOperator::kRSqr,
          basis::TwoBodySpeciesPN::kPP
        );
      std::string filename = fmt::format("quadrupole_test_Nmax{}_{}_rel.dat",Nmax,basis::kTwoBodySpeciesPNCodeChar[int(operator_species)]);
      basis::WriteRelativeOperatorLSJT(
          filename,
          relative_space,
          operator_labels,relative_component_sectors,relative_component_matrices,
          true  // verbose
        );
    }
}


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
  basis::ConstructZeroOperatorRelativeLSJT(
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
          basis::TwoBodySpeciesPN::kPP,
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

  TestQuadrupole();
  // TestCoulombSU11();
  // TestCoulombSpline();

  // termination
  return 0;
}
