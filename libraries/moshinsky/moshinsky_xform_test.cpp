/****************************************************************
  moshinsky_xform_test.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include <iomanip>

#include "moshinsky/moshinsky_xform.h"

void test_moshinsky_matrix()
{

  std::cout << "Inspecting a Moshinsky matrix" << std::endl;

  // set up sectors
  int L=2;
  int S=1;
  int J=2;
  int T=0;
  int g=0;
  int N=4;

  std::cout << "  relative-cm" << std::endl;
  basis::RelativeCMSubspaceNLSJT relative_cm_subspace(L,S,J,T,g,N);
  std::cout << relative_cm_subspace.DebugStr();

  std::cout << "  two-body" << std::endl;
  basis::TwoBodySubspaceNLSJT two_body_subspace (L,S,J,T,g,N);
  std::cout << two_body_subspace.DebugStr();

  std::cout << "  Moshinsky matrix" << std::endl;
  Eigen::MatrixXd matrix = moshinsky::MoshinskyMatrixNLSJT(relative_cm_subspace,two_body_subspace);
  std::cout << matrix << std::endl;

  std::cout << "  Orthogonality test" << std::endl;
  std::cout << matrix.transpose()*matrix << std::endl;

}

void test_transform_isoscalar()
{

  ////////////////////////////////////////////////////////////////
  // define relative identity operator
  ////////////////////////////////////////////////////////////////

  std::cout << "Setup" << std::endl;

  // set up space
  int Nmax = 2;
  int Jmax = Nmax+1;
  basis::RelativeSpaceLSJT space(Nmax,Jmax);

  // set up operator containers
  //
  // These are vectors to store information for T0=0/1/2 components.
  std::vector<basis::RelativeSectorsLSJT> component_sectors(3);
  std::vector<basis::MatrixVector> component_matrices(3);

  // populate operator containers
  int J0 = 0;
  int g0 = 0;
  for (int T0=0; T0<=2; ++T0)
    // for each isospin component
    {

      // enumerate sectors
      component_sectors[T0] = basis::RelativeSectorsLSJT(space,J0,T0,g0);
      std::cout << " T0 " << T0 << " size " << component_sectors[T0].size() << std::endl;
          
      // populate matrices
      if (T0==0)
        basis::SetOperatorToIdentity(component_sectors[T0],component_matrices[T0]);
      else
        basis::SetOperatorToZero(component_sectors[T0],component_matrices[T0]);
    }

  ////////////////////////////////////////////////////////////////
  // augment to relative-cm
  ////////////////////////////////////////////////////////////////

  // // configuration parameters
  // const int Nmax_relative = 2;
  // const int J0 = 0;
  // const int g0 = 0;
  // const int Nmax_two_body = 2;

  // // set up relative space
  // std::cout << std::endl;
  // std::cout << "Setting up relative space..." << std::endl;
  // std::cout << "Nmax " << Nmax_relative << std::endl;
  // RelativeSpaceLSJT relative_space(Nmax_relative);
  // // WriteRelativeSubspaces(std::cout,relative_space);

  // // set up relative operator
  // std::cout << std::endl;
  // std::cout << "Setting up relative sectors..." << std::endl;
  // std::cout << "J0 " << J0 << " g0 " << g0 << std::endl;
  // RelativeSectorsLSJT relative_sectors(relative_space,J0,g0);  // for scalar operator
  // SectorMatrices Ar_matrices;
  // std::cout << std::endl;

  // // input relative operator
  // std::cout << "Reading relative operator..." << std::endl;
  // // ReadRelativeOperator(std::cin,relative_space,relative_sectors,Ar_matrices);
  // SetOperatorToIdentityReduced(relative_space,relative_sectors,Ar_matrices);

  // // output relative operator
  // std::cout << std::endl;
  // std::cout << "Relative operator" << std::endl;
  // // caution: affects future output precision
  // WriteRelativeOperatorMatrices(std::cout,relative_space,relative_sectors,Ar_matrices,10,7);

  // // set up two-body space
  // std::cout << std::endl;
  // std::cout << "Setting up two-body space..." << std::endl;
  // TwoBodySpaceLSJT two_body_space(Nmax_two_body);
  // std::cout << "Nmax " << Nmax_two_body << std::endl;
  // // WriteRelativeSubspaces(std::cout,relative_space);

  // // construct two-body operator
  // std::cout << std::endl;
  // std::cout << "Transforming..." << std::endl;
  // TwoBodySectorsLSJT two_body_sectors(two_body_space,J0,g0);  // for scalar operator
  // SectorMatrices A_matrices;
  // for (int two_body_sector_index = 0; two_body_sector_index < two_body_sectors.size(); ++two_body_sector_index)
  //   {
  //     Sector two_body_sector = two_body_sectors.GetSector(two_body_sector_index);
  //     const TwoBodySubspaceLSJT& two_body_subspace2 = two_body_space.GetSubspace(two_body_sector.index2());
  //     const TwoBodySubspaceLSJT& two_body_subspace1 = two_body_space.GetSubspace(two_body_sector.index1());
  //     Eigen::MatrixXd two_body_sector_matrix = TransformedSector(
  //         two_body_subspace2,two_body_subspace1,
  //         relative_space,relative_sectors,Ar_matrices,
  //         J0 // for scalar operator
  //       );

  //     A_matrices.push_back(two_body_sector_matrix);
  //   }

  // // dump two_body operator
  // std::cout << std::endl;
  // std::cout << "Two-body operator" << std::endl;
  // WriteTwoBodyOperator(std::cout,two_body_space,two_body_sectors,A_matrices);
  // WriteTwoBodyOperatorMatrices(std::cout,two_body_space,two_body_sectors,A_matrices,10,7);

  
}

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  test_moshinsky_matrix();
  // test_transform_isoscalar();

  // termination
  return 0;
}
