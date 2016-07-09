/****************************************************************
  moshinsky.cpp

  Perform Moshinsky transformation of general relative operator in
  LSJT scheme.

  Language: C++11

  Mark A. Caprio
  University of Notre Dame

  7/8/16 (mac): Created.

****************************************************************/

#include <iomanip>

#include "am/wigner_gsl.h"
#include "basis/lsjt_scheme.h"
#include "basis/lsjt_operator.h"
// #include "basis/indexing_jjjt.h"
#include "moshinsky/moshinsky_bracket.h"




////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

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
  WriteTwoBodyOperatorMatrices(std::cout,two_body_space,two_body_sectors,A_matrices,10,7);

  

  // termination
  return 0;
}
