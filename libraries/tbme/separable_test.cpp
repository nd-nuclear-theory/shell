/****************************************************************
  separable_test.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "tbme/separable.h"

#include "cppformat/format.h"

void TestIdentity()
{

  std::cout << "TestIdentity" << std::endl;

  // source space
  const int Nmax = 2;
  const basis::OrbitalSpacePN orbital_space(Nmax);
  const basis::TwoBodySpaceJJJPN space(orbital_space,basis::WeightMax(basis::Rank::kTwoBody,Nmax));
  const basis::TwoBodySectorsJJJPN sectors(space,0,0,0);

  for (int sector_index=0; sector_index<sectors.size(); ++sector_index)
    {

      // extract sector
      const typename basis::TwoBodySectorsJJJPN::SectorType& sector = sectors.GetSector(sector_index);

      // populate sector
      Eigen::MatrixXd matrix = shell::IdentityOperatorSectorMatrixJJJPN(sector);

      std::cout
        << fmt::format("Sector {}: {}",sector_index,sector.ket_subspace().LabelStr()) << std::endl
        << sector.ket_subspace().DebugStr()
        << std::endl
        << matrix << std::endl
        << std::endl;
        
    }
}

void TestAngularMomentum()
{

  std::cout << "TestAngularMomentum" << std::endl;

  // source space
  const int Nmax = 2;
  const basis::OrbitalSpacePN orbital_space(Nmax);
  const basis::TwoBodySpaceJJJPN space(orbital_space,basis::WeightMax(basis::Rank::kTwoBody,Nmax));
  const basis::TwoBodySectorsJJJPN sectors(space,0,0,0);

  // operator choice
  shell::AngularMomentumOperatorFamily operator_family = shell::AngularMomentumOperatorFamily::kTotal;
  shell::AngularMomentumOperatorSpecies operator_species = shell::AngularMomentumOperatorSpecies::kTotal;
  int A = 2;


  for (int sector_index=0; sector_index<sectors.size(); ++sector_index)
    {

      // extract sector
      const typename basis::TwoBodySectorsJJJPN::SectorType& sector = sectors.GetSector(sector_index);

      // populate sector
      Eigen::MatrixXd matrix = shell::AngularMomentumSectorMatrixJJJPN(operator_family,operator_species,A,sector);

      std::cout
        << fmt::format("Sector {}: {}",sector_index,sector.ket_subspace().LabelStr()) << std::endl
        << sector.ket_subspace().DebugStr()
        << std::endl
        << matrix << std::endl
        << std::endl;
        
    }
}


////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  TestIdentity();
  TestAngularMomentum();

  // termination
  return 0;
}
