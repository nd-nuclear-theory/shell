/****************************************************************
  separable_test.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "mcutils/eigen.h"
#include "obme/obme.h"
#include "tbme/tbme_separable.h"

#include "fmt/format.h"

void TestIdentity()
{

  std::cout << "TestIdentity" << std::endl;

  // source space
  const int Nmax = 2;
  const basis::OrbitalSpacePN orbital_space(Nmax);
  const basis::TwoBodySpaceJJJPN space(orbital_space,basis::WeightMax(basis::Rank::kTwoBody,Nmax));
  const basis::TwoBodySectorsJJJPN sectors(space,0,0,0);

  for (std::size_t sector_index=0; sector_index<sectors.size(); ++sector_index)
    {

      // extract sector
      const typename basis::TwoBodySectorsJJJPN::SectorType& sector = sectors.GetSector(sector_index);

      // populate sector
      int A = 2;
      Eigen::MatrixXd matrix = shell::IdentityOperatorMatrixJJJPN(sector,A);

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
  am::AngularMomentumOperatorType operator_type = am::AngularMomentumOperatorType::kTotal;
  int A = 2;

  // set up one-body operators
  const basis::OrbitalSpaceLJPN orbital_space_ljpn(orbital_space);
  const basis::OrbitalSectorsLJPN am_sectors(orbital_space_ljpn,1,0,0);
  const basis::OrbitalSectorsLJPN am_squared_sectors(orbital_space_ljpn,0,0,0);
  basis::OperatorBlocks<double> am_blocks, am_squared_blocks;
  shell::AngularMomentumOneBodyOperator(
      operator_type,
      orbital_space_ljpn,
      am_sectors,
      am_blocks
    );
  shell::AngularMomentumSquaredOneBodyOperator(
      operator_type,
      orbital_space_ljpn,
      am_squared_sectors,
      am_squared_blocks
    );


  for (std::size_t sector_index=0; sector_index<sectors.size(); ++sector_index)
    {

      // extract sector
      const typename basis::TwoBodySectorsJJJPN::SectorType& sector = sectors.GetSector(sector_index);
      assert(sector.IsDiagonal());
      const basis::TwoBodySubspaceJJJPN& subspace = sector.ket_subspace();

      // populate sector
      Eigen::MatrixXd matrix =
        shell::UpgradeOneBodyOperatorJJJPN(
            orbital_space_ljpn,
            am_squared_sectors, am_squared_blocks,
            sector, A
          )
        + 2*(-std::sqrt(3)) * shell::RacahReduceTensorProductJJJPN(
            orbital_space_ljpn,
            am_sectors, am_blocks, am_sectors, am_blocks,
            sector, /*J0=*/0
          );

      int J = subspace.J();
      std::cout
        << fmt::format("Sector {}: {}",sector_index,subspace.LabelStr()) << " expect " << J*(J+1) << std::endl
        << sector.ket_subspace().DebugStr()
        << std::endl
        << mcutils::ChopMatrix(matrix, 1e-12) << std::endl
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
