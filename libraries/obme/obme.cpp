/****************************************************************
  obme.cpp

  Patrick J. Fasano
  University of Notre Dame

****************************************************************/

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "eigen3/Eigen/Core"

#include "am/halfint.h"
#include "am/rme.h"
#include "basis/operator.h"
#include "obme/radial.h"

#include "obme/obme.h"

namespace shell
{
void SolidHarmonicOneBodyOperator(
    shell::RadialBasisType basis_type,
    shell::RadialOperatorType operator_type,
    int order,
    const basis::OrbitalSpaceLJPN& space,
    const basis::OrbitalSectorsLJPN& sectors,
    basis::OperatorBlocks<double>& matrices)
{
  GenerateRadialOperator(basis_type, operator_type, order, space, sectors, matrices);

  for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
  {
    // get sector and reference to block
    const basis::OrbitalSectorsLJPN::SectorType sector =
        sectors.GetSector(sector_index);

    double angular_factor = am::LJCoupledSphericalHarmonicCRME(
        sector.bra_subspace().l(), sector.bra_subspace().j(),
        sector.ket_subspace().l(), sector.ket_subspace().j(), sectors.j0());
    matrices[sector_index] *= angular_factor;
  }
}

void OneBodyOperatorTensorProduct(
    const basis::OrbitalSpaceLJPN& space,
    const basis::OrbitalSectorsLJPN& sectors_a,
    const basis::OperatorBlocks<double>& matrices_a,
    const basis::OrbitalSectorsLJPN& sectors_b,
    const basis::OperatorBlocks<double>& matrices_b,
    const basis::OrbitalSectorsLJPN& sectors,
    basis::OperatorBlocks<double>& matrices)
{
  // convenience variables
  const int& j0_a = sectors_a.j0();
  const int& g0_a = sectors_a.g0();
  const int& Tz0_a = sectors_a.Tz0();

  const int& j0_b = sectors_b.j0();
  const int& g0_b = sectors_b.g0();
  const int& Tz0_b = sectors_b.Tz0();

  const int& j0 = sectors.j0();
  const int& g0 = sectors.g0();
  const int& Tz0 = sectors.Tz0();

  // check that this is a valid tensor product
  assert(am::AllowedTriangle(sectors_a.j0(), sectors_b.j0(), j0));
  assert((sectors_a.g0() + sectors_b.g0() + g0) % 2 == 0);
  assert(sectors_a.Tz0() + sectors_b.Tz0() == Tz0);

  // initialize output matrices
  basis::SetOperatorToZero(sectors, matrices);

  for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
  {
    // get sector and reference to block
    const basis::OrbitalSectorsLJPN::SectorType sector =
        sectors.GetSector(sector_index);
    const int bra_subspace_index = sector.bra_subspace_index();
    const int ket_subspace_index = sector.ket_subspace_index();
    // get angular momentum quantum numbers
    const HalfInt bra_j = sector.bra_subspace().j();
    const HalfInt ket_j = sector.ket_subspace().j();

    // loop over intermediate subspaces
    for (int inner_subspace_index = 0; inner_subspace_index < space.size();
         ++inner_subspace_index)
    {
      if (sectors_a.ContainsSector(bra_subspace_index, inner_subspace_index)
          && sectors_b.ContainsSector(inner_subspace_index, ket_subspace_index))
      {
        const int sector_index_a =
            sectors_a.LookUpSectorIndex(bra_subspace_index, inner_subspace_index);
        const int sector_index_b =
            sectors_b.LookUpSectorIndex(inner_subspace_index, ket_subspace_index);
        const HalfInt inner_j = space.GetSubspace(inner_subspace_index).j();
        matrices[sector_index] +=
            am::RacahReductionFactorRose(bra_j, ket_j, inner_j, j0_a, j0_b, j0)
            * matrices_a[sector_index_a] * matrices_b[sector_index_b];
      }
    }
  }
}

}  // namespace shell
