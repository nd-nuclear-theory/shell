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
#include "analytic/radial_oscillator_me.h"
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
    const basis::OrbitalSectorsLJPN::SectorType& sector =
        sectors.GetSector(sector_index);

    double angular_factor = am::LJCoupledSphericalHarmonicCRME(
        sector.bra_subspace().l(), sector.bra_subspace().j(),
        sector.ket_subspace().l(), sector.ket_subspace().j(), sectors.j0());
    matrices[sector_index] *= angular_factor;
  }
}

void LadderOneBodyOperator(
    shell::RadialBasisType basis_type,
    shell::LadderOperatorType operator_type,
    const basis::OrbitalSpaceLJPN& space,
    const basis::OrbitalSectorsLJPN& sectors,
    basis::OperatorBlocks<double>& matrices
  )
{
  // convenience variables
  const int j0 = sectors.j0();
  const int g0 = sectors.g0();
  const int Tz0 = sectors.Tz0();

  // initialize output matrices
  basis::SetOperatorToZero(sectors, matrices);

  // TODO(pjf): generalize for non-oscillator bases
  assert(basis_type==shell::RadialBasisType::kOscillator);

  // loop over sectors
  for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
  {
    // get sector and reference to block
    const basis::OrbitalSectorsLJPN::SectorType sector =
        sectors.GetSector(sector_index);
    basis::OperatorBlock<double>& sector_matrix = matrices[sector_index];

    // get sizes
    const int bra_subspace_size = sector.bra_subspace().size();
    const int ket_subspace_size = sector.ket_subspace().size();

    // angular momentum factor
    double angular_factor = am::LJCoupledSphericalHarmonicCRME(
        sector.bra_subspace().l(), sector.bra_subspace().j(),
        sector.ket_subspace().l(), sector.ket_subspace().j(), sectors.j0());

    // main loop
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < bra_subspace_size; ++j)
    {
      for (int k = 0; k < ket_subspace_size; ++k)
      {
        // get states
        basis::OrbitalStateLJPN bra_state(sector.bra_subspace(), j);
        basis::OrbitalStateLJPN ket_state(sector.ket_subspace(), k);

        // convenience variables
        const int bra_n = bra_state.n();
        const int bra_l = bra_state.l();
        const int bra_N = 2*bra_n + bra_l;
        const int ket_n = ket_state.n();
        const int ket_l = ket_state.l();
        const int ket_N = 2*ket_n + ket_l;

        double matrix_element = 0;
        if (operator_type == shell::LadderOperatorType::kRaising)
        {
          matrix_element = analytic::CDaggerOscillatorMatrixElement(
              bra_N, bra_l, ket_N, ket_l
            );
        } else if (operator_type == shell::LadderOperatorType::kLowering) {
          matrix_element = analytic::COscillatorMatrixElement(
              bra_N, bra_l, ket_N, ket_l
            );
        } else { assert(0); }
        sector_matrix(j, k) = angular_factor * matrix_element;
      }
    }
  }

}

void AngularMomentumOneBodyOperator(
    am::AngularMomentumOperatorType operator_type,
    const basis::OrbitalSpaceLJPN& space,
    const basis::OrbitalSectorsLJPN& sectors,
    basis::OperatorBlocks<double>& matrices
  )
{
  // must have valid quantum numbers
  assert(sectors.j0() == 1);
  assert(sectors.g0() == 0);
  assert(sectors.Tz0() == 0);
  // initialize output matrices
  basis::SetOperatorToZero(sectors, matrices);

  for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
  {
    double factor;
    // get sector and reference to block
    const basis::OrbitalSectorsLJPN::SectorType& sector =
        sectors.GetSector(sector_index);
        // extract sector subspaces
    const auto& bra_subspace = sector.bra_subspace();
    const auto& ket_subspace = sector.ket_subspace();

    if (operator_type == am::AngularMomentumOperatorType::kOrbital)
      {
        factor = am::jjJCoupledAngularMomentumJ1RME(
            bra_subspace.l(), HalfInt(1,2), bra_subspace.j(),
            ket_subspace.l(), HalfInt(1,2), ket_subspace.j()
          );
      }
    else if (operator_type == am::AngularMomentumOperatorType::kSpin)
      {
        factor = am::jjJCoupledAngularMomentumJ2RME(
            bra_subspace.l(), HalfInt(1,2), bra_subspace.j(),
            ket_subspace.l(), HalfInt(1,2), ket_subspace.j()
          );
      }
    else // if (operator_type == am::AngularMomentumOperatorType::kTotal)
      {
        factor = am::jjJCoupledAngularMomentumJRME(
            bra_subspace.l(), HalfInt(1,2), bra_subspace.j(),
            ket_subspace.l(), HalfInt(1,2), ket_subspace.j()
          );
      }

    matrices[sector_index] = factor *
      basis::OperatorBlock<double>::Identity(
          bra_subspace.size(), ket_subspace.size()
        );
  }
}

void AngularMomentumSquaredOneBodyOperator(
    am::AngularMomentumOperatorType operator_type,
    const basis::OrbitalSpaceLJPN& space,
    const basis::OrbitalSectorsLJPN& sectors,
    basis::OperatorBlocks<double>& matrices
  )
{
  assert(sectors.j0() == 0);
  assert(sectors.g0() == 0);
  assert(sectors.Tz0() == 0);
  // initialize output matrices
  basis::SetOperatorToIdentity(sectors, matrices);

  for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
  {
    double factor;
    // get sector and reference to block
    const basis::OrbitalSectorsLJPN::SectorType& sector =
        sectors.GetSector(sector_index);
    const auto& bra_subspace = sector.bra_subspace();
    const auto& ket_subspace = sector.ket_subspace();

    assert(sector.IsDiagonal());
    if (operator_type == am::AngularMomentumOperatorType::kOrbital)
      {
        factor = bra_subspace.l()*(bra_subspace.l()+1);
      }
    else if (operator_type == am::AngularMomentumOperatorType::kSpin)
      {
        factor = 3./4.;
      }
    else // if (operator_type == am::AngularMomentumOperatorType::kTotal)
      {
        factor = double(bra_subspace.j())*double(bra_subspace.j()+1);
      }
    matrices[sector_index] *= factor;
  }
}

void IsospinOneBodyOperator(
    const basis::OrbitalSpaceLJPN& space,
    const basis::OrbitalSectorsLJPN& sectors,
    basis::OperatorBlocks<double>& matrices
  )
{
  // convenience variables
  const int& j0 = sectors.j0();
  const int& g0 = sectors.g0();
  const int& Tz0 = sectors.Tz0();

  // check for valid tensor character
  assert(j0==0);
  assert(g0==0);

  // initialize output matrices
  matrices.resize(sectors.size());

  for (int sector_index=0; sector_index < sectors.size(); ++sector_index)
  {
    const auto& sector = sectors.GetSector(sector_index);

    double matrix_element = ((Tz0==0) ? double(sector.ket_subspace().Tz()) : 1.);
    matrices[sector_index] = matrix_element
      * basis::OperatorBlock<double>::Identity(
            sector.bra_subspace().size(), sector.ket_subspace().size()
          );
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
    const basis::OrbitalSectorsLJPN::SectorType& sector =
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
