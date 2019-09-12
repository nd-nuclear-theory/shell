/****************************************************************
  radial.cpp

  Patrick J. Fasano
  University of Notre Dame

****************************************************************/

#include "am/am.h"
#include "am/halfint.h"
#include "analytic/radial_oscillator_me.h"
#include "basis/operator.h"
#include "spline/spline_me.h"

#include "obme/radial.h"

namespace shell
{
void GenerateRadialOperator(
    shell::RadialBasisType basis_type,
    shell::RadialOperatorType operator_type,
    int order,
    const basis::OrbitalSpaceLJPN& space,
    const basis::OrbitalSectorsLJPN& sectors,
    basis::OperatorBlocks<double>& matrices)
{
  // convenience variables
  const int j0 = sectors.j0();
  const int g0 = sectors.g0();
  const int Tz0 = sectors.Tz0();

  // initialize output matrices
  basis::SetOperatorToZero(sectors, matrices);

  // loop over sectors
  for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
  {
    // get sector and reference to block
    const basis::OrbitalSectorsLJPN::SectorType sector =
        sectors.GetSector(sector_index);
    basis::OperatorBlock<double>& sector_matrix = matrices[sector_index];

    // get sizes
    const std::size_t bra_subspace_size = sector.bra_subspace().size();
    const std::size_t ket_subspace_size = sector.ket_subspace().size();

    // main loop
    #pragma omp parallel for collapse(2)
    for (std::size_t j = 0; j < bra_subspace_size; ++j)
    {
      for (std::size_t k = 0; k < ket_subspace_size; ++k)
      {
        // get states
        basis::OrbitalStateLJPN bra_state(sector.bra_subspace(), j);
        basis::OrbitalStateLJPN ket_state(sector.ket_subspace(), k);

        // convenience variables
        const int bra_n = bra_state.n();
        const int bra_l = bra_state.l();
        const int ket_n = ket_state.n();
        const int ket_l = ket_state.l();

        double matrix_element = 0.;
        // if applicable, use analytic formulae
        if (order == 0)
        {
          if ((bra_n == ket_n) && (bra_l == ket_l))
          {
            matrix_element = 1.;
          }
          else
          {
            matrix_element = 0.;
          }
        }
        else if (basis_type == shell::RadialBasisType::kOscillator)
        {
          const int bra_N = 2 * bra_n + bra_l;
          const int ket_N = 2 * ket_n + ket_l;
          const int operator_sign = (operator_type == shell::RadialOperatorType::kK) ? -1 : +1;
          if ((order == 1) && (j0 == 1) && (g0 == 1))
          {
            matrix_element = analytic::CoordinateOscillatorMatrixElement(
                bra_N, bra_l, ket_N, ket_l, operator_sign);
          }
          else if (order == 2 && ((j0 == 0) || (j0 == 2)) && (g0 == 0))
          {
            matrix_element = analytic::CoordinateSqrOscillatorMatrixElement(
                bra_N, bra_l, ket_N, ket_l, operator_sign);
          }
          else
          {
            matrix_element = spline::RadialMatrixElement(
                bra_n, bra_l, 1., spline::BasisType::kOscillator,
                ket_n, ket_l, 1., spline::BasisType::kOscillator,
                static_cast<spline::OperatorType>(operator_type), order);
          }
        }
        else if (basis_type == shell::RadialBasisType::kLaguerre)
        {
          matrix_element = spline::RadialMatrixElement(
              bra_n, bra_l, 1., spline::BasisType::kLaguerre,
              ket_n, ket_l, 1., spline::BasisType::kLaguerre,
              static_cast<spline::OperatorType>(operator_type), order);
        }

        sector_matrix(j, k) = matrix_element;
      }
    }
  }
}

void GenerateRadialOverlaps(
    shell::RadialBasisType bra_basis_type,
    shell::RadialBasisType ket_basis_type,
    double scale_factor,
    const basis::OrbitalSpaceLJPN& bra_space,
    const basis::OrbitalSpaceLJPN& ket_space,
    const basis::OrbitalSectorsLJPN& sectors,
    basis::OperatorBlocks<double>& matrices)
{
  assert(sectors.j0() == 0);
  assert(sectors.g0() == 0);
  assert(sectors.Tz0() == 0);

  // initialize output matrices
  basis::SetOperatorToZero(sectors, matrices);

  // loop over sectors
  for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
  {
    // get sector and reference to block
    const basis::OrbitalSectorsLJPN::SectorType sector =
        sectors.GetSector(sector_index);
    basis::OperatorBlock<double>& sector_matrix = matrices[sector_index];

    // get sizes
    const std::size_t bra_subspace_size = sector.bra_subspace().size();
    const std::size_t ket_subspace_size = sector.ket_subspace().size();

    // main loop
    #pragma omp parallel for collapse(2)
    for (std::size_t j = 0; j < bra_subspace_size; ++j)
    {
      for (std::size_t k = 0; k < ket_subspace_size; ++k)
      {
        // get states
        basis::OrbitalStateLJPN bra_state(sector.bra_subspace(), j);
        basis::OrbitalStateLJPN ket_state(sector.ket_subspace(), k);

        // convenience variables
        const int bra_n = bra_state.n();
        const int bra_l = bra_state.l();
        const int ket_n = ket_state.n();
        const int ket_l = ket_state.l();

        double matrix_element = 0.;
        if (scale_factor == 1. && bra_basis_type == ket_basis_type)
        {
          if ((bra_n == ket_n) && (bra_l == ket_l))
          {
            matrix_element = 1.;
          }
          else
          {
            matrix_element = 0.;
          }
        }
        else
        {
          matrix_element = spline::RadialMatrixElement(
              bra_n, bra_l, 1., static_cast<spline::BasisType>(bra_basis_type),
              ket_n, ket_l, scale_factor,
              static_cast<spline::BasisType>(ket_basis_type),
              spline::OperatorType::kR, 0);
        }

        sector_matrix(j, k) = matrix_element;
      }
    }
  }
}

void ComposeRadialOperators(
    const basis::OrbitalSpaceLJPN& bra_space_a,
    const basis::OrbitalSpaceLJPN& ket_space_a,
    const basis::OrbitalSectorsLJPN& sectors_a,
    const basis::OperatorBlocks<double>& matrices_a,
    const basis::OrbitalSpaceLJPN& bra_space_b,
    const basis::OrbitalSpaceLJPN& ket_space_b,
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

  // check that spaces are valid
  assert(ket_space_a.OrbitalInfo() == bra_space_b.OrbitalInfo());
  const auto& inner_space = ket_space_a;

  // initialize output matrices
  basis::SetOperatorToZero(sectors, matrices);

  for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
  {
    // get sector and reference to block
    const basis::OrbitalSectorsLJPN::SectorType sector =
        sectors.GetSector(sector_index);
    const std::size_t bra_subspace_index = sector.bra_subspace_index();
    const std::size_t ket_subspace_index = sector.ket_subspace_index();

    // loop over intermediate subspaces
    for (std::size_t inner_subspace_index = 0; inner_subspace_index < inner_space.size();
         ++inner_subspace_index)
    {
      if (sectors_a.ContainsSector(bra_subspace_index, inner_subspace_index)
          && sectors_b.ContainsSector(inner_subspace_index, ket_subspace_index))
      {
        const std::size_t sector_index_a =
            sectors_a.LookUpSectorIndex(bra_subspace_index, inner_subspace_index);
        const std::size_t sector_index_b =
            sectors_b.LookUpSectorIndex(inner_subspace_index, ket_subspace_index);
        matrices[sector_index] +=
            matrices_a[sector_index_a] * matrices_b[sector_index_b];
      }
    }
  }
}

void SimilarityTransformOperator(
    const basis::OrbitalSpaceLJPN& source_space,
    const basis::OrbitalSpaceLJPN& target_space,
    const basis::OrbitalSectorsLJPN& xform_sectors,
    const basis::OperatorBlocks<double>& xform_matrices,
    const basis::OrbitalSectorsLJPN& operator_sectors,
    const basis::OperatorBlocks<double>& operator_matrices,
    const basis::OrbitalSectorsLJPN& output_sectors,
    basis::OperatorBlocks<double>& output_matrices)
{
  // ensure that sectors for xform are valid for an xform
  assert(xform_sectors.j0() == 0);
  assert(xform_sectors.g0() == 0);
  assert(xform_sectors.Tz0() == 0);

  // initialize output matrices
  basis::SetOperatorToZero(output_sectors, output_matrices);

  for (std::size_t sector_index = 0; sector_index < output_sectors.size(); ++sector_index)
  {
    // get sector and reference to block
    const auto& sector = output_sectors.GetSector(sector_index);
    const std::size_t target_bra_subspace_index = sector.bra_subspace_index();
    const auto& target_bra_subspace_labels = sector.bra_subspace().labels();
    const std::size_t target_ket_subspace_index = sector.ket_subspace_index();
    const auto& target_ket_subspace_labels = sector.ket_subspace().labels();

    // find the bra and ket subspace indices in the source space
    // if either doesn't exist, this sector will be left/padded with zeros
    std::size_t source_bra_subspace_index =
        source_space.LookUpSubspaceIndex(target_bra_subspace_labels);
    if (source_bra_subspace_index == basis::kNone) { continue; }
    std::size_t source_ket_subspace_index =
        source_space.LookUpSubspaceIndex(target_ket_subspace_labels);
    if (source_ket_subspace_index == basis::kNone) { continue; }

    // find the sector index in the input operator sectors
    std::size_t operator_sector_index = operator_sectors.LookUpSectorIndex(
        source_bra_subspace_index, source_ket_subspace_index);
    assert(operator_sector_index != basis::kNone);

    // find xform sectors for transforming the bra and ket subspaces
    std::size_t left_xform_sector_index = xform_sectors.LookUpSectorIndex(
        source_bra_subspace_index, target_bra_subspace_index);
    assert(left_xform_sector_index != basis::kNone);
    std::size_t right_xform_sector_index = xform_sectors.LookUpSectorIndex(
        source_ket_subspace_index, target_ket_subspace_index);
    assert(right_xform_sector_index != basis::kNone);

    // transform sector matrix using U^T O U, and store in output matrix
    output_matrices[sector_index] =
        xform_matrices[left_xform_sector_index].transpose()
        * operator_matrices[operator_sector_index]
        * xform_matrices[right_xform_sector_index];
  }
}

}  // namespace shell
