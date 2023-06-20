/****************************************************************

  relative_xform.cpp

  Patrick J. Fasano, University of Notre Dame.
  SPDX-License-Identifier: MIT

****************************************************************/

#include "Daejeon16_wrapper.h"

extern "C" {
    double v_daejeon16_(
        const int* nrel, const int* nprel,
        const int* lrel, const int* lprel,
        const int* JJ, const int* IS,
        const int* keyhw, const int* keyph
    );
}

// wrapper to convert references to pointers
double v_daejeon16(
    const int& nrel, const int& nprel,
    const int& lrel, const int& lprel,
    const int& JJ, const int& IS,
    const int& keyhw, const int& keyph
  )
{ return v_daejeon16_(&nrel, &nprel, &lrel, &lprel, &JJ, &IS, &keyhw, &keyph); }

namespace contrib {
  void ConstructDaejeon16Operator(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_blocks
    )
  {
    // validate operator parameters
    assert(operator_labels.J0==0);
    assert(operator_labels.g0==0);
    assert((operator_labels.T0_min==0)&&(operator_labels.T0_max>=0));

    // zero initialize operator
    basis::ConstructZeroOperatorRelativeLSJT(
        operator_labels,relative_space,relative_component_sectors,relative_component_blocks
      );

    // select T0=0 component
    const basis::RelativeSectorsLSJT& sectors = relative_component_sectors[/*T0=*/0];
    basis::OperatorBlocks<double>& matrices = relative_component_blocks[/*T0=*/0];

    // iterate over sectors
    for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
    {
      // extract sector and subspaces
      const basis::RelativeSectorsLSJT::SectorType& sector = sectors.GetSector(sector_index);
      const basis::RelativeSubspaceLSJT& bra_subspace = sector.bra_subspace();
      const basis::RelativeSubspaceLSJT& ket_subspace = sector.ket_subspace();

      // alias to matrix
      basis::OperatorBlock<double>& matrix = matrices[sector_index];

      // extract labels
      const int bra_Nmax = bra_subspace.Nmax();
      const int ket_Nmax = ket_subspace.Nmax();
      const int bra_L = bra_subspace.L();
      const int ket_L = ket_subspace.L();
      const int bra_S = bra_subspace.S();
      const int ket_S = ket_subspace.S();
      const int bra_J = bra_subspace.J();
      const int ket_J = ket_subspace.J();
      const int bra_T = bra_subspace.T();
      const int ket_T = ket_subspace.T();
      const int bra_nmax = (bra_Nmax-bra_L)/2;   // max radial quantum number
      const int ket_nmax = (ket_Nmax-ket_L)/2;

      // Daejeon16 doesn't flip spins (?) -- must be related to being isoscalar
      if (bra_S != ket_S)
        continue;

      // populate matrix
      for (int ket_n=0; ket_n<= ket_nmax; ++ket_n)
      {
        for (int bra_n=0; bra_n<= bra_nmax; ++bra_n)
        {
          constexpr int keyhw = 1;  // interaction in MeV
          constexpr int keyph = 0;  // positive at origin convention for H.O. functions
          matrix(bra_n,ket_n) = v_daejeon16(
              bra_n, ket_n, bra_L, ket_L, bra_J, bra_S, keyhw, keyph
            );
        }
      }
    }
  }
}  // namespace contrib

