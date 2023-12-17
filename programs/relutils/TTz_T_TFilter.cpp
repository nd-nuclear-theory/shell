#include <fstream>

#include "basis/jjjt_operator.h"
#include "basis/jjjttz_operator.h"
#include "mcutils/eigen.h"
#include "mcutils/parsing.h"
#include "mcutils/profiling.h"
#include "tbme/me2j_io.h"


void TransformOperatorTwoBodyJJJTToTwoBodyJJJTTz( // assumptions: J, T, g, Tz are the same between bra and ket for a given matrix element
    const basis::TwoBodySpaceJJJT& two_body_jjjt_space,
    const std::array<basis::TwoBodySectorsJJJT,3>& two_body_jjjt_component_sectors,
    const std::array<basis::OperatorBlocks<double>,3>& two_body_jjjt_component_matrices,
    const basis::TwoBodySpaceJJJTTz& two_body_jjjttz_space,
    basis::TwoBodySectorsJJJTTz& two_body_jjjttz_sectors,
    basis::OperatorBlocks<double>& two_body_jjjttz_matrices
  ) {
    // enumerate target sectors
    two_body_jjjttz_sectors = basis::TwoBodySectorsJJJTTz(two_body_jjjttz_space,0,0,0);

    // populate matrices
    two_body_jjjttz_matrices.resize(two_body_jjjttz_sectors.size());
    for (std::size_t two_body_jjjttz_sector_index=0; two_body_jjjttz_sector_index<two_body_jjjttz_sectors.size(); two_body_jjjttz_sector_index++)
      {
        // make reference to target sector
        const basis::TwoBodySectorsJJJTTz::SectorType& two_body_jjjttz_sector
          = two_body_jjjttz_sectors.GetSector(two_body_jjjttz_sector_index);

        int J = two_body_jjjttz_sector.bra_subspace().J();
        int T = two_body_jjjttz_sector.bra_subspace().T();
        int g = two_body_jjjttz_sector.bra_subspace().g();
        int Tz = two_body_jjjttz_sector.bra_subspace().Tz();

        Eigen::MatrixXd& matrix = two_body_jjjttz_matrices[two_body_jjjttz_sector_index];

        // find sector index for corresponding source
        std::size_t two_body_jjjt_bra_subspace_index = two_body_jjjt_space.LookUpSubspaceIndex(basis::TwoBodySubspaceJJJT::SubspaceLabelsType(J,T,g));

        std::size_t two_body_jjjt_sector_index_0 = two_body_jjjt_component_sectors[0].LookUpSectorIndex(two_body_jjjt_bra_subspace_index,two_body_jjjt_bra_subspace_index);
        basis::TwoBodySectorsJJJT::SectorType two_body_jjjt_sector_0 = two_body_jjjt_component_sectors[0].GetSector(two_body_jjjt_sector_index_0);
        std::size_t two_body_jjjt_sector_index_1 = two_body_jjjt_component_sectors[1].LookUpSectorIndex(two_body_jjjt_bra_subspace_index,two_body_jjjt_bra_subspace_index);
        basis::TwoBodySectorsJJJT::SectorType two_body_jjjt_sector_1 = two_body_jjjt_component_sectors[1].GetSector(two_body_jjjt_sector_index_1);
        std::size_t two_body_jjjt_sector_index_2 = two_body_jjjt_component_sectors[2].LookUpSectorIndex(two_body_jjjt_bra_subspace_index,two_body_jjjt_bra_subspace_index);
        basis::TwoBodySectorsJJJT::SectorType two_body_jjjt_sector_2 = two_body_jjjt_component_sectors[2].GetSector(two_body_jjjt_sector_index_2);

        for (std::size_t two_body_jjjttz_bra_state_index=0; two_body_jjjttz_bra_state_index<two_body_jjjttz_sector.bra_subspace().size(); two_body_jjjttz_bra_state_index++)
          {
            for (std::size_t two_body_jjjttz_ket_state_index=two_body_jjjttz_bra_state_index; two_body_jjjttz_ket_state_index<two_body_jjjttz_sector.ket_subspace().size(); two_body_jjjttz_ket_state_index++)
              {
                basis::TwoBodySubspaceJJJTTz::StateLabelsType two_body_jjjttz_bra_state_labels = two_body_jjjttz_sector.bra_subspace().GetStateLabels(two_body_jjjttz_bra_state_index);
                basis::TwoBodySubspaceJJJTTz::StateLabelsType two_body_jjjttz_ket_state_labels = two_body_jjjttz_sector.ket_subspace().GetStateLabels(two_body_jjjttz_ket_state_index);
                std::size_t two_body_jjjt_bra_state_index_0 = two_body_jjjt_sector_0.bra_subspace().LookUpStateIndex(two_body_jjjttz_bra_state_labels);
                std::size_t two_body_jjjt_ket_state_index_0 = two_body_jjjt_sector_0.ket_subspace().LookUpStateIndex(two_body_jjjttz_ket_state_labels);
                std::size_t two_body_jjjt_bra_state_index_1 = two_body_jjjt_sector_1.bra_subspace().LookUpStateIndex(two_body_jjjttz_bra_state_labels);
                std::size_t two_body_jjjt_ket_state_index_1 = two_body_jjjt_sector_1.ket_subspace().LookUpStateIndex(two_body_jjjttz_ket_state_labels);
                std::size_t two_body_jjjt_bra_state_index_2 = two_body_jjjt_sector_2.bra_subspace().LookUpStateIndex(two_body_jjjttz_bra_state_labels);
                std::size_t two_body_jjjt_ket_state_index_2 = two_body_jjjt_sector_2.ket_subspace().LookUpStateIndex(two_body_jjjttz_ket_state_labels);
                if (T==0) {
                  matrix(two_body_jjjttz_bra_state_index,two_body_jjjttz_ket_state_index) = two_body_jjjt_component_matrices[0][two_body_jjjt_sector_index_0](two_body_jjjt_bra_state_index_0,two_body_jjjt_ket_state_index_0);
                } else if (T==1) {
                  if (Tz==1) {
                    matrix(two_body_jjjttz_bra_state_index,two_body_jjjttz_ket_state_index)
                      = two_body_jjjt_component_matrices[0][two_body_jjjt_sector_index_0](two_body_jjjt_bra_state_index_0,two_body_jjjt_ket_state_index_0)
                      + 1.0/std::sqrt(2.0)*two_body_jjjt_component_matrices[1][two_body_jjjt_sector_index_1](two_body_jjjt_bra_state_index_1,two_body_jjjt_ket_state_index_1)
                      + 1.0/std::sqrt(10.0)*two_body_jjjt_component_matrices[2][two_body_jjjt_sector_index_2](two_body_jjjt_bra_state_index_2,two_body_jjjt_ket_state_index_2);
                  } else if (Tz==0) {
                    matrix(two_body_jjjttz_bra_state_index,two_body_jjjttz_ket_state_index)
                      = two_body_jjjt_component_matrices[0][two_body_jjjt_sector_index_0](two_body_jjjt_bra_state_index_0,two_body_jjjt_ket_state_index_0)
                      + std::sqrt(2.0/5.0)*two_body_jjjt_component_matrices[2][two_body_jjjt_sector_index_2](two_body_jjjt_bra_state_index_2,two_body_jjjt_ket_state_index_2);
                  } else if (Tz==-1) {
                    matrix(two_body_jjjttz_bra_state_index,two_body_jjjttz_ket_state_index)
                      = two_body_jjjt_component_matrices[0][two_body_jjjt_sector_index_0](two_body_jjjt_bra_state_index_0,two_body_jjjt_ket_state_index_0)
                      - 1.0/std::sqrt(2.0)*two_body_jjjt_component_matrices[1][two_body_jjjt_sector_index_1](two_body_jjjt_bra_state_index_1,two_body_jjjt_ket_state_index_1)
                      + 1.0/std::sqrt(10.0)*two_body_jjjt_component_matrices[2][two_body_jjjt_sector_index_2](two_body_jjjt_bra_state_index_2,two_body_jjjt_ket_state_index_2);
                  }
                }
              }
          }
      }
  }

void TransformOperatorTwoBodyJJJTTzToTwoBodyJJJT(
    const basis::TwoBodySpaceJJJTTz& two_body_jjjttz_space,
    const basis::TwoBodySectorsJJJTTz& two_body_jjjttz_sectors,
    const basis::OperatorBlocks<double>& two_body_jjjttz_matrices,
    const basis::TwoBodySpaceJJJT& two_body_jjjt_space,
    std::array<basis::TwoBodySectorsJJJT,3>& two_body_jjjt_component_sectors,
    std::array<basis::OperatorBlocks<double>,3>& two_body_jjjt_component_matrices
  ){
    for (int T0 = 0; T0 <= 2; T0++) {
      // enumerate target sectors
      two_body_jjjt_component_sectors[T0] = basis::TwoBodySectorsJJJT(two_body_jjjt_space,0,0,0);

      // populate matrices
      two_body_jjjt_component_matrices[T0].resize(two_body_jjjt_component_sectors[T0].size());
      for (std::size_t two_body_jjjt_sector_index=0; two_body_jjjt_sector_index<two_body_jjjt_component_sectors[T0].size(); two_body_jjjt_sector_index++)
        {
          // make reference to target sector
          const basis::TwoBodySectorsJJJT::SectorType& two_body_jjjt_sector
            = two_body_jjjt_component_sectors[0].GetSector(two_body_jjjt_sector_index);

          int J = two_body_jjjt_sector.bra_subspace().J();
          int T = two_body_jjjt_sector.bra_subspace().T();
          int g = two_body_jjjt_sector.bra_subspace().g();

          Eigen::MatrixXd& matrix = two_body_jjjt_component_matrices[T0][two_body_jjjt_sector_index];

          if (T==0) {
            std::size_t two_body_jjjttz_bra_subspace_index = two_body_jjjttz_space.LookUpSubspaceIndex(basis::TwoBodySubspaceJJJTTz::SubspaceLabelsType(J,T,g,0));
            std::size_t two_body_jjjttz_sector_index = two_body_jjjttz_sectors.LookUpSectorIndex(two_body_jjjttz_bra_subspace_index,two_body_jjjttz_bra_subspace_index);
            basis::TwoBodySectorsJJJTTz::SectorType two_body_jjjttz_sector = two_body_jjjttz_sectors.GetSector(two_body_jjjttz_sector_index);
            for (std::size_t two_body_jjjt_bra_state_index=0; two_body_jjjt_bra_state_index<two_body_jjjt_sector.bra_subspace().size(); two_body_jjjt_bra_state_index++)
              {
                for (std::size_t two_body_jjjt_ket_state_index=two_body_jjjt_bra_state_index; two_body_jjjt_ket_state_index<two_body_jjjt_sector.ket_subspace().size(); two_body_jjjt_ket_state_index++)
                  {
                    basis::TwoBodySubspaceJJJT::StateLabelsType two_body_jjjt_bra_state_labels = two_body_jjjt_sector.bra_subspace().GetStateLabels(two_body_jjjt_bra_state_index);
                    basis::TwoBodySubspaceJJJT::StateLabelsType two_body_jjjt_ket_state_labels = two_body_jjjt_sector.ket_subspace().GetStateLabels(two_body_jjjt_ket_state_index);
                    std::size_t two_body_jjjttz_bra_state_index = two_body_jjjttz_sector.bra_subspace().LookUpStateIndex(two_body_jjjt_bra_state_labels);
                    std::size_t two_body_jjjttz_ket_state_index = two_body_jjjttz_sector.ket_subspace().LookUpStateIndex(two_body_jjjt_ket_state_labels);
                    matrix(two_body_jjjt_bra_state_index,two_body_jjjt_ket_state_index) = two_body_jjjttz_matrices[two_body_jjjttz_sector_index](two_body_jjjttz_bra_state_index,two_body_jjjttz_ket_state_index);
                  }
              }
          } else if (T==1) {
            std::size_t two_body_jjjttz_bra_subspace_index_Tz_p1 = two_body_jjjttz_space.LookUpSubspaceIndex(basis::TwoBodySubspaceJJJTTz::SubspaceLabelsType(J,T,g,1)); // for Tz=+1
            std::size_t two_body_jjjttz_bra_subspace_index_Tz_0 = two_body_jjjttz_space.LookUpSubspaceIndex(basis::TwoBodySubspaceJJJTTz::SubspaceLabelsType(J,T,g,0)); // for Tz=0
            std::size_t two_body_jjjttz_bra_subspace_index_Tz_n1 = two_body_jjjttz_space.LookUpSubspaceIndex(basis::TwoBodySubspaceJJJTTz::SubspaceLabelsType(J,T,g,-1)); // for Tz=-1
            std::size_t two_body_jjjttz_sector_index_Tz_p1 = two_body_jjjttz_sectors.LookUpSectorIndex(two_body_jjjttz_bra_subspace_index_Tz_p1,two_body_jjjttz_bra_subspace_index_Tz_p1);
            std::size_t two_body_jjjttz_sector_index_Tz_0 = two_body_jjjttz_sectors.LookUpSectorIndex(two_body_jjjttz_bra_subspace_index_Tz_0,two_body_jjjttz_bra_subspace_index_Tz_0);
            std::size_t two_body_jjjttz_sector_index_Tz_n1 = two_body_jjjttz_sectors.LookUpSectorIndex(two_body_jjjttz_bra_subspace_index_Tz_n1,two_body_jjjttz_bra_subspace_index_Tz_n1);

            basis::TwoBodySectorsJJJTTz::SectorType two_body_jjjttz_sector_Tz_p1 = two_body_jjjttz_sectors.GetSector(two_body_jjjttz_sector_index_Tz_p1);
            basis::TwoBodySectorsJJJTTz::SectorType two_body_jjjttz_sector_Tz_0 = two_body_jjjttz_sectors.GetSector(two_body_jjjttz_sector_index_Tz_0);
            basis::TwoBodySectorsJJJTTz::SectorType two_body_jjjttz_sector_Tz_n1 = two_body_jjjttz_sectors.GetSector(two_body_jjjttz_sector_index_Tz_n1);
            for (std::size_t two_body_jjjt_bra_state_index=0; two_body_jjjt_bra_state_index<two_body_jjjt_sector.bra_subspace().size(); two_body_jjjt_bra_state_index++)
              {
                for (std::size_t two_body_jjjt_ket_state_index=two_body_jjjt_bra_state_index; two_body_jjjt_ket_state_index<two_body_jjjt_sector.ket_subspace().size(); two_body_jjjt_ket_state_index++)
                  {
                    basis::TwoBodySubspaceJJJT::StateLabelsType two_body_jjjt_bra_state_labels = two_body_jjjt_sector.bra_subspace().GetStateLabels(two_body_jjjt_bra_state_index);
                    basis::TwoBodySubspaceJJJT::StateLabelsType two_body_jjjt_ket_state_labels = two_body_jjjt_sector.ket_subspace().GetStateLabels(two_body_jjjt_ket_state_index);
                    std::size_t two_body_jjjttz_bra_state_index_Tz_p1 = two_body_jjjttz_sector_Tz_p1.bra_subspace().LookUpStateIndex(two_body_jjjt_bra_state_labels);
                    std::size_t two_body_jjjttz_ket_state_index_Tz_p1 = two_body_jjjttz_sector_Tz_p1.ket_subspace().LookUpStateIndex(two_body_jjjt_ket_state_labels);
                    std::size_t two_body_jjjttz_bra_state_index_Tz_0 = two_body_jjjttz_sector_Tz_0.bra_subspace().LookUpStateIndex(two_body_jjjt_bra_state_labels);
                    std::size_t two_body_jjjttz_ket_state_index_Tz_0 = two_body_jjjttz_sector_Tz_0.ket_subspace().LookUpStateIndex(two_body_jjjt_ket_state_labels);
                    std::size_t two_body_jjjttz_bra_state_index_Tz_n1 = two_body_jjjttz_sector_Tz_n1.bra_subspace().LookUpStateIndex(two_body_jjjt_bra_state_labels);
                    std::size_t two_body_jjjttz_ket_state_index_Tz_n1 = two_body_jjjttz_sector_Tz_n1.ket_subspace().LookUpStateIndex(two_body_jjjt_ket_state_labels);
                    if (T0==0) {
                      matrix(two_body_jjjt_bra_state_index,two_body_jjjt_ket_state_index)
                       = 1/3.0*(two_body_jjjttz_matrices[two_body_jjjttz_sector_index_Tz_p1](two_body_jjjttz_bra_state_index_Tz_p1,two_body_jjjttz_ket_state_index_Tz_p1)
                       + two_body_jjjttz_matrices[two_body_jjjttz_sector_index_Tz_0](two_body_jjjttz_bra_state_index_Tz_0,two_body_jjjttz_ket_state_index_Tz_0)
                       + two_body_jjjttz_matrices[two_body_jjjttz_sector_index_Tz_n1](two_body_jjjttz_bra_state_index_Tz_n1,two_body_jjjttz_ket_state_index_Tz_n1));
                    } else if (T0==1) {
                      matrix(two_body_jjjt_bra_state_index,two_body_jjjt_ket_state_index)
                       = 1/3.0*(std::sqrt(9.0/2.0)*two_body_jjjttz_matrices[two_body_jjjttz_sector_index_Tz_p1](two_body_jjjttz_bra_state_index_Tz_p1,two_body_jjjttz_ket_state_index_Tz_p1)
                       + std::sqrt(-9.0/2.0)*two_body_jjjttz_matrices[two_body_jjjttz_sector_index_Tz_n1](two_body_jjjttz_bra_state_index_Tz_n1,two_body_jjjttz_ket_state_index_Tz_n1));
                    } else if (T0==2) {
                      matrix(two_body_jjjt_bra_state_index,two_body_jjjt_ket_state_index)
                       = 1/3.0*(std::sqrt(5.0/2.0)*two_body_jjjttz_matrices[two_body_jjjttz_sector_index_Tz_p1](two_body_jjjttz_bra_state_index_Tz_p1,two_body_jjjttz_ket_state_index_Tz_p1)
                       - std::sqrt(10.0)*two_body_jjjttz_matrices[two_body_jjjttz_sector_index_Tz_0](two_body_jjjttz_bra_state_index_Tz_0,two_body_jjjttz_ket_state_index_Tz_0)
                       + std::sqrt(5.0/2.0)*two_body_jjjttz_matrices[two_body_jjjttz_sector_index_Tz_n1](two_body_jjjttz_bra_state_index_Tz_n1,two_body_jjjttz_ket_state_index_Tz_n1));
                    }
                  }
              }
          }
        }
    }

  }

void TFilter(
    size_t mode, // 0 for isoscalar, 1 for non-isoscalar
    const std::array<basis::TwoBodySectorsJJJT,3>& two_body_jjjt_component_sectors,
    std::array<basis::OperatorBlocks<double>,3>& two_body_jjjt_component_matrices
  ) {
    if (mode == 0) {
      basis::SetOperatorToZero(two_body_jjjt_component_sectors[1],two_body_jjjt_component_matrices[1]);
      basis::SetOperatorToZero(two_body_jjjt_component_sectors[2],two_body_jjjt_component_matrices[2]);
    } else if (mode == 1) {
      basis::SetOperatorToZero(two_body_jjjt_component_sectors[0],two_body_jjjt_component_matrices[0]);
    }
  }

int main(int argc, char **argv)
{
  int Nmax = 4;
  int J0 = 0;
  int g0 = 0;
  int Tz0 = 0;
  const basis::TwoBodySpaceJJJT& two_body_jjjt_space(basis::Rank::kTwoBody,Nmax),
  const std::array<basis::TwoBodySectorsJJJT,3>& two_body_jjjt_component_sectors,
  const std::array<basis::OperatorBlocks<double>,3>& two_body_jjjt_component_matrices,
  const basis::TwoBodySpaceJJJTTz& two_body_jjjttz_space,
  basis::TwoBodySectorsJJJTTz& two_body_jjjttz_sectors,
  basis::OperatorBlocks<double>& two_body_jjjttz_matrices
  TransformOperatorTwoBodyJJJTToTwoBodyJJJTTz( // assumptions: J, T, g, Tz are the same between bra and ket for a given matrix element
      const basis::TwoBodySpaceJJJT& two_body_jjjt_space,
      const std::array<basis::TwoBodySectorsJJJT,3>& two_body_jjjt_component_sectors,
      const std::array<basis::OperatorBlocks<double>,3>& two_body_jjjt_component_matrices,
      const basis::TwoBodySpaceJJJTTz& two_body_jjjttz_space,
      basis::TwoBodySectorsJJJTTz& two_body_jjjttz_sectors,
      basis::OperatorBlocks<double>& two_body_jjjttz_matrices
    )
}
