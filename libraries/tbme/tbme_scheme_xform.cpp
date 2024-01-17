/****************************************************************
  tbme_scheme_xform.cpp

  Zhou Zhou
  University of Notre Dame

****************************************************************/

#include "tbme/tbme_scheme_xform.h"

namespace shell {
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  void TransformOperatorTwoBodyJJJTToTwoBodyJJJTTz( // assumptions: J, T, g, Tz are the same between bra and ket for a given matrix element
      const basis::TwoBodySpaceJJJT& two_body_jjjt_space,
      const std::array<basis::TwoBodySectorsJJJT,3>& two_body_jjjt_component_sectors,
      const std::array<basis::OperatorBlocks<double>,3>& two_body_jjjt_component_matrices,
      const basis::TwoBodySpaceJJJTTz& two_body_jjjttz_space,
      basis::TwoBodySectorsJJJTTz& two_body_jjjttz_sectors,
      basis::OperatorBlocks<double>& two_body_jjjttz_matrices
    ) {
      int J0=0;
      int g0=0;
      int Tz0=0;
      // enumerate target sectors
      two_body_jjjttz_sectors = basis::TwoBodySectorsJJJTTz(two_body_jjjttz_space,J0,g0,Tz0);

      // populate matrices
      basis::SetOperatorToZero(two_body_jjjttz_sectors,two_body_jjjttz_matrices);
      for (std::size_t two_body_jjjttz_sector_index=0; two_body_jjjttz_sector_index<two_body_jjjttz_sectors.size(); two_body_jjjttz_sector_index++)
        {
          // make reference to target sector
          const basis::TwoBodySectorsJJJTTz::SectorType& two_body_jjjttz_sector
            = two_body_jjjttz_sectors.GetSector(two_body_jjjttz_sector_index);

          // matrix elements between states with different T are assumed to vanish and do not exist in me2j files (or jjjttz format)
          if (two_body_jjjttz_sector.bra_subspace().T()!=two_body_jjjttz_sector.ket_subspace().T()) {
            continue;
          }

          int J = two_body_jjjttz_sector.bra_subspace().J();
          int T = two_body_jjjttz_sector.bra_subspace().T();
          int g = two_body_jjjttz_sector.bra_subspace().g();
          int Tz = two_body_jjjttz_sector.bra_subspace().Tz();

          Eigen::MatrixXd& matrix = two_body_jjjttz_matrices[two_body_jjjttz_sector_index];
          // find sector index for corresponding source
          std::size_t two_body_jjjt_bra_subspace_index = two_body_jjjt_space.LookUpSubspaceIndex(basis::TwoBodySubspaceJJJT::SubspaceLabelsType(J,T,g));
          std::size_t two_body_jjjt_ket_subspace_index = two_body_jjjt_space.LookUpSubspaceIndex(basis::TwoBodySubspaceJJJT::SubspaceLabelsType(J,T,g));
          // std::cout << "two_body_jjjt_bra_subspace_index" << two_body_jjjt_bra_subspace_index << std::endl;

          for (std::size_t two_body_jjjttz_bra_state_index=0; two_body_jjjttz_bra_state_index<two_body_jjjttz_sector.bra_subspace().size(); two_body_jjjttz_bra_state_index++)
            {
              for (std::size_t two_body_jjjttz_ket_state_index=two_body_jjjttz_bra_state_index; two_body_jjjttz_ket_state_index<two_body_jjjttz_sector.ket_subspace().size(); two_body_jjjttz_ket_state_index++)
                {
                  // find bra and ket state labels for jjjttz and pass to jjjt to find bra and ket states indices for jjjt
                  // going from target indices to source indices
                  basis::TwoBodySubspaceJJJTTz::StateLabelsType two_body_jjjttz_bra_state_labels = two_body_jjjttz_sector.bra_subspace().GetStateLabels(two_body_jjjttz_bra_state_index);
                  basis::TwoBodySubspaceJJJTTz::StateLabelsType two_body_jjjttz_ket_state_labels = two_body_jjjttz_sector.ket_subspace().GetStateLabels(two_body_jjjttz_ket_state_index);
                  std::array<double,3> matrix_element_by_T0;
                  // find source sector indices from subspace indices
                  for (int T0 = 0; T0 <= 2; T0++) {
                    std::size_t two_body_jjjt_sector_index = two_body_jjjt_component_sectors[T0].LookUpSectorIndex(two_body_jjjt_bra_subspace_index,two_body_jjjt_ket_subspace_index);
                    if (two_body_jjjt_sector_index == basis::kNone) {
                      matrix_element_by_T0[T0] = 0;
                    } else {
                      basis::TwoBodySectorsJJJT::SectorType two_body_jjjt_sector = two_body_jjjt_component_sectors[T0].GetSector(two_body_jjjt_sector_index);
                      std::size_t two_body_jjjt_bra_state_index = two_body_jjjt_sector.bra_subspace().LookUpStateIndex(two_body_jjjttz_bra_state_labels);
                      std::size_t two_body_jjjt_ket_state_index = two_body_jjjt_sector.ket_subspace().LookUpStateIndex(two_body_jjjttz_ket_state_labels);
                      // std::cout << "two_body_jjjt_bra/ket_state_index" << two_body_jjjt_bra_state_index << " " << two_body_jjjt_ket_state_index << std::endl;
                      matrix_element_by_T0[T0] = two_body_jjjt_component_matrices[T0][two_body_jjjt_sector_index](two_body_jjjt_bra_state_index,two_body_jjjt_ket_state_index);
                    }
                  }
                  if (T==0) {
                    matrix(two_body_jjjttz_bra_state_index,two_body_jjjttz_ket_state_index)
                      = matrix_element_by_T0[0];
                  } else if (T==1) {
                    if (Tz==1) {
                      matrix(two_body_jjjttz_bra_state_index,two_body_jjjttz_ket_state_index)
                        = matrix_element_by_T0[0] + 1.0/std::sqrt(2.0)*matrix_element_by_T0[1] + 1.0/std::sqrt(10.0)*matrix_element_by_T0[2];
                    } else if (Tz==0) {
                      matrix(two_body_jjjttz_bra_state_index,two_body_jjjttz_ket_state_index)
                        = matrix_element_by_T0[0] - std::sqrt(2.0/5.0)*matrix_element_by_T0[2];
                    } else if (Tz==-1) {
                      matrix(two_body_jjjttz_bra_state_index,two_body_jjjttz_ket_state_index)
                        = matrix_element_by_T0[0] - 1.0/std::sqrt(2.0)*matrix_element_by_T0[1] + 1.0/std::sqrt(10.0)*matrix_element_by_T0[2];
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
        int J0=0;
        int g0=0;
        // enumerate target sectors
        two_body_jjjt_component_sectors[T0] = basis::TwoBodySectorsJJJT(two_body_jjjt_space,J0,T0,g0);
        // populate matrices
        basis::SetOperatorToZero(two_body_jjjt_component_sectors[T0],two_body_jjjt_component_matrices[T0]);
        // two_body_jjjt_component_matrices[T0].resize(two_body_jjjt_component_sectors[T0].size());
        for (std::size_t two_body_jjjt_sector_index=0; two_body_jjjt_sector_index<two_body_jjjt_component_sectors[T0].size(); two_body_jjjt_sector_index++)
          {
            // make reference to target sector
            const basis::TwoBodySectorsJJJT::SectorType& two_body_jjjt_sector
              = two_body_jjjt_component_sectors[T0].GetSector(two_body_jjjt_sector_index);

            if (two_body_jjjt_sector.bra_subspace().T()!=two_body_jjjt_sector.ket_subspace().T()) {
              continue;
            }

            int J = two_body_jjjt_sector.bra_subspace().J();
            int T = two_body_jjjt_sector.bra_subspace().T();
            int g = two_body_jjjt_sector.bra_subspace().g();

            Eigen::MatrixXd& matrix = two_body_jjjt_component_matrices[T0][two_body_jjjt_sector_index];
            for (std::size_t two_body_jjjt_bra_state_index=0; two_body_jjjt_bra_state_index<two_body_jjjt_sector.bra_subspace().size(); two_body_jjjt_bra_state_index++)
              {
                for (std::size_t two_body_jjjt_ket_state_index=two_body_jjjt_bra_state_index; two_body_jjjt_ket_state_index<two_body_jjjt_sector.ket_subspace().size(); two_body_jjjt_ket_state_index++)
                  {
                    basis::TwoBodySubspaceJJJT::StateLabelsType two_body_jjjt_bra_state_labels = two_body_jjjt_sector.bra_subspace().GetStateLabels(two_body_jjjt_bra_state_index);
                    basis::TwoBodySubspaceJJJT::StateLabelsType two_body_jjjt_ket_state_labels = two_body_jjjt_sector.ket_subspace().GetStateLabels(two_body_jjjt_ket_state_index);
                    std::array<double,3> matrix_element_by_Tz; // for Tz = -1, 0, 1
                    for (int Tz = -T; Tz <= T; Tz++) {
                      std::size_t two_body_jjjttz_bra_subspace_index_Tz = two_body_jjjttz_space.LookUpSubspaceIndex(basis::TwoBodySubspaceJJJTTz::SubspaceLabelsType(J,T,g,Tz));
                      std::size_t two_body_jjjttz_sector_index_Tz = two_body_jjjttz_sectors.LookUpSectorIndex(two_body_jjjttz_bra_subspace_index_Tz,two_body_jjjttz_bra_subspace_index_Tz);
                      if (two_body_jjjttz_sector_index_Tz == basis::kNone) {
                        matrix_element_by_Tz[Tz+1] = 0;
                      } else {
                        basis::TwoBodySectorsJJJTTz::SectorType two_body_jjjttz_sector_Tz = two_body_jjjttz_sectors.GetSector(two_body_jjjttz_sector_index_Tz);
                        std::size_t two_body_jjjttz_bra_state_index_Tz = two_body_jjjttz_sector_Tz.bra_subspace().LookUpStateIndex(two_body_jjjt_bra_state_labels);
                        std::size_t two_body_jjjttz_ket_state_index_Tz = two_body_jjjttz_sector_Tz.ket_subspace().LookUpStateIndex(two_body_jjjt_ket_state_labels);
                        matrix_element_by_Tz[Tz+1] = two_body_jjjttz_matrices[two_body_jjjttz_sector_index_Tz](two_body_jjjttz_bra_state_index_Tz,two_body_jjjttz_ket_state_index_Tz);
                      }
                    }
                    if (T==0) {
                      matrix(two_body_jjjt_bra_state_index,two_body_jjjt_ket_state_index)
                        = matrix_element_by_Tz[1];
                    } else if (T==1) {
                      if (T0==0) {
                        matrix(two_body_jjjt_bra_state_index,two_body_jjjt_ket_state_index)
                          = 1/3.0*(matrix_element_by_Tz[0]+matrix_element_by_Tz[1]+matrix_element_by_Tz[2]);
                      } else if (T0==1) {
                        matrix(two_body_jjjt_bra_state_index,two_body_jjjt_ket_state_index)
                          = 1/3.0*(-std::sqrt(9.0/2.0)*matrix_element_by_Tz[0]+std::sqrt(9.0/2.0)*matrix_element_by_Tz[2]);
                      } else if (T0==2) {
                        matrix(two_body_jjjt_bra_state_index,two_body_jjjt_ket_state_index)
                          = 1/3.0*(std::sqrt(5.0/2.0)*matrix_element_by_Tz[0]-std::sqrt(10.0)*matrix_element_by_Tz[1]+std::sqrt(5.0/2.0)*matrix_element_by_Tz[2]);
                      }
                    }
                  }
              }
          }
      }

    }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
