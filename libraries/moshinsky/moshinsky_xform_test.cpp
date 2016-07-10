/****************************************************************
  moshinsky_xform_test.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include <fstream>
#include <iomanip>

#include "moshinsky/moshinsky_xform.h"

void test_relative_cm()
{

  std::cout << "Obtaining the relative-cm matrix" << std::endl;

  // set up relative space
  std::cout << "  relative" << std::endl;
  int Nmax = 2;
  int Jmax = Nmax+1;
  basis::RelativeSpaceLSJT relative_space(Nmax,Jmax);
  std::cout << relative_space.DebugStr();

  // set up relative identity operator
  int J0 = 0;
  int T0 = 0;
  int g0 = 0;
  basis::SymmetryPhaseMode symmetry_phase_mode = basis::SymmetryPhaseMode::kHermitian;

  basis::RelativeSectorsLSJT relative_sectors(relative_space,J0,T0,g0);
  basis::MatrixVector relative_matrices;
  basis::SetOperatorToIdentity(relative_sectors,relative_matrices);

  // inspect relative identity operator (sanity check)
  for (int sector_index=0; sector_index<relative_sectors.size(); ++sector_index)
    {
      const basis::RelativeSectorsLSJT::SectorType& relative_sector = relative_sectors.GetSector(sector_index);
      std::cout << " sector " << sector_index << " diagonal " << relative_sector.IsDiagonal() << std::endl;

      // if(relative_sector.IsDiagonal())
      //   {
      //     std::cout << relative_sector.bra_subspace().LabelStr() << std::endl;
      //     std::cout << relative_sector.bra_subspace().DebugStr() << std::endl;
      //   }

      std::cout << relative_matrices[sector_index] << std::endl;
    }

  // set up target space
  std::cout << "  relative-cm" << std::endl;
  basis::RelativeCMSpaceNLSJT relative_cm_space(Nmax);
  basis::RelativeCMSectorsNLSJT relative_cm_sectors(relative_cm_space,J0,T0,g0);
  basis::MatrixVector relative_cm_matrices;
  relative_cm_matrices.resize(relative_cm_sectors.size());

  for (int sector_index=0; sector_index<relative_cm_sectors.size(); ++sector_index)
    {
      const basis::RelativeCMSectorsNLSJT::SectorType& relative_cm_sector = relative_cm_sectors.GetSector(sector_index);

      std::cout << " sector " << sector_index << " diagonal " << relative_cm_sector.IsDiagonal() << std::endl;

      // if(relative_cm_sector.IsDiagonal())
      //   {
      //     std::cout << relative_cm_sector.bra_subspace().LabelStr() << std::endl;
      //     std::cout << relative_cm_sector.bra_subspace().DebugStr() << std::endl;
      //   }

      relative_cm_matrices[sector_index] = moshinsky::RelativeCMMatrixNLSJT(
          relative_space,relative_sectors,relative_matrices,
          relative_cm_sector,
          J0, T0, g0,
          symmetry_phase_mode
        );

      std::cout << relative_cm_matrices[sector_index] << std::endl;
    }
}


void test_moshinsky_matrix()
{

  std::cout << "Inspecting a Moshinsky matrix" << std::endl;

  // set up target sector
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
  basis::TwoBodySubspaceNLSJT two_body_subspace(L,S,J,T,g,N);
  std::cout << two_body_subspace.DebugStr();

  std::cout << "  Moshinsky transform matrix (includes antisymmetry factor)" << std::endl;
  Eigen::MatrixXd matrix = moshinsky::TransformationMatrixRelativeCMTwoBodyNLSJT(relative_cm_subspace,two_body_subspace);
  std::cout << matrix << std::endl;

  std::cout << "  Orthogonality test (expansion of two-body AS basis in antisymmetric relative-cm basis)" << std::endl;

  // Discussion
  //
  // If one uses bare Moshinsky brackets, normality is not expected,
  // and orthogonality is not immediately obvious.  (The resulting
  // matrix turns out to be diagonal with entries 0.5 except for 1 on
  // the like-orbital state.)  The Moshinsky brackets are unitary on
  // the *distinguishable* particle space.  However, our indexing does
  // not run over the full range of distinguishable particle indices.
  // The relative-cm states have the antisymmetry constraint lr+S+T~1
  // (or, equivalentsly, Nr+S+T~1), while the two-body states have the
  // constraint of canonical ordering on the single-particle orbitals.
  // We thus expect to "lose half the norm", roughly speaking, going
  // either way.  As far as the Moshinsky brackets are concerned, in
  // terms of distinguishable particle states, we are expanding each
  // relative-cm state in terms of an incomplete set of two-body
  // states (about half missing) or vice versa.
  //
  // However, a factor of sqrt(2) brings us to the transformation
  // brackets between antisymmetry-allowed relative-cm states to
  // *antisymmetrized* (AS) two-body states.  We thus expect
  // orthogonality, and normality for the expansion of two-body states
  // in terms of relative-cm states, except on the like-orbital
  // states, which will have a norm of 2.

  std::cout << matrix.transpose()*matrix << std::endl;

}

void test_transform_identity()
{

  std::cout << "test_transform_identity" << std::endl;

  ////////////////////////////////////////////////////////////////
  // define relative identity operator
  ////////////////////////////////////////////////////////////////

  // define operator properties
  int Nmax = 2;
  int J0 = 0;
  int T0 = 0;
  int g0 = 0;
  basis::SymmetryPhaseMode symmetry_phase_mode = basis::SymmetryPhaseMode::kHermitian;

  // construct relative identity operator
  int Jmax = Nmax+1;
  basis::RelativeSpaceLSJT relative_space(Nmax,Jmax);
  std::vector<basis::RelativeSectorsLSJT> relative_component_sectors(3);
  std::vector<basis::MatrixVector> relative_component_matrices(3);
  for (int T0=0; T0<=2; ++T0)
    // for each isospin component
    {

      // enumerate sectors
      relative_component_sectors[T0] = basis::RelativeSectorsLSJT(relative_space,J0,T0,g0);
         
      // populate matrices
      relative_component_matrices[T0].resize(relative_component_sectors[T0].size());
      if (T0==0)
        basis::SetOperatorToIdentity(relative_component_sectors[T0],relative_component_matrices[T0]);
      else
        basis::SetOperatorToZero(relative_component_sectors[T0],relative_component_matrices[T0]);
    }

  ////////////////////////////////////////////////////////////////
  // augment to relative-cm NLSJT
  ////////////////////////////////////////////////////////////////

  std::cout << "  relative-cm NLSJT" << std::endl;

  // construct augmented operator
  basis::RelativeCMSpaceNLSJT relative_cm_space(Nmax);
  std::vector<basis::RelativeCMSectorsNLSJT> relative_cm_component_sectors(3);
  std::vector<basis::MatrixVector> relative_cm_component_matrices(3);
  for (int T0=0; T0<=2; ++T0)
    // for each isospin component
    {

      // enumerate sectors
      relative_cm_component_sectors[T0]
        = basis::RelativeCMSectorsNLSJT(relative_cm_space,J0,T0,g0);

      // std::cout << " T0 " << T0
      //           << " sectors " << relative_cm_component_sectors[T0].size()
      //           << std::endl;

      // populate matrices
      relative_cm_component_matrices[T0].resize(relative_cm_component_sectors[T0].size());
      for (int sector_index=0; sector_index<relative_cm_component_sectors[T0].size(); ++sector_index)
        {
          const basis::RelativeCMSectorsNLSJT::SectorType& relative_cm_sector
            = relative_cm_component_sectors[T0].GetSector(sector_index);
          relative_cm_component_matrices[T0][sector_index] = moshinsky::RelativeCMMatrixNLSJT(
              relative_space,relative_component_sectors[T0],relative_component_matrices[T0],
              relative_cm_sector,
              J0, T0, g0,
              symmetry_phase_mode
            );
        }
    }

  ////////////////////////////////////////////////////////////////
  // transform to two-body NLSJT
  ////////////////////////////////////////////////////////////////

  std::cout << "  two-body NLSJT" << std::endl;

  // construct transformed operator
  basis::TwoBodySpaceNLSJT two_body_nlsjt_space(Nmax);
  std::vector<basis::TwoBodySectorsNLSJT> two_body_nlsjt_component_sectors(3);
  std::vector<basis::MatrixVector> two_body_nlsjt_component_matrices(3);
  for (int T0=0; T0<=2; ++T0)
    // for each isospin component
    {

      // enumerate sectors
      two_body_nlsjt_component_sectors[T0]
        = basis::TwoBodySectorsNLSJT(two_body_nlsjt_space,J0,T0,g0);
         
      // std::cout << " T0 " << T0
      //           << " sectors " << two_body_nlsjt_component_sectors[T0].size()
      //           << std::endl;

      // populate matrices
      two_body_nlsjt_component_matrices[T0].resize(two_body_nlsjt_component_sectors[T0].size());
      for (int sector_index=0; sector_index<two_body_nlsjt_component_sectors[T0].size(); ++sector_index)
        {
          // target sector
          const basis::TwoBodySectorsNLSJT::SectorType& two_body_nlsjt_sector
            = two_body_nlsjt_component_sectors[T0].GetSector(sector_index);

          // look up source sector
          int relative_cm_bra_subspace_index = relative_cm_space.LookUpSubspaceIndex(
              basis::RelativeCMSubspaceNLSJTLabels(
                  two_body_nlsjt_sector.bra_subspace().L(),
                  two_body_nlsjt_sector.bra_subspace().S(),
                  two_body_nlsjt_sector.bra_subspace().J(),
                  two_body_nlsjt_sector.bra_subspace().T(),
                  two_body_nlsjt_sector.bra_subspace().g(),
                  two_body_nlsjt_sector.bra_subspace().N()
                )
            );
          int relative_cm_ket_subspace_index = relative_cm_space.LookUpSubspaceIndex(
              basis::RelativeCMSubspaceNLSJTLabels(
                  two_body_nlsjt_sector.ket_subspace().L(),
                  two_body_nlsjt_sector.ket_subspace().S(),
                  two_body_nlsjt_sector.ket_subspace().J(),
                  two_body_nlsjt_sector.ket_subspace().T(),
                  two_body_nlsjt_sector.ket_subspace().g(),
                  two_body_nlsjt_sector.ket_subspace().N()
                )
            );
          int relative_cm_sector_index = relative_cm_component_sectors[T0].LookUpSectorIndex(
                relative_cm_bra_subspace_index,
                relative_cm_ket_subspace_index
              );
          const basis::RelativeCMSectorsNLSJT::SectorType& relative_cm_sector
            = relative_cm_component_sectors[T0].GetSector(relative_cm_sector_index);
          const Eigen::MatrixXd& relative_cm_matrix
            = relative_cm_component_matrices[T0][relative_cm_sector_index];

          // transform
          two_body_nlsjt_component_matrices[T0][sector_index] = moshinsky::TwoBodyMatrixNLSJT(
              relative_cm_sector,
              two_body_nlsjt_sector,
              relative_cm_matrix
            );
        }
    }

  // write sector matrices for inspection
  //
  // Discussion: The T0=0 diagonal sectors should be identity-like,
  // but with 2's on the diagonal if the state is a like-orbital
  // state.  To aid in checking this, we display the subspace contents
  // for the diagonal sectors.  All other sectors should be vanishing.

  std::cout << "two-body NLSJT matrices" << std::endl;
  for (int T0=0; T0<=2; ++T0)
    for (int sector_index=0; sector_index<two_body_nlsjt_component_sectors[T0].size(); ++sector_index)
      {
      const basis::TwoBodySectorsNLSJT::SectorType& two_body_nlsjt_sector
        = two_body_nlsjt_component_sectors[T0].GetSector(sector_index);

      std::cout << " T0 " << T0
                << " sector " << sector_index
                << " diagonal " << two_body_nlsjt_sector.IsDiagonal() << std::endl;
      std::cout << two_body_nlsjt_component_matrices[T0][sector_index] << std::endl;
      if (two_body_nlsjt_sector.IsDiagonal())
        {
          // std::cout << two_body_nlsjt_sector.ket_subspace().LabelStr() << std::endl;
          std::cout << two_body_nlsjt_sector.ket_subspace().DebugStr();
          std::cout << std::endl;
        }
      std::cout << std::endl;
  
    }

  ////////////////////////////////////////////////////////////////
  // reassemble to two-body LSJT as intermediate output
  ////////////////////////////////////////////////////////////////

  std::cout << "  two-body LSJT" << std::endl;

  // construct reassembled operator
  basis::TwoBodySpaceLSJT two_body_lsjt_space(Nmax);
  std::vector<basis::TwoBodySectorsLSJT> two_body_lsjt_component_sectors(3);
  std::vector<basis::MatrixVector> two_body_lsjt_component_matrices(3);

  for (int T0=0; T0<=2; ++T0)
    // for each isospin component
    {

      // enumerate sectors
      two_body_lsjt_component_sectors[T0]
        = basis::TwoBodySectorsLSJT(two_body_lsjt_space,J0,T0,g0);

      // populate matrices
      two_body_lsjt_component_matrices[T0].resize(two_body_lsjt_component_sectors[T0].size());
      for (int sector_index=0; sector_index<two_body_lsjt_component_sectors[T0].size(); ++sector_index)
        {
          // retrieve target sector
          const basis::TwoBodySectorsLSJT::SectorType& two_body_lsjt_sector
            = two_body_lsjt_component_sectors[T0].GetSector(sector_index);

          // initialize matrix
          Eigen::MatrixXd& two_body_lsjt_matrix = two_body_lsjt_component_matrices[T0][sector_index];
          two_body_lsjt_matrix = Eigen::MatrixXd::Zero(
              two_body_lsjt_sector.bra_subspace().size(),
              two_body_lsjt_sector.ket_subspace().size()
            );

          // populate matrix elements
          for (int bra_index = 0; bra_index < two_body_lsjt_sector.bra_subspace().size(); ++bra_index)
            for (int ket_index = 0; ket_index < two_body_lsjt_sector.ket_subspace().size(); ++ket_index)
              // for each target matrix element
              {
                
                // retrieve target states
                basis::TwoBodyStateLSJT two_body_lsjt_bra(two_body_lsjt_sector.bra_subspace(),bra_index);
                basis::TwoBodyStateLSJT two_body_lsjt_ket(two_body_lsjt_sector.ket_subspace(),ket_index);

                // look up source subspace indices
                int two_body_nlsjt_bra_subspace_index = two_body_nlsjt_space.LookUpSubspaceIndex(
                    basis::TwoBodySubspaceNLSJTLabels(
                        two_body_lsjt_bra.L(),
                        two_body_lsjt_bra.S(),
                        two_body_lsjt_bra.J(),
                        two_body_lsjt_bra.T(),
                        two_body_lsjt_bra.g(),
                        two_body_lsjt_bra.N()
                      )
                  );
                int two_body_nlsjt_ket_subspace_index = two_body_nlsjt_space.LookUpSubspaceIndex(
                    basis::TwoBodySubspaceNLSJTLabels(
                        two_body_lsjt_ket.L(),
                        two_body_lsjt_ket.S(),
                        two_body_lsjt_ket.J(),
                        two_body_lsjt_ket.T(),
                        two_body_lsjt_ket.g(),
                        two_body_lsjt_ket.N()
                      )
                  );

                // look up source matrix element indices
                const basis::TwoBodySubspaceNLSJT& two_body_nlsjt_bra_subspace
                  = two_body_nlsjt_space.GetSubspace(two_body_nlsjt_bra_subspace_index);
                const basis::TwoBodySubspaceNLSJT& two_body_nlsjt_ket_subspace
                  = two_body_nlsjt_space.GetSubspace(two_body_nlsjt_ket_subspace_index);
                int two_body_nlsjt_bra_index = two_body_nlsjt_bra_subspace.LookUpStateIndex(
                    basis::TwoBodyStateNLSJT::StateLabelsType(
                        two_body_lsjt_bra.N1(),
                        two_body_lsjt_bra.l1(),
                        two_body_lsjt_bra.N2(),
                        two_body_lsjt_bra.l2()
                      )
                  );
                int two_body_nlsjt_ket_index = two_body_nlsjt_ket_subspace.LookUpStateIndex(
                    basis::TwoBodyStateNLSJT::StateLabelsType(
                        two_body_lsjt_ket.N1(),
                        two_body_lsjt_ket.l1(),
                        two_body_lsjt_ket.N2(),
                        two_body_lsjt_ket.l2()
                      )
                  );

                // canonicalize indices for matrix element lookup
                //
                // We must ensure that we look up a canonical (upper
                // triangular) NLSJT sector.  Looking up a canonical
                // (upper triangular) matrix element within a diagonal
                // sector is not essential, since the Moshinsky
                // transformation machinery up until this point has
                // actually been populating the full (square) matrices
                // for the diagonal sectors.
                //
                // Note that no canonicalization factor is needed.
                // Since N is the trailing entry in the subspace label
                // tuple, used in the canonical ordering, canonical
                // swaps will never entail swapping subspace LSJT
                // labels, just the N labels.

                // std::cout << " pre-lookup "
                //           << " " << two_body_nlsjt_bra_subspace_index
                //           << " " << two_body_nlsjt_ket_subspace_index
                //           << " " << ";"
                //           << " " << two_body_nlsjt_bra_index
                //           << " " << two_body_nlsjt_ket_index
                //           << " " << ";"
                //           << " " << two_body_nlsjt_bra_subspace.size()
                //           << " " << two_body_nlsjt_ket_subspace.size()
                //           << std::endl;

                bool swapped_subspaces, swapped_states;
                basis::CanonicalizeIndices(
                    two_body_nlsjt_bra_subspace_index,two_body_nlsjt_ket_subspace_index,
                    swapped_subspaces,
                    two_body_nlsjt_bra_index,two_body_nlsjt_ket_index,
                    swapped_states
                  );

                // look up matrix element
                int two_body_nlsjt_sector_index
                  = two_body_nlsjt_component_sectors[T0].LookUpSectorIndex(
                      two_body_nlsjt_bra_subspace_index,
                      two_body_nlsjt_ket_subspace_index
                    );

                Eigen::MatrixXd& two_body_nlsjt_matrix
                  = two_body_nlsjt_component_matrices[T0][two_body_nlsjt_sector_index];
                // std::cout << " lookup "
                //           << " " << two_body_nlsjt_bra_subspace_index
                //           << " " << two_body_nlsjt_ket_subspace_index
                //           << " " << ";"
                //           << " " << two_body_nlsjt_bra_index
                //           << " " << two_body_nlsjt_ket_index
                //           << " " << ";"
                //           << " " << two_body_nlsjt_matrix.rows()
                //           << " " << two_body_nlsjt_matrix.cols()
                //           << std::endl;
                double two_body_nlsjt_matrix_element = two_body_nlsjt_matrix(
                    two_body_nlsjt_bra_index,two_body_nlsjt_ket_index
                  );

                // re-save matrix element
                // std::cout << " write "
                //           << " " << bra_index
                //           << " " << ket_index
                //           << " " << ";"
                //           << " " << two_body_lsjt_matrix.rows()
                //           << " " << two_body_lsjt_matrix.cols()
                //           << std::endl;

                two_body_lsjt_matrix(bra_index,ket_index) = two_body_nlsjt_matrix_element;

              }
        }
    }


  // write operator
  std::string relative_cm_filename("test/moshinsky_xform_test_two_body_identity_AS.dat");
  std::ostringstream os;
  for (int T0=0; T0<=2; ++T0)
    {
      basis::WriteTwoBodyOperatorComponentLSJT(
          os,
          T0,
          two_body_lsjt_component_sectors[T0],two_body_lsjt_component_matrices[T0],
          basis::NormalizationConversion::kNone
        );
    }
  std::ofstream ofile(relative_cm_filename.c_str());
  ofile << os.str();



  
}

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  test_relative_cm();
  test_moshinsky_matrix();
  test_transform_identity();

  // termination
  return 0;
}
