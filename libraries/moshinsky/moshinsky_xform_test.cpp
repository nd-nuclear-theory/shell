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
  // define operator labels
  ////////////////////////////////////////////////////////////////

  // define operator properties
  int Nmax_relative = 4;
  int Jmax_relative = Nmax_relative+1;
  int Nmax = 2;  // target operator

  basis::OperatorLabelsJT operator_labels;
  operator_labels.J0 = 0;
  operator_labels.g0 = 0;
  operator_labels.T0_min = 0;
  operator_labels.T0_max = 2;
  operator_labels.symmetry_phase_mode = basis::SymmetryPhaseMode::kHermitian;

  ////////////////////////////////////////////////////////////////
  // construct relative identity operator
  ////////////////////////////////////////////////////////////////

  // define space and operator containers
  basis::RelativeSpaceLSJT relative_space(Nmax_relative,Jmax_relative);
  std::array<basis::RelativeSectorsLSJT,3> relative_component_sectors;
  std::array<basis::MatrixVector,3> relative_component_matrices;

  // do construction
  ConstructIdentityOperatorRelativeLSJT(
      operator_labels,
      relative_space,relative_component_sectors,relative_component_matrices
    );

  ////////////////////////////////////////////////////////////////
  // augment to relative-cm NLSJT
  ////////////////////////////////////////////////////////////////

  std::cout << "  relative-cm NLSJT" << std::endl;

  // define space and operator containers
  basis::RelativeCMSpaceNLSJT relative_cm_nlsjt_space(Nmax);
  std::array<basis::RelativeCMSectorsNLSJT,3> relative_cm_nlsjt_component_sectors;
  std::array<basis::MatrixVector,3> relative_cm_nlsjt_component_matrices;

  // do transformation
  moshinsky::TransformOperatorRelativeLSJTToRelativeCMNLSJT(
      operator_labels,
      relative_space,relative_component_sectors,relative_component_matrices,
      relative_cm_nlsjt_space,relative_cm_nlsjt_component_sectors,relative_cm_nlsjt_component_matrices
    );

  ////////////////////////////////////////////////////////////////
  // transform to two-body NLSJT
  ////////////////////////////////////////////////////////////////

  std::cout << "  two-body NLSJT" << std::endl;

  // define space and operator containers
  basis::TwoBodySpaceNLSJT two_body_nlsjt_space(Nmax);
  std::array<basis::TwoBodySectorsNLSJT,3> two_body_nlsjt_component_sectors;
  std::array<basis::MatrixVector,3> two_body_nlsjt_component_matrices;

  // do transformation
  moshinsky::TransformOperatorRelativeCMNLSJTToTwoBodyNLSJT(
      operator_labels,
      relative_cm_nlsjt_space,relative_cm_nlsjt_component_sectors,relative_cm_nlsjt_component_matrices,
      two_body_nlsjt_space,two_body_nlsjt_component_sectors,two_body_nlsjt_component_matrices
  );

  // write sector matrices for inspection
  //
  // Discussion: The T0=0 diagonal sectors should be identity-like,
  // but with 2's on the diagonal if the state is a like-orbital
  // state.  To aid in checking this, we display the subspace contents
  // for the diagonal sectors.  All other sectors should be vanishing.

  std::cout << "two-body NLSJT matrices" << std::endl;
  for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
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
  // reassemble to two-body LSJT as intermediate diagnostic output
  ////////////////////////////////////////////////////////////////

  std::cout << "  two-body LSJT" << std::endl;

  // define space and operator containers
  basis::TwoBodySpaceLSJT two_body_lsjt_space(Nmax);
  std::array<basis::TwoBodySectorsLSJT,3> two_body_lsjt_component_sectors;
  std::array<basis::MatrixVector,3> two_body_lsjt_component_matrices;

  // construct gathered operator
  basis::GatherBlocksTwoBodyNLSJTToTwoBodyLSJT(
      operator_labels,
      two_body_nlsjt_space,two_body_nlsjt_component_sectors,two_body_nlsjt_component_matrices,
      two_body_lsjt_space,two_body_lsjt_component_sectors,two_body_lsjt_component_matrices
    );

  // write operator
  std::string two_body_lsjt_filename("test/moshinsky_xform_test_two_body_identity_AS.dat");
  std::ostringstream os;
  for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
    {
      basis::WriteTwoBodyOperatorComponentLSJT(
          os,
          T0,
          two_body_lsjt_component_sectors[T0],two_body_lsjt_component_matrices[T0],
          basis::NormalizationConversion::kNone
        );
    }
  std::ofstream ofile(two_body_lsjt_filename.c_str());
  ofile << os.str();

  ////////////////////////////////////////////////////////////////
  // recouple to two-body NJJJT
  ////////////////////////////////////////////////////////////////

  std::cout << "  two-body NJJJT" << std::endl;

  // define space and operator containers
  basis::TwoBodySpaceNJJJT two_body_njjjt_space(Nmax);
  std::array<basis::TwoBodySectorsNJJJT,3> two_body_njjjt_component_sectors;
  std::array<basis::MatrixVector,3> two_body_njjjt_component_matrices;

  // do recoupling
  // TODO

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
