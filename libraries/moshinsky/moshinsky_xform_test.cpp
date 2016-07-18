/****************************************************************
  moshinsky_xform_test.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include <fstream>
#include <iomanip>

#include "basis/jjjt_operator.h"

#include "mcpp/profiling.h"
#include "moshinsky/construct_relative.h"
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
  basis::RelativeCMSpaceLSJTN relative_cm_space(Nmax);
  basis::RelativeCMSectorsLSJTN relative_cm_sectors(relative_cm_space,J0,T0,g0);
  basis::MatrixVector relative_cm_matrices;
  relative_cm_matrices.resize(relative_cm_sectors.size());

  for (int sector_index=0; sector_index<relative_cm_sectors.size(); ++sector_index)
    {
      const basis::RelativeCMSectorsLSJTN::SectorType& relative_cm_sector = relative_cm_sectors.GetSector(sector_index);

      std::cout << " sector " << sector_index << " diagonal " << relative_cm_sector.IsDiagonal() << std::endl;

      // if(relative_cm_sector.IsDiagonal())
      //   {
      //     std::cout << relative_cm_sector.bra_subspace().LabelStr() << std::endl;
      //     std::cout << relative_cm_sector.bra_subspace().DebugStr() << std::endl;
      //   }

      relative_cm_matrices[sector_index] = moshinsky::RelativeCMMatrixLSJTN(
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
  basis::RelativeCMSubspaceLSJTN relative_cm_subspace(L,S,J,T,g,N);
  std::cout << relative_cm_subspace.DebugStr();

  std::cout << "  two-body" << std::endl;
  basis::TwoBodySubspaceLSJTN two_body_subspace(L,S,J,T,g,N);
  std::cout << two_body_subspace.DebugStr();

  std::cout << "  Moshinsky transform matrix (includes antisymmetry factor)" << std::endl;
  Eigen::MatrixXd matrix = moshinsky::TransformationMatrixRelativeCMTwoBodyLSJTN(relative_cm_subspace,two_body_subspace);
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

////////////////////////////////////////////////////////////////
void test_transform_simple(
    std::string lsjt_filename,
    std::string jjjt_filename,
    char operator_code
  )
{

  std::cout << "test_transform_simple" << " " << operator_code << std::endl;

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
  operator_labels.T0_max = 0;
  operator_labels.symmetry_phase_mode = basis::SymmetryPhaseMode::kHermitian;

  ////////////////////////////////////////////////////////////////
  // construct relative identity operator
  ////////////////////////////////////////////////////////////////

  // define space and operator containers
  basis::RelativeSpaceLSJT relative_space(Nmax_relative,Jmax_relative);
  std::array<basis::RelativeSectorsLSJT,3> relative_component_sectors;
  std::array<basis::MatrixVector,3> relative_component_matrices;

  // do construction
  if (operator_code=='I')
    ConstructIdentityOperatorRelativeLSJT(
        operator_labels,
        relative_space,relative_component_sectors,relative_component_matrices
      );
  else if (operator_code=='K')
     ConstructKinematicOperator(
         operator_labels,
         relative_space,relative_component_sectors,relative_component_matrices,
         relative::KinematicOperator::kKSqr
       );

  ////////////////////////////////////////////////////////////////
  // augment to relative-cm LSJTN
  ////////////////////////////////////////////////////////////////

  std::cout << "  relative-cm LSJTN" << std::endl;

  // define space and operator containers
  basis::RelativeCMSpaceLSJTN relative_cm_lsjtn_space(Nmax);
  std::array<basis::RelativeCMSectorsLSJTN,3> relative_cm_lsjtn_component_sectors;
  std::array<basis::MatrixVector,3> relative_cm_lsjtn_component_matrices;

  // do transformation
  moshinsky::TransformOperatorRelativeLSJTToRelativeCMLSJTN(
      operator_labels,
      relative_space,relative_component_sectors,relative_component_matrices,
      relative_cm_lsjtn_space,relative_cm_lsjtn_component_sectors,relative_cm_lsjtn_component_matrices
    );

  ////////////////////////////////////////////////////////////////
  // transform to two-body LSJTN
  ////////////////////////////////////////////////////////////////

  std::cout << "  two-body LSJTN" << std::endl;

  // define space and operator containers
  basis::TwoBodySpaceLSJTN two_body_lsjtn_space(Nmax);
  std::array<basis::TwoBodySectorsLSJTN,3> two_body_lsjtn_component_sectors;
  std::array<basis::MatrixVector,3> two_body_lsjtn_component_matrices;

  // do transformation
  Timer two_body_lsjtn_timer;
  two_body_lsjtn_timer.Start();
  moshinsky::TransformOperatorRelativeCMLSJTNToTwoBodyLSJTN(
      operator_labels,
      relative_cm_lsjtn_space,relative_cm_lsjtn_component_sectors,relative_cm_lsjtn_component_matrices,
      two_body_lsjtn_space,two_body_lsjtn_component_sectors,two_body_lsjtn_component_matrices
  );
  two_body_lsjtn_timer.Stop();
  std::cout << "Time: " << two_body_lsjtn_timer.ElapsedTime() << std::endl;

  // write sector matrices for inspection
  //
  // Discussion: The T0=0 diagonal sectors should be identity-like,
  // but with 2's on the diagonal if the state is a like-orbital
  // state.  To aid in checking this, we display the subspace contents
  // for the diagonal sectors.  All other sectors should be vanishing.

  std::cout << "writing two-body LSJTN matrices" << std::endl;
  for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
    for (int sector_index=0; sector_index<two_body_lsjtn_component_sectors[T0].size(); ++sector_index)
      {
      const basis::TwoBodySectorsLSJTN::SectorType& two_body_lsjtn_sector
        = two_body_lsjtn_component_sectors[T0].GetSector(sector_index);
      const Eigen::MatrixXd& two_body_lsjtn_matrix
        = two_body_lsjtn_component_matrices[T0][sector_index];


      std::cout << " T0 " << T0
                << " sector " << sector_index
                << " diagonal " << two_body_lsjtn_sector.IsDiagonal() << std::endl;
      std::cout << two_body_lsjtn_matrix << std::endl;
      if (two_body_lsjtn_sector.IsDiagonal())
        {
          // std::cout << two_body_lsjtn_sector.ket_subspace().LabelStr() << std::endl;
          std::cout << two_body_lsjtn_sector.ket_subspace().DebugStr();
          std::cout << std::endl;
        }
      std::cout << std::endl;
  
    }

  ////////////////////////////////////////////////////////////////
  // gather to two-body LSJT as intermediate diagnostic output
  ////////////////////////////////////////////////////////////////

  std::cout << "  two-body LSJT" << std::endl;

  // define space and operator containers
  basis::TwoBodySpaceLSJT two_body_lsjt_space(Nmax);
  std::array<basis::TwoBodySectorsLSJT,3> two_body_lsjt_component_sectors;
  std::array<basis::MatrixVector,3> two_body_lsjt_component_matrices;

  // construct gathered operator
  basis::GatherOperatorTwoBodyLSJTNToTwoBodyLSJT(
      operator_labels,
      two_body_lsjtn_space,two_body_lsjtn_component_sectors,two_body_lsjtn_component_matrices,
      two_body_lsjt_space,two_body_lsjt_component_sectors,two_body_lsjt_component_matrices
    );

  ////////////////////////////////////////////////////////////////
  // write as two-body LSJT
  ////////////////////////////////////////////////////////////////

  std::string two_body_lsjt_filename(lsjt_filename);
  std::ostringstream lsjt_sstream;
  for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
    {
      basis::WriteTwoBodyOperatorComponentLSJT(
          lsjt_sstream,
          T0,
          two_body_lsjt_component_sectors[T0],two_body_lsjt_component_matrices[T0],
          basis::NormalizationConversion::kNone
        );
    }
  std::ofstream lsjt_stream(two_body_lsjt_filename.c_str());
  lsjt_stream << lsjt_sstream.str();

  ////////////////////////////////////////////////////////////////
  // recouple to two-body JJJTN
  ////////////////////////////////////////////////////////////////

  std::cout << "  two-body JJJTN" << std::endl;

  // define space and operator containers
  basis::TwoBodySpaceJJJTN two_body_jjjtn_space(Nmax);
  std::array<basis::TwoBodySectorsJJJTN,3> two_body_jjjtn_component_sectors;
  std::array<basis::MatrixVector,3> two_body_jjjtn_component_matrices;

  // do recoupling
  Timer two_body_jjjtn_timer;
  two_body_jjjtn_timer.Start();
  moshinsky::TransformOperatorTwoBodyLSJTNToTwoBodyJJJTN(
      operator_labels,
      two_body_lsjtn_space,two_body_lsjtn_component_sectors,two_body_lsjtn_component_matrices,
      two_body_jjjtn_space,two_body_jjjtn_component_sectors,two_body_jjjtn_component_matrices
    );
  two_body_jjjtn_timer.Stop();
  std::cout << "Time: " << two_body_jjjtn_timer.ElapsedTime() << std::endl;

  // write sector matrices for inspection
  //
  // Discussion: The T0=0 diagonal sectors should be identity-like,
  // but with 2's on the diagonal if the state is a like-orbital
  // state.  To aid in checking this, we display the subspace contents
  // for the diagonal sectors.  All other sectors should be vanishing.

  std::cout << "writing two-body JJJTN matrices" << std::endl;
  for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
    for (int sector_index=0; sector_index<two_body_jjjtn_component_sectors[T0].size(); ++sector_index)
      {
      const basis::TwoBodySectorsJJJTN::SectorType& two_body_jjjtn_sector
        = two_body_jjjtn_component_sectors[T0].GetSector(sector_index);
      const Eigen::MatrixXd& two_body_jjjtn_matrix
        = two_body_jjjtn_component_matrices[T0][sector_index];

      std::cout << " T0 " << T0
                << " sector " << sector_index
                << " diagonal " << two_body_jjjtn_sector.IsDiagonal() << std::endl;
      std::cout << two_body_jjjtn_matrix << std::endl;
      if (two_body_jjjtn_sector.IsDiagonal())
        {
          // std::cout << two_body_lsjtn_sector.ket_subspace().LabelStr() << std::endl;
          std::cout << two_body_jjjtn_sector.ket_subspace().DebugStr();
          std::cout << std::endl;
        }
      std::cout << std::endl;
  
    }

  ////////////////////////////////////////////////////////////////
  // gather to two-body JJJT
  ////////////////////////////////////////////////////////////////

  std::cout << "  two-body JJJT" << std::endl;

  // define space and operator containers
  basis::TwoBodySpaceJJJT two_body_jjjt_space(Nmax);
  std::array<basis::TwoBodySectorsJJJT,3> two_body_jjjt_component_sectors;
  std::array<basis::MatrixVector,3> two_body_jjjt_component_matrices;

  // construct gathered operator
  basis::GatherOperatorTwoBodyJJJTNToTwoBodyJJJT(
      operator_labels,
      two_body_jjjtn_space,two_body_jjjtn_component_sectors,two_body_jjjtn_component_matrices,
      two_body_jjjt_space,two_body_jjjt_component_sectors,two_body_jjjt_component_matrices
    );

  ////////////////////////////////////////////////////////////////
  // write as two-body JJJT
  ////////////////////////////////////////////////////////////////

  std::string two_body_jjjt_filename(jjjt_filename);
  std::ostringstream jjjt_sstream;
  for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
    {
      basis::WriteTwoBodyOperatorComponentJJJT(
          jjjt_sstream,
          T0,
          two_body_jjjt_component_sectors[T0],two_body_jjjt_component_matrices[T0],
          basis::NormalizationConversion::kNone
        );
    }
  std::ofstream jjjt_stream(two_body_jjjt_filename.c_str());
  jjjt_stream << jjjt_sstream.str();

}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void test_transform_timing(
  )
{

  std::cout << "test_transform_timing" << std::endl;

  ////////////////////////////////////////////////////////////////
  // define operator labels
  ////////////////////////////////////////////////////////////////

  // define operator properties
  int Nmax_relative = 26;
  int Jmax_relative = Nmax_relative+1;
  int Nmax = 10;  // target operator

  basis::OperatorLabelsJT operator_labels;
  operator_labels.J0 = 0;
  operator_labels.g0 = 0;
  operator_labels.T0_min = 0;
  operator_labels.T0_max = 0;
  operator_labels.symmetry_phase_mode = basis::SymmetryPhaseMode::kHermitian;

  std::cout << " Nmax_relative " << Nmax_relative << " Nmax " << Nmax << std::endl;

  ////////////////////////////////////////////////////////////////
  // construct relative identity operator
  ////////////////////////////////////////////////////////////////

  std::cout << "  relative LSJT" << std::endl;

  // define space and operator containers
  basis::RelativeSpaceLSJT relative_space(Nmax_relative,Jmax_relative);
  std::array<basis::RelativeSectorsLSJT,3> relative_component_sectors;
  std::array<basis::MatrixVector,3> relative_component_matrices;

  // do construction
  Timer relative_lsjt_timer;
  relative_lsjt_timer.Start();
  ConstructIdentityOperatorRelativeLSJT(
      operator_labels,
      relative_space,relative_component_sectors,relative_component_matrices
    );
  relative_lsjt_timer.Stop();
  std::cout << "Time: " << relative_lsjt_timer.ElapsedTime() << std::endl;

  ////////////////////////////////////////////////////////////////
  // augment to relative-cm LSJTN
  ////////////////////////////////////////////////////////////////

  std::cout << "  relative-cm LSJTN" << std::endl;

  // define space and operator containers
  basis::RelativeCMSpaceLSJTN relative_cm_lsjtn_space(Nmax);
  std::array<basis::RelativeCMSectorsLSJTN,3> relative_cm_lsjtn_component_sectors;
  std::array<basis::MatrixVector,3> relative_cm_lsjtn_component_matrices;

  // do transformation
  Timer relative_cm_lsjtn_timer;
  relative_cm_lsjtn_timer.Start();
  moshinsky::TransformOperatorRelativeLSJTToRelativeCMLSJTN(
      operator_labels,
      relative_space,relative_component_sectors,relative_component_matrices,
      relative_cm_lsjtn_space,relative_cm_lsjtn_component_sectors,relative_cm_lsjtn_component_matrices
    );
  relative_cm_lsjtn_timer.Stop();
  std::cout << "Time: " << relative_cm_lsjtn_timer.ElapsedTime() << std::endl;

  ////////////////////////////////////////////////////////////////
  // transform to two-body LSJTN
  ////////////////////////////////////////////////////////////////

  std::cout << "  two-body LSJTN" << std::endl;

  // define space and operator containers
  basis::TwoBodySpaceLSJTN two_body_lsjtn_space(Nmax);
  std::array<basis::TwoBodySectorsLSJTN,3> two_body_lsjtn_component_sectors;
  std::array<basis::MatrixVector,3> two_body_lsjtn_component_matrices;

  // do transformation
  Timer two_body_lsjtn_timer;
  two_body_lsjtn_timer.Start();
  moshinsky::TransformOperatorRelativeCMLSJTNToTwoBodyLSJTN(
      operator_labels,
      relative_cm_lsjtn_space,relative_cm_lsjtn_component_sectors,relative_cm_lsjtn_component_matrices,
      two_body_lsjtn_space,two_body_lsjtn_component_sectors,two_body_lsjtn_component_matrices
  );
  two_body_lsjtn_timer.Stop();
  std::cout << "Time: " << two_body_lsjtn_timer.ElapsedTime() << std::endl;

  ////////////////////////////////////////////////////////////////
  // recouple to two-body JJJTN
  ////////////////////////////////////////////////////////////////

  std::cout << "  two-body JJJTN" << std::endl;

  // define space and operator containers
  basis::TwoBodySpaceJJJTN two_body_jjjtn_space(Nmax);
  std::array<basis::TwoBodySectorsJJJTN,3> two_body_jjjtn_component_sectors;
  std::array<basis::MatrixVector,3> two_body_jjjtn_component_matrices;

  // do recoupling
  Timer two_body_jjjtn_timer;
  two_body_jjjtn_timer.Start();
  moshinsky::TransformOperatorTwoBodyLSJTNToTwoBodyJJJTN(
      operator_labels,
      two_body_lsjtn_space,two_body_lsjtn_component_sectors,two_body_lsjtn_component_matrices,
      two_body_jjjtn_space,two_body_jjjtn_component_sectors,two_body_jjjtn_component_matrices
    );
  two_body_jjjtn_timer.Stop();
  std::cout << "Time: " << two_body_jjjtn_timer.ElapsedTime() << std::endl;

  ////////////////////////////////////////////////////////////////
  // gather to two-body JJJT
  ////////////////////////////////////////////////////////////////

  std::cout << "  two-body JJJT" << std::endl;

  // define space and operator containers
  basis::TwoBodySpaceJJJT two_body_jjjt_space(Nmax);
  std::array<basis::TwoBodySectorsJJJT,3> two_body_jjjt_component_sectors;
  std::array<basis::MatrixVector,3> two_body_jjjt_component_matrices;

  // construct gathered operator
  Timer two_body_jjjt_timer;
  two_body_jjjt_timer.Start();
  basis::GatherOperatorTwoBodyJJJTNToTwoBodyJJJT(
      operator_labels,
      two_body_jjjtn_space,two_body_jjjtn_component_sectors,two_body_jjjtn_component_matrices,
      two_body_jjjt_space,two_body_jjjt_component_sectors,two_body_jjjt_component_matrices
    );
  two_body_jjjt_timer.Stop();
  std::cout << "Time: " << two_body_jjjt_timer.ElapsedTime() << std::endl;

}

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  // test_relative_cm();
  // test_moshinsky_matrix();
  // test_transform_simple(
  //     "test/moshinsky_xform_test_two_body_lsjt_identity_AS.dat",
  //     "test/moshinsky_xform_test_two_body_jjjt_identity_AS.dat",
  //     'I'
  //   );
  // test_transform_simple(
  //     "test/moshinsky_xform_test_two_body_lsjt_kinetic_AS.dat",
  //     "test/moshinsky_xform_test_two_body_jjjt_kinetic_AS.dat",
  //     'K'
  //   );
  test_transform_timing();

  // termination
  return 0;
}
