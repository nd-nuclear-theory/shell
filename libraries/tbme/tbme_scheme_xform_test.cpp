/****************************************************************
  tbme_scheme_xform_test.cpp

  Zhou Zhou
  University of Notre Dame
****************************************************************/

#include <fstream>

#include "mcutils/eigen.h"
#include "mcutils/parsing.h"
#include "mcutils/profiling.h"
#include "tbme/me2j_io.h"
#include "tbme/tbme_scheme_xform.h"

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

void TestTToTz () {
  int Nmax = 4;
  int J0 = 0;
  int g0 = 0;
  int Tz0 = 0;
  const basis::TwoBodySpaceJJJT two_body_jjjt_space(basis::Rank::kTwoBody,Nmax);
  std::array<basis::TwoBodySectorsJJJT,3> two_body_jjjt_component_sectors;
  std::array<basis::OperatorBlocks<double>,3> two_body_jjjt_component_matrices;
  for (int T0 = 0; T0 <= 2; T0++) {
    two_body_jjjt_component_sectors[T0] = basis::TwoBodySectorsJJJT(two_body_jjjt_space,J0,T0,g0);
    if (T0 == 0) { // close but not identity operator
      basis::SetOperatorToIdentity(two_body_jjjt_component_sectors[T0],two_body_jjjt_component_matrices[T0]);
    } else {
      basis::SetOperatorToZero(two_body_jjjt_component_sectors[T0],two_body_jjjt_component_matrices[T0]);
    }
  }
  const basis::TwoBodySpaceJJJTTz two_body_jjjttz_space(basis::Rank::kTwoBody,Nmax);
  basis::TwoBodySectorsJJJTTz two_body_jjjttz_sectors;
  basis::OperatorBlocks<double> two_body_jjjttz_matrices;
  shell::TransformOperatorTwoBodyJJJTToTwoBodyJJJTTz(
    two_body_jjjt_space,two_body_jjjt_component_sectors,two_body_jjjt_component_matrices,
    two_body_jjjttz_space,two_body_jjjttz_sectors,two_body_jjjttz_matrices
  );
  // write matrices
  std::ios_base::openmode mode_argument = std::ios_base::out;
  std::ofstream os("TtoTztest.txt", mode_argument);
  basis::WriteTwoBodyOperatorJJJTTz(
    os,
    two_body_jjjttz_sectors,two_body_jjjttz_matrices,
    basis::NormalizationConversion::kNone
  );
}

void TestTzToT () {
  int Nmax = 4;
  int J0 = 0;
  int g0 = 0;
  int Tz0 = 0;
  const basis::TwoBodySpaceJJJTTz two_body_jjjttz_space(basis::Rank::kTwoBody,Nmax);
  const basis::TwoBodySectorsJJJTTz two_body_jjjttz_sectors(two_body_jjjttz_space,J0,g0,Tz0);
  basis::OperatorBlocks<double> two_body_jjjttz_matrices;
  basis::SetOperatorToIdentity(two_body_jjjttz_sectors,two_body_jjjttz_matrices);
  const basis::TwoBodySpaceJJJT two_body_jjjt_space(basis::Rank::kTwoBody,Nmax);
  std::array<basis::TwoBodySectorsJJJT,3> two_body_jjjt_component_sectors;
  std::array<basis::OperatorBlocks<double>,3> two_body_jjjt_component_matrices;
  shell::TransformOperatorTwoBodyJJJTTzToTwoBodyJJJT(two_body_jjjttz_space,two_body_jjjttz_sectors,two_body_jjjttz_matrices,
    two_body_jjjt_space,two_body_jjjt_component_sectors,two_body_jjjt_component_matrices
    );
  // write matrices
  std::ios_base::openmode mode_argument = std::ios_base::out;
  std::ofstream os("TztoTtest.txt", mode_argument);
  int T0 = 0;
  WriteTwoBodyOperatorComponentJJJT(
      os,
      T0,
      two_body_jjjt_component_sectors[T0],
      two_body_jjjt_component_matrices[T0],
      basis::NormalizationConversion::kNone
    );
}

template <typename tSectorsType, typename tFloat>
  void SetOperatorToRandom(
      const tSectorsType& sectors,
      basis::OperatorBlocks<tFloat>& matrices
    )
    // Set operator blocks to zero.
    //
    // Arguments:
    //   sectors (input): the set of sectors on which the operator is
    //     defined
    //   matrices (output): matrices to hold blocks
  {

    // clear vector of matrices
    matrices.clear();
    matrices.resize(sectors.size());

    // iterate over sectors
    for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
      {

        // extract sector
        const typename tSectorsType::SectorType& sector = sectors.GetSector(sector_index);

        // extract sector subspaces
        const typename tSectorsType::SubspaceType& bra_subspace = sector.bra_subspace();
        const typename tSectorsType::SubspaceType& ket_subspace = sector.ket_subspace();

        // generate matrix for sector
        matrices[sector_index] = basis::OperatorBlock<tFloat>::Random(bra_subspace.dimension(),ket_subspace.dimension());

      }
  }

void TestTzToTToTz () {
  int Nmax = 4;
  int J0 = 0;
  int g0 = 0;
  int Tz0 = 0;
  basis::TwoBodySpaceJJJTTz two_body_jjjttz_space(basis::Rank::kTwoBody,Nmax);
  basis::TwoBodySectorsJJJTTz two_body_jjjttz_sectors(two_body_jjjttz_space,J0,g0,Tz0);
  basis::OperatorBlocks<double> two_body_jjjttz_matrices;
  basis::SetOperatorToZero(two_body_jjjttz_sectors,two_body_jjjttz_matrices);
  // std::string input_filename = "vrelccm_me2j.txt";
  // read me2j file
  // shell::ReadMe2jFile(two_body_jjjttz_space,two_body_jjjttz_sectors,two_body_jjjttz_matrices,input_filename);
  SetOperatorToRandom(two_body_jjjttz_sectors,two_body_jjjttz_matrices);
  std::string random_input_filename = "random_me2j.txt";
  shell::WriteMe2jFile(two_body_jjjttz_space,two_body_jjjttz_sectors,two_body_jjjttz_matrices,random_input_filename);
  basis::TwoBodySpaceJJJT two_body_jjjt_space(basis::Rank::kTwoBody,Nmax);
  std::array<basis::TwoBodySectorsJJJT,3> two_body_jjjt_component_sectors;
  std::array<basis::OperatorBlocks<double>,3> two_body_jjjt_component_matrices;
  // write matrices
  std::ios_base::openmode mode_argument = std::ios_base::out;
  std::ofstream os1("TestTzToTToTz1.txt", mode_argument);
  basis::WriteTwoBodyOperatorJJJTTz(
    os1,
    two_body_jjjttz_sectors,two_body_jjjttz_matrices,
    basis::NormalizationConversion::kNone
  );
  shell::TransformOperatorTwoBodyJJJTTzToTwoBodyJJJT(two_body_jjjttz_space,two_body_jjjttz_sectors,two_body_jjjttz_matrices,
    two_body_jjjt_space,two_body_jjjt_component_sectors,two_body_jjjt_component_matrices
    );
  basis::SetOperatorToZero(two_body_jjjttz_sectors,two_body_jjjttz_matrices);
  shell::TransformOperatorTwoBodyJJJTToTwoBodyJJJTTz(
    two_body_jjjt_space,two_body_jjjt_component_sectors,two_body_jjjt_component_matrices,
    two_body_jjjttz_space,two_body_jjjttz_sectors,two_body_jjjttz_matrices
  );
  // write matrices
  std::ofstream os2("TestTzToTToTz2.txt", mode_argument);
  basis::WriteTwoBodyOperatorJJJTTz(
    os2,
    two_body_jjjttz_sectors,two_body_jjjttz_matrices,
    basis::NormalizationConversion::kNone
  );
  // std::string output_filename = "vrelccm_Tz2T2Tz_me2j.txt";
  std::string output_filename = "random_Tz2T2Tz_me2j.txt";
  shell::WriteMe2jFile(two_body_jjjttz_space,two_body_jjjttz_sectors,two_body_jjjttz_matrices,output_filename);
}

void TestTToTzToT () {

  int Nmax = 4;
  int J0 = 0;
  int g0 = 0;
  int Tz0 = 0;

  basis::TwoBodySpaceJJJT two_body_jjjt_space(basis::Rank::kTwoBody,Nmax);
  std::array<basis::TwoBodySectorsJJJT,3> two_body_jjjt_component_sectors;
  std::array<basis::OperatorBlocks<double>,3> two_body_jjjt_component_matrices;

  for (size_t T0 = 0; T0 < 3; T0++) {
    two_body_jjjt_component_sectors[T0] = basis::TwoBodySectorsJJJT(two_body_jjjt_space,J0,T0,g0);
    SetOperatorToRandom(two_body_jjjt_component_sectors[T0],two_body_jjjt_component_matrices[T0]);
  }

  // write matrices
  std::ios_base::openmode mode_argument = std::ios_base::out;
  std::ofstream os1("TestTToTzToT1.txt", mode_argument);
  for (size_t T0 = 0; T0 < 3; T0++) {
    basis::WriteTwoBodyOperatorComponentJJJT(
      os1, T0,
      two_body_jjjt_component_sectors[T0],two_body_jjjt_component_matrices[T0],
      basis::NormalizationConversion::kNone
    );
  }

  const basis::TwoBodySpaceJJJTTz two_body_jjjttz_space(basis::Rank::kTwoBody,Nmax);
  basis::TwoBodySectorsJJJTTz two_body_jjjttz_sectors;
  basis::OperatorBlocks<double> two_body_jjjttz_matrices;

  basis::SetOperatorToZero(two_body_jjjttz_sectors,two_body_jjjttz_matrices);

  shell::TransformOperatorTwoBodyJJJTToTwoBodyJJJTTz(
    two_body_jjjt_space,two_body_jjjt_component_sectors,two_body_jjjt_component_matrices,
    two_body_jjjttz_space,two_body_jjjttz_sectors,two_body_jjjttz_matrices
  );
  for (size_t T0 = 0; T0 < 3; T0++) {
    basis::SetOperatorToZero(two_body_jjjt_component_sectors[T0],two_body_jjjt_component_matrices[T0]);
  }
  shell::TransformOperatorTwoBodyJJJTTzToTwoBodyJJJT(two_body_jjjttz_space,two_body_jjjttz_sectors,two_body_jjjttz_matrices,
    two_body_jjjt_space,two_body_jjjt_component_sectors,two_body_jjjt_component_matrices
    );

  // write matrices
  std::ofstream os2("TestTToTzToT2.txt", mode_argument);
  for (size_t T0 = 0; T0 < 3; T0++) {
    basis::WriteTwoBodyOperatorComponentJJJT(
      os2, T0,
      two_body_jjjt_component_sectors[T0],two_body_jjjt_component_matrices[T0],
      basis::NormalizationConversion::kNone
    );
  }
}

int main(int argc, char **argv)
{
  TestTToTz();
  TestTzToT();
  TestTzToTToTz(); // see if the output is the same as the input
  TestTToTzToT();
}
