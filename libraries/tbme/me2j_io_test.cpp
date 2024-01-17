/****************************************************************
  me2j_io_test.cpp

  Zhou Zhou
  University of Notre Dame

****************************************************************/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

#include "fmt/format.h"
#include "mcutils/eigen.h"
#include "mcutils/profiling.h"
#include "tbme/me2j_io.h"
#include "basis/jjjttz_operator.h"

void TestMe2jio () {
  ////////////////////////////////////////////////////////////////
  // initialization
  ////////////////////////////////////////////////////////////////

  // header
  std::cout << std::endl;
  std::cout << "me2jstat  -- MFDn ME2J io tests" << std::endl;
  std::cout << std::endl;

  // set up parameters
  std::string input_filename = "TBME_TUD.int";
  int Nmax = 4;
  int J0 = 0;
  int g0 = 0;
  int Tz0 = 0;

  // initialize
  const basis::TwoBodySpaceJJJTTz space(basis::Rank::kTwoBody,Nmax);
  const basis::TwoBodySectorsJJJTTz sectors(space,J0,g0,Tz0);
  basis::OperatorBlocks<double> matrices;
  SetOperatorToZero(sectors,matrices);

  // read me2j file
  shell::ReadMe2jFile(space,sectors,matrices,input_filename);

  // write jjjttz format matrices
  std::ios_base::openmode mode_argument = std::ios_base::out;
  std::ofstream os("TBME_TUD.int_jjjttz.txt", mode_argument);
  basis::WriteTwoBodyOperatorJJJTTz(
    os,
    sectors,matrices,
    basis::NormalizationConversion::kNone
  );

  // write me2j file
  std::string output_filename = "TBME_TUD.int.txt";
  shell::WriteMe2jFile(space,sectors,matrices,output_filename);
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

void TestTzToMe2jToTz () {
  int Nmax = 4;
  int J0 = 0;
  int g0 = 0;
  int Tz0 = 0;
  basis::TwoBodySpaceJJJTTz two_body_jjjttz_space(basis::Rank::kTwoBody,Nmax);
  basis::TwoBodySectorsJJJTTz two_body_jjjttz_sectors(two_body_jjjttz_space,J0,g0,Tz0);
  basis::OperatorBlocks<double> two_body_jjjttz_matrices;

  SetOperatorToRandom(two_body_jjjttz_sectors,two_body_jjjttz_matrices);

  std::ios_base::openmode mode_argument = std::ios_base::out;
  std::ofstream os1("TestTzToMe2jToTz1.txt", mode_argument);
  basis::WriteTwoBodyOperatorJJJTTz(
    os1,
    two_body_jjjttz_sectors,two_body_jjjttz_matrices,
    basis::NormalizationConversion::kNone
  );

  std::string random_input_filename = "random_me2j.txt";
  shell::WriteMe2jFile(two_body_jjjttz_space,two_body_jjjttz_sectors,two_body_jjjttz_matrices,random_input_filename);

  basis::SetOperatorToZero(two_body_jjjttz_sectors,two_body_jjjttz_matrices);

  shell::ReadMe2jFile(two_body_jjjttz_space,two_body_jjjttz_sectors,two_body_jjjttz_matrices,random_input_filename);

  // write matrices
  std::ofstream os2("TestTzToMe2jToTz2.txt", mode_argument);
  basis::WriteTwoBodyOperatorJJJTTz(
    os2,
    two_body_jjjttz_sectors,two_body_jjjttz_matrices,
    basis::NormalizationConversion::kNone
  );

}

////////////////////////////////////////////////////////////////
// main program
/////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  // TestMe2jio();

  TestTzToMe2jToTz();

  // exit
  return EXIT_SUCCESS;
}
