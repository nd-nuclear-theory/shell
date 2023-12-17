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

////////////////////////////////////////////////////////////////
// main program
/////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

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

  const basis::TwoBodySpaceJJJTTz space(basis::Rank::kTwoBody,Nmax);
  const basis::TwoBodySectorsJJJTTz sectors(space,J0,g0,Tz0);
  basis::OperatorBlocks<double> matrices;
  SetOperatorToZero(sectors,matrices);
  shell::ReadMe2jFile(space,sectors,matrices,input_filename);
  // write matrices
  std::ios_base::openmode mode_argument = std::ios_base::out;
  std::ofstream os("me2j_io_test_input.txt", mode_argument);
  basis::WriteTwoBodyOperatorJJJTTz(
    os,
    sectors,matrices,
    basis::NormalizationConversion::kNone
  );
  std::string output_filename = "output.txt";
  shell::WriteMe2jFile(space,sectors,matrices,output_filename);

  // exit
  return EXIT_SUCCESS;
}
