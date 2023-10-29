/****************************************************************
  pn2rel.cpp

  Populate Hamiltonian-like operator from Petr Navratil ("PN") format
  relative operator file.

  See lsjt_operator.h for documentation of operator storage and the
  relative operator file format.

  commandline arguments: Nmax Jmax source_filename output_filename

  Language: C++11

  Anna E. McCoy
  TRIUMF

  5/2/18 (aem): Created, based upon pn2rel.cpp in SpNCCI.
****************************************************************/

#include "fmt/format.h"
#include "mcutils/parsing.h"
#include "relative/pn_io.h"

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////




int main(int argc, char **argv)
{
  // header
  if(argc<3)
    {
      std::cout<<"Syntax: Nmax Jmax source_filename output_filename"<<std::endl;
      std::exit(EXIT_FAILURE);
    }

  int Nmax=std::stoi(argv[1]);
  int Jmax=std::stoi(argv[2]);
  std::string source_filename=argv[3];
  std::string output_filename=argv[4];

  std::cout << std::endl;
  std::cout << "pn2rel -- PN to relative file conversion" << std::endl;
  std::cout << std::endl;

  // set up zero operator
  basis::RelativeSpaceLSJT relative_space(Nmax,Jmax);
  basis::OperatorLabelsJT operator_labels(0,0,0,2,basis::SymmetryPhaseMode::kHermitian);
  basis::RelativeOperatorParametersLSJT operator_parameters(operator_labels,Nmax,Jmax);
  std::array<basis::RelativeSectorsLSJT,3> relative_component_sectors;
  std::array<basis::OperatorBlocks<double>,3> relative_component_matrices;
  basis::ConstructZeroOperatorRelativeLSJT(
      operator_parameters,
      relative_space,relative_component_sectors,relative_component_matrices
    );

  // Read in Petr Navratil formatted matrix elements 
  relative::ReadPNOperatorPN(
    source_filename,
    relative_space,
    operator_parameters,
    relative_component_sectors,
    relative_component_matrices,
    false  // verbose
    );

  // write operator
  basis::WriteRelativeOperatorLSJT(
      output_filename,
      relative_space,
      operator_labels,
      relative_component_sectors,
      relative_component_matrices,
      true  // verbose
    );

  // termination
  return 0;
}
