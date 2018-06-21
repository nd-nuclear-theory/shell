/******************************************************************************

  jpvstat.cpp -- scan sector headers for JPV relative file

  LIMITATION: currently subject to hard-coded maximum Nmax

  Syntax:
    jpvstat input_filename

  Mark A. Caprio
  University of Notre Dame

  3/28/17 (mac): Created, from structure of h2stat and code from jpv2rel.
  11/28/17 (pjf): Include version in header.

******************************************************************************/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

#include "fmt/format.h"
#include "mcutils/parsing.h"
#include "relative/jpv_io.h"

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

struct RunParameters
// Stores simple parameters for run
{
  // filenames
  std::string input_filename;
};

void ProcessArguments(int argc, char **argv, RunParameters& run_parameters)
{
  // usage message
  if ((argc-1 < 1) || (std::string(argv[1])=="--help"))
    {
      std::cout
        << "Syntax: jpvstat input_filename" << std::endl
        << std::endl
        << std::endl;

      std::exit(EXIT_SUCCESS);
    }

  if (argc-1 == 1)
    // just filename
    {
      run_parameters.input_filename = argv[1];
    }
}


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
  std::cout << "jpvstat  -- JPV relative file sector summary" << std::endl;
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  // read parameters
  RunParameters run_parameters;
  ProcessArguments(argc,argv,run_parameters);

  ////////////////////////////////////////////////////////////////
  // process file
  ////////////////////////////////////////////////////////////////

  // set up dummy empty operator
  //
  // ReadJPVOperator requires an output container into which to
  // nominally put matrix elements, though this one will be of null
  // size, and thus no matrix elements will be stored.

  // Problem: JPV input routine currently requires dimension info from
  // target indexing.  So Nmax must be chosen large enough to yield
  // all LSJT subspaces found in the input file.

  int Nmax=99;  // ad hoc
  int Jmax=99;  // to avoid failing assertion that sector J < Jmax
  basis::RelativeSpaceLSJT relative_space(Nmax,Jmax);
  basis::OperatorLabelsJT operator_labels(0,0,0,0,basis::SymmetryPhaseMode::kHermitian);
  std::array<basis::RelativeSectorsLSJT,3> relative_component_sectors;
  std::array<basis::OperatorBlocks<double>,3> relative_component_matrices;
  basis::ConstructZeroOperatorRelativeLSJT(
      basis::RelativeOperatorParametersLSJT(operator_labels,Nmax,Jmax),
      relative_space,relative_component_sectors,relative_component_matrices
    );

  // parse JPV file
  const basis::OperatorLabelsJT operator_labels_isoscalar(0,0,0,0,basis::SymmetryPhaseMode::kHermitian);
  relative::ReadJPVOperator(
      run_parameters.input_filename,
      relative_space,
      operator_labels_isoscalar,
      relative_component_sectors,
      relative_component_matrices,
      true  // verbose
    );


  ////////////////////////////////////////////////////////////////
  // termination
  ////////////////////////////////////////////////////////////////

  // exit
  return EXIT_SUCCESS;
}
