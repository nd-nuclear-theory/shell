/****************************************************************
  jpv2rel.cpp

  Populate Hamiltonian-like operator from Iowa State ("JPV") format
  relative operator file.

  Initial support is for isoscalar operators, though the format has
  also been extended to non-isoscalar operators.

  See lsjt_operator.h for documentation of operator storage and the
  relative operator file format.

  Standard input:
    Nmax Jmax
    source_filename
    target_filename

  Language: C++11
                                 
  Mark A. Caprio
  University of Notre Dame

  7/25/16 (mac): Created, based upon writerel.cpp.

****************************************************************/

#include <fstream>

#include "basis/lsjt_operator.h"
#include "mcpp/parsing.h"
#include "relative/construct_relative.h"

////////////////////////////////////////////////////////////////
// parameter input
////////////////////////////////////////////////////////////////

struct Parameters
// Container for run input parameters.
{
  int Nmax, Jmax;
  std::string source_filename;
  std::string target_filename;
};

void ReadParameters(Parameters& parameters)
// Read run parameters from stdin.
//
// Arguments:
//   parameters (Parameters, output) :
//     container for input parameters
{

  // set up line counter for use in error messages
  std::string line;
  int line_count = 0;


  // line 1: truncation properties
  {
    ++line_count;
    std::getline(std::cin,line);
    std::istringstream line_stream(line);
    line_stream >> parameters.Nmax
                >> parameters.Jmax;
    ParsingCheck(line_stream,line_count,line);
  }

  // line 2: source filename
  {
    ++line_count;
    std::getline(std::cin,line);
    std::istringstream line_stream(line);
    line_stream >> parameters.source_filename;
    ParsingCheck(line_stream,line_count,line);
  }

  // line 3: target filename
  {
    ++line_count;
    std::getline(std::cin,line);
    std::istringstream line_stream(line);
    line_stream >> parameters.target_filename;
    ParsingCheck(line_stream,line_count,line);
  }

}

////////////////////////////////////////////////////////////////
// operator work
////////////////////////////////////////////////////////////////

void ReadJPVOperator(
    const std::string& source_filename,
    const basis::RelativeSpaceLSJT& relative_space,
    const basis::OperatorLabelsJT& operator_labels,
    const std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
    std::array<basis::MatrixVector,3>& relative_component_matrices
  )
// Define operator.
//
// Arguments:
//   source_filename (std::string) : input filename
//   relative_space (...) : target space
//   relative_component_sectors (...) : target sectors
//   relative_component_matrices (..., output) : target matrices to be populated
{

}

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  // parameter input
  Parameters parameters;
  ReadParameters(parameters);

  // set up zero operator
  std::cout << "Operator setup..." << std::endl;
  basis::RelativeSpaceLSJT relative_space(parameters.Nmax,parameters.Jmax);
  basis::OperatorLabelsJT operator_labels;
  operator_labels.J0 = 0;
  operator_labels.g0 = 0;
  operator_labels.T0_min = 0;
  operator_labels.T0_max = 0;
  operator_labels.symmetry_phase_mode = basis::SymmetryPhaseMode::kHermitian;
  std::array<basis::RelativeSectorsLSJT,3> relative_component_sectors;
  std::array<basis::MatrixVector,3> relative_component_matrices;
  relative::ConstructDiagonalConstantOperator(
      basis::RelativeOperatorParametersLSJT(operator_labels,parameters.Nmax,parameters.Jmax),
      relative_space,relative_component_sectors,relative_component_matrices,
      0.
    );

  // operator diagnostics
  std::cout << "  Truncation:"
            << " Nmax " << parameters.Nmax
            << " Jmax " << parameters.Jmax
            << std::endl;
  std::cout << "  Matrix elements:";
  for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
    std::cout << " " << basis::UpperTriangularEntries(relative_component_sectors[T0]);
  std::cout << std::endl;
  std::cout << "  Allocated:";
  for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
    std::cout << " " << basis::AllocatedEntries(relative_component_matrices[T0]);
  std::cout << std::endl;
        
  // populate matrix elements
  ReadJPVOperator(
      parameters.source_filename,
      relative_space,operator_labels,relative_component_sectors,relative_component_matrices
    );

  // write operator
  basis::WriteRelativeOperatorLSJT(
      parameters.target_filename,
      relative_space,
      operator_labels,  // only need operator labels
      relative_component_sectors,
      relative_component_matrices,
      true  // verbose
    );

  // termination
  return 0;
}
