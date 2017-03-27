/****************************************************************
  jpv2rel.cpp

  Populate Hamiltonian-like operator from Iowa State ("JPV") format
  relative operator file.

  See lsjt_operator.h for documentation of operator storage and the
  relative operator file format.

  Notes on the interpretation of JPV relative files are provided in
  comments below.

  The basic format as considered here provides support only for
  isoscalar operators, though the format has also been extended to
  non-isoscalar operators through the use of separate "pp", "nn", and
  "pn" files.

  Standard input:
    Nmax Jmax
    source_filename
    target_filename

  Language: C++11
                                 
  Mark A. Caprio
  University of Notre Dame

  7/25/16 (mac): Created, based upon writerel.cpp.
  8/10/16 (mac): Fix input indexing.
  10/9/16 (pjf): Rename mcpp -> mcutils.
  10/19/16 (mac): Remove superflous debugging options.

****************************************************************/

#include <fstream>

#include "basis/lsjt_operator.h"
#include "cppformat/format.h"
#include "mcutils/parsing.h"
#include "relative/construct_relative.h"
#include "relative/jpv_io.h"

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
  relative::ReadJPVOperator(
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
