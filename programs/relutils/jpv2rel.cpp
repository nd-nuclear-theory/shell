/****************************************************************
  jpv2rel.cpp

  Populate Hamiltonian-like operator from Iowa State ("JPV") format
  relative operator file.

  See lsjt_operator.h for documentation of operator storage and the
  relative operator file format.

  Notes on the interpretation of JPV relative files are provided in
  comments below.

  The basic format provides support only for isoscalar operators, but
  the format has been extended to non-isoscalar operators through the
  use of separate "pp", "nn", and "pn" files.

  Standard input:
    Nmax Jmax T0_max
    source_filename  OR  source_filename_pp source_filename_nn source_filename_pn
    target_filename

    T0_max: 0 for isoscalar, 2 for pp/nn/pn matrix elements

  Language: C++11

  Mark A. Caprio
  University of Notre Dame

  7/25/16 (mac): Created, based upon writerel.cpp.
  8/10/16 (mac): Fix input indexing.
  10/9/16 (pjf): Rename mcpp -> mcutils.
  10/19/16 (mac): Remove superflous debugging options.
  3/28/17 (mac): Add support for non-isoscalar operators.
  11/28/17 (pjf): Print header with version.

****************************************************************/

#include "fmt/format.h"
#include "mcutils/parsing.h"
#include "relative/jpv_io.h"

////////////////////////////////////////////////////////////////
// parameter input
////////////////////////////////////////////////////////////////

struct Parameters
// Container for run input parameters.
{
  int Nmax, Jmax;
  int T0_max;  // 0 for isoscalar; 2 for pn
  std::string source_filename;  // for isoscalar
  std::array<std::string,3> source_filenames;  // for pn
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
                >> parameters.Jmax
                >> parameters.T0_max;
    mcutils::ParsingCheck(line_stream,line_count,line);
    if (!((parameters.T0_max==0)||(parameters.T0_max==2)))
        mcutils::ParsingError(line_count,line,"Invalid T0_max");
  }

  // line 2: source filename(s)
  {
    ++line_count;
    std::getline(std::cin,line);
    std::istringstream line_stream(line);
    if (parameters.T0_max==0)
      line_stream >> parameters.source_filename;
    else
      line_stream
        >> parameters.source_filenames[0]
        >> parameters.source_filenames[1]
        >> parameters.source_filenames[2];
    mcutils::ParsingCheck(line_stream,line_count,line);
  }

  // line 3: target filename
  {
    ++line_count;
    std::getline(std::cin,line);
    std::istringstream line_stream(line);
    line_stream >> parameters.target_filename;
    mcutils::ParsingCheck(line_stream,line_count,line);
  }

}

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  // header
  std::cout << std::endl;
  std::cout << "jpv2rel -- JPV to relative file conversion" << std::endl;
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  // parameter input
  Parameters parameters;
  ReadParameters(parameters);

  // set up zero operator
  std::cout << "Operator setup..." << std::endl;
  basis::RelativeSpaceLSJT relative_space(parameters.Nmax,parameters.Jmax);
  basis::OperatorLabelsJT operator_labels(0,0,0,parameters.T0_max,basis::SymmetryPhaseMode::kHermitian);
  basis::RelativeOperatorParametersLSJT operator_parameters(operator_labels,parameters.Nmax,parameters.Jmax);
  std::array<basis::RelativeSectorsLSJT,3> relative_component_sectors;
  std::array<basis::OperatorBlocks<double>,3> relative_component_blocks;
  basis::ConstructZeroOperatorRelativeLSJT(
      operator_parameters,
      relative_space,relative_component_sectors,relative_component_blocks
    );

  // operator diagnostics
  std::cout << "  Truncation:"
            << " Nmax " << parameters.Nmax
            << " Jmax " << parameters.Jmax
            << " T0_max " << parameters.T0_max
            << std::endl;
  std::cout << "  Matrix elements:";
  for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
    std::cout << " " << basis::UpperTriangularEntries(relative_component_sectors[T0]);
  std::cout << std::endl;
  std::cout << "  Allocated:";
  for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
    std::cout << " " << basis::AllocatedEntries(relative_component_blocks[T0]);
  std::cout << std::endl;

  // populate matrix elements
  if (parameters.T0_max==0)
    {
      relative::ReadJPVOperator(
          parameters.source_filename,
          relative_space,operator_labels,relative_component_sectors,relative_component_blocks,
          false  // verbose
        );
    }
  else
    {
      relative::ReadJPVOperatorPN(
          parameters.source_filenames,
          relative_space,operator_parameters,relative_component_sectors,relative_component_blocks,
          false  // verbose
        );
    }

  // write operator
  basis::WriteRelativeOperatorLSJT(
      parameters.target_filename,
      relative_space,
      operator_labels,
      relative_component_sectors,
      relative_component_blocks,
      true  // verbose
    );

  // termination
  return 0;
}
