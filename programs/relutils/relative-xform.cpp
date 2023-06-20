/****************************************************************
  relative-xform.cpp

  Transform (dilate) relative matrix elements.

  See lsjt_operator.h for documentation of operator storage and the
  relative operator file format.

  Standard input:
    source_filename
    target_Nmax b_ratio num_steps
    target_filename

  The parameters are

    source_filename: input relative operator file

    b_ratio: ratio of target oscillator length parameter to source oscillator
      length parameter

    num_steps: number of steps for numerical integration of overlaps

    target_Nmax: target relative basis trunctation

  Other operator parameter values are taken from the source file.

  Example:
    Daejeon16_Nmax40_hw25.0_rel.dat
    40 1.5811388300841898 4000
    Daejeon16_Nmax40_hw10.0_rel.dat

  Language: C++11

  Patrick J. Fasano
  University of Notre Dame

  + 10/13/20 (mac): Created.

****************************************************************/

#include "basis/lsjt_operator.h"
#include "mcutils/parsing.h"
#include "mcutils/profiling.h"
#include "relative/relative_xform.h"

////////////////////////////////////////////////////////////////
// parameter input
////////////////////////////////////////////////////////////////

struct Parameters
// Container for run input parameters.
{
  std::string source_filename;
  int Nmax;
  double b_ratio;
  int num_steps;
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


  // line 1: source filename
  {
    mcutils::GetLine(std::cin, line, line_count);
    std::istringstream line_stream(line);
    line_stream >> parameters.source_filename;
    mcutils::ParsingCheck(line_stream,line_count,line);
    mcutils::FileExistCheck(parameters.source_filename, true, false);
  }

  // line 2: transformation parameters
  {
    mcutils::GetLine(std::cin, line, line_count);
    std::istringstream line_stream(line);
    line_stream >> parameters.Nmax >> parameters.b_ratio >> parameters.num_steps;
    mcutils::ParsingCheck(line_stream,line_count,line);
  }

  // line 3: target filename
  {
    mcutils::GetLine(std::cin, line, line_count);
    std::istringstream line_stream(line);
    line_stream >> parameters.target_filename;
    mcutils::ParsingCheck(line_stream,line_count,line);
    mcutils::FileExistCheck(parameters.target_filename, false, true);
  }

}

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  // header
  std::cout << std::endl;
  std::cout << "relative-xform -- relative operator matrix transformation" << std::endl;
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  // parameter input
  Parameters parameters;
  ReadParameters(parameters);

  // set up source operator
  basis::RelativeSpaceLSJT source_relative_space;
  basis::RelativeOperatorParametersLSJT source_operator_parameters;
  std::array<basis::RelativeSectorsLSJT,3> source_relative_component_sectors;
  std::array<basis::OperatorBlocks<double>,3> source_relative_component_blocks;

  mcutils::SteadyTimer read_relative_lsjt_timer;
  read_relative_lsjt_timer.Start();
  basis::ReadRelativeOperatorLSJT(
      parameters.source_filename,
      source_relative_space,
      source_operator_parameters,
      source_relative_component_sectors,
      source_relative_component_blocks,
      true  // verbose
    );
  read_relative_lsjt_timer.Stop();
  std::cout << "  Time: " << read_relative_lsjt_timer.ElapsedTime() << std::endl;

  // set up target operator
  basis::RelativeSpaceLSJT target_relative_space;
  basis::RelativeOperatorParametersLSJT target_operator_parameters;
  std::array<basis::RelativeSectorsLSJT,3> target_relative_component_sectors;
  std::array<basis::OperatorBlocks<double>,3> target_relative_component_blocks;

  // copy operator parameters and override Nmax
  target_operator_parameters = source_operator_parameters;
  target_operator_parameters.Nmax = parameters.Nmax;

  // construct target space
  target_relative_space = basis::RelativeSpaceLSJT(
      target_operator_parameters.Nmax, target_operator_parameters.Jmax
    );

  // transform source operator to target operator
  mcutils::SteadyTimer transform_relative_lsjt_timer;
  transform_relative_lsjt_timer.Start();
  relative::TransformRelativeOperatorLSJT(
      source_operator_parameters,
      source_relative_space,
      source_relative_component_sectors,
      source_relative_component_blocks,
      target_operator_parameters,
      target_relative_space,
      target_relative_component_sectors,
      target_relative_component_blocks,
      parameters.b_ratio, parameters.num_steps,
      true  // verbose
    );
  transform_relative_lsjt_timer.Stop();
  std::cout << "  Time: " << transform_relative_lsjt_timer.ElapsedTime() << std::endl;

  // write operator
  mcutils::SteadyTimer write_relative_lsjt_timer;
  write_relative_lsjt_timer.Start();
  basis::WriteRelativeOperatorLSJT(
      parameters.target_filename,
      target_relative_space,
      target_operator_parameters,  // only need operator labels
      target_relative_component_sectors,
      target_relative_component_blocks,
      true  // verbose
    );
  write_relative_lsjt_timer.Stop();
  std::cout << "  Time: " << write_relative_lsjt_timer.ElapsedTime() << std::endl;

  // termination
  return 0;
}
