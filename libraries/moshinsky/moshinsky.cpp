/****************************************************************
  moshinsky.cpp

  Perform Moshinsky transformation of general relative operator.

  See lsjt_operator.h for documentation of operator storage and the
  relative operator file format.

  Standard input:
    ob|tb truncation_cutoff coupling
    relative_filename
    two_body_filename

    N1max, N2max : one-body and two-body basis truncations

    coupling : target coupling scheme (lsjt, ...)

  Language: C++11

  Mark A. Caprio
  University of Notre Dame

  7/16/16 (mac): Created.

****************************************************************/

#include <fstream>

#include "am/wigner_gsl.h"
#include "basis/lsjt_operator.h"
#include "basis/jjjt_operator.h"
#include "mcpp/parsing.h"
#include "mcpp/profiling.h"

#include "moshinsky/moshinsky_xform.h"


////////////////////////////////////////////////////////////////
// parameter input
////////////////////////////////////////////////////////////////

struct Parameters
// Container for run input parameters.
{
  int truncation_rank, truncation_cutoff;
  std::string coupling;
  std::string relative_filename;
  std::string two_body_filename;
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


  // line 1: target basis truncation and coupling
  {
    ++line_count;
    std::getline(std::cin,line);
    std::istringstream line_stream(line);
    std::string truncation_rank_code;
    line_stream >> truncation_rank_code
                >> parameters.truncation_cutoff
                >> parameters.coupling;
    ParsingCheck(line_stream,line_count,line);

    // process truncation rank code
    if (truncation_rank_code == "ob")
      parameters.truncation_rank = 1;
    else if (truncation_rank_code == "tb")
      parameters.truncation_rank = 2;
    else
    {
      std::cerr << "ERROR: unrecognized truncation rank code" << std::endl;
      std::exit(EXIT_FAILURE);
    }

  }

  // line 2: relative (source) filename
  {
    ++line_count;
    std::getline(std::cin,line);
    std::istringstream line_stream(line);
    line_stream >> parameters.relative_filename;
    ParsingCheck(line_stream,line_count,line);
  }

  // line 3: two-body (target) filename
  {
    ++line_count;
    std::getline(std::cin,line);
    std::istringstream line_stream(line);
    line_stream >> parameters.two_body_filename;
    ParsingCheck(line_stream,line_count,line);
  }
}

////////////////////////////////////////////////////////////////
// relative operator input
////////////////////////////////////////////////////////////////

void ReadRelative(
    const std::string& relative_filename,
    basis::RelativeSpaceLSJT& relative_space,
    basis::OperatorLabelsJT& operator_labels,
    std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
    std::array<basis::MatrixVector,3>& relative_component_matrices
  )
// Set up and read relative operator.
//
// Arguments:
//   parameters (Parameters) : includes tensorial properties of operator
//      choice of operator to use 
//   relative_space (..., output) : target space
//   relative_component_sectors (..., output) : target sectors
//   relative_component_matrices (..., output) : target matrices
{
  std::cout << "Read relative operator..." << std::endl;

  // set up stream for readback
  std::ifstream is(relative_filename.c_str());

  // read header parameters
  basis::RelativeOperatorParametersLSJT operator_parameters;
  basis::ReadRelativeOperatorParametersLSJT(is,operator_parameters);
  operator_labels = static_cast<basis::OperatorLabelsJT>(operator_parameters);
  std::cout
    << " J0 " << operator_parameters.J0 << " g0 " << operator_parameters.g0
    << " T0 " << operator_parameters.T0_min << " .. " << operator_parameters.T0_max
    << " symmetry " << int(operator_parameters.symmetry_phase_mode) << std::endl
    << " Nmax " << operator_parameters.Nmax << " Jmax " << operator_parameters.Jmax << std::endl;

  // set up relative space
  relative_space = basis::RelativeSpaceLSJT(
      operator_parameters.Nmax,operator_parameters.Jmax
    );
  
  // populate sectors and matrices
  for (int T0=operator_parameters.T0_min; T0<=operator_parameters.T0_max; ++T0)
    // for each isospin component
    {
      // enumerate sectors
      relative_component_sectors[T0]
        = basis::RelativeSectorsLSJT(relative_space,operator_labels.J0,T0,operator_labels.g0);

      // read matrices
      basis::ReadRelativeOperatorComponentLSJT(
          is,
          T0,
          relative_component_sectors[T0],relative_component_matrices[T0]
        );
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

  // set up operator
  basis::RelativeSpaceLSJT relative_space;
  basis::OperatorLabelsJT operator_labels;
  std::array<basis::RelativeSectorsLSJT,3> relative_component_sectors;
  std::array<basis::MatrixVector,3> relative_component_matrices;
  ReadRelative(
      parameters.relative_filename,
      relative_space,
      operator_labels,relative_component_sectors,relative_component_matrices
    );

  // AD HOC: following is initial working code for two-body lsjt, to
  // be refactored and extended to jjjt and jjjpn

  // int N1max = parameters.N1max; // TODO
  int Nmax = parameters.truncation_cutoff;

  ////////////////////////////////////////////////////////////////
  // augment to relative-cm NLSJT
  ////////////////////////////////////////////////////////////////

  std::cout << "Augment to relative-cm NLSJT..." << std::endl;

  // define space and operator containers
  basis::RelativeCMSpaceNLSJT relative_cm_nlsjt_space(Nmax);
  std::array<basis::RelativeCMSectorsNLSJT,3> relative_cm_nlsjt_component_sectors;
  std::array<basis::MatrixVector,3> relative_cm_nlsjt_component_matrices;

  // do transformation
  Timer relative_cm_nlsjt_timer;
  relative_cm_nlsjt_timer.Start();
  moshinsky::TransformOperatorRelativeLSJTToRelativeCMNLSJT(
      operator_labels,
      relative_space,relative_component_sectors,relative_component_matrices,
      relative_cm_nlsjt_space,relative_cm_nlsjt_component_sectors,relative_cm_nlsjt_component_matrices
    );
  relative_cm_nlsjt_timer.Stop();
  std::cout << "  Time: " << relative_cm_nlsjt_timer.ElapsedTime() << std::endl;

  ////////////////////////////////////////////////////////////////
  // transform to two-body NLSJT
  ////////////////////////////////////////////////////////////////

  std::cout << "Transform to two-body NLSJT..." << std::endl;

  // define space and operator containers
  basis::TwoBodySpaceNLSJT two_body_nlsjt_space(Nmax);
  std::array<basis::TwoBodySectorsNLSJT,3> two_body_nlsjt_component_sectors;
  std::array<basis::MatrixVector,3> two_body_nlsjt_component_matrices;

  // do transformation
  Timer two_body_nlsjt_timer;
  two_body_nlsjt_timer.Start();
  moshinsky::TransformOperatorRelativeCMNLSJTToTwoBodyNLSJT(
      operator_labels,
      relative_cm_nlsjt_space,relative_cm_nlsjt_component_sectors,relative_cm_nlsjt_component_matrices,
      two_body_nlsjt_space,two_body_nlsjt_component_sectors,two_body_nlsjt_component_matrices
  );
  two_body_nlsjt_timer.Stop();
  std::cout << "  Time: " << two_body_nlsjt_timer.ElapsedTime() << std::endl;

  ////////////////////////////////////////////////////////////////
  // gather to two-body LSJT
  ////////////////////////////////////////////////////////////////

  std::cout << "Gather two-body LSJT..." << std::endl;

  // define space and operator containers
  basis::TwoBodySpaceLSJT two_body_lsjt_space(Nmax);
  std::array<basis::TwoBodySectorsLSJT,3> two_body_lsjt_component_sectors;
  std::array<basis::MatrixVector,3> two_body_lsjt_component_matrices;

  // construct gathered operator
  Timer two_body_lsjt_timer;
  two_body_lsjt_timer.Start();
  basis::GatherOperatorTwoBodyNLSJTToTwoBodyLSJT(
      operator_labels,
      two_body_nlsjt_space,two_body_nlsjt_component_sectors,two_body_nlsjt_component_matrices,
      two_body_lsjt_space,two_body_lsjt_component_sectors,two_body_lsjt_component_matrices
    );
  two_body_lsjt_timer.Stop();
  std::cout << "  Time: " << two_body_lsjt_timer.ElapsedTime() << std::endl;

  ////////////////////////////////////////////////////////////////
  // write as two-body LSJT
  ////////////////////////////////////////////////////////////////

  std::cout << "Write two-body NLSJT..." << std::endl;
  std::string two_body_lsjt_filename(parameters.two_body_filename);
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
  // recouple to two-body NJJJT
  ////////////////////////////////////////////////////////////////

  // std::cout << "Recouple to two-body NJJJT..." << std::endl;
  // 
  // // define space and operator containers
  // basis::TwoBodySpaceNJJJT two_body_njjjt_space(Nmax);
  // std::array<basis::TwoBodySectorsNJJJT,3> two_body_njjjt_component_sectors;
  // std::array<basis::MatrixVector,3> two_body_njjjt_component_matrices;
  // 
  // // do recoupling
  // Timer two_body_njjjt_timer;
  // two_body_njjjt_timer.Start();
  // moshinsky::TransformOperatorTwoBodyNLSJTToTwoBodyNJJJT(
  //     operator_labels,
  //     two_body_nlsjt_space,two_body_nlsjt_component_sectors,two_body_nlsjt_component_matrices,
  //     two_body_njjjt_space,two_body_njjjt_component_sectors,two_body_njjjt_component_matrices
  //   );
  // two_body_njjjt_timer.Stop();
  // std::cout << "  Time: " << two_body_njjjt_timer.ElapsedTime() << std::endl;




  // termination
  return 0;
}
