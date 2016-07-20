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
#include "basis/jjjpnorb_operator.h"
#include "mcpp/parsing.h"
#include "mcpp/profiling.h"

#include "moshinsky/moshinsky_xform.h"


////////////////////////////////////////////////////////////////
// parameter input
////////////////////////////////////////////////////////////////

enum class Coupling {kLSJT, kJJJT, kJJJPN};

struct Parameters
// Container for run input parameters.
{
  basis::Rank truncation_rank;
  int truncation_cutoff;
  Coupling coupling;
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
    std::string truncation_rank_code, coupling_code;
    line_stream >> truncation_rank_code
                >> parameters.truncation_cutoff
                >> coupling_code;
    ParsingCheck(line_stream,line_count,line);

    // process truncation rank code
    if (truncation_rank_code == "ob")
      parameters.truncation_rank = basis::Rank::kOneBody;
    else if (truncation_rank_code == "tb")
      parameters.truncation_rank = basis::Rank::kTwoBody;
    else
    {
      std::cerr << "ERROR: unrecognized truncation rank code" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    // process coupling
    if (coupling_code == "lsjt")
      parameters.coupling = Coupling::kLSJT;
    else if (coupling_code == "jjjt")
      parameters.coupling = Coupling::kJJJT;
    else if (coupling_code == "jjjpn")
      parameters.coupling = Coupling::kJJJPN;
    else
    {
      std::cerr << "ERROR: unrecognized coupling scheme code" << std::endl;
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
    << " "
    << " J0 " << operator_parameters.J0
    << " g0 " << operator_parameters.g0
    << " T0_min " << operator_parameters.T0_min
    << " T0_max " << operator_parameters.T0_max
    << " symmetry " << int(operator_parameters.symmetry_phase_mode)
    << std::endl
    << " "
    << " Nmax " << operator_parameters.Nmax
    << " Jmax " << operator_parameters.Jmax
    << std::endl;

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

  // write diagnostics
  std::cout << "  Allocated matrix elements:";
  for (int T0=operator_parameters.T0_min; T0<=operator_parameters.T0_max; ++T0)
    std::cout << " " << basis::AllocatedEntries(relative_component_matrices[T0]);
  std::cout << std::endl;

}



////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  // parameter input
  Parameters parameters;
  ReadParameters(parameters);
  std::cout << "Parameters"
            << std::endl
            << " "
            << " rank " << int(parameters.truncation_rank)
            << " cutoff " << parameters.truncation_cutoff
            << std::endl;


  // process truncation cutoff
  int N1max, N2max;
  std::tie(N1max,N2max) = basis::TwoBodyCutoffs(parameters.truncation_rank,parameters.truncation_cutoff);

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

  
  // AD HOC: following is initial working code, to be refactored and
  // extended to jjjt and jjjpn

  ////////////////////////////////////////////////////////////////
  // augment to relative-cm LSJTN
  ////////////////////////////////////////////////////////////////

  std::cout << "Augment to relative-cm LSJTN..." << std::endl;

  // define space and operator containers
  basis::RelativeCMSpaceLSJTN relative_cm_lsjtn_space(N2max);
  std::array<basis::RelativeCMSectorsLSJTN,3> relative_cm_lsjtn_component_sectors;
  std::array<basis::MatrixVector,3> relative_cm_lsjtn_component_matrices;

  // do transformation
  Timer relative_cm_lsjtn_timer;
  relative_cm_lsjtn_timer.Start();
  moshinsky::TransformOperatorRelativeLSJTToRelativeCMLSJTN(
      operator_labels,
      relative_space,relative_component_sectors,relative_component_matrices,
      relative_cm_lsjtn_space,relative_cm_lsjtn_component_sectors,relative_cm_lsjtn_component_matrices
    );
  relative_cm_lsjtn_timer.Stop();
  std::cout << "  Time: " << relative_cm_lsjtn_timer.ElapsedTime() << std::endl;

  // write diagnostics
  std::cout << "  Allocated matrix elements:";
  for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
    std::cout << " " << basis::AllocatedEntries(relative_cm_lsjtn_component_matrices[T0]);
  std::cout << std::endl;

  ////////////////////////////////////////////////////////////////
  // transform to two-body LSJTN
  ////////////////////////////////////////////////////////////////

  std::cout << "Transform to two-body LSJTN..." << std::endl;

  // define space and operator containers
  basis::TwoBodySpaceLSJTN two_body_lsjtn_space(parameters.truncation_rank,parameters.truncation_cutoff);
  std::array<basis::TwoBodySectorsLSJTN,3> two_body_lsjtn_component_sectors;
  std::array<basis::MatrixVector,3> two_body_lsjtn_component_matrices;

  // do transformation
  Timer two_body_lsjtn_timer;
  two_body_lsjtn_timer.Start();
  moshinsky::TransformOperatorRelativeCMLSJTNToTwoBodyLSJTN(
      operator_labels,
      relative_cm_lsjtn_space,relative_cm_lsjtn_component_sectors,relative_cm_lsjtn_component_matrices,
      two_body_lsjtn_space,two_body_lsjtn_component_sectors,two_body_lsjtn_component_matrices
  );
  two_body_lsjtn_timer.Stop();
  std::cout << "  Time: " << two_body_lsjtn_timer.ElapsedTime() << std::endl;

  // write diagnostics
  std::cout << "  Allocated matrix elements:";
  for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
    std::cout << " " << basis::AllocatedEntries(two_body_lsjtn_component_matrices[T0]);
  std::cout << std::endl;

  ////////////////////////////////////////////////////////////////
  // BEGIN: LSJT coupling only
  ////////////////////////////////////////////////////////////////

  if (parameters.coupling == Coupling::kLSJT)
    {

      ////////////////////////////////////////////////////////////////
      // gather to two-body LSJT
      ////////////////////////////////////////////////////////////////

      std::cout << "Gather two-body LSJT..." << std::endl;

      // define space and operator containers
      basis::TwoBodySpaceLSJT two_body_lsjt_space(parameters.truncation_rank,parameters.truncation_cutoff);
      std::array<basis::TwoBodySectorsLSJT,3> two_body_lsjt_component_sectors;
      std::array<basis::MatrixVector,3> two_body_lsjt_component_matrices;

      // construct gathered operator
      Timer two_body_lsjt_timer;
      two_body_lsjt_timer.Start();
      basis::GatherOperatorTwoBodyLSJTNToTwoBodyLSJT(
          operator_labels,
          two_body_lsjtn_space,two_body_lsjtn_component_sectors,two_body_lsjtn_component_matrices,
          two_body_lsjt_space,two_body_lsjt_component_sectors,two_body_lsjt_component_matrices
        );
      two_body_lsjt_timer.Stop();
      std::cout << "  Time: " << two_body_lsjt_timer.ElapsedTime() << std::endl;

      ////////////////////////////////////////////////////////////////
      // write as two-body LSJT
      ////////////////////////////////////////////////////////////////

      std::cout << "Write two-body LSJT..." << std::endl;
      std::ostringstream lsjt_sstream;
      for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
        {
          basis::WriteTwoBodyOperatorComponentLSJT(
              lsjt_sstream,
              T0,
              two_body_lsjt_component_sectors[T0],two_body_lsjt_component_matrices[T0],
              basis::NormalizationConversion::kASToNAS
            );
        }
      std::ofstream lsjt_stream(parameters.two_body_filename.c_str());
      lsjt_stream << lsjt_sstream.str();

      ////////////////////////////////////////////////////////////////
      // termination
      ////////////////////////////////////////////////////////////////

      std::exit(0);
    }

  ////////////////////////////////////////////////////////////////
  // END: LSJT coupling only
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  // recouple to two-body JJJTN
  ////////////////////////////////////////////////////////////////

  std::cout << "Recouple to two-body JJJTN..." << std::endl;
  
  // define space and operator containers
  basis::TwoBodySpaceJJJTN two_body_jjjtn_space(parameters.truncation_rank,parameters.truncation_cutoff);
  std::array<basis::TwoBodySectorsJJJTN,3> two_body_jjjtn_component_sectors;
  std::array<basis::MatrixVector,3> two_body_jjjtn_component_matrices;
  
  // do recoupling
  Timer two_body_jjjtn_timer;
  two_body_jjjtn_timer.Start();
  moshinsky::TransformOperatorTwoBodyLSJTNToTwoBodyJJJTN(
      operator_labels,
      two_body_lsjtn_space,two_body_lsjtn_component_sectors,two_body_lsjtn_component_matrices,
      two_body_jjjtn_space,two_body_jjjtn_component_sectors,two_body_jjjtn_component_matrices
    );
  two_body_jjjtn_timer.Stop();
  std::cout << "  Time: " << two_body_jjjtn_timer.ElapsedTime() << std::endl;

  ////////////////////////////////////////////////////////////////
  // gather to two-body JJJT
  ////////////////////////////////////////////////////////////////

  std::cout << "Gather to two-body JJJT..." << std::endl;

  // define space and operator containers
  basis::TwoBodySpaceJJJT two_body_jjjt_space(parameters.truncation_rank,parameters.truncation_cutoff);
  std::array<basis::TwoBodySectorsJJJT,3> two_body_jjjt_component_sectors;
  std::array<basis::MatrixVector,3> two_body_jjjt_component_matrices;

  // construct gathered operator
  Timer two_body_jjjt_timer;
  two_body_jjjt_timer.Start();
  basis::GatherOperatorTwoBodyJJJTNToTwoBodyJJJT(
      operator_labels,
      two_body_jjjtn_space,two_body_jjjtn_component_sectors,two_body_jjjtn_component_matrices,
      two_body_jjjt_space,two_body_jjjt_component_sectors,two_body_jjjt_component_matrices
    );
  two_body_jjjt_timer.Stop();
  std::cout << "Time: " << two_body_jjjt_timer.ElapsedTime() << std::endl;

  ////////////////////////////////////////////////////////////////
  // BEGIN: JJJT coupling only
  ////////////////////////////////////////////////////////////////

  if (parameters.coupling == Coupling::kJJJT)
    {

      ////////////////////////////////////////////////////////////////
      // write as two-body JJJT
      ////////////////////////////////////////////////////////////////

      std::cout << "Write two-body JJJT..." << std::endl;
      std::string two_body_jjjt_filename(parameters.two_body_filename);
      std::ostringstream jjjt_sstream;
      for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
        {
          basis::WriteTwoBodyOperatorComponentJJJT(
              jjjt_sstream,
              T0,
              two_body_jjjt_component_sectors[T0],two_body_jjjt_component_matrices[T0],
              basis::NormalizationConversion::kASToNAS
            );
        }
      std::ofstream jjjt_stream(parameters.two_body_filename.c_str());
      jjjt_stream << jjjt_sstream.str();

      ////////////////////////////////////////////////////////////////
      // termination
      ////////////////////////////////////////////////////////////////

      std::exit(0);
    }

  ////////////////////////////////////////////////////////////////
  // END: JJJT coupling only
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  // BEGIN: JJJPN coupling only
  ////////////////////////////////////////////////////////////////

  if (parameters.coupling == Coupling::kJJJPN)
    {

      ////////////////////////////////////////////////////////////////
      // branch to two-body JJJPN
      ////////////////////////////////////////////////////////////////

      std::cout << "Branch to two-body JJJPN..." << std::endl;

      // define space and operator containers
      basis::OrbitalSpacePN orbital_space(N1max);
      basis::TwoBodySpaceJJJPN two_body_jjjpn_space(
          orbital_space,
          basis::WeightMax(parameters.truncation_rank,parameters.truncation_cutoff)
        );
      basis::TwoBodySectorsJJJPN two_body_jjjpn_sectors;
      basis::MatrixVector two_body_jjjpn_matrices;

      // do branching
      Timer two_body_jjjpn_timer;
      two_body_jjjpn_timer.Start();
      moshinsky::TransformOperatorTwoBodyJJJTToTwoBodyJJJPN(
          operator_labels,
          two_body_jjjt_space,two_body_jjjt_component_sectors,two_body_jjjt_component_matrices,
          two_body_jjjpn_space,two_body_jjjpn_sectors,two_body_jjjpn_matrices
        );
      two_body_jjjpn_timer.Stop();
      std::cout << "Time: " << two_body_jjjpn_timer.ElapsedTime() << std::endl;

      ////////////////////////////////////////////////////////////////
      // write as two-body JJJPN
      ////////////////////////////////////////////////////////////////

      // output as NAS and with 1-based indexing for comparison with MFDn
      // h2 format

      std::ostringstream jjjpn_sstream;
      basis::WriteTwoBodyOperatorJJJPN(
          jjjpn_sstream,
          two_body_jjjpn_sectors,two_body_jjjpn_matrices,
          basis::NormalizationConversion::kASToNAS,
          1
        );
      std::ofstream jjjpn_stream(parameters.two_body_filename.c_str());
      jjjpn_stream << jjjpn_sstream.str();

      ////////////////////////////////////////////////////////////////
      // termination
      ////////////////////////////////////////////////////////////////

      std::exit(0);
    }

  // We should never get here...

  // termination
  return 0;
}
