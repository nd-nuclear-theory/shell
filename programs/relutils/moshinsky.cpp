/****************************************************************
  moshinsky.cpp

  Perform Moshinsky transformation of general relative operator.

  See lsjt_operator.h for documentation of operator storage and the
  relative operator file format.

  Standard input:
    ob|tb truncation_cutoff
    source_filename coupling
    target_filename coupling [parameters...]

    coupling : source/target coupling scheme or h2 version
      "rel" -- relative matrix elements
      "relcm" -- relative-cm matrix elements
      "lsjt" -- two-body LSJT-coupled matrix elements
      "jjjt" -- two-body jjJT-coupled matrix elements
      "jjjpn" -- two-body jjJpn-coupled matrix elements
        parameters: h2_format Tz0
          h2_format -- h2 version number
          Tz0 -- operator branch from T->Tz

  Note: Currently supported source couplings are "rel" and "relcm".

  Language: C++11

  Mark A. Caprio
  University of Notre Dame

  + 07/16/16 (mac): Created.
  + 08/16/16 (mac): Add diagnostic output of relative-cm matrix elements.
  + 10/09/16 (pjf): Rename mcpp -> mcutils.
  + 10/25/16 (mac): Update use of mcutils::ParsingError.
  + 11/28/17 (pjf): Print header with version.
  + 07/31/18 (pjf): Call ChopMatrix before final write, for simple comparison.
  + 12/21/18 (pjf): Add support for reading relative-cm files.
  + 02/21/19 (pjf): Add H2 Version15200 support.
  + 03/28/19 (pjf): Allow input of relative-cm (rcmlsjt) matrix elements.
  + 04/07/19 (pjf): Change keyword for relative-cm from rcmlsjt to relcm.
  + 05/09/19 (pjf): Use std::size_t for basis indices and sizes.
  + 06/20/19 (pjf):
    - Allow final branching to Tz0 > 0.
    - Branch to jjJpn and write sector-by-sector.
  + 10/06/19 (pjf): Fix input file for target coupling other than jjjpn.

****************************************************************/

#include <fstream>

#include "am/wigner_gsl.h"
#include "basis/lsjt_operator.h"
#include "basis/jjjt_operator.h"
#include "basis/jjjpn_scheme.h"
#include "basis/jjjpn_operator.h"
#include "mcutils/eigen.h"
#include "mcutils/parsing.h"
#include "mcutils/profiling.h"
#include "tbme/h2_io.h"

#include "moshinsky/moshinsky_xform.h"


////////////////////////////////////////////////////////////////
// parameter input
////////////////////////////////////////////////////////////////

enum class Coupling : int {
    kRelative = 0, kRelativeCMLSJT = 1, kLSJT = 2, kJJJT = 3, kJJJPN = 4
  };

std::unordered_map<std::string,Coupling>
kCouplingDefinitions =
  {
    {"rel",      Coupling::kRelative},
    {"relcm",    Coupling::kRelativeCMLSJT},
    {"lsjt",     Coupling::kLSJT},
    {"jjjt",     Coupling::kJJJT},
    {"jjjpn",    Coupling::kJJJPN}
  };

struct Parameters
// Container for run input parameters.
{

  // input
  std::string input_filename;

  // transformation
  basis::Rank truncation_rank;
  int truncation_cutoff;
  Coupling source_coupling;
  Coupling target_coupling;

  // output
  shell::H2Format output_h2_format;
  int Tz0;
  std::string output_filename;

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
                >> parameters.truncation_cutoff;
    mcutils::ParsingCheck(line_stream,line_count,line);

    // process truncation rank code
    if (truncation_rank_code == "ob")
      parameters.truncation_rank = basis::Rank::kOneBody;
    else if (truncation_rank_code == "tb")
      parameters.truncation_rank = basis::Rank::kTwoBody;
    else
      mcutils::ParsingError(line_count,line,"unrecognized truncation rank code");


  }

  // line 2: relative (source) filename
  {
    ++line_count;
    std::getline(std::cin,line);
    std::istringstream line_stream(line);
    std::string source_coupling_code;
    line_stream >> parameters.input_filename >> source_coupling_code;
    mcutils::ParsingCheck(line_stream,line_count,line);
    if (kCouplingDefinitions.count(source_coupling_code))
      parameters.source_coupling = kCouplingDefinitions.at(source_coupling_code);
    else
      mcutils::ParsingError(line_count,line,"unrecognized coupling scheme code");
  }

  // line 3: two-body (target) filename
  {
    ++line_count;
    std::getline(std::cin,line);
    std::istringstream line_stream(line);
    std::string target_coupling_code;

    line_stream >> parameters.output_filename
                >> target_coupling_code;
    mcutils::ParsingCheck(line_stream,line_count,line);
    if (kCouplingDefinitions.count(target_coupling_code))
      parameters.target_coupling = kCouplingDefinitions.at(target_coupling_code);
    else
      mcutils::ParsingError(line_count,line,"unrecognized coupling scheme code");

    if (parameters.target_coupling == Coupling::kJJJPN)
    {
      line_stream >> parameters.output_h2_format
                  >> parameters.Tz0;
      mcutils::ParsingCheck(line_stream,line_count,line);
    }
    else
    {
      parameters.output_h2_format = basis::kNone;
      parameters.Tz0 = 0;
    }
  }
}

////////////////////////////////////////////////////////////////
// gathering and writing intermediate couplings
////////////////////////////////////////////////////////////////

void GatherAndWriteRelativeCM(
    int N2max,
    const basis::OperatorLabelsJT& operator_labels,
    const basis::RelativeCMSpaceLSJTN& relative_cm_lsjtn_space,
    const std::array<basis::RelativeCMSectorsLSJTN,3>& relative_cm_lsjtn_component_sectors,
    const std::array<basis::OperatorBlocks<double>,3>& relative_cm_lsjtn_component_blocks,
    const std::string& output_filename
  )
{
  std::cout << "Gather relative-cm LSJT..." << std::endl;

  // define space and operator containers
  basis::RelativeCMSpaceLSJT relative_cm_lsjt_space(N2max);
  std::array<basis::RelativeCMSectorsLSJT,3> relative_cm_lsjt_component_sectors;
  std::array<basis::OperatorBlocks<double>,3> relative_cm_lsjt_component_blocks;

  // construct gathered operator
  mcutils::SteadyTimer relative_cm_lsjt_timer;
  relative_cm_lsjt_timer.Start();
  basis::GatherOperatorRelativeCMLSJTNToRelativeCMLSJT(
      operator_labels,
      relative_cm_lsjtn_space,relative_cm_lsjtn_component_sectors,relative_cm_lsjtn_component_blocks,
      relative_cm_lsjt_space,relative_cm_lsjt_component_sectors,relative_cm_lsjt_component_blocks
    );
  relative_cm_lsjt_timer.Stop();
  std::cout << "  Time: " << relative_cm_lsjt_timer.ElapsedTime() << std::endl;

  // write

  std::cout << "Write relative-cm LSJT..." << std::endl;
  mcutils::SteadyTimer write_relative_cm_lsjt_timer;
  write_relative_cm_lsjt_timer.Start();
  basis::WriteRelativeCMOperatorLSJT(
      output_filename,
      relative_cm_lsjt_space,
      operator_labels,relative_cm_lsjt_component_sectors,relative_cm_lsjt_component_blocks,
      true
    );
  write_relative_cm_lsjt_timer.Stop();
  std::cout << "  Time: " << write_relative_cm_lsjt_timer.ElapsedTime() << std::endl;
}

void GatherAndWriteTwoBodyLSJT(
    const Parameters& parameters,
    const basis::OperatorLabelsJT& operator_labels,
    const basis::TwoBodySpaceLSJTN& two_body_lsjtn_space,
    const std::array<basis::TwoBodySectorsLSJTN,3>& two_body_lsjtn_component_sectors,
    const std::array<basis::OperatorBlocks<double>,3>& two_body_lsjtn_component_blocks,
    const std::string& output_filename
  )
{
  std::cout << "Gather two-body LSJT..." << std::endl;

  // define space and operator containers
  basis::TwoBodySpaceLSJT two_body_lsjt_space(parameters.truncation_rank,parameters.truncation_cutoff);
  std::array<basis::TwoBodySectorsLSJT,3> two_body_lsjt_component_sectors;
  std::array<basis::OperatorBlocks<double>,3> two_body_lsjt_component_blocks;

  // construct gathered operator
  mcutils::SteadyTimer two_body_lsjt_timer;
  two_body_lsjt_timer.Start();
  basis::GatherOperatorTwoBodyLSJTNToTwoBodyLSJT(
      operator_labels,
      two_body_lsjtn_space,two_body_lsjtn_component_sectors,two_body_lsjtn_component_blocks,
      two_body_lsjt_space,two_body_lsjt_component_sectors,two_body_lsjt_component_blocks
    );
  two_body_lsjt_timer.Stop();
  std::cout << "  Time: " << two_body_lsjt_timer.ElapsedTime() << std::endl;

  // write

  std::cout << "Write two-body LSJT..." << std::endl;
  std::ostringstream lsjt_sstream;
  for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
    {
      basis::WriteTwoBodyOperatorComponentLSJT(
          lsjt_sstream,
          T0,
          two_body_lsjt_component_sectors[T0],two_body_lsjt_component_blocks[T0],
          basis::NormalizationConversion::kASToNAS
        );
    }
  std::ofstream lsjt_stream(output_filename);
  lsjt_stream << lsjt_sstream.str();

}

void WriteTwoBodyJJJT(
    const basis::OperatorLabelsJT& operator_labels,
    const basis::TwoBodySpaceJJJT& two_body_jjjt_space,
    const std::array<basis::TwoBodySectorsJJJT,3>& two_body_jjjt_component_sectors,
    const std::array<basis::OperatorBlocks<double>,3>& two_body_jjjt_component_blocks,
    const std::string& output_filename
  )
{
  std::cout << "Write two-body JJJT..." << std::endl;
  std::ostringstream jjjt_sstream;
  for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
    {
      basis::WriteTwoBodyOperatorComponentJJJT(
          jjjt_sstream,
          T0,
          two_body_jjjt_component_sectors[T0],two_body_jjjt_component_blocks[T0],
          basis::NormalizationConversion::kASToNAS
        );
    }
  std::ofstream jjjt_stream(output_filename);
  jjjt_stream << jjjt_sstream.str();
}

void WriteTwoBodyJJJPN(
    const basis::OperatorLabelsJT& operator_labels,
    const basis::TwoBodySpaceJJJPN& two_body_jjjpn_space,
    const basis::TwoBodySectorsJJJPN& two_body_jjjpn_sectors,
    const basis::OperatorBlocks<double> two_body_jjjpn_blocks,
    const std::string& output_filename
  )
{
  // output as NAS and with 1-based indexing for comparison with MFDn
  // h2 format

  std::ostringstream jjjpn_sstream;
  basis::WriteTwoBodyOperatorJJJPN(
      jjjpn_sstream,
      two_body_jjjpn_sectors,two_body_jjjpn_blocks,
      basis::NormalizationConversion::kASToNAS,
      1
    );
  std::ofstream jjjpn_stream(output_filename);
  jjjpn_stream << jjjpn_sstream.str();
}

void WriteTwoBodyH2(
    const Parameters& parameters,
    const basis::OperatorLabelsJT& operator_labels,
    const basis::OrbitalSpacePN& orbital_space,
    const basis::TwoBodySpaceJJJPN& two_body_jjjpn_space,
    const basis::TwoBodySectorsJJJPN& two_body_jjjpn_sectors,
    const basis::OperatorBlocks<double> two_body_jjjpn_blocks
  )
{
  // stream initialization
  std::cout << "Output stream" << std::endl;
  shell::OutH2Stream output_stream(
      parameters.output_filename,
      orbital_space,two_body_jjjpn_space,two_body_jjjpn_sectors,
      parameters.output_h2_format
    );
  std::cout << output_stream.DiagnosticStr();
  std::cout << std::endl;

  // iterate over sectors and write out
  for (std::size_t sector_index = 0; sector_index < output_stream.num_sectors(); ++sector_index)
    {
      // "neaten" output by eliminating near-zero values
      basis::OperatorBlock<double> matrix = two_body_jjjpn_blocks[sector_index];
      mcutils::ChopMatrix(matrix);
      output_stream.WriteSector(
          sector_index,
          matrix,
          basis::NormalizationConversion::kASToNAS
        );
      std::cout << "." << std::flush;
    }
  std::cout << std::endl;
}

////////////////////////////////////////////////////////////////
// reading and scattering intermediate couplings
////////////////////////////////////////////////////////////////

void ReadAndScatterRelativeCM(
    const std::string& input_filename,
    basis::OperatorLabelsJT& operator_labels,
    basis::RelativeCMSpaceLSJTN& relative_cm_lsjtn_space,
    std::array<basis::RelativeCMSectorsLSJTN,3>& relative_cm_lsjtn_component_sectors,
    std::array<basis::OperatorBlocks<double>,3>& relative_cm_lsjtn_component_blocks
  )
{
  // define space and operator containers
  basis::RelativeCMOperatorParametersLSJT operator_parameters;
  basis::RelativeCMSpaceLSJT relative_cm_lsjt_space;
  std::array<basis::RelativeCMSectorsLSJT,3> relative_cm_lsjt_component_sectors;
  std::array<basis::OperatorBlocks<double>,3> relative_cm_lsjt_component_blocks;

  // read
  mcutils::SteadyTimer read_relative_cm_lsjt_timer;
  read_relative_cm_lsjt_timer.Start();
  basis::ReadRelativeCMOperatorLSJT(
      input_filename,
      relative_cm_lsjt_space,
      operator_parameters,relative_cm_lsjt_component_sectors,relative_cm_lsjt_component_blocks,
      true
    );
  read_relative_cm_lsjt_timer.Stop();
  std::cout << "  Time: " << read_relative_cm_lsjt_timer.ElapsedTime() << std::endl;

  // construct scattered operator
  std::cout << "Scatter relative-cm LSJTN..." << std::endl;
  mcutils::SteadyTimer relative_cm_lsjtn_timer;
  relative_cm_lsjtn_timer.Start();
  relative_cm_lsjtn_space = basis::RelativeCMSpaceLSJTN(operator_parameters.Nmax);
  operator_labels = operator_parameters;
  basis::ScatterOperatorRelativeCMLSJTToRelativeCMLSJTN(
      operator_labels,
      relative_cm_lsjt_space,relative_cm_lsjt_component_sectors,relative_cm_lsjt_component_blocks,
      relative_cm_lsjtn_space,relative_cm_lsjtn_component_sectors,relative_cm_lsjtn_component_blocks
    );
  relative_cm_lsjtn_timer.Stop();
  std::cout << "  Time: " << relative_cm_lsjtn_timer.ElapsedTime() << std::endl;

}

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  // header
  std::cout << std::endl;
  std::cout << "moshinsky -- Moshinsky transformation of general relative operator" << std::endl;
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  // parameter input
  Parameters parameters;
  ReadParameters(parameters);
  std::cout << "Parameters"
            << std::endl
            << " "
            << " rank " << int(parameters.truncation_rank)
            << " cutoff " << parameters.truncation_cutoff
            << std::endl;


  // check sanity of requested recoupling
  if (parameters.source_coupling == parameters.target_coupling)
    {
      std::cout << "Source coupling and target coupling are the same. "
                << "Nothing to do, exiting."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  if (parameters.source_coupling > parameters.target_coupling)
    {
      std::cout << "Inverse Moshinsky transformation not supported, exiting."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }

  // process truncation cutoff
  int N1max, N2max;
  std::tie(N1max,N2max) = basis::TwoBodyCutoffs(parameters.truncation_rank,parameters.truncation_cutoff);

  // operator tensorial character
  basis::OperatorLabelsJT operator_labels;


  ////////////////////////////////////////////////////////////////
  // total run timer
  ////////////////////////////////////////////////////////////////
  mcutils::SteadyTimer total_run_timer;
  total_run_timer.Start();

  ////////////////////////////////////////////////////////////////
  // relative LSJT
  ////////////////////////////////////////////////////////////////

  // declare relative space and operator containers
  basis::RelativeSpaceLSJT relative_space;
  std::array<basis::RelativeSectorsLSJT,3> relative_component_sectors;
  std::array<basis::OperatorBlocks<double>,3> relative_component_blocks;

  // read space and operator
  if (parameters.source_coupling == Coupling::kRelative)
    {
      mcutils::SteadyTimer read_relative_lsjt_timer;
      read_relative_lsjt_timer.Start();
      basis::RelativeOperatorParametersLSJT operator_parameters;
      basis::ReadRelativeOperatorLSJT(
          parameters.input_filename,
          relative_space,
          operator_parameters,relative_component_sectors,relative_component_blocks,
          true
        );
      read_relative_lsjt_timer.Stop();
      std::cout << "  Time: " << read_relative_lsjt_timer.ElapsedTime() << std::endl;
      operator_labels = operator_parameters;
    }


  ////////////////////////////////////////////////////////////////
  // relative-cm LSJTN
  ////////////////////////////////////////////////////////////////

  // declare relative-cm space and operator containers
  basis::RelativeCMSpaceLSJTN relative_cm_lsjtn_space;
  std::array<basis::RelativeCMSectorsLSJTN,3> relative_cm_lsjtn_component_sectors;
  std::array<basis::OperatorBlocks<double>,3> relative_cm_lsjtn_component_blocks;

  if (parameters.source_coupling < Coupling::kRelativeCMLSJT)
    {
      std::cout << "Augment to relative-cm LSJTN..." << std::endl;

      // do transformation
      mcutils::SteadyTimer relative_cm_lsjtn_timer;
      relative_cm_lsjtn_timer.Start();
      relative_cm_lsjtn_space = basis::RelativeCMSpaceLSJTN(N2max);
      moshinsky::TransformOperatorRelativeLSJTToRelativeCMLSJTN(
          operator_labels,
          relative_space,relative_component_sectors,relative_component_blocks,
          relative_cm_lsjtn_space,relative_cm_lsjtn_component_sectors,relative_cm_lsjtn_component_blocks,
          true
        );
      relative_cm_lsjtn_timer.Stop();
      std::cout << "  Time: " << relative_cm_lsjtn_timer.ElapsedTime() << std::endl;
    }
  else if (parameters.source_coupling == Coupling::kRelativeCMLSJT)
    {
      ReadAndScatterRelativeCM(
          parameters.input_filename,
          operator_labels,
          relative_cm_lsjtn_space, relative_cm_lsjtn_component_sectors, relative_cm_lsjtn_component_blocks
        );
    }

  if (parameters.source_coupling <= Coupling::kRelativeCMLSJT)
    {
      // write diagnostics
      std::cout << "  Allocated matrix elements:";
      for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
        std::cout << " " << basis::AllocatedEntries(relative_cm_lsjtn_component_blocks[T0]);
      std::cout << std::endl;
    }

  // gather and write
  if (parameters.target_coupling == Coupling::kRelativeCMLSJT)
    {
      GatherAndWriteRelativeCM(
          N2max, operator_labels,
          relative_cm_lsjtn_space,relative_cm_lsjtn_component_sectors,relative_cm_lsjtn_component_blocks,
          parameters.output_filename
        );
      std::exit(EXIT_SUCCESS);
    }

  ////////////////////////////////////////////////////////////////
  // Moshinsky transform to two-body LSJTN
  ////////////////////////////////////////////////////////////////

  // declare space and operator containers
  basis::TwoBodySpaceLSJTN two_body_lsjtn_space;
  std::array<basis::TwoBodySectorsLSJTN,3> two_body_lsjtn_component_sectors;
  std::array<basis::OperatorBlocks<double>,3> two_body_lsjtn_component_blocks;

  if (parameters.source_coupling < Coupling::kLSJT)
    {
      std::cout << "Transform to two-body LSJTN..." << std::endl;

      // do transformation
      mcutils::SteadyTimer two_body_lsjtn_timer;
      two_body_lsjtn_timer.Start();
      two_body_lsjtn_space = basis::TwoBodySpaceLSJTN(parameters.truncation_rank,parameters.truncation_cutoff);
      moshinsky::TransformOperatorRelativeCMLSJTNToTwoBodyLSJTN(
          operator_labels,
          relative_cm_lsjtn_space,relative_cm_lsjtn_component_sectors,relative_cm_lsjtn_component_blocks,
          two_body_lsjtn_space,two_body_lsjtn_component_sectors,two_body_lsjtn_component_blocks
        );
      two_body_lsjtn_timer.Stop();
      std::cout << "  Time: " << two_body_lsjtn_timer.ElapsedTime() << std::endl;
    }
  else if (parameters.source_coupling == Coupling::kLSJT)
    {
      std::cerr << "ERROR: lsjt input not yet implemented" << std::endl;
      std::exit(EXIT_FAILURE);
    }

  // write diagnostics
  if (parameters.source_coupling <= Coupling::kLSJT)
    {
      std::cout << "  Allocated matrix elements:";
      for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
        std::cout << " " << basis::AllocatedEntries(two_body_lsjtn_component_blocks[T0]);
      std::cout << std::endl;
    }

  if (parameters.target_coupling == Coupling::kLSJT)
    {
      GatherAndWriteTwoBodyLSJT(
          parameters,operator_labels,
          two_body_lsjtn_space,
          two_body_lsjtn_component_sectors,
          two_body_lsjtn_component_blocks,
          parameters.output_filename
        );
      std::exit(EXIT_SUCCESS);
    }

  ////////////////////////////////////////////////////////////////
  // recouple to two-body JJJTN
  ////////////////////////////////////////////////////////////////

  // declare JJJTN space and operator containers
  basis::TwoBodySpaceJJJTN two_body_jjjtn_space;
  std::array<basis::TwoBodySectorsJJJTN,3> two_body_jjjtn_component_sectors;
  std::array<basis::OperatorBlocks<double>,3> two_body_jjjtn_component_blocks;

  if (parameters.source_coupling < Coupling::kJJJT)
    {
      std::cout << "Recouple to two-body JJJTN..." << std::endl;

      // do recoupling
      mcutils::SteadyTimer two_body_jjjtn_timer;
      two_body_jjjtn_timer.Start();
      two_body_jjjtn_space = basis::TwoBodySpaceJJJTN(parameters.truncation_rank,parameters.truncation_cutoff);
      moshinsky::TransformOperatorTwoBodyLSJTNToTwoBodyJJJTN(
          operator_labels,
          two_body_lsjtn_space,two_body_lsjtn_component_sectors,two_body_lsjtn_component_blocks,
          two_body_jjjtn_space,two_body_jjjtn_component_sectors,two_body_jjjtn_component_blocks
        );
      two_body_jjjtn_timer.Stop();
      std::cout << "  Time: " << two_body_jjjtn_timer.ElapsedTime() << std::endl;
    }

  ////////////////////////////////////////////////////////////////
  // gather to two-body JJJT
  ////////////////////////////////////////////////////////////////

  // declare JJJT space and operator containers
  basis::TwoBodySpaceJJJT two_body_jjjt_space;
  std::array<basis::TwoBodySectorsJJJT,3> two_body_jjjt_component_sectors;
  std::array<basis::OperatorBlocks<double>,3> two_body_jjjt_component_blocks;

  if (parameters.source_coupling < Coupling::kJJJT)
    {
      std::cout << "Gather to two-body JJJT..." << std::endl;

      // construct gathered operator
      mcutils::SteadyTimer two_body_jjjt_timer;
      two_body_jjjt_timer.Start();
      two_body_jjjt_space = basis::TwoBodySpaceJJJT(parameters.truncation_rank,parameters.truncation_cutoff);
      basis::GatherOperatorTwoBodyJJJTNToTwoBodyJJJT(
          operator_labels,
          two_body_jjjtn_space,two_body_jjjtn_component_sectors,two_body_jjjtn_component_blocks,
          two_body_jjjt_space,two_body_jjjt_component_sectors,two_body_jjjt_component_blocks
        );
      two_body_jjjt_timer.Stop();
      std::cout << "  Time: " << two_body_jjjt_timer.ElapsedTime() << std::endl;
    }
  else if (parameters.source_coupling == Coupling::kJJJT)
    {
      std::cerr << "ERROR: jjjt input not yet implemented" << std::endl;
      std::exit(EXIT_FAILURE);
    }

  if (parameters.target_coupling == Coupling::kJJJT)
    {
      WriteTwoBodyJJJT(
          operator_labels,
          two_body_jjjt_space,
          two_body_jjjt_component_sectors,
          two_body_jjjt_component_blocks,
          parameters.output_filename
        );
      std::exit(EXIT_SUCCESS);
    }

  ////////////////////////////////////////////////////////////////
  // branch to two-body JJJPN
  ////////////////////////////////////////////////////////////////

  std::cout << "Branch to two-body JJJPN and write to file..." << std::endl;

  // define space and operator containers
  basis::OrbitalSpacePN orbital_space(N1max);
  basis::TwoBodySpaceJJJPNOrdering two_body_jjjpn_space_ordering =
    shell::kH2SpaceOrdering.at(parameters.output_h2_format);
  basis::TwoBodySpaceJJJPN two_body_jjjpn_space(
      orbital_space,
      basis::WeightMax(parameters.truncation_rank,parameters.truncation_cutoff),
      two_body_jjjpn_space_ordering
    );
  basis::TwoBodySectorsJJJPN two_body_jjjpn_sectors(
      two_body_jjjpn_space,
      operator_labels.J0, operator_labels.g0, parameters.Tz0
    );

  // stream initialization
  mcutils::SteadyTimer two_body_jjjpn_timer;
  two_body_jjjpn_timer.Start();
  shell::OutH2Stream output_stream(
      parameters.output_filename,
      orbital_space, two_body_jjjpn_space, two_body_jjjpn_sectors,
      parameters.output_h2_format
    );
  std::cout << output_stream.DiagnosticStr();

  // do branching
  for (std::size_t sector_index=0; sector_index<two_body_jjjpn_sectors.size(); ++sector_index)
    {
      // make reference to target sector
      const basis::TwoBodySectorsJJJPN::SectorType& two_body_jjjpn_sector
        = two_body_jjjpn_sectors.GetSector(sector_index);

      // transform
      auto matrix = moshinsky::TwoBodyMatrixJJJPN(
          operator_labels,
          two_body_jjjt_space,
          two_body_jjjt_component_sectors,
          two_body_jjjt_component_blocks,
          two_body_jjjpn_sector
        );
      mcutils::ChopMatrix(matrix);
      output_stream.WriteSector(
          sector_index,
          matrix,
          basis::NormalizationConversion::kASToNAS
        );
      std::cout << "." << std::flush;
    }
  output_stream.Close();
  two_body_jjjpn_timer.Stop();
  std::cout << std::endl;
  std::cout << "  Time: " << two_body_jjjpn_timer.ElapsedTime() << std::endl;

  ////////////////////////////////////////////////////////////////
  // termination
  ////////////////////////////////////////////////////////////////

  total_run_timer.Stop();
  std::cout << "Total run time: " << total_run_timer.ElapsedTime() << std::endl;

  std::exit(EXIT_SUCCESS);
}
