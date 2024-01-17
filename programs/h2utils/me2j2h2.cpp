#include <fstream>

#include "basis/jjjt_operator.h"
#include "basis/jjjpn_scheme.h"
#include "basis/jjjpn_operator.h"
#include "mcutils/eigen.h"
#include "mcutils/parsing.h"
#include "mcutils/profiling.h"
#include "tbme/h2_io.h"
#include "tbme/me2j_io.h"
#include "tbme/tbme_scheme_xform.h"

#include "moshinsky/moshinsky_xform.h"

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

struct RunParameters
// Stores simple parameters for run
{
  // filenames
  std::string input_filename;
  std::string output_filename;
  // mode
  int Nmax;
};

void ProcessArguments(int argc, char **argv, RunParameters& run_parameters)
{
  // usage message
  if (argc-1 < 3)
    {
      std::cout << "Syntax: me2j2h2 Nmax input_filename output_filename" << std::endl;
      std::exit(EXIT_SUCCESS);
    }

  // format
  std::istringstream parameter_stream(argv[1]);
  parameter_stream >> run_parameters.Nmax;
  if (!parameter_stream)
    {
      std::cerr << "Expecting numeric value for Nmax argument" << std::endl;
      std::exit(EXIT_FAILURE);
    }

  // input filename
  run_parameters.input_filename = argv[2];

  // output filename
  run_parameters.output_filename = argv[3];
}

int main(int argc, char **argv)
{
  // header
  std::cout << std::endl;
  std::cout << "me2j2h2  -- ME2J TO H2" << std::endl;
  std::cout << std::endl;

  // read parameters
  RunParameters run_parameters;
  ProcessArguments(argc,argv,run_parameters);

  // set up parameters
  // std::string input_filename = "vrelccm_me2j.txt"; // for testing
  // std::string input_filename = "chi2bSMSI2B_srg0800_eMax16_EMax16_hwHO020.me2j.bin";
  int Nmax = run_parameters.Nmax;
  int J0 = 0;
  int g0 = 0;
  int Tz0 = 0;

  // initialize
  basis::TwoBodySpaceJJJTTz two_body_jjjttz_space(basis::Rank::kTwoBody,Nmax);
  basis::TwoBodySectorsJJJTTz two_body_jjjttz_sectors(two_body_jjjttz_space,J0,g0,Tz0);
  basis::OperatorBlocks<double> two_body_jjjttz_matrices;
  SetOperatorToZero(two_body_jjjttz_sectors,two_body_jjjttz_matrices);

  // read me2j file
  shell::ReadMe2jFile(two_body_jjjttz_space,two_body_jjjttz_sectors,two_body_jjjttz_matrices,run_parameters.input_filename);

  // transform
  basis::TwoBodySpaceJJJT two_body_jjjt_space(basis::Rank::kTwoBody,Nmax);
  std::array<basis::TwoBodySectorsJJJT,3> two_body_jjjt_component_sectors;
  std::array<basis::OperatorBlocks<double>,3> two_body_jjjt_component_matrices;
  shell::TransformOperatorTwoBodyJJJTTzToTwoBodyJJJT(two_body_jjjttz_space,two_body_jjjttz_sectors,two_body_jjjttz_matrices,
    two_body_jjjt_space,two_body_jjjt_component_sectors,two_body_jjjt_component_matrices
    );
  // TFilter(0,two_body_jjjt_component_sectors,two_body_jjjt_component_matrices);
  ////////////////////////////////////////////////////////////////
  // branch to two-body JJJPN
  ////////////////////////////////////////////////////////////////

  std::cout << "Branch to two-body JJJPN and write to file..." << std::endl;

  // define space and operator containers
  basis::OrbitalSpacePN orbital_space(Nmax);
  basis::TwoBodySpaceJJJPNOrdering two_body_jjjpn_space_ordering =
    shell::kH2SpaceOrdering.at(0);
  basis::TwoBodySpaceJJJPN two_body_jjjpn_space(
      orbital_space,
      basis::WeightMax(basis::Rank::kTwoBody,Nmax),
      two_body_jjjpn_space_ordering
    );
  basis::TwoBodySectorsJJJPN two_body_jjjpn_sectors(
      two_body_jjjpn_space,
      J0, g0, Tz0
    );

  // stream initialization
  mcutils::SteadyTimer two_body_jjjpn_timer;
  two_body_jjjpn_timer.Start();
  shell::OutH2Stream output_stream(
      // "vrel_ccm_h2fromme2j_Nmax4.dat", // for testing
      run_parameters.output_filename,
      orbital_space, two_body_jjjpn_space, two_body_jjjpn_sectors,
      shell::kVersion0
    );
  std::cout << output_stream.DiagnosticStr();

  basis::OperatorLabelsJT operator_labels(J0,g0,0,2,basis::SymmetryPhaseMode::kHermitian);

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
          two_body_jjjt_component_matrices,
          two_body_jjjpn_sector
        );
      // mcutils::ChopMatrix(matrix);
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

  std::exit(EXIT_SUCCESS);
}
