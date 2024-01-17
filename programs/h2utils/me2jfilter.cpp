#include <fstream>

#include "tbme/me2j_io.h"
#include "tbme/tbme_scheme_xform.h"

void TFilter(
    size_t mode, // 0 for isoscalar, 1 for non-isoscalar
    const std::array<basis::TwoBodySectorsJJJT,3>& two_body_jjjt_component_sectors,
    std::array<basis::OperatorBlocks<double>,3>& two_body_jjjt_component_matrices
  ) {
    if (mode == 0) {
      basis::SetOperatorToZero(two_body_jjjt_component_sectors[1],two_body_jjjt_component_matrices[1]);
      basis::SetOperatorToZero(two_body_jjjt_component_sectors[2],two_body_jjjt_component_matrices[2]);
    } else if (mode == 1) {
      basis::SetOperatorToZero(two_body_jjjt_component_sectors[0],two_body_jjjt_component_matrices[0]);
    }
  }

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
      std::cout << "Syntax: me2jfilter Nmax input_filename output_filename" << std::endl;
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
  std::cout << "me2jfilter  -- ME2J T Filtering" << std::endl;
  std::cout << std::endl;

  // read parameters
  RunParameters run_parameters;
  ProcessArguments(argc,argv,run_parameters);

  // set up parameters
  // std::string input_filename = "chi2bSMSI2C_srg0400_eMax16_EMax16_hwHO014.me2j.bin";
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
  TFilter(0,two_body_jjjt_component_sectors,two_body_jjjt_component_matrices);
  shell::TransformOperatorTwoBodyJJJTToTwoBodyJJJTTz(
    two_body_jjjt_space,two_body_jjjt_component_sectors,two_body_jjjt_component_matrices,
    two_body_jjjttz_space,two_body_jjjttz_sectors,two_body_jjjttz_matrices
  );

  // write me2j file
  // std::string output_filename = "chi2bSMSI2C_srg0400is_eMax16_EMax16_hwHO014.me2j.bin";
  shell::WriteMe2jFile(two_body_jjjttz_space,two_body_jjjttz_sectors,two_body_jjjttz_matrices,run_parameters.output_filename);

  // exit
  return EXIT_SUCCESS;
}
