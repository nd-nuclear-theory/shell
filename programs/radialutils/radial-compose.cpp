/**************************************************************************/
/**
  @file radial-compose.cpp

  compose radial transformations

  Syntax:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    radial-compose in_olap_file1 in_olap_file2 out_olap_file
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  @author Patrick J. Fasano
  University of Notre Dame

  + 1/16/17 (pjf): Created, based on radial-xform.cpp.
  + 10/12/17 (pjf): Update for changes to radial_io.
  + 11/28/17 (pjf): Print header with version.
  + 10/18/18 (pjf): Update to use new obme library.
  + 08/16/19 (pjf): Remove radial operator type and power from OutOBMEStream.

******************************************************************************/

#include <sys/stat.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "basis/nlj_operator.h"
#include "basis/nlj_orbital.h"
#include "mcutils/profiling.h"
#include "obme/obme_io.h"
#include "obme/radial.h"

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

// Stores simple parameters for run
struct RunParameters
{
  // filenames
  std::string orbital_filename;
  std::string olap1_filename;
  std::string olap2_filename;
  std::string output_filename;
};

void PrintUsage(const char* argv[])
{
  std::cout << "Usage: " << argv[0]
            << " in_olap_file1 in_olap_file2 out_olap_file" << std::endl;
}

void ProcessArguments(const int argc, const char* argv[], RunParameters& run_parameters)
{
  // header
  std::cout << std::endl;
  std::cout << "radial-compose -- radial matrix element composition" << std::endl;
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  // argument counter
  int arg = 0;

  // usage message
  if (argc - 1 != 3)
  {
    PrintUsage(argv);
    std::cerr << "ERROR: Wrong number of arguments." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // olap file 1
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.olap1_filename;
    if (!parameter_stream)
    {
      PrintUsage(argv);
      std::cerr << "ERROR: Invalid olap1 filename." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    mcutils::FileExistCheck(run_parameters.olap1_filename, true, false);
  }

  // olap file 2
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.olap2_filename;
    if (!parameter_stream)
    {
      PrintUsage(argv);
      std::cerr << "ERROR: Invalid olap2 filename." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    mcutils::FileExistCheck(run_parameters.olap2_filename, true, false);
  }

  // output olap file
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.output_filename;
    if (!parameter_stream)
    {
      PrintUsage(argv);
      std::cerr << "ERROR: Invalid output filename." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    mcutils::FileExistCheck(run_parameters.output_filename, false, true);
  }
}

int main(int argc, const char* argv[])
{
  RunParameters run_parameters;
  ProcessArguments(argc, argv, run_parameters);

  // Read input
  shell::InOBMEStream olaps1(run_parameters.olap1_filename);
  shell::InOBMEStream olaps2(run_parameters.olap2_filename);

  // get indexing
  basis::OrbitalSpaceLJPN olap1_bra_space, olap1_ket_space;
  basis::OrbitalSectorsLJPN olap1_sectors;
  olaps1.SetToIndexing(olap1_bra_space, olap1_ket_space, olap1_sectors);

  basis::OrbitalSpaceLJPN olap2_bra_space, olap2_ket_space;
  basis::OrbitalSectorsLJPN olap2_sectors;
  olaps2.SetToIndexing(olap2_bra_space, olap2_ket_space, olap2_sectors);

  // check that overlaps are valid
  if (olaps1.operator_type() != basis::OneBodyOperatorType::kRadial)
  {
    std::cerr << "ERROR: Invalid olap1." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (olaps1.operator_type() != basis::OneBodyOperatorType::kRadial)
  {
    std::cerr << "ERROR: Invalid olap2." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if (olap1_ket_space.OrbitalInfo() != olap2_bra_space.OrbitalInfo())
  {
    std::cerr << "ERROR: Incompatible overlaps." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // construct new indexing
  const basis::OneBodyOperatorType& out_operator_type =
      basis::OneBodyOperatorType::kRadial;
  const basis::OrbitalSpaceLJPN& out_bra_space = olap1_bra_space;
  const basis::OrbitalSpaceLJPN& out_ket_space = olap2_ket_space;
  basis::OrbitalSectorsLJPN out_sectors(out_bra_space, out_ket_space, 0, 0, 0);

  // Eigen initialization
  basis::OperatorBlocks<double> olap1_matrices, olap2_matrices, output_matrices;
  olaps1.Read(olap1_matrices);
  olaps1.Close();
  olaps2.Read(olap2_matrices);
  olaps2.Close();

  shell::ComposeRadialOperators(
      olap1_bra_space, olap1_ket_space, olap1_sectors, olap1_matrices,
      olap2_bra_space, olap2_ket_space, olap2_sectors, olap2_matrices,
      out_sectors, output_matrices);

  // write out to file
  std::cout << "INFO: Writing to file " << run_parameters.output_filename
            << std::endl;
  shell::OutOBMEStream os(
      run_parameters.output_filename, out_bra_space, out_ket_space, out_sectors,
      out_operator_type);
  os.Write(output_matrices);
  os.Close();

  return EXIT_SUCCESS;
}
