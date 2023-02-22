/**************************************************************************/
/**
  @file radial-xform.cpp

  transform single-particle matrix elements

  Syntax:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    radial-xform orbital_file input_xform_file input_me_file output_me_file
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  @author Patrick J. Fasano
  University of Notre Dame

  + 11/7/16 (pjf): Created, based on radial-scale.cpp.
  + 11/22/16 (pjf): Implemented similarity transforms.
  + 09/19/17 (pjf): Improved transformation to accept more types of transformations.
  + 10/12/17 (pjf): Update for changes to radial_io.
  + 10/18/17 (pjf): Correctly set target_operator_order.
  + 11/28/17 (pjf): Print header with version.
  + 10/18/18 (pjf): Update to use new obme library.
  + 08/16/19 (pjf): Remove radial operator type and power from OutOBMEStream.

******************************************************************************/

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
  std::string olap_filename;
  std::string input_filename;
  std::string output_filename;
};

void PrintUsage(const char* argv[])
{
  std::cout << "Usage: " << argv[0]
            << " orbital_file input_olap_file input_me_file output_me_file"
            << std::endl;
}

void ProcessArguments(const int argc, const char* argv[], RunParameters& run_parameters)
{
  // usage message
  if (argc - 1 != 4)
  {
    PrintUsage(argv);
    std::cerr << "ERROR: Wrong number of arguments." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // orbital file
  {
    std::istringstream parameter_stream(argv[1]);
    parameter_stream >> run_parameters.orbital_filename;
    if (!parameter_stream)
    {
      PrintUsage(argv);
      std::cerr << "ERROR: Invalid input filename." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    /** @todo check for file existence */
  }

  // olap file
  {
    std::istringstream parameter_stream(argv[2]);
    parameter_stream >> run_parameters.olap_filename;
    if (!parameter_stream)
    {
      PrintUsage(argv);
      std::cerr << "ERROR: Invalid input filename." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    /** @todo check for file existence */
  }

  // input file
  {
    std::istringstream parameter_stream(argv[3]);
    parameter_stream >> run_parameters.input_filename;
    if (!parameter_stream)
    {
      PrintUsage(argv);
      std::cerr << "ERROR: Invalid input filename." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    /** @todo check for file existence */
  }

  // output file
  {
    std::istringstream parameter_stream(argv[4]);
    parameter_stream >> run_parameters.output_filename;
    if (!parameter_stream)
    {
      PrintUsage(argv);
      std::cerr << "ERROR: Invalid output filename." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    /** @todo check and warn if overwriting existing file */
  }
}

int main(int argc, const char* argv[])
{
  // header
  std::cout << std::endl;
  std::cout << "radial-xform -- radial matrix element transformation" << std::endl;
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  RunParameters run_parameters;
  ProcessArguments(argc, argv, run_parameters);

  // Read input
  shell::InOBMEStream is(run_parameters.input_filename);
  shell::InOBMEStream xforms(run_parameters.olap_filename);

  // Read orbitals
  std::ifstream orbitals(run_parameters.orbital_filename);
  std::vector<basis::OrbitalPNInfo> input_orbitals =
      basis::ParseOrbitalPNStream(orbitals, true);
  basis::OrbitalSpaceLJPN space(input_orbitals);

  // get indexing
  basis::OrbitalSpaceLJPN source_bra_space, source_ket_space;
  basis::OrbitalSectorsLJPN source_sectors;
  basis::OneBodyOperatorType source_operator_type = is.operator_type();
  is.SetToIndexing(source_bra_space, source_ket_space, source_sectors);

  basis::OrbitalSpaceLJPN xform_bra_space, xform_ket_space;
  basis::OrbitalSectorsLJPN xform_sectors;
  basis::OneBodyOperatorType xform_operator_type = xforms.operator_type();
  xforms.SetToIndexing(xform_bra_space, xform_ket_space, xform_sectors);

  // check that operator is transformable (bra and ket spaces are the same)
  if (source_bra_space.OrbitalInfo() != source_ket_space.OrbitalInfo())
  {
    std::cerr << "ERROR: Bra and ket spaces of this operator are not the same. "
              << "Cannot transform." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // check that operator matches provided orbital file
  if (space.OrbitalInfo() != source_ket_space.OrbitalInfo())
  {
    std::cerr << "ERROR: Operator space does not match orbitals provided."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // check that overlaps are valid
  if ((xform_operator_type != basis::OneBodyOperatorType::kRadial)
      || (source_ket_space.OrbitalInfo() != xform_bra_space.OrbitalInfo())
      || (xform_sectors.Tz0() != 0))
  {
    std::cerr << "ERROR: Invalid olap for this input." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // construct new indexing
  const basis::OneBodyOperatorType& target_operator_type = source_operator_type;
  const basis::OrbitalSpaceLJPN& target_space = xform_ket_space;
  basis::OrbitalSectorsLJPN target_sectors;
  target_sectors = basis::OrbitalSectorsLJPN(
      target_space, target_space,
      source_sectors.j0(), source_sectors.g0(), source_sectors.Tz0());

  // Eigen initialization
  basis::OperatorBlocks<double> xform_matrices, source_matrices, target_matrices;
  is.Read(source_matrices);
  xforms.Read(xform_matrices);
  is.Close();
  xforms.Close();

  shell::SimilarityTransformOperator(
      space, target_space,
      xform_sectors, xform_matrices,
      source_sectors, source_matrices,
      target_sectors, target_matrices);

  // write out to file
  std::cout << "INFO: Writing to file " << run_parameters.output_filename
            << std::endl;
  shell::OutOBMEStream os(
      run_parameters.output_filename, target_space, target_space, target_sectors,
      target_operator_type);
  os.Write(target_matrices);

  os.Close();

  return EXIT_SUCCESS;
}
