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

******************************************************************************/

#include <sys/stat.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "mcutils/profiling.h"
#include "basis/nlj_orbital.h"
#include "basis/nlj_operator.h"
#include "radial/radial_io.h"

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

// Stores simple parameters for run
struct RunParameters {
  // filenames
  std::string orbital_filename;
  std::string olap1_filename;
  std::string olap2_filename;
  std::string output_filename;
};

void PrintUsage(const char *argv[]) {
  std::cout << "Usage: "
            << argv[0] << " in_olap_file1 in_olap_file2 out_olap_file"
            << std::endl;
}

void ProcessArguments(const int argc, const char *argv[], RunParameters& run_parameters) {
  // argument counter
  int arg = 0;

  // usage message
  if (argc-1 != 3) {
    PrintUsage(argv);
    std::cerr << "ERROR: Wrong number of arguments." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // olap file 1
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.olap1_filename;
    if (!parameter_stream) {
      PrintUsage(argv);
      std::cerr << "ERROR: Invalid olap1 filename." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    struct stat st;
    if (stat(run_parameters.olap1_filename.c_str(), &st) != 0) {
      std::cerr << "ERROR: file " << run_parameters.olap1_filename << " does not exist!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // olap file 2
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.olap2_filename;
    if (!parameter_stream) {
      PrintUsage(argv);
      std::cerr << "ERROR: Invalid olap2 filename." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    struct stat st;
    if (stat(run_parameters.olap2_filename.c_str(), &st) != 0) {
      std::cerr << "ERROR: file " << run_parameters.olap2_filename << " does not exist!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // output olap file
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.output_filename;
    if (!parameter_stream) {
      PrintUsage(argv);
      std::cerr << "ERROR: Invalid output filename." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    struct stat st;
    if (stat(run_parameters.output_filename.c_str(), &st) == 0) {
      std::cerr << "WARN: overwriting file " << run_parameters.output_filename << std::endl;
    }
  }
}

int main(int argc, const char *argv[]) {
  RunParameters run_parameters;
  ProcessArguments(argc, argv, run_parameters);

  // Read input
  shell::InRadialStream olaps1(run_parameters.olap1_filename);
  shell::InRadialStream olaps2(run_parameters.olap2_filename);

  // get indexing
  basis::OrbitalSpaceLJPN olap1_bra_space, olap1_ket_space;
  basis::OrbitalSectorsLJPN olap1_sectors;
  shell::RadialOperatorType olap1_operator_type = olaps1.radial_operator_type();
  olaps1.SetToIndexing(olap1_bra_space, olap1_ket_space, olap1_sectors);

  basis::OrbitalSpaceLJPN olap2_bra_space, olap2_ket_space;
  basis::OrbitalSectorsLJPN olap2_sectors;
  shell::RadialOperatorType olap2_operator_type = olaps2.radial_operator_type();
  olaps2.SetToIndexing(olap2_bra_space, olap2_ket_space, olap2_sectors);

  // check that overlaps are valid
  if (olap1_operator_type != shell::RadialOperatorType::kO)
  {
    std::cerr << "ERROR: Invalid olap1." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (olap2_operator_type != shell::RadialOperatorType::kO)
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
  const shell::RadialOperatorType& out_operator_type = shell::RadialOperatorType::kO;
  const basis::OrbitalSpaceLJPN& out_bra_space = olap1_bra_space;
  const basis::OrbitalSpaceLJPN& out_ket_space = olap2_ket_space;
  basis::OrbitalSectorsLJPN out_sectors(out_bra_space, out_ket_space, 0, 0);

  // Eigen initialization
  basis::OperatorBlocks<double> olap1_matrices, olap2_matrices, output_matrices;
  olaps1.Read(olap1_matrices);
  olaps1.Close();
  olaps2.Read(olap2_matrices);
  olaps2.Close();

  // main loop
  for (int sector_index=0; sector_index < out_sectors.size(); ++sector_index) {
    // get indexing
    const auto& out_sector = out_sectors.GetSector(sector_index);
    const auto& out_bra_subspace = out_sector.bra_subspace();
    const int out_bra_subspace_index = out_sector.bra_subspace_index();
    const auto& out_ket_subspace = out_sector.ket_subspace();
    const int out_ket_subspace_index = out_sector.ket_subspace_index();

    // make sure we can look up subspace index for shared subspace by labels
    assert(out_bra_subspace.labels() == out_ket_subspace.labels());
    assert(olap1_ket_space.ContainsSubspace(out_bra_subspace.labels()));

    // get index of shared subspace (olap1_ket_space, olap2_bra_space)
    const int shared_subspace_index = olap1_ket_space.LookUpSubspaceIndex(out_bra_subspace.labels());

    // Sanity check on olap sector
    // This is only true in general for similarity transforms.
    assert(olap1_sectors.ContainsSector(out_bra_subspace_index, shared_subspace_index));
    assert(olap2_sectors.ContainsSector(shared_subspace_index, out_ket_subspace_index));

    // Get olap sector indices
    int left_olap_sector_index  = olap1_sectors.LookUpSectorIndex(out_bra_subspace_index, shared_subspace_index);
    int right_olap_sector_index = olap2_sectors.LookUpSectorIndex(shared_subspace_index, out_ket_subspace_index);

    // do matrix multipication and add to the output sectors
    output_matrices.push_back(
      olap1_matrices[left_olap_sector_index] *  olap2_matrices[right_olap_sector_index]
    );
  }

  // write out to file
  std::cout << "INFO: Writing to file " << run_parameters.output_filename << std::endl;
  shell::OutRadialStream os(run_parameters.output_filename,
                            out_bra_space, out_ket_space, out_sectors,
                            out_operator_type);
  os.Write(output_matrices);
  os.Close();


  return EXIT_SUCCESS;
}
