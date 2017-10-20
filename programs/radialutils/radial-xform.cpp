/**************************************************************************/
/**
  @file radial-xform.cpp

  transform single-particle matrix elements

  Syntax:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    radial-xform orbital_file input_olap_file input_me_file output_me_file
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  @author Patrick J. Fasano
  University of Notre Dame

  + 11/7/16 (pjf): Created, based on radial-scale.cpp.
  + 11/22/16 (pjf): Implemented similarity transforms.
  + 09/19/17 (pjf): Improved transformation to accept more types of transformations.
  + 10/12/17 (pjf): Update for changes to radial_io.
  + 10/18/17 (pjf): Correctly set target_operator_order.

******************************************************************************/

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
  std::string olap_filename;
  std::string input_filename;
  std::string output_filename;
};

void PrintUsage(const char *argv[]) {
  std::cout << "Usage: "
            << argv[0] << " orbital_file input_olap_file input_me_file output_me_file"
            << std::endl;
}

void ProcessArguments(const int argc, const char *argv[], RunParameters& run_parameters) {
  // usage message
  if (argc-1 != 4) {
    PrintUsage(argv);
    std::cerr << "ERROR: Wrong number of arguments." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // orbital file
  {
    std::istringstream parameter_stream(argv[1]);
    parameter_stream >> run_parameters.orbital_filename;
    if (!parameter_stream) {
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
    if (!parameter_stream) {
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
    if (!parameter_stream) {
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
    if (!parameter_stream) {
      PrintUsage(argv);
      std::cerr << "ERROR: Invalid output filename." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    /** @todo check and warn if overwriting existing file */
  }
}

int main(int argc, const char *argv[]) {
  RunParameters run_parameters;
  ProcessArguments(argc, argv, run_parameters);

  // Read input
  shell::InRadialStream is(run_parameters.input_filename);
  shell::InRadialStream olaps(run_parameters.olap_filename);

  // Read orbitals
  std::ifstream orbitals(run_parameters.orbital_filename);
  std::vector<basis::OrbitalPNInfo> input_orbitals = basis::ParseOrbitalPNStream(orbitals, true);
  basis::OrbitalSpaceLJPN space(input_orbitals);

  // get indexing
  basis::OrbitalSpaceLJPN source_bra_space, source_ket_space;
  basis::OrbitalSectorsLJPN source_sectors;
  shell::RadialOperatorType source_operator_type = is.radial_operator_type();
  int source_operator_order = is.radial_operator_power();
  is.SetToIndexing(source_bra_space, source_ket_space, source_sectors);

  basis::OrbitalSpaceLJPN olap_bra_space, olap_ket_space;
  basis::OrbitalSectorsLJPN olap_sectors;
  shell::RadialOperatorType olap_type = olaps.radial_operator_type();
  olaps.SetToIndexing(olap_bra_space, olap_ket_space, olap_sectors);

  // check that operator is transformable (bra and ket spaces are the same)
  if (source_bra_space.OrbitalInfo() != source_ket_space.OrbitalInfo()) {
    std::cerr << "ERROR: Bra and ket spaces of this operator are not the same. "
              << "Cannot transform." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // check that operator matches provided orbital file
  if (space.OrbitalInfo() != source_ket_space.OrbitalInfo()) {
    std::cerr << "ERROR: Operator space does not match orbitals provided." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // check that overlaps are valid
  if ((olap_type != shell::RadialOperatorType::kO)
      || (source_ket_space.OrbitalInfo() != olap_bra_space.OrbitalInfo())
      || (olap_sectors.Tz0() != 0)
    )
  {
    std::cerr << "ERROR: Invalid olap for this input." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // construct new indexing
  const shell::RadialOperatorType& target_operator_type = source_operator_type;
  const int target_operator_order = source_operator_order;
  const basis::OrbitalSpaceLJPN& target_space = olap_ket_space;
  basis::OrbitalSectorsLJPN target_sectors;
  if (source_sectors.mode() == basis::SectorsConstraintMode::kAll) {
    target_sectors = basis::OrbitalSectorsLJPN(target_space, target_space);
  } else if (source_sectors.mode() == basis::SectorsConstraintMode::kRadial) {
    target_sectors = basis::OrbitalSectorsLJPN(target_space, target_space,
      source_sectors.l0max(), source_sectors.Tz0());
  } else if (source_sectors.mode() == basis::SectorsConstraintMode::kSpherical) {
    target_sectors = basis::OrbitalSectorsLJPN(target_space, target_space,
      source_sectors.j0(), source_sectors.g0(), source_sectors.Tz0());
  }

  // Eigen initialization
  basis::OperatorBlocks<double> olap_matrices, input_matrices, output_matrices;
  is.Read(input_matrices);
  olaps.Read(olap_matrices);

  // main loop
  for (int target_sector_index=0; target_sector_index < target_sectors.size(); ++target_sector_index) {
    const auto& target_sector = target_sectors.GetSector(target_sector_index);

    // get subspace indices and labels for target sector
    int target_bra_subspace_index = target_sector.bra_subspace_index();
    const auto& target_bra_subspace_labels = target_sector.bra_subspace().labels();
    int target_ket_subspace_index = target_sector.ket_subspace_index();
    const auto& target_ket_subspace_labels = target_sector.ket_subspace().labels();

    // Get subspace indices in source space
    int source_bra_subspace_index
      = source_bra_space.LookUpSubspaceIndex(target_bra_subspace_labels);
    if (source_bra_subspace_index == basis::kNone) {
      std::cerr << "ERROR: required subspace missing from input space:" << std::endl;
      std::cerr << target_sector.bra_subspace().LabelStr() << std::endl;
      std::exit(EXIT_FAILURE);
    }
    int source_ket_subspace_index
      = source_ket_space.LookUpSubspaceIndex(target_ket_subspace_labels);
    if (source_ket_subspace_index == basis::kNone) {
      std::cerr << "ERROR: required subspace missing from input space:" << std::endl;
      std::cerr << target_sector.ket_subspace().LabelStr() << std::endl;
      std::exit(EXIT_FAILURE);
    }

    // Get source operator sector index
    int source_sector_index =
      source_sectors.LookUpSectorIndex(source_bra_subspace_index, source_ket_subspace_index);
    if (source_sector_index == basis::kNone) {
      std::cerr << "ERROR: required input sector missing:" << std::endl;
      std::cerr << "  bra subspace: " << target_sector.bra_subspace().LabelStr() << std::endl;
      std::cerr << "  ket subspace: " << target_sector.ket_subspace().LabelStr() << std::endl;
      std::exit(EXIT_FAILURE);
    }

    // Get olap sector index
    int left_olap_sector_index
      = olap_sectors.LookUpSectorIndex(source_bra_subspace_index, target_bra_subspace_index);
    if (left_olap_sector_index == basis::kNone) {
      std::cerr << "ERROR: required overlap sector missing:" << std::endl;
      std::cerr << "  subspace: " << target_sector.bra_subspace().LabelStr() << std::endl;
      std::exit(EXIT_FAILURE);
    }
    int right_olap_sector_index
      = olap_sectors.LookUpSectorIndex(source_ket_subspace_index, target_ket_subspace_index);
    if (right_olap_sector_index == basis::kNone) {
      std::cerr << "ERROR: required overlap sector missing:" << std::endl;
      std::cerr << "  subspace: " << target_sector.ket_subspace().LabelStr() << std::endl;
      std::exit(EXIT_FAILURE);
    }

    // do matrix multipication and add to the output sectors
    output_matrices.push_back(
      olap_matrices[left_olap_sector_index].transpose()
      * input_matrices[source_sector_index]
      * olap_matrices[right_olap_sector_index]
    );
  }

  // write out to file
  std::cout << "INFO: Writing to file " << run_parameters.output_filename << std::endl;
  shell::OutRadialStream os(run_parameters.output_filename,
                            target_space, target_space, target_sectors,
                            target_operator_type, target_operator_order);
  os.Write(output_matrices);

  is.Close();
  os.Close();
  olaps.Close();

  return EXIT_SUCCESS;
}
