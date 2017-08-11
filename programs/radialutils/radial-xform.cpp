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
  basis::OrbitalSpaceLJPN in_bra_space, in_ket_space;
  basis::OrbitalSectorsLJPN in_sectors;
  shell::RadialOperatorType in_operator_type = is.radial_operator_type();
  is.SetToIndexing(in_bra_space, in_ket_space, in_sectors);

  basis::OrbitalSpaceLJPN olap_bra_space, olap_ket_space;
  basis::OrbitalSectorsLJPN olap_sectors;
  shell::RadialOperatorType olap_type = olaps.radial_operator_type();
  olaps.SetToIndexing(olap_bra_space, olap_ket_space, olap_sectors);

  // check that operator is transformable (bra and ket spaces are the same)
  if (in_bra_space.OrbitalInfo() != in_ket_space.OrbitalInfo()) {
    std::cerr << "ERROR: Bra and ket spaces of this operator are not the same. "
              << "Cannot transform." << std::endl;
    ////////////// CHECKS TURNED OFF FOR NATURAL ORBITAL TESTING! //////////////
    // std::exit(EXIT_FAILURE);
  }

  // check that operator matches provided orbital file
  if (space.OrbitalInfo() != in_ket_space.OrbitalInfo()) {
    std::cerr << "ERROR: Operator space does not match orbitals provided." << std::endl;
    ////////////// CHECKS TURNED OFF FOR NATURAL ORBITAL TESTING! //////////////
    // std::exit(EXIT_FAILURE);
  }

  // check that overlaps are valid
  if ((olap_type != shell::RadialOperatorType::kO) ||
      (in_ket_space.OrbitalInfo() != olap_bra_space.OrbitalInfo()))
  {
    std::cerr << "ERROR: Invalid olap for this input." << std::endl;
    ////////////// CHECKS TURNED OFF FOR NATURAL ORBITAL TESTING! //////////////
    // std::exit(EXIT_FAILURE);
  }

  // check that we know how to handle this type of transformation
  /// @note Only similarity transforms are currently supported.
  if (olap_bra_space.OrbitalInfo() != olap_ket_space.OrbitalInfo()) {
    std::cerr << "ERROR: Overlap bra and ket spaces differ. Only similarity "
              << "transforms are currently supported." << std::endl;
    ////////////// CHECKS TURNED OFF FOR NATURAL ORBITAL TESTING! //////////////
    // std::exit(EXIT_FAILURE);
  }

  // check that overlaps match provided orbital file
  if (space.OrbitalInfo() != olap_ket_space.OrbitalInfo()) {
    std::cerr << "ERROR: Overlaps space does not match orbitals provided." << std::endl;
    ////////////// CHECKS TURNED OFF FOR NATURAL ORBITAL TESTING! //////////////
    // std::exit(EXIT_FAILURE);
  }

  // construct new indexing
  const shell::RadialOperatorType& out_operator_type = in_operator_type;
  const basis::OrbitalSpaceLJPN& out_bra_space = olap_ket_space;
  const basis::OrbitalSpaceLJPN& out_ket_space = olap_ket_space;
  basis::OrbitalSectorsLJPN out_sectors(out_ket_space, out_ket_space, in_sectors.l0max(), in_sectors.Tz0());

  // Eigen initialization
  basis::OperatorBlocks<double> olap_matrices, input_matrices, output_matrices;
  is.Read(input_matrices);
  olaps.Read(olap_matrices);

  // main loop
  for (int sector_index=0; sector_index < in_sectors.size(); ++sector_index) {
    const auto& operator_sector = in_sectors.GetSector(sector_index);
    int bra_subspace_index = operator_sector.bra_subspace_index();
    int ket_subspace_index = operator_sector.ket_subspace_index();

    // Sanity check on olap sector
    // This is only true in general for similarity transforms.
    ////////////// CHECKS TURNED OFF FOR NATURAL ORBITAL TESTING! //////////////
    // assert(olap_sectors.ContainsSector(bra_subspace_index, bra_subspace_index));
    // assert(olap_sectors.ContainsSector(ket_subspace_index, ket_subspace_index));

    // Get olap sector index
    int left_olap_sector_index  = olap_sectors.LookUpSectorIndex(bra_subspace_index, bra_subspace_index);
    int right_olap_sector_index = olap_sectors.LookUpSectorIndex(ket_subspace_index, ket_subspace_index);

    // do matrix multipication and add to the output sectors
    output_matrices.push_back(
      olap_matrices[left_olap_sector_index].transpose() * input_matrices[sector_index] * olap_matrices[right_olap_sector_index]
    );
  }

  // write out to file
  std::cout << "INFO: Writing to file " << run_parameters.output_filename << std::endl;
  shell::OutRadialStream os(run_parameters.output_filename,
                            out_ket_space, out_ket_space, out_sectors,
                            in_operator_type);
  os.Write(output_matrices);

  is.Close();
  os.Close();
  olaps.Close();

  return EXIT_SUCCESS;
}
