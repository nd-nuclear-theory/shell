/******************************************************************************/
/**
  @file natorb-gen.cpp

  generate natural orbitals

  Syntax:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obdme-compare_test tolerance input_orbitals robdme_info robdme robdme_info robdme
    obdme-compare_test tolerance input_orbitals robdme_info robdme robdme
    obdme-compare_test tolerance input_orbitals robdme robdme
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Patrick J. Fasano
  University of Notre Dame

  + 07/25/18 (pjf): Created, based on natorb-gen.
  + 04/07/19 (pjf): Modify to try to determine overall phase in addition to
      checking consistency (to provided tolerance) between matrix elements.
  + 05/09/19 (pjf): Use std::size_t for basis indices and sizes.

******************************************************************************/

#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>

#include "basis/nlj_orbital.h"
#include "density/obdme_io.h"
#include "basis/operator.h"
#include "mcutils/eigen.h"
#include "mcutils/parsing.h"

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

// Stores simple parameters for run
struct RunParameters {
  double tolerance;
  // filenames
  std::string input_orbital_file;
  std::string robdme_info1;
  std::string robdme1;
  std::string robdme_info2;
  std::string robdme2;
};

void PrintUsage(const char **argv) {
  std::cout << "Usage: " << argv[0]
                         << " tolerance input_orbitals"
                         << " [robdme_info] robdme [robdme_info] robdme"
                         << std::endl;
}

void ProcessArguments(int argc, const char *argv[], RunParameters& run_parameters) {
  // argument counter
  int arg = 0;

  // usage message
  if ((argc-1 < 4) || (argc-1 > 6)) {
    PrintUsage(argv);
    std::cerr << "ERROR: Incorrect number of arguments." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // tolerance
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.tolerance;
  }

  // input orbital file
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.input_orbital_file;
    mcutils::FileExistCheck(run_parameters.input_orbital_file, true, false);
  }

  // input ROBDME info file
  if (argc-1 >= 5)
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.robdme_info1;
    mcutils::FileExistCheck(run_parameters.robdme_info1, true, false);
  }
  else
  {
    run_parameters.robdme_info1 = "";
  }

  // input ROBDME file
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.robdme1;
    mcutils::FileExistCheck(run_parameters.robdme1, true, false);
  }

  // input ROBDME info file
  if (argc-1 == 6)
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.robdme_info2;
    mcutils::FileExistCheck(run_parameters.robdme_info2, true, false);
  }
  else
  {
    run_parameters.robdme_info2 = "";
  }

  // input static ROBDME file
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.robdme2;
    mcutils::FileExistCheck(run_parameters.robdme2, true, false);
  }
}


////////////////////////////////////////////////////////////////
// main program
////////////////////////////////////////////////////////////////
int main(int argc, const char *argv[]) {
  // // header
  // std::cout << std::endl;
  // std::cout << "obdme-compare_test -- natural orbital generation" << std::endl;
  // std::cout << "version: " VCS_REVISION << std::endl;
  // std::cout << std::endl;

  RunParameters run_parameters;
  ProcessArguments(argc, argv, run_parameters);

  // Read input orbitals
  std::ifstream input_orbital_s(run_parameters.input_orbital_file);
  basis::OrbitalSpaceLJPN input_space(basis::ParseOrbitalPNStream(input_orbital_s, true));
  input_orbital_s.close();

  // Get OBDMEs
  shell::InOBDMEStream in_obdme1, in_obdme2;
  if (run_parameters.robdme_info1 == "") {
    in_obdme1 = shell::InOBDMEStreamSingle(run_parameters.robdme1, input_space);
  } else {
    in_obdme1 = shell::InOBDMEStreamMulti(
        run_parameters.robdme_info1, run_parameters.robdme1, input_space
      );
  }
  if (run_parameters.robdme_info2 == "") {
    in_obdme2 = shell::InOBDMEStreamSingle(run_parameters.robdme2, input_space);
  } else {
    in_obdme2 = shell::InOBDMEStreamMulti(
        run_parameters.robdme_info2, run_parameters.robdme2, input_space
      );
  }

  assert(in_obdme1.g0() == in_obdme2.g0());
  assert(in_obdme1.Tz0() == in_obdme2.Tz0());

  // overall phase freedom on OBDMEs -- zero means undetermined
  int overall_phase=0;

  // convenience definitions
  const double& tolerance = run_parameters.tolerance;
  int min_K = std::max(in_obdme1.min_K(), in_obdme2.min_K());
  int max_K = std::min(in_obdme1.max_K(), in_obdme2.max_K());
  for (int K = min_K; K <= max_K; K++) {
    // std::cout << "Comparing K=" << K << std::endl;
    basis::OrbitalSectorsLJPN sectors;
    basis::OperatorBlocks<double> matrices1, matrices2;
    in_obdme1.GetMultipole(K, sectors, matrices1);
    in_obdme2.GetMultipole(K, sectors, matrices2);

    for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
    {
      // get next sector
      const basis::OrbitalSectorsLJPN::SectorType sector = sectors.GetSector(sector_index);
      // get sizes
      const std::size_t bra_subspace_size = sector.bra_subspace().size();
      const std::size_t ket_subspace_size = sector.ket_subspace().size();
      const basis::OperatorBlock<double>& matrix1 = matrices1[sector_index];
      const basis::OperatorBlock<double>& matrix2 = matrices2[sector_index];

      // main loop
      for (std::size_t j=0; j < bra_subspace_size; ++j)
      {
        for (std::size_t k=0; k < ket_subspace_size; ++k)
        {
          double difference;
          // numerical zeros -- skip
          if ((matrix1(j,k) < tolerance) && (matrix2(j,k) < tolerance))
            continue;

          if (overall_phase == 0)
          {
            overall_phase = +1;
            difference = matrix1(j,k) - overall_phase*matrix2(j,k);
            if (std::abs(difference) < tolerance) continue;  // matrix elements agree
            overall_phase = -1;
            difference = matrix1(j,k) - overall_phase*matrix2(j,k);
            if (std::abs(difference) < tolerance) continue;  // matrix elements agree
          }
          else
          {
            difference = matrix1(j,k) - overall_phase*matrix2(j,k);
            if (std::abs(difference) < tolerance) continue;  // matrix elements agree
          }

          basis::OrbitalStateLJPN bra_state(sector.bra_subspace(), j);
          basis::OrbitalStateLJPN ket_state(sector.ket_subspace(), k);
          std::cerr << "Matrix elements disagree!" << std::endl;
          std::cerr << " " << bra_state.LabelStr()
                    << " " << bra_state.LabelStr()
                    << " K=" << K
                    << "  (" << matrix1(j,k)
                    << ")-(" << overall_phase
                    << ")*(" << matrix2(j,k)
                    << ") = " << difference
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
      }
    }
  }

  /* return code */
  return EXIT_SUCCESS;
}
