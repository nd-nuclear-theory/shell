/******************************************************************************/
/**
  @file natorb-gen.cpp

  generate natural orbitals

  Syntax:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    natorb-gen
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Input format:

    set-indexing orbital_filename
    define-densities J g n robdme_filename [robdme_info_filename]
    set-output-xform output_orbital_filename output_xform_filename [output_info_filename]

  Patrick J. Fasano
  University of Notre Dame

  + 11/13/16 (pjf): Created, based on orbital-gen.
  + 10/12/17 (pjf): Update for changes to radial_io.
  + 11/28/17 (pjf): Print header with version.
  + 12/29/17 (pjf): Use input orbital indexing for output orbitals.
  + 12/30/17 (pjf): Ensure orbital file is properly sorted.
  + 07/27/18 (pjf): Update for new OBDME input routines.
  + 05/09/19 (pjf): Use std::size_t for basis indices and sizes.
  + 08/16/19 (pjf): Remove radial operator type and power from OutOBMEStream.
  + 04/09/20 (pjf):
    - Rewrite to take input on stdin.
    - Require quantum numbers of many-body state.
    - Use Eigen's reverse() rather than custom sorting function.

******************************************************************************/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include "basis/nlj_orbital.h"
#include "density/obdme_io.h"
#include "obme/obme_io.h"
#include "basis/operator.h"
#include "mcutils/eigen.h"
#include "mcutils/parsing.h"
#include "am/halfint_fmt.h"

////////////////////////////////////////////////////////////////
// define a function mapping eigenvalue to weight
/////////////////////////////////////////////////////////////////

// traditional oscillator weight
void ComputeOrbitalWeights(std::vector<basis::OrbitalPNInfo>& orbitals) {
  for (auto& orbital : orbitals) {
    orbital.weight = 2 * orbital.n + orbital.l;
  }
}

// // weight = 0.25 * |log(eigenvalue/max_eigenvalue)|
// //   -- max_eigenvalue specific to species
// void ComputeOrbitalWeights(std::vector<basis::OrbitalPNInfo>& orbitals) {
//   double max_eigenvalue[2];
//   for (auto& orbital : orbitals) {
//     if (orbital.weight > max_eigenvalue[static_cast<int>(orbital.orbital_species)]) {
//       max_eigenvalue[static_cast<int>(orbital.orbital_species)] = orbital.weight;
//     }
//   }
//
//   for (auto& orbital : orbitals) {
//     if (orbital.weight > 0 ) {
//       orbital.weight = 0.25 *
//         std::abs(std::log(orbital.weight / max_eigenvalue[static_cast<int>(orbital.orbital_species)]));
//     } else {
//       orbital.weight = 999;
//     }
//   }
// }

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

// Stores simple parameters for run
struct RunParameters {
  // filenames
  basis::OrbitalSpaceLJPN input_space;
  std::unique_ptr<shell::InOBDMEStream> density_stream;
  std::string output_xform_filename;
  std::string output_orbital_filename;
  std::string output_info_filename;
};

void PrintUsage(char **argv) { std::cout << "Usage: " << argv[0] << std::endl; }

void ReadParameters(RunParameters& run_parameters) {
  std::string line;
  int line_count = 0;

  while (mcutils::GetLine(std::cin, line, line_count)) {
    // set up for line parsing
    std::istringstream line_stream(line);
    std::string keyword;
    line_stream >> keyword;

    // select action based on keyword
    if (keyword == "set-output-xform") {
      line_stream >> run_parameters.output_orbital_filename
                  >> run_parameters.output_xform_filename;
      mcutils::ParsingCheck(line_stream, line_count, line);
      mcutils::FileExistCheck(run_parameters.output_orbital_filename, false, true);
      mcutils::FileExistCheck(run_parameters.output_xform_filename, false, true);
      if (!line_stream.eof()) {
        line_stream >> run_parameters.output_info_filename;
        mcutils::ParsingCheck(line_stream, line_count, line);
        mcutils::FileExistCheck(run_parameters.output_info_filename,false,true);
      } else {
        run_parameters.output_info_filename = "";
      }
    }
    else if (keyword == "set-indexing")
    {
      std::string orbital_filename;
      line_stream >> orbital_filename;
      mcutils::ParsingCheck(line_stream, line_count, line);
      mcutils::FileExistCheck(orbital_filename, true, false);

      std::ifstream orbital_stream(orbital_filename);
      std::vector<basis::OrbitalPNInfo> input_orbitals =
          basis::ParseOrbitalPNStream(orbital_stream, true);
      run_parameters.input_space = basis::OrbitalSpaceLJPN(input_orbitals);
    }
    else if (keyword == "define-densities")
    {
      std::unique_ptr<shell::InOBDMEStream> density_stream;
      float J;
      int g, n;
      std::string robdme_filename, robdme_info_filename="";
      line_stream >> J >> g >> n >> robdme_filename;
      mcutils::ParsingCheck(line_stream, line_count, line);
      mcutils::FileExistCheck(robdme_filename, true, false);
      if (!line_stream.eof()) {
        line_stream >> robdme_info_filename;
        mcutils::ParsingCheck(line_stream, line_count, line);
        mcutils::FileExistCheck(robdme_info_filename, true, false);

        // construct multi-file stream
        run_parameters.density_stream =
          std::unique_ptr<shell::InOBDMEStream>(new shell::InOBDMEStreamMulti(
            robdme_info_filename, robdme_filename, run_parameters.input_space,
            HalfInt(2*J,2), g, n, HalfInt(2*J, 2), g, n
          ));
      } else {
        // construct single-file stream
        run_parameters.density_stream =
          std::unique_ptr<shell::InOBDMEStream>(new shell::InOBDMEStreamSingle(
            robdme_filename, run_parameters.input_space,
            HalfInt(2*J,2), g, n, HalfInt(2*J, 2), g, n
          ));
      }
    }
  }
}


////////////////////////////////////////////////////////////////
// main program
////////////////////////////////////////////////////////////////
int main(int argc, const char *argv[]) {
  // header
  std::cout << std::endl;
  std::cout << "natorb-gen -- natural orbital generation" << std::endl;
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  RunParameters run_parameters;
  ReadParameters(run_parameters);

  // Get sectors and blocks
  basis::OrbitalSectorsLJPN sectors;
  basis::OperatorBlocks<double> density_matrices;
  run_parameters.density_stream->GetMultipole(0, sectors, density_matrices);

  // Here we loop through the density matrices and diagonalize each sector.
  basis::OperatorBlocks<double> xform_matrices;
  // TODO(pjf) generalize to use eigenvalues for orbital weights
  std::vector<basis::OrbitalPNInfo> orbital_occupation_info;
  for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index) {
    // get next sector
    const auto& sector = sectors.GetSector(sector_index);
    const auto& species = sector.bra_subspace().orbital_species();
    const auto& l = sector.bra_subspace().l();
    const auto& j = sector.bra_subspace().j();

    // the one-body static density matrix should be real and symmetric, so
    // we can use the SelfAdjointEigenSolver from Eigen
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(density_matrices[sector_index]);

    // SelfAdjointEigenSolver returns sorted by eigenvalue in ascending order.
    // We need to be in descending order, so we reverse
    Eigen::VectorXd eigenvalues = eigensolver.eigenvalues().reverse();
    Eigen::MatrixXd eigenvectors = eigensolver.eigenvectors().rowwise().reverse();

    // add to output
    xform_matrices.push_back(eigenvectors);
    // TODO(pjf) generalize to use eigenvalues for orbital weights
    for (int i=0; i < eigenvalues.size(); ++i) {
      orbital_occupation_info.emplace_back(
          species, i, l, j,
          eigenvalues(i)
        );
    }
  }

  // TODO(pjf) generalize to use eigenvalues for orbital weights
  // ComputeOrbitalWeights(output_orbitals);
  basis::OrbitalPNList output_orbitals = run_parameters.input_space.OrbitalInfo();

  // sort orbitals and write out
  std::sort(output_orbitals.begin(), output_orbitals.end(), basis::OrbitalSortCmpWeight);
  std::ofstream output_orbital_s(run_parameters.output_orbital_filename);
  output_orbital_s << basis::OrbitalDefinitionStr(output_orbitals, true);
  output_orbital_s.close();

  if (run_parameters.output_info_filename != "") {
    std::sort(
        orbital_occupation_info.begin(), orbital_occupation_info.end(),
        basis::OrbitalSortCmpWeight
      );
    std::reverse(orbital_occupation_info.begin(), orbital_occupation_info.end());
    std::ofstream output_info_s(run_parameters.output_info_filename);
    output_info_s << fmt::format("# {:>4s}  {:>4s}  {:>4s}  {:>4s}  {:>15s}\n",
        "n", "l", "j", "tz", "occupation"
      );
    for (const auto& orbital : orbital_occupation_info) {
      output_info_s << fmt::format(
          "  {:>4d}  {:>4d}  {:>4.1f}  {:>+4.1f}  {:>15.8e}\n",
          orbital.n, orbital.l, float(orbital.j),
          float(basis::kOrbitalSpeciesPNCodeTz[static_cast<int>(orbital.orbital_species)]),
          orbital.weight
        );
    }
    output_info_s.close();
  }

  // write xform out to file
  const basis::OrbitalSpaceLJPN output_space(output_orbitals);
  basis::OrbitalSectorsLJPN output_sectors(
      run_parameters.input_space, output_space, 0, 0, 0
      // note hard-coded J0=0, g0, Tz0=0
    );
  shell::OutOBMEStream xs(
      run_parameters.output_xform_filename,
      run_parameters.input_space, output_space, output_sectors,
      basis::OneBodyOperatorType::kRadial
    );
  xs.Write(xform_matrices);

  /* return code */
  return EXIT_SUCCESS;
}
