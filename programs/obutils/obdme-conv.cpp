/******************************************************************************/
/*
  @file obdme-conv.cpp

  convert OBDME to OBME format

  Syntax:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obdme-conv
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Input format:

    set-indexing <orbital_filename>
    define-densities <J_bra> <g_bra> <n_bra> <J_ket> <g_ket> <n_ket> <robdme_filename> [robdme_info_filename]
    set-output-filename <output_filename> <J0>

  Patrick J. Fasano
  Argonne National Laboratory

  + 01/08/24 (pjf): Created, based on espi-gen.
*/
/******************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <optional>
#include <vector>
#include <cmath>

#include <Eigen/Core>
#include <fmt/format.h>
#include <fmt/os.h>

#include "am/halfint.h"
#include "basis/nlj_orbital.h"
#include "basis/operator.h"
#include "mcutils/eigen.h"
#include "mcutils/io.h"
#include "mcutils/parsing.h"
#include "obme/obme_io.h"
#include "density/obdme_io.h"

// Stores simple parameters for run
struct RunParameters {
  // filenames
  std::string orbital_filename;
  std::string input_filename;
  std::optional<std::string> input_info_filename;
  HalfInt J_bra, J_ket;
  int g_bra, g_ket;
  int n_bra, n_ket;
  std::vector<std::tuple<std::string,int>> multipole_filenames;
};

RunParameters ReadParameters() {
  RunParameters run_parameters;
  std::string line;
  int line_count = 0;

  while (mcutils::GetLine(std::cin, line, line_count)) {
    // set up for line parsing
    std::istringstream line_stream(line);
    std::string keyword;
    line_stream >> keyword;

    // select action based on keyword
    if (keyword == "set-output-filename") {
      std::string filename;
      int J0;
      line_stream >> filename >> J0;
      mcutils::ParsingCheck(line_stream, line_count, line);
      mcutils::FileExistCheck(filename, false, true);
      run_parameters.multipole_filenames.emplace_back(filename, J0);
    }
    else if (keyword == "set-indexing")
    {
      line_stream >> run_parameters.orbital_filename;
      mcutils::ParsingCheck(line_stream, line_count, line);
      mcutils::FileExistCheck(run_parameters.orbital_filename, true, false);
    }
    else if (keyword == "define-densities")
    {
      float J_bra, J_ket;
      int g_bra, n_bra, g_ket, n_ket;
      std::string robdme_filename, robdme_info_filename="";
      line_stream >> J_bra >> g_bra >> n_bra >> J_ket >> g_ket >> n_ket >> robdme_filename;
      mcutils::ParsingCheck(line_stream, line_count, line);
      mcutils::FileExistCheck(robdme_filename, true, false);

      run_parameters.J_bra = HalfInt(2*J_bra,2);
      run_parameters.g_bra = g_bra;
      run_parameters.n_bra = n_bra;
      run_parameters.J_ket = HalfInt(2*J_ket,2);
      run_parameters.g_ket = g_ket;
      run_parameters.n_ket = n_ket;
      run_parameters.input_filename = std::move(robdme_filename);

      if (!line_stream.eof()) {
        line_stream >> robdme_info_filename;
        mcutils::ParsingCheck(line_stream, line_count, line);
        mcutils::FileExistCheck(robdme_info_filename, true, false);
        run_parameters.input_info_filename = std::move(robdme_info_filename);
      }
    }
  }

  return run_parameters;
}

int main(int argc, const char* argv[])
{
  // header
  fmt::print("\n");
  fmt::print("obdme-conv -- one-body density matrix conversion\n");
  fmt::print("version: {}\n", VCS_REVISION);
  fmt::print("\n");

  auto run_parameters = ReadParameters();

  // initialize orbitals
  std::ifstream orbital_stream(run_parameters.orbital_filename);
  auto orbital_space = basis::OrbitalSpaceLJPN(
      basis::ParseOrbitalPNStream(orbital_stream, true)
    );
  orbital_stream.close();

  // initialize density I/O
  std::unique_ptr<shell::InOBDMEStream> density_stream_ptr;
  if (run_parameters.input_info_filename)
  {
    density_stream_ptr = std::make_unique<shell::InOBDMEStreamMulti>(
        run_parameters.input_info_filename.value(), run_parameters.input_filename,
        orbital_space,
        run_parameters.J_bra, run_parameters.g_bra, run_parameters.n_bra,
        run_parameters.J_ket, run_parameters.g_ket, run_parameters.n_ket
      );
  }
  else
  {
    density_stream_ptr = std::make_unique<shell::InOBDMEStreamSingle>(
        run_parameters.input_filename,
        orbital_space,
        run_parameters.J_bra, run_parameters.g_bra, run_parameters.n_bra,
        run_parameters.J_ket, run_parameters.g_ket, run_parameters.n_ket
      );
  }

  // loop over multipoles and write out
  for (auto&& [output_filename, J0] : run_parameters.multipole_filenames)
  {
    basis::OrbitalSectorsLJPN density_sectors;
    basis::OperatorBlocks<double> density_matrices;
    density_stream_ptr->GetMultipole(
        J0, density_sectors, density_matrices
      );
    shell::OutOBMEStream dos(
        output_filename,
        orbital_space, orbital_space, density_sectors,
        basis::OneBodyOperatorType::kSpherical
      );
    dos.Write(density_matrices);
    dos.Close();
  }

  return EXIT_SUCCESS;
}
