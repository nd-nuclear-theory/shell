/******************************************************************************/
/**
  @file orbital-extract.cpp

  extract orbital info from OBME file

  Syntax:
    + orbital-extract --weight <obme_file> <output_file>
    + orbital-extract --diagonal <obme_file> <output_file>

  Patrick J. Fasano
  Argonne National Laboratory

  + 01/09/24 (pjf): Created, based on orbital-gen.
*/
/******************************************************************************/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>

#include <fmt/format.h>
#include <fmt/os.h>

#include "basis/nlj_orbital.h"
#include "mcutils/io.h"
#include "mcutils/parsing.h"
#include "obme/obme_io.h"

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

enum class ExtractionMode {kWeight, kDiagonal};

// Stores simple parameters for run
struct RunParameters {
  // filename
  std::string input_filename;
  std::string output_filename;
  // mode
  ExtractionMode mode;
};

void PrintUsage(const int argc, const char **argv) {
  fmt::print(
      stderr, "Usage: {} <--weight | --diagonal> <obme_file> <output_file>\n",
      (argc > 0) ? argv[0] : "orbital-extract"
    );
}

RunParameters ProcessArguments(int argc, const char *argv[]) {
  RunParameters run_parameters;

  // argument counter
  int arg = 1;

  // usage message
  if (argc-1 != 3) {
    PrintUsage(argc, argv);
    std::cerr << "Incorrect number of arguments." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // operation mode
  {
    std::istringstream parameter_stream(argv[arg++]);
    if (parameter_stream.str() == "--weight") {
      run_parameters.mode = ExtractionMode::kWeight;
    } else if (parameter_stream.str() == "--diagonal") {
      run_parameters.mode = ExtractionMode::kDiagonal;
    } else {
      PrintUsage(argc, argv);
      fmt::print(stderr, "Invalid operation mode.\n");
      std::exit(EXIT_FAILURE);
    }
  }

  // input filename
  {
    std::istringstream parameter_stream(argv[arg++]);
    parameter_stream >> run_parameters.input_filename;
    if (!parameter_stream) {
      PrintUsage(argc, argv);
      fmt::print(stderr, "Invalid input filename\n");
      std::exit(EXIT_FAILURE);
    }
    mcutils::FileExistCheck(run_parameters.input_filename, true, false);
  }

  // output filename
  {
    std::istringstream parameter_stream(argv[arg++]);
    parameter_stream >> run_parameters.output_filename;
    if (!parameter_stream) {
      PrintUsage(argc, argv);
      fmt::print(stderr, "Invalid output filename.\n");
      std::exit(EXIT_FAILURE);
    }
    mcutils::FileExistCheck(run_parameters.output_filename, false, true);
  }

  return run_parameters;
}

int main(int argc, const char *argv[]) {
  // header
  fmt::print("\n");
  fmt::print("orbital-extract -- orbital info extraction\n");
  fmt::print("version: {}\n", VCS_REVISION);
  fmt::print("\n");

  auto run_parameters = ProcessArguments(argc, argv);
  std::vector<basis::OrbitalPNInfo> orbitals_info;

  // read obme file
  auto input_stream = shell::InOBMEStream(run_parameters.input_filename);
  basis::OrbitalSpaceLJPN space, ket_space;
  basis::OrbitalSectorsLJPN sectors;
  basis::OperatorBlocks<double> matrices;
  input_stream.SetToIndexing(space, ket_space, sectors);
  if (space.OrbitalInfo() != ket_space.OrbitalInfo())
  {
    fmt::print(
        "ERROR: Bra and ket orbitals are not equivalent. "
        "Extraction is ambiguous for xform files.\n"
      );
    std::exit(EXIT_FAILURE);
  }

  if (run_parameters.mode == ExtractionMode::kWeight)
  {
    orbitals_info = space.OrbitalInfo();
  }
  else if (run_parameters.mode == ExtractionMode::kDiagonal)
  {
    if (!sectors.IsScalar())
    {
      fmt::print("ERROR: Diagonal extraction is ambiguous for nonscalar operator.\n");
      std::exit(EXIT_FAILURE);
    }

    input_stream.Read(matrices);
    for (int sector_index=0; sector_index<sectors.size(); ++sector_index)
    {
      const auto& sector = sectors.GetSector(sector_index);
      assert(sector.IsDiagonal());
      const auto& subspace = sector.bra_subspace();
      for (int i=0; i<subspace.size(); ++i)
      {
        auto info = subspace.GetState(i).OrbitalInfo();
        info.weight = matrices[sector_index](i,i);
        orbitals_info.push_back(info);
      }
    }
  }

  if (orbitals_info.size() == 0) {
    PrintUsage(argc, argv);
    fmt::print(stderr, "ERROR: No orbitals extracted.\n");
    std::exit(EXIT_FAILURE);
  }

  auto output_file = fmt::output_file(run_parameters.output_filename);
  for (const auto& info : orbitals_info)
  {
    output_file.print(
        "  {:>4d}  {:>4d}  {:>4.1f}  {:>+4.1f}  {:>15.8e}\n",
        info.n, info.l, float(info.j),
        float(basis::kOrbitalSpeciesPNCodeTz[static_cast<int>(info.orbital_species)]),
        info.weight
      );
  }
  output_file.close();

  /* return code */
  return EXIT_SUCCESS;
}
