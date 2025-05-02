/******************************************************************************/
/*
  @file obsolve

  Solve eigenproblem for one-body Hamiltonian

  Syntax:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    obsolve <input_filename> <res_filename> <xform_filename>
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Patrick J. Fasano
    Argonne National Laboratory

    + 01/07/24 (pjf): Created.
*/
/******************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <optional>
#include <vector>
#include <cmath>
#include <ranges>

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

// Stores simple parameters for run
struct RunParameters {
  // filenames
  std::string input_filename;
  std::string res_filename;
  std::string xform_filename;
};

void PrintUsage(const int argc, const char* argv[])
{
  fmt::print(
      stderr,
      "Usage: {} input_filename res_filename xform_filename\n",
      (argc > 0) ? argv[0] : "obsolve"
    );
}

RunParameters ProcessArguments(const int argc, const char* argv[])
{
  RunParameters run_parameters;

  // usage message
  if (argc-1 != 3)
  {
      PrintUsage(argc, argv);
      fmt::print(stderr, "ERROR: Wrong number of arguments\n");
      std::exit(EXIT_FAILURE);
  }

  int count = 1;

  // input operator
  {
    run_parameters.input_filename = argv[count++];
    mcutils::FileExistCheck(run_parameters.input_filename, true, false);
  }

  // result file
  {
    run_parameters.res_filename = argv[count++];
    mcutils::FileExistCheck(run_parameters.res_filename, false, true);
  }

  // xform (eigenvector) file
  {
    run_parameters.xform_filename = argv[count++];
    mcutils::FileExistCheck(run_parameters.xform_filename, false, true);
  }

  return run_parameters;
}

auto DiagonalizeOperator(
    const basis::OrbitalSpaceLJPN& space,
    const basis::OrbitalSectorsLJPN& sectors,
    const basis::OperatorBlocks<double>& matrices
  )
{
  struct return_t {std::vector<basis::OrbitalPNInfo> eigenvalue_info; basis::OperatorBlocks<double> eigenvector_matrices;};
  return_t retval{};
  auto& [eigenvalue_info, eigenvector_matrices] = retval;

  for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
  {
    const auto& matrix = matrices[sector_index];
    const auto& sector = sectors.GetSector(sector_index);
    const auto& species = sector.bra_subspace().orbital_species();
    const auto& l = sector.bra_subspace().l();
    const auto& j = sector.bra_subspace().j();

    if (((matrix-matrix.transpose()).array().abs() > 1e-6).any())
    {
      fmt::print(stderr, "WARN: sector {} may not be symmetric", sector_index);
    }

    auto eigensolver = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(matrix);

    eigenvector_matrices.push_back(eigensolver.eigenvectors());
    for (int i=0; i < eigensolver.eigenvalues().size(); ++i)
    {
      eigenvalue_info.emplace_back(
          species, i, l, j, eigensolver.eigenvalues()(i)
        );
    }
  }

  return retval;
}

int main(int argc, const char* argv[])
{
  // header
  fmt::print("\n");
  fmt::print("obsolve -- one-body eigensolver\n");
  fmt::print("version: {}\n", VCS_REVISION);
  fmt::print("\n");

  auto run_parameters = ProcessArguments(argc, argv);

  // read operator
  auto input_stream = shell::InOBMEStream(run_parameters.input_filename);

  basis::OrbitalSpaceLJPN space{};
  basis::OrbitalSectorsLJPN sectors{};

  {
    basis::OrbitalSpaceLJPN ket_space{};
    input_stream.SetToIndexing(space, ket_space, sectors);
    if (space.OrbitalInfo() != ket_space.OrbitalInfo())
    {
      fmt::print(stderr, "ERROR: Bra and ket space of this operator are not the same.");
      std::exit(EXIT_FAILURE);
    }
  }

  if ((sectors.J0() != 0) || (sectors.g0() != 0) || (sectors.Tz0() != 0))
  {
    fmt::print(
      stderr,
      "ERROR: Cannot diagonalize nonscalar operator with J0={} g0={} Tz0{}",
      sectors.J0(), sectors.g0(), sectors.Tz0()
    );
    std::exit(EXIT_FAILURE);
  }

  basis::OperatorBlocks<double> operator_matrices;
  input_stream.Read(operator_matrices);
  input_stream.Close();

  const auto [eigenvalue_info, eigenvector_matrices] = DiagonalizeOperator(
      space, sectors, operator_matrices
    );

  auto res_file = fmt::output_file(run_parameters.res_filename);
  for (const auto& eigenvalue : eigenvalue_info)
  {
    res_file.print(
        "  {:>4d}  {:>4d}  {:>4.1f}  {:>+4.1f}  {:>15.8e}\n",
        eigenvalue.n, eigenvalue.l, float(eigenvalue.j),
        float(basis::kOrbitalSpeciesPNCodeTz[static_cast<int>(eigenvalue.orbital_species)]),
        eigenvalue.weight
      );
  }
  res_file.close();

  auto xform_stream = shell::OutOBMEStream(
      run_parameters.xform_filename,
      space, space, sectors,
      basis::OneBodyOperatorType::kRadial
    );
  xform_stream.Write(eigenvector_matrices);
  xform_stream.Close();

  return EXIT_SUCCESS;
}
