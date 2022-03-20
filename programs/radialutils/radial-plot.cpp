/******************************************************************************/
/**
  @file radial-plot.cpp

  generate tables of radial functions, suitable for plotting

  Syntax:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    radial-plot analytic_basis_type scale l j species num_points r_max xform_file output_file
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  analytic_basis_type: oscillator|laguerre

  Patrick J. Fasano
  University of Notre Dame

  + 09/25/18 (pjf): Created, based on natorb-gen.
  + 05/09/19 (pjf): Use std::size_t for basis indices and sizes.

******************************************************************************/

#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>

#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Eigenvalues"

#include "am/halfint.h"
#include "basis/nlj_orbital.h"
#include "basis/operator.h"
#include "density/obdme_io.h"
#include "mcutils/eigen.h"
#include "mcutils/parsing.h"
#include "obme/obme_io.h"
#include "spline/wavefunction_basis.h"

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

std::map<std::string, shell::RadialBasisType> kBasisTypeDefinitions(
    {{"oscillator", shell::RadialBasisType::kOscillator},
     {"laguerre", shell::RadialBasisType::kLaguerre}});

// Stores simple parameters for run
struct RunParameters {
  // filenames
  shell::RadialBasisType basis_type;
  double scale;
  int l;
  HalfInt j;
  basis::OrbitalSpeciesPN species;
  int num_points;
  double r_max;
  std::string xform_file;
  std::string output_file;
};

void PrintUsage(const char **argv) {
  std::cout << "Usage: " << argv[0]
            << " analytic_basis_type scale l j species num_points r_max xform_file output_file"
            << std::endl;
}

void ProcessArguments(int argc, const char *argv[], RunParameters& run_parameters) {
  // argument counter
  int arg = 0;

  // usage message
  if (argc-1 != 9) {
    PrintUsage(argv);
    std::cerr << "ERROR: Incorrect number of arguments." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // analytic basis type
  {
    std::string basis_type_str;
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> basis_type_str;
    // convert basis
    if (kBasisTypeDefinitions.count(basis_type_str) == 0) {
      std::cerr <<  "Valid analytic basis types: oscillator|laguerre" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    run_parameters.basis_type = kBasisTypeDefinitions.at(basis_type_str);
  }

  // scale
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.scale;
  }

  // l
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.l;
  }

  // j
  {
    std::istringstream parameter_stream(argv[++arg]);
    float j;
    parameter_stream >> j;
    run_parameters.j = HalfInt(2*j, 2);
  }

  // species
  {
    std::istringstream parameter_stream(argv[++arg]);
    int species;
    parameter_stream >> species;
    run_parameters.species = basis::OrbitalSpeciesPN(species);
  }

  // num points
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.num_points;
  }

  // r max
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.r_max;
  }

  // xform file
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.xform_file;
    mcutils::FileExistCheck(run_parameters.xform_file, true, false);
  }

  // output file
  {
    std::istringstream parameter_stream(argv[++arg]);
    parameter_stream >> run_parameters.output_file;
    mcutils::FileExistCheck(run_parameters.output_file, false, true);
  }
}

////////////////////////////////////////////////////////////////
// main program
////////////////////////////////////////////////////////////////
int main(int argc, const char *argv[]) {
  // header
  std::cout << std::endl;
  std::cout << "radial-plot -- radial function plotting" << std::endl;
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  RunParameters run_parameters;
  ProcessArguments(argc, argv, run_parameters);

  // get indexing and blocks
  basis::OrbitalSpaceLJPN bra_space, ket_space;
  basis::OrbitalSectorsLJPN sectors;
  basis::OperatorBlocks<double> xform_matrices;
  shell::InOBMEStream xform_s(run_parameters.xform_file);
  xform_s.SetToIndexing(bra_space, ket_space, sectors);
  xform_s.Read(xform_matrices);
  xform_s.Close();

  // check that xform is valid
  assert(xform_s.operator_type() == basis::OneBodyOperatorType::kRadial);
  assert(sectors.j0() == 0);
  assert(sectors.g0() == 0);
  assert(sectors.Tz0() == 0);

  // get subspace for desired quantum numbers
  const basis::OrbitalSubspaceLJPNLabels subspace_labels =
    {run_parameters.species, run_parameters.l, run_parameters.j};
  const std::size_t subspace_index = bra_space.LookUpSubspaceIndex(subspace_labels);
  const auto subspace = bra_space.GetSubspace(subspace_index);

  // get xform matrix for desired quantum numbers
  const std::size_t sector_index = sectors.LookUpSectorIndex(subspace_index, subspace_index);

  // get spline wave function type
  spline::Basis basis_type;
  if (run_parameters.basis_type == shell::RadialBasisType::kOscillator) {
    basis_type = spline::Basis::HC;
  } else if (run_parameters.basis_type == shell::RadialBasisType::kLaguerre) {
    basis_type = spline::Basis::LC;
  }

  // construct matrix representation of analytic basis
  auto r_values = Eigen::VectorXd::LinSpaced(run_parameters.num_points, 0, run_parameters.r_max);
  Eigen::MatrixXd wf_values(run_parameters.num_points, subspace.size());
  for (std::size_t state_index = 0; state_index < subspace.size(); ++state_index)
  {
    basis::OrbitalStateLJPN state(subspace, state_index);
    spline::WaveFunction wavefunction(state.n(), state.l(), run_parameters.scale, basis_type);
    #pragma omp parallel for
    for (int i = 0; i < run_parameters.num_points; ++i)
    {
      const auto& r = r_values(i);
      wf_values(i, state_index) = spline::basis::WaveFunctionValue(r, wavefunction);
    }
  }

  // apply xform matrix
  wf_values = wf_values * xform_matrices[sector_index];

  // generate output matrix
  Eigen::MatrixXd output_matrix(run_parameters.num_points, subspace.size()+1);
  output_matrix.block(0, 0, run_parameters.num_points, 1) = r_values;
  output_matrix.block(0, 1, run_parameters.num_points, subspace.size()) = wf_values;

  // open output file for writing
  std::ofstream out_s(run_parameters.output_file);
  out_s << output_matrix << std::endl;
  out_s.close();

  /* return code */
  return EXIT_SUCCESS;
}
