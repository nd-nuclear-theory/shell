/******************************************************************************
  em-gen.cpp

  generate EM matrix element files

  Syntax:
    + em-gen

  Input format:

    set-indexing <orbital_filename>
    set-basis-scale-factor <scale_factor>
    define-radial-source <type> <power> <radial_filename>
      type = r|k
    define-am-source <type> <filename>
      type = l|s
    define-target <type> <lambda> <species> <output_filename>
      type = E|Dl|Ds
      species = p|n|total

  Patrick J. Fasano
  University of Notre Dame

  + 09/29/17 (pjf): Created, based on orbital-gen.
  + 10/24/17 (pjf): Rewritten, generating E/Dl/Ds, with multiple operators
      generated in one invocation.
  + 11/28/17 (pjf): Print header with version.
  + 04/04/19 (pjf): Write matrix elements in Rose convention, for consistency.
  + 05/09/19 (pjf): Use std::size_t for basis indices and sizes.
  + 08/17/19 (pjf): Refactor to use only spherical tensor operators, such as
      solid harmonics, l, and s. (Fixes #4.)
  + 08/27/19 (pjf):
    - Change prefactors to use solid harmonic convention.
    - Use shell::kCharCodeRadialOperatorType instead of static_cast.

******************************************************************************/

#include <sys/stat.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "am/am.h"
#include "am/halfint.h"
#include "am/wigner_gsl.h"
#include "am/rme.h"
#include "basis/nlj_orbital.h"
#include "basis/proton_neutron.h"
#include "fmt/format.h"
#include "obme/obme.h"
#include "obme/obme_operator.h"
#include "obme/obme_io.h"

const double kPi = 3.1415926535897;

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

enum class OperatorType : char { kE, kDl, kDs };
std::map<std::string, OperatorType> operator_map = {
    {"E", OperatorType::kE}, {"Dl", OperatorType::kDl}, {"Ds", OperatorType::kDs}};

// Label for radial operator as (type,power).
typedef std::tuple<shell::RadialOperatorType, int> RadialOperatorLabels;
typedef std::tuple<am::AngularMomentumOperatorType> AngularMomentumOperatorLabels;

// Indexing and matrix elements for a radial operator.
//
// Caveat: Just "copying" indexing in here is not good enough, since
// then sectors can be left pointing to deleted temporaries for the
// orbital subspaces.
struct OperatorData {
  OperatorData() = default;

  OperatorData(const std::string& filename_) : filename(filename_)
  {
    // open radial operator file
    basis::OrbitalSpaceLJPN ket_space;
    shell::InOBMEStream stream(filename_);
    stream.SetToIndexing(space, ket_space, sectors);

    // sanity check on input file
    if (space.OrbitalInfo() != ket_space.OrbitalInfo()) {
      std::cerr << "ERROR: Bra and ket spaces of radial matrix elements are "
                   "not the same."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (stream.operator_type() != basis::OneBodyOperatorType::kSpherical) {
      std::cerr << "ERROR: Matrix elements must be of a spherical tensor operator." << std::endl;
      std::exit(EXIT_FAILURE);
    }

    // read matrices
    stream.Read(matrices);

    // close stream
    stream.Close();
  }

  // input filename
  std::string filename;

  // operator indexing and storage
  basis::OrbitalSpaceLJPN space;
  basis::OrbitalSectorsLJPN sectors;
  basis::OperatorBlocks<double> matrices;
};

// Map to hold all loaded radial and angular momentum operators
typedef std::map<RadialOperatorLabels, OperatorData> RadialOperatorMap;
typedef std::map<AngularMomentumOperatorLabels, OperatorData> AngularMomentumOperatorMap;

struct TargetChannel {
  OperatorType operator_type;
  int lambda;
  basis::OperatorTypePN species;
  std::string output_filename;

  TargetChannel(OperatorType op_, int lambda_, basis::OperatorTypePN species_,
                const std::string& output_filename_)
      : operator_type(op_),
        lambda(lambda_),
        species(species_),
        output_filename(output_filename_) {}
};

// Stores simple parameters for run
struct RunParameters {
  basis::OrbitalSpaceLJPN space;
  double scale_factor;

  RadialOperatorMap radial_operators;
  AngularMomentumOperatorMap am_operators;
  std::vector<TargetChannel> targets;
};

void PrintUsage(const char** argv) {
  std::cout << "Usage: " << argv[0] << std::endl;
}

void ReadParameters(RunParameters& run_parameters) {
  std::string line;
  int line_count = 0;

  while (std::getline(std::cin, line)) {
    ++line_count;

    // set up for line parsing
    std::istringstream line_stream(line);
    std::string keyword;
    line_stream >> keyword;

    // skip blank line or hash comment line
    if ((keyword == "") || (keyword == "#")) continue;

    // select action based on keyword
    if (keyword == "set-indexing") {
      std::string orbital_filename;
      line_stream >> orbital_filename;
      mcutils::ParsingCheck(line_stream, line_count, line);
      struct stat st;
      if (stat(orbital_filename.c_str(), &st) != 0) {
        std::cerr << "ERROR: file " << orbital_filename << " does not exist!"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }

      std::ifstream orbital_stream(orbital_filename);
      auto input_orbitals = basis::ParseOrbitalPNStream(orbital_stream, true);
      run_parameters.space = basis::OrbitalSpaceLJPN(input_orbitals);
    } else if (keyword == "set-basis-scale-factor") {
      line_stream >> run_parameters.scale_factor;
      mcutils::ParsingCheck(line_stream, line_count, line);
    } else if (keyword == "define-radial-source") {
      std::string radial_me_filename;
      std::string radial_type;
      int radial_power;
      line_stream >> radial_type >> radial_power >> radial_me_filename;
      mcutils::ParsingCheck(line_stream, line_count, line);
      RadialOperatorLabels labels
        = {shell::kCharCodeRadialOperatorType.at(radial_type), radial_power};

      mcutils::FileExistCheck(radial_me_filename, true, false);
      OperatorData operator_data(radial_me_filename);

      run_parameters.radial_operators.insert({labels, operator_data});
    } else if (keyword == "define-am-source") {
      std::string am_me_filename;
      char am_type;
      line_stream >> am_type >> am_me_filename;
      mcutils::ParsingCheck(line_stream, line_count, line);
      AngularMomentumOperatorLabels labels
        = {static_cast<am::AngularMomentumOperatorType>(am_type)};

      mcutils::FileExistCheck(am_me_filename, true, false);
      OperatorData operator_data(am_me_filename);

      run_parameters.am_operators.insert({labels, operator_data});
    } else if (keyword == "define-target") {
      // type lambda species output_filename
      std::string operator_string, output_filename;
      int lambda;
      std::string species_code;

      line_stream >> operator_string >> lambda >> species_code >> output_filename;
      mcutils::ParsingCheck(line_stream, line_count, line);

      // operator_type
      auto it = operator_map.find(operator_string);
      if (it == operator_map.end()) {
        std::cerr << "Valid transition types: E, Dl, Ds" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      OperatorType operator_type = it->second;

      // species
      basis::OperatorTypePN species = basis::kCharCodeOperatorTypePN.at(species_code);

      // check output file
      mcutils::FileExistCheck(output_filename, false, true);

      run_parameters.targets.emplace_back(
          operator_type, lambda, species, output_filename
        );
    }
  }
}

double CalculatePrefactorE(
    double scale_factor,
    const TargetChannel& target,
    const basis::OrbitalSubspaceLJPN& a,
    const basis::OrbitalSubspaceLJPN& b
  )
{
  assert(a.orbital_species() == b.orbital_species());
  if (target.species != basis::OperatorTypePN::kTotal)
  {
    if (static_cast<int>(a.orbital_species()) != static_cast<int>(target.species)) return 0.;
  }

  const int& lambda = target.lambda;

  double factor = std::pow(scale_factor, lambda);
  return factor;
}

void GenerateTargetE(
    const RunParameters& run_parameters,
    const TargetChannel& target
  )
{
  int radial_order = target.lambda;
  int j0 = target.lambda;
  int g0 = radial_order % 2;
  int Tz0 = 0;

  // construct output sectors
  basis::OrbitalSectorsLJPN output_sectors(
      run_parameters.space, run_parameters.space,
      j0, g0, Tz0
    );
  // declare output storage
  basis::OperatorBlocks<double> output_matrices;

  // get solid harmonic operator
  RadialOperatorLabels radial_key{shell::RadialOperatorType::kR, radial_order};
  if (run_parameters.radial_operators.count(radial_key) == 0) {
    std::cerr << fmt::format(
                    "ERROR: Missing radial matrix elements: {}^{:d} j0={:d} "
                    "g0={:d} Tz0={:d}",
                    static_cast<char>(shell::RadialOperatorType::kR),
                    radial_order, radial_order, g0, Tz0
                    )
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const OperatorData& radial_operator_data = run_parameters.radial_operators.at(radial_key);
  assert(radial_operator_data.space.OrbitalInfo() == run_parameters.space.OrbitalInfo());

  // populate matrices with unscaled matrix elements
  output_matrices = radial_operator_data.matrices;

  // scale operator and apply species constraint
  for (std::size_t sector_index = 0; sector_index < output_sectors.size(); ++sector_index) {
    const auto& sector = output_sectors.GetSector(sector_index);
    output_matrices[sector_index] *= CalculatePrefactorE(
        run_parameters.scale_factor, target,
        sector.bra_subspace(), sector.ket_subspace()
      );
  }

  // write operator to file
  shell::OutOBMEStream os(
      target.output_filename,
      run_parameters.space, run_parameters.space, output_sectors,
      basis::OneBodyOperatorType::kSpherical
    );
  os.Write(output_matrices);
}

double CalculatePrefactorM(
    double scale_factor,
    const TargetChannel& target,
    const basis::OrbitalSubspaceLJPN& a,
    const basis::OrbitalSubspaceLJPN& b
  )
{
  assert(a.orbital_species() == b.orbital_species());
  if (target.species != basis::OperatorTypePN::kTotal)
  {
    if (static_cast<int>(a.orbital_species()) != static_cast<int>(target.species)) return 0.;
  }

  const int& lambda = target.lambda;

  // radial scaling
  double factor = std::pow(scale_factor, lambda - 1);
  // gradient of solid harmonic -- Heyde (1990), eq. 4.12
  factor *= Hat(lambda) * std::sqrt(lambda);

  if (target.operator_type == OperatorType::kDl) {
    factor *= (2 / (lambda + 1));
  }
  return factor;
}

void GenerateTargetM(
    const RunParameters& run_parameters,
    const TargetChannel& target
  )
{
  int radial_order = target.lambda - 1;
  int j0 = target.lambda;
  int g0 = radial_order % 2;
  int Tz0 = 0;

  // construct output sectors
  basis::OrbitalSectorsLJPN output_sectors(
      run_parameters.space, run_parameters.space,
      j0, g0, Tz0
    );
  // declare output storage
  basis::OperatorBlocks<double> output_matrices;

  // get angular momentum operator
  AngularMomentumOperatorLabels am_key;
  if (target.operator_type == OperatorType::kDl)
    am_key = {am::AngularMomentumOperatorType::kOrbital};
  else if (target.operator_type == OperatorType::kDs)
    am_key = {am::AngularMomentumOperatorType::kSpin};

  if (run_parameters.am_operators.count(am_key) == 0) {
    std::cerr << fmt::format(
                    "ERROR: Missing angular-momentum matrix elements: {}",
                    static_cast<char>(am::AngularMomentumOperatorType::kOrbital)
                  )
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const OperatorData& am_operator_data = run_parameters.am_operators.at(am_key);

  // populate matrices with unscaled matrix elements
  if (radial_order == 0)
  {
    output_matrices = am_operator_data.matrices;
  }
  else
  {
    RadialOperatorLabels radial_key{shell::RadialOperatorType::kR, radial_order};
    if (run_parameters.radial_operators.count(radial_key) == 0) {
      std::cerr << fmt::format(
                      "ERROR: Missing radial matrix elements: {}^{:d} j0={:d} "
                      "g0={:d} Tz0={:d}",
                      static_cast<char>(shell::RadialOperatorType::kR),
                      radial_order, radial_order, g0, Tz0
                      )
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    const OperatorData& radial_operator_data = run_parameters.radial_operators.at(radial_key);
    assert(radial_operator_data.space.OrbitalInfo() == run_parameters.space.OrbitalInfo());

    // calculate product operator
    shell::OneBodyOperatorTensorProduct(
        run_parameters.space,
        am_operator_data.sectors, am_operator_data.matrices,
        radial_operator_data.sectors, radial_operator_data.matrices,
        output_sectors, output_matrices
      );
  }

  // scale operator and apply species constraint
  for (std::size_t sector_index = 0; sector_index < output_sectors.size(); ++sector_index) {
    const auto& sector = output_sectors.GetSector(sector_index);
    output_matrices[sector_index] *= CalculatePrefactorM(
        run_parameters.scale_factor, target,
        sector.bra_subspace(), sector.ket_subspace()
      );
  }

  // write operator to file
  shell::OutOBMEStream os(
      target.output_filename,
      run_parameters.space, run_parameters.space, output_sectors,
      basis::OneBodyOperatorType::kSpherical
    );
  os.Write(output_matrices);
}

void GenerateTarget(
    const RunParameters& run_parameters,
    const TargetChannel& target
  )
{
  if (target.operator_type == OperatorType::kE)
  {
    GenerateTargetE(run_parameters, target);
  }
  else if (
      (target.operator_type == OperatorType::kDl)
      ||(target.operator_type == OperatorType::kDs)
    )
  {
    GenerateTargetM(run_parameters, target);
  }
}

int main(int argc, const char* argv[]) {
  if (argc != 1) PrintUsage(argv);

  // header
  std::cout << std::endl;
  std::cout << "em-gen -- electromagnetic matrix element generation" << std::endl;
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  RunParameters run_parameters;
  ReadParameters(run_parameters);

  for (const auto& target : run_parameters.targets) {
    GenerateTarget(run_parameters, target);
  }

  /* return code */
  return EXIT_SUCCESS;
}
