/******************************************************************************
  em-gen.cpp

  generate EM matrix element files

  Syntax:
    + em-gen

  Input format:

    set-indexing orbital_filename
    set-basis-scale-factor scale_factor
    define-radial-source radial_filename
    define-target type lambda species output_filename
      type = E|Dl|Ds
      species = p|n

  Patrick J. Fasano
  University of Notre Dame

  + 09/29/17 (pjf): Created, based on orbital-gen.
  + 10/24/17 (pjf): Rewritten, generating E/Dl/Ds, with multiple operators
      generated in one invocation.
  + 11/28/17 (pjf): Print header with version.
  + 04/04/19 (pjf): Write matrix elements in Rose convention, for consistency.

******************************************************************************/

#include <sys/stat.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "am/am.h"
#include "am/halfint.h"
#include "am/wigner_gsl.h"
#include "basis/nlj_orbital.h"
#include "basis/proton_neutron.h"
#include "fmt/format.h"
#include "obme/obme_operator.h"
#include "obme/obme_io.h"

const double kPi = 3.1415926535897;

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

enum class OperatorType : char { kE, kDl, kDs };
std::map<std::string, OperatorType> operator_map = {
    {"E", OperatorType::kE}, {"Dl", OperatorType::kDl}, {"Ds", OperatorType::kDs}};

// Label for radial operator as (type,power,j0,g0,Tz0) of r, k, or o.
typedef std::tuple<shell::RadialOperatorType, int, int, int, int> RadialOperatorLabels;

// Indexing and matrix elements for a radial operator.
//
// Caveat: Just "copying" indexing in here is not good enough, since
// then sectors can be left pointing to deleted temporaries for the
// orbital subspaces.
struct RadialOperatorData {
  RadialOperatorData() = default;

  RadialOperatorData(const std::string& filename_) : filename(filename_)
  {
    // open radial operator file
    basis::OrbitalSpaceLJPN ket_space;
    shell::InOBMEStream stream(filename_);
    stream.SetToIndexing(space, ket_space, sectors);
    auto radial_operator_type = stream.radial_operator_type();
    auto radial_operator_power = stream.radial_operator_power();

    // sanity check on input file
    if (space.OrbitalInfo() != ket_space.OrbitalInfo()) {
      std::cerr << "ERROR: Bra and ket spaces of radial matrix elements are "
                   "not the same."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (stream.operator_type() != basis::OneBodyOperatorType::kRadial) {
      std::cerr << "ERROR: Matrix elements must be radial matrix elements." << std::endl;
      std::exit(EXIT_FAILURE);
    }

    // construct labels for input
    labels = RadialOperatorLabels{radial_operator_type,
                                  radial_operator_power,
                                  sectors.j0(), sectors.g0(), sectors.Tz0()
                                 };

    // read matrices
    stream.Read(matrices);

    // close stream
    stream.Close();
  }

  // operator identification
  RadialOperatorLabels labels;

  // input filename
  std::string filename;

  // operator indexing and storage
  basis::OrbitalSpaceLJPN space;
  basis::OrbitalSectorsLJPN sectors;
  basis::OperatorBlocks<double> matrices;
};

// Map to hold all loaded radial operators.
typedef std::map<RadialOperatorLabels, RadialOperatorData> RadialOperatorMap;

struct TargetChannel {
  OperatorType operator_type;
  int lambda;
  basis::OrbitalSpeciesPN species;
  std::string output_filename;

  TargetChannel(OperatorType op_, int lambda_, basis::OrbitalSpeciesPN species_,
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
      ParsingCheck(line_stream, line_count, line);
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
      ParsingCheck(line_stream, line_count, line);
    } else if (keyword == "define-radial-source") {
      std::string radial_me_filename;
      line_stream >> radial_me_filename;
      ParsingCheck(line_stream, line_count, line);
      struct stat st;
      if (stat(radial_me_filename.c_str(), &st) != 0) {
        std::cerr << "ERROR: file " << radial_me_filename << " does not exist!"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }

      RadialOperatorData radial_operator_data(radial_me_filename);
      run_parameters.radial_operators.insert(
          {radial_operator_data.labels, radial_operator_data});
    } else if (keyword == "define-target") {
      // type lambda species output_filename
      std::string operator_string, output_filename;
      int lambda;
      char species_code;

      line_stream >> operator_string >> lambda >> species_code >> output_filename;
      ParsingCheck(line_stream, line_count, line);

      // operator_type
      auto it = operator_map.find(operator_string);
      if (it == operator_map.end()) {
        std::cerr << "Valid transition types: E, Dl, Ds" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      OperatorType operator_type = it->second;

      // species
      basis::OrbitalSpeciesPN species;
      if (species_code == 'p') {
        species = basis::OrbitalSpeciesPN::kP;
      } else if (species_code == 'n') {
        species = basis::OrbitalSpeciesPN::kN;
      } else {
        std::cerr << "Invalid species code: " << species_code << std::endl;
        std::exit(EXIT_FAILURE);
      }

      // check output file
      struct stat st;
      if (stat(output_filename.c_str(), &st) == 0) {
        std::cerr << "WARN: overwriting file " << output_filename << std::endl;
      }

      run_parameters.targets.emplace_back(operator_type, lambda, species,
                                          output_filename);
    }
  }
}

double CalculatePrefactorE(double scale_factor,
                           const TargetChannel& target,
                           const basis::OrbitalSubspaceLJPN& a,
                           const basis::OrbitalSubspaceLJPN& b) {
  assert(a.orbital_species() == b.orbital_species());
  if (a.orbital_species() != target.species) return 0;

  const int& lambda = target.lambda;

  double factor = std::pow(scale_factor, lambda);
  factor *= 1. / std::sqrt(4. * kPi);
  factor *= ParitySign(b.j() + lambda - HalfInt(1, 2));

  // parity conservation -- enforced by sector enumeration
  // factor *= ((a.l() + b.l() + lambda) % 2 == 0);

  // Racah reduction formula
  factor *= Hat(a.j()) * Hat(lambda)
    * am::Wigner3J(a.j(), b.j(), lambda, HalfInt(1, 2), -HalfInt(1, 2), 0);

  return factor;
}

double CalculatePrefactorM(double scale_factor,
                           const TargetChannel& target,
                           const basis::OrbitalSubspaceLJPN& a,
                           const basis::OrbitalSubspaceLJPN& b) {
  assert(a.orbital_species() == b.orbital_species());
  if (a.orbital_species() != target.species) return 0;

  const int& lambda = target.lambda;

  double factor = std::pow(scale_factor, lambda - 1);
  factor *= 1 / std::sqrt(4. * kPi);
  factor *= ParitySign(b.j() + lambda - HalfInt(1, 2));

  // parity conservation -- enforced by sector enumeration
  // factor *= ((a.l() + b.l() + lambda + 1) % 2 == 0);

  // Racah reduction formula
  factor *= Hat(a.j()) * Hat(lambda)
    * am::Wigner3J(a.j(), b.j(), lambda, HalfInt(1, 2), -HalfInt(1, 2), 0);

  double kappa = double(
      ParitySign(a.l() + a.j() + HalfInt(1, 2)) * (a.j() + HalfInt(1, 2))
      + ParitySign(b.l() + b.j() + HalfInt(1, 2)) * (b.j() + HalfInt(1, 2)));

  factor *= (lambda - kappa);
  if (target.operator_type == OperatorType::kDl) {
    factor *= (1 + kappa / (lambda + 1));
  } else if (target.operator_type == OperatorType::kDs) {
    factor *= -0.5;
  }
  return factor;
}

void GenerateTarget(const RunParameters& run_parameters,
                    const TargetChannel& target) {
  // get radial matrix elements
  int radial_order =
      target.lambda - ((target.operator_type == OperatorType::kDl)
                       || (target.operator_type == OperatorType::kDs));
  int g0 = radial_order % 2;
  int Tz0 = 0;
  RadialOperatorLabels key{shell::RadialOperatorType::kR, radial_order,
                           target.lambda, g0, 0};
  if (run_parameters.radial_operators.count(key) == 0) {
    std::cerr << fmt::format(
                     "ERROR: Missing radial matrix elements: {}^{:d} j0={:d} "
                     "g0={:d} Tz0={:d}",
                     static_cast<char>(shell::RadialOperatorType::kR),
                     radial_order, target.lambda, g0, Tz0)
              << std::endl;
  }
  const RadialOperatorData& radial_operator_data =
      run_parameters.radial_operators.at(key);
  assert(radial_operator_data.space.OrbitalInfo()
         == run_parameters.space.OrbitalInfo());

  // construct new sectors
  basis::OrbitalSectorsLJPN output_sectors(
      run_parameters.space, run_parameters.space, target.lambda, g0, Tz0);
  // default initialization
  basis::OperatorBlocks<double> output_matrices;

  // loop over blocks and apply factor
  for (int sector_index = 0; sector_index < output_sectors.size(); ++sector_index) {
    const auto& output_sector = output_sectors.GetSector(sector_index);
    const auto& bra_subspace = output_sector.bra_subspace();
    const auto& ket_subspace = output_sector.ket_subspace();

    const int bra_subspace_index = output_sector.bra_subspace_index();
    const int ket_subspace_index = output_sector.ket_subspace_index();
    int radial_sector_index = radial_operator_data.sectors.LookUpSectorIndex(
        bra_subspace_index, ket_subspace_index);
    assert(radial_sector_index != basis::kNone);

    double prefactor = 0.;
    if (target.operator_type == OperatorType::kE) {
      prefactor = CalculatePrefactorE(run_parameters.scale_factor, target,
                                      bra_subspace, ket_subspace);
    } else if ((target.operator_type == OperatorType::kDl)
               || (target.operator_type == OperatorType::kDs)) {
      prefactor = CalculatePrefactorM(run_parameters.scale_factor, target,
                                      bra_subspace, ket_subspace);
    }
    output_matrices.push_back(
        prefactor * radial_operator_data.matrices[radial_sector_index]);
  }

  // write operator to file
  shell::OutOBMEStream os(target.output_filename, run_parameters.space,
                          run_parameters.space, output_sectors,
                          basis::OneBodyOperatorType::kSpherical, shell::RadialOperatorType::kGeneric, radial_order);
  os.Write(output_matrices);
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
