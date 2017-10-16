/******************************************************************************
  em-gen.cpp

  generate EM matrix element files

  Syntax:
    + em-gen

  Input format:

    set-transition-type type
      type = E|M
    set-transition-rank lambda
    set-basis-scale-factor scale_factor
    define-radial-source radial_filename
    define-output output_filename
    set-charge charge species value (optional)
      charge = q|gl|gs
      species = p|n

  Patrick J. Fasano
  University of Notre Dame

  + 09/29/17 (pjf): Created, based on orbital-gen.

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
#include "radial/radial_io.h"

const double kPi = 3.1415926535897;
const double kQp = 1.;
const double kQn = 0.;
const double kGp = 5.585694702;  // NIST CODATA 2014
const double kGn = -3.82608545;  // NIST CODATA 2014
const double kGlp = 1.;
const double kGln = 0.;

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

enum class TransitionType : char { kE = 'E', kM = 'M' };

// Stores simple parameters for run
struct RunParameters {
  // filename
  std::string orbital_filename;
  std::string radial_me_filename;
  std::string output_filename;
  // mode
  TransitionType transition_type;
  int lambda;
  double scale_factor;
  std::array<double, 2> electric_charge;
  std::array<double, 2> g;
  std::array<double, 2> orbital_g;

  /** default constructor */
  RunParameters()
      : radial_me_filename(""),
        output_filename(""),
        transition_type(TransitionType::kE),
        lambda(0),
        scale_factor(1.0),
        electric_charge{kQp, kQn},
        g{kGp, kGn},
        orbital_g{kGlp, kGln} {}
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
    if (keyword == "set-transition-type") {
      char transition_type;
      line_stream >> transition_type;
      if ((transition_type != 'E') && (transition_type != 'M')) {
        std::cerr << "Valid transition types: E, M" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      run_parameters.transition_type = static_cast<TransitionType>(transition_type);
    } else if (keyword == "set-transition-rank") {
      line_stream >> run_parameters.lambda;
      ParsingCheck(line_stream, line_count, line);
    } else if (keyword == "set-basis-scale-factor") {
      line_stream >> run_parameters.scale_factor;
      ParsingCheck(line_stream, line_count, line);
    } else if (keyword == "define-radial-source") {
      line_stream >> run_parameters.radial_me_filename;
      ParsingCheck(line_stream, line_count, line);
      struct stat st;
      if (stat(run_parameters.radial_me_filename.c_str(), &st) != 0) {
        std::cerr << "ERROR: file " << run_parameters.radial_me_filename
                  << " does not exist!" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    } else if (keyword == "define-output") {
      line_stream >> run_parameters.output_filename;
      ParsingCheck(line_stream, line_count, line);
      struct stat st;
      if (stat(run_parameters.output_filename.c_str(), &st) == 0) {
        std::cerr << "WARN: overwriting file " << run_parameters.output_filename
                  << std::endl;
      }
    } else if (keyword == "set-charge") {
      std::string charge, species_code;
      double value;
      line_stream >> charge >> species_code >> value;
      ParsingCheck(line_stream, line_count, line);

      int species_index;
      if (species_code == "p" || species_code == "P") {
        species_index = 0;
      } else if (species_code == "n" || species_code == "N") {
        species_index = 1;
      }

      if (charge == "q") {
        run_parameters.electric_charge[species_index] = value;
      } else if (charge == "gl") {
        run_parameters.orbital_g[species_index] = value;
      } else if (charge == "gs") {
        run_parameters.g[species_index] = value;
      }
    }
  }
}

double CalculatePrefactorE(const RunParameters& run_parameters,
                           const basis::OrbitalSubspaceLJPN& a,
                           const basis::OrbitalSubspaceLJPN& b) {
  assert(a.orbital_species() == b.orbital_species());
  int lambda = run_parameters.lambda;
  double charge =
      run_parameters.electric_charge[static_cast<int>(a.orbital_species())];
  double factor = std::pow(run_parameters.scale_factor, lambda);

  factor *= charge / std::sqrt(4. * kPi);
  factor *= ParitySign(b.j() + lambda - HalfInt(1, 2));

  // parity conservation -- enforced by sector enumeration
  // factor *= ((a.l() + b.l() + lambda) % 2 == 0);

  // Racah reduction formula
  //
  // Note: Suhonen 6.23 has an extra factor of Hat(lambda), due to difference
  // from group theory convention
  factor *= Hat(a.j()) * Hat(b.j())
            * am::Wigner3J(a.j(), b.j(), lambda, HalfInt(1, 2), -HalfInt(1, 2), 0);

  return factor;
}

double CalculatePrefactorM(const RunParameters& run_parameters,
                           const basis::OrbitalSubspaceLJPN& a,
                           const basis::OrbitalSubspaceLJPN& b) {
  assert(a.orbital_species() == b.orbital_species());
  int lambda = run_parameters.lambda;
  double gl = run_parameters.orbital_g[static_cast<int>(a.orbital_species())];
  double gs = run_parameters.g[static_cast<int>(a.orbital_species())];

  double factor = std::pow(run_parameters.scale_factor, lambda - 1);

  factor *= 1 / std::sqrt(4. * kPi);
  factor *= ParitySign(b.j() + lambda - HalfInt(1, 2));

  // parity conservation -- enforced by sector enumeration
  // factor *= ((a.l() + b.l() + lambda + 1) % 2 == 0);

  // Racah reduction formula
  //
  // Note: Suhonen 6.23 has an extra factor of Hat(lambda), due to difference
  // from group theory convention
  factor *= Hat(a.j()) * Hat(b.j())
            * am::Wigner3J(a.j(), b.j(), lambda, HalfInt(1, 2), -HalfInt(1, 2), 0);

  double kappa = double(
      ParitySign(a.l() + a.j() + HalfInt(1, 2)) * (a.j() + HalfInt(1, 2))
      + ParitySign(b.l() + b.j() + HalfInt(1, 2)) * (b.j() + HalfInt(1, 2)));

  factor *= (lambda - kappa);
  factor *= (gl * (1 + kappa / (lambda + 1)) - (gs / 2.));
  return factor;
}

int main(int argc, const char* argv[]) {
  if (argc != 1) PrintUsage(argv);

  RunParameters run_parameters;
  ReadParameters(run_parameters);

  // Read radial matrix elements
  shell::InRadialStream radial_me_s(run_parameters.radial_me_filename);

  // get indexing
  basis::OrbitalSpaceLJPN bra_space, ket_space;
  basis::OrbitalSectorsLJPN radial_sectors;
  shell::RadialOperatorType radial_operator_type =
      radial_me_s.radial_operator_type();
  int radial_operator_power = radial_me_s.radial_operator_power();
  radial_me_s.SetToIndexing(bra_space, ket_space, radial_sectors);

  if (bra_space.OrbitalInfo() != ket_space.OrbitalInfo()) {
    std::cerr << "ERROR: Bra and ket spaces of radial matrix elements are not "
                 "the same. Cannot transform."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  int radial_order = run_parameters.lambda
                     - (run_parameters.transition_type == TransitionType::kM);
  assert(radial_operator_type == shell::RadialOperatorType::kR);
  assert(radial_operator_power == radial_order);
  assert(radial_sectors.Tz0() == 0);

  // get radial matrix elements
  basis::OperatorBlocks<double> radial_matrices;
  radial_me_s.Read(radial_matrices);
  radial_me_s.Close();

  // construct new sectors
  basis::OperatorBlocks<double> operator_matrices;
  int g0 = radial_order % 2;
  basis::OrbitalSectorsLJPN operator_sectors(bra_space, ket_space,
                                             run_parameters.lambda, g0, 0);

  for (int sector_index = 0; sector_index < operator_sectors.size(); ++sector_index) {
    const auto& operator_sector = operator_sectors.GetSector(sector_index);
    const auto& bra_subspace = operator_sector.bra_subspace();
    const auto& ket_subspace = operator_sector.ket_subspace();

    const int bra_subspace_index = operator_sector.bra_subspace_index();
    const int ket_subspace_index = operator_sector.ket_subspace_index();
    int radial_sector_index =
        radial_sectors.LookUpSectorIndex(bra_subspace_index, ket_subspace_index);
    assert(radial_sector_index != basis::kNone);

    double prefactor = 0.;
    if (run_parameters.transition_type == TransitionType::kE) {
      prefactor = CalculatePrefactorE(run_parameters, bra_subspace, ket_subspace);
    } else {
      prefactor = CalculatePrefactorM(run_parameters, bra_subspace, ket_subspace);
    }

    operator_matrices.push_back(prefactor * radial_matrices[radial_sector_index]);
  }

  // write operator to file
  shell::OutRadialStream os(run_parameters.output_filename, bra_space,
                            ket_space, operator_sectors,
                            shell::RadialOperatorType::kGeneric, radial_order);
  os.Write(operator_matrices);

  /* return code */
  return EXIT_SUCCESS;
}
