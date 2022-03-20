/****************************************************************

  relcm_lenpic_me.cpp

  Patrick J. Fasano
  University of Notre Dame

****************************************************************/

#include "relcm_lenpic_me.h"

#include <Eigen/Dense>
#include <stdexcept>
#include <vector>
#if __has_include(<numbers>)
#  include <numbers>
#endif

#include "am/racah_reduction.h"
#include "am/rme.h"
#include "am/wigner_gsl.h"
#include "analytic/radial_oscillator_me.h"
#include "fmt/format.h"
#include "quadrature/mesh.h"
#include "quadrature/oscillator_wf.h"
#include "spline/spline.h"

namespace relative::lenpic
{

////////////////////////////////////////////////////////////////
// physical constants
//
// constants modified from chime/programs/constants.h (S. Pal, MIT license)
//
// Note: some alternative values of constants are included as comments below
//
// updated with data from:
// [1] E. Tiesinga, P. J. Mohr, D. B. Newell, and B. N. Taylor (2020), "The 2018
//     CODATA Recommended Values of the Fundamental Physical Constants" (Web
//     Version 8.1). Database developed by J. Baker, M. Douma, and S.
//     Kotochigova. Available at http://physics.nist.gov/constants, National
//     Institute of Standards and Technology, Gaithersburg, MD
//     20899.
// [2] P.A. Zyla et al. (Particle Data Group), Prog. Theor. Exp. Phys. 2020,
//     083C01 (2020).
// [3] Epelbaum, Krebs, Meissner, Eur. Phys. J. A 51, 53 (2015).
////////////////////////////////////////////////////////////////
namespace constants
{

// mathematical constants
#ifdef __cpp_lib_math_constants
using std::numbers::pi;
#else
inline constexpr double pi = 3.14159265358979323846;
#endif

// Constants copied from S. Pal:
// constexpr double hbarc = 197.326'960'2;  // (in MeV fm)
// constexpr double charged_pion_mass_MeV = 139.570'61;
// constexpr double neutral_pion_mass_MeV = 134.977'0;
// constexpr double pion_decay_constant_MeV = 92.4;
// constexpr double pion_mass_MeV =
//     (2 * charged_pion_mass_MeV + neutral_pion_mass_MeV) / 3.;

// fundamental constants
constexpr double hbarc = 197.326'980'4;     // (in MeV fm) [1]
constexpr double hbarc_GeV = 1000 * hbarc;  // (in GeV fm)

// masses in MeV
// constexpr double charged_pion_mass_MeV = 139.57039;  // [2]
// constexpr double neutral_pion_mass_MeV = 134.9768;   // [2]
constexpr double charged_pion_mass_MeV = 139.57;  // [3]
constexpr double neutral_pion_mass_MeV = 134.98;  // [3]
constexpr double pion_mass_MeV = 138.03;          // [3]
// constexpr double proton_mass_MeV = 938.272'088'16;  // [1]
constexpr double proton_mass_MeV = 938.272;  // [3]
// constexpr double neutron_mass_MeV = 939.565'420'52;  // [1]
constexpr double neutron_mass_MeV = 939.565;  // [3]

// masses fm^{-1}
constexpr double pion_mass_fm = pion_mass_MeV / hbarc;
constexpr double proton_mass_fm = proton_mass_MeV / hbarc;
constexpr double neutron_mass_fm = neutron_mass_MeV / hbarc;
constexpr double nucleon_mass_MeV = ((proton_mass_MeV + neutron_mass_MeV) / 2);
constexpr double nucleon_mass_fm = nucleon_mass_MeV / hbarc;

// electroweak properties
constexpr double nuclear_magneton_MeV =
    1.0 / (proton_mass_MeV * 2);  // (e = 1), note m_p in denominator [1]
constexpr double nuclear_magneton_fm =
    1.0 / (proton_mass_fm * 2);  // (e = 1), note m_p in denominator [1]
// constexpr double pion_decay_constant_MeV = 92.3;  // [2]
constexpr double pion_decay_constant_MeV = 92.4;  // [3]
constexpr double pion_decay_constant_fm = pion_decay_constant_MeV / hbarc;
// constexpr double gA = 1.2754;  // [2]
constexpr double gA = 1.29;  // [3], see note about Goldberger-Treiman discrepancy
}  // namespace constants

////////////////////////////////////////////////////////////////
// N2LO M1 operator
////////////////////////////////////////////////////////////////

void ConstructNLOM1Operator(
    const basis::OperatorLabelsJT& operator_labels,
    const basis::RelativeCMSpaceLSJT& relative_cm_space,
    std::array<basis::RelativeCMSectorsLSJT, 3>& relative_cm_component_sectors,
    std::array<basis::OperatorBlocks<double>, 3>& relative_cm_component_matrices,
    double regulator_R,
    double single_particle_b,
    std::size_t num_points
  )
{
  // diagnostic output
  fmt::print("  Generating LENPIC SCS NLO M1 two-body operator");

  // (iso)spin coefficients for [sigma1 x sigma2]_1
  //
  // evaluated using an immediately-invoked lambda function so kSpinFactors can
  // be const
  static const auto kSpinFactors{
    []() {
      std::array<std::array<std::array<double, 3>, 2>, 2> temp{};
      for (int bra_S : {0, 1})
        for (int ket_S : {0, 1})
          for (int S0 : {0, 1, 2})
            if (am::AllowedTriangle(bra_S, ket_S, S0))
              temp[bra_S][ket_S][S0] = am::RacahReductionFactor12Rose(
                  0.5_hi, 0.5_hi, bra_S,
                  0.5_hi, 0.5_hi, ket_S,
                  1, 1, S0
                )
                * 2 * am::AngularMomentumJRME(0.5_hi, 0.5_hi)
                * 2 * am::AngularMomentumJRME(0.5_hi, 0.5_hi);
      return temp;
    }()
  };

  // check quantum numbers
  if (
      (operator_labels.J0 != 1) || (operator_labels.g0 != 0)
      || (operator_labels.T0_min > 1) || (operator_labels.T0_max < 1)
    )
  {
    throw std::invalid_argument(fmt::format(
        "invalid operator labels: J0={} g0={} T0_min={} T0_max={}",
        operator_labels.J0,
        operator_labels.g0,
        operator_labels.T0_min,
        operator_labels.T0_max
      ));
  }

  // zero initialize operator
  basis::ConstructZeroOperatorRelativeCMLSJT(
      operator_labels,
      relative_cm_space,
      relative_cm_component_sectors,
      relative_cm_component_matrices
    );

  // only T0=1 component exists
  const int T0 = 1;
  const auto& sectors = relative_cm_component_sectors[T0];
  auto& matrices = relative_cm_component_matrices[T0];

  // common setup
  const double b_rel = single_particle_b * std::sqrt(2);
  const double b_cm = single_particle_b / std::sqrt(2);
  const Eigen::ArrayXd z_grid = shell::UniformMesh(num_points);
  const Eigen::ArrayXd r_grid = shell::TransformedCoordinate<+1>(b_rel, z_grid);
  const Eigen::ArrayXd jacobian_grid =
      shell::TransformedCoordinateJacobian<+1>(b_rel, z_grid);
  const Eigen::ArrayXd exp_mpi_r_grid = exp(-constants::pion_mass_fm * r_grid);
  const Eigen::ArrayXd regulator_grid =
      (regulator_R == 0)
          ? Eigen::ArrayXd::Ones(num_points).eval()
          : (pow(1. - exp(-r_grid * r_grid / (regulator_R * regulator_R)), 6)).eval();

  // here we precompute the oscillator functions, but we do it three times so at
  // to absorb the factors of r^-1 and r^-2 into the oscillator function
  // recurrences; for instance for the integral exp(-m_pi*r)/(m_pi*r), we use
  // absorb r^-0.5 into the each of bra and ket radial wave functions
  std::vector<std::vector<Eigen::ArrayXd>> oscillator_functions_grid,
      oscillator_functions_rm1_grid, oscillator_functions_rm2_grid;
  oscillator_functions_grid.reserve(relative_cm_space.Nmax() + 1);
  oscillator_functions_rm1_grid.reserve(relative_cm_space.Nmax() + 1);
  oscillator_functions_rm2_grid.reserve(relative_cm_space.Nmax() + 1);
  for (int l = 0; l <= relative_cm_space.Nmax(); ++l)
  {
    oscillator_functions_grid.push_back(shell::OscillatorWaveFunctions<+1>(
        (relative_cm_space.Nmax() - l) / 2, l, b_rel, r_grid, /*power_shift=*/0.
      ));
    oscillator_functions_rm1_grid.push_back(shell::OscillatorWaveFunctions<+1>(
        (relative_cm_space.Nmax() - l) / 2, l, b_rel, r_grid, /*power_shift=*/-0.5
      ));
    oscillator_functions_rm2_grid.push_back(shell::OscillatorWaveFunctions<+1>(
        (relative_cm_space.Nmax() - l) / 2, l, b_rel, r_grid, /*power_shift=*/-1.
      ));
  }

// iterate over sectors
#pragma omp parallel for default(none) schedule(dynamic) shared( \
    std::cout,                                                   \
    b_cm,                                                        \
    sectors,                                                     \
    matrices,                                                    \
    kSpinFactors,                                                \
    z_grid,                                                      \
    jacobian_grid,                                               \
    regulator_grid,                                              \
    oscillator_functions_grid,                                   \
    oscillator_functions_rm1_grid,                               \
    oscillator_functions_rm2_grid,                               \
    exp_mpi_r_grid                                               \
  )
  for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
  {
    // diagnostic output
    #pragma omp critical
    if (sector_index % 10 == 0)
      std::cout << "." << std::flush;

    // extract sector, subspace, and labels
    const auto& sector = sectors.GetSector(sector_index);
    auto& matrix = matrices.at(sector_index);
    const auto& bra_subspace = sector.bra_subspace();
    const auto& ket_subspace = sector.ket_subspace();
    const auto& [bra_L, bra_S, bra_J, bra_T, bra_g] = bra_subspace.labels();
    const auto& [ket_L, ket_S, ket_J, ket_T, ket_g] = ket_subspace.labels();

    // only |delta_T| = 1 sectors contribute
    if (bra_T == ket_T)
      continue;

    // extract spin and isospin factors/RMEs
    const auto spin_rmes = kSpinFactors[bra_S][ket_S];
    const double isospin_rme = kSpinFactors[bra_T][ket_T][1];

    for (std::size_t bra_index = 0; bra_index < bra_subspace.size(); ++bra_index)
    {
      for (std::size_t ket_index = 0; ket_index < ket_subspace.size(); ++ket_index)
      {
        const auto bra_state = bra_subspace.GetState(bra_index);
        const auto ket_state = ket_subspace.GetState(ket_index);
        const auto& [bra_Nr, bra_lr, bra_Nc, bra_lc] = bra_state.labels();
        const auto& [ket_Nr, ket_lr, ket_Nc, ket_lc] = ket_state.labels();
        const int bra_nr = (bra_Nr - bra_lr) / 2;
        const int ket_nr = (ket_Nr - ket_lr) / 2;

        bool cm_changing =
            (std::abs(bra_Nc - ket_Nc) == 1) && (std::abs(bra_lc - ket_lc) == 1);
        bool cm_preserving = ((bra_Nc == ket_Nc) && (bra_lc == ket_lc));
        if (!cm_changing && !cm_preserving)
          continue;

        // short circuit on orbital angular momentum selection rules
        if (!((bra_lc == ket_lc) || am::AllowedTriangle(bra_lc, ket_lc, 1)))
          continue;
        if (!((bra_lr == ket_lr) || am::AllowedTriangle(bra_lr, ket_lr, 1)
              || am::AllowedTriangle(bra_lr, ket_lr, 2)
              || am::AllowedTriangle(bra_lr, ket_lr, 3)))
          continue;


        // labeled by {rel/cm,l0}
        const std::array<std::array<double, 4>, 2> spherical_harmonic_rmes{
            {{(bra_lr == ket_lr)
                  ? 1. /*am::SphericalHarmonicCRME(bra_lr, ket_lr, 0)*/
                  : 0,
              am::AllowedTriangle(bra_lr, ket_lr, 1)
                  ? am::SphericalHarmonicCRME(bra_lr, ket_lr, 1)
                  : 0,
              am::AllowedTriangle(bra_lr, ket_lr, 2)
                  ? am::SphericalHarmonicCRME(bra_lr, ket_lr, 2)
                  : 0,
              am::AllowedTriangle(bra_lr, ket_lr, 3)
                  ? am::SphericalHarmonicCRME(bra_lr, ket_lr, 3)
                  : 0},
             {(bra_lc == ket_lc)
                  ? 1. /*am::SphericalHarmonicCRME(bra_lc, ket_lc, 0)*/
                  : 0,
              am::AllowedTriangle(bra_lc, ket_lc, 1)
                  ? am::SphericalHarmonicCRME(bra_lc, ket_lc, 1)
                  : 0,
              0. /*am::SphericalHarmonicCRME(bra_lc, ket_lc, 2)*/,
              0. /*am::SphericalHarmonicCRME(bra_lc, ket_lc, 3)*/}}
          };

        // initialize angular factors
        std::array<double, 7> tensor_product_rmes{};
        // aggregate initialization -> all zero

        ////////////////////////////////////////////////////////////////
        // angular momentum factors
        //
        // tensor_product_rmes[0] is [sigma1 x sigma2]_1 in the LSJ basis
        // tensor_product_rmes[1] through tensor_product_rmes[6] are the
        // angular factors A_1 through A_6 as defined in the DNP19 slides
        // from S. Sarker
        ////////////////////////////////////////////////////////////////

        // this is just S_1, coupled with the identity on orbital angular momenta
        if (cm_preserving)
        {
          if (am::AllowedTriangle(bra_S, ket_S, 1) && (bra_lr == ket_lr) && (bra_L == ket_L))
            tensor_product_rmes[0] =
                am::RacahReductionFactor2Rose(
                    bra_L, bra_S, bra_J, ket_L, ket_S, ket_J, 1
                  )
                * spin_rmes[1];
        }
        // A_1 through A_5 are only needed for the center-of-mass changing part
        // of the operator
        if (cm_changing)
        {
          if (
              am::AllowedTriangle(bra_S, ket_S, 0)
              && am::AllowedTriangle(bra_lr, ket_lr, 1)
              && am::AllowedTriangle(bra_lc, ket_lc, 1)
              && am::AllowedTriangle(bra_L, ket_L, 1)
            )
            // "A_1"
            tensor_product_rmes[1] =
                -std::sqrt(3.)
                * am::RacahReductionFactor12Rose(
                    bra_lr, bra_lc, bra_L, ket_lr, ket_lc, ket_L, 1, 1, 1
                  )
                * spherical_harmonic_rmes[0][1] * spherical_harmonic_rmes[1][1]
                * am::RacahReductionFactor12Rose(
                    bra_L, bra_S, bra_J, ket_L, ket_S, ket_J, 1, 0, 1
                  )
                * spin_rmes[0];

          if (
              am::AllowedTriangle(bra_S, ket_S, 2)
              && am::AllowedTriangle(bra_lr, ket_lr, 1)
              && am::AllowedTriangle(bra_lc, ket_lc, 1)
              && am::AllowedTriangle(bra_L, ket_L, 1)
            )
            // "A_2"
            tensor_product_rmes[2] =
                std::sqrt(3. / 5.)
                * am::RacahReductionFactor12Rose(
                    bra_lr, bra_lc, bra_L, ket_lr, ket_lc, ket_L, 1, 1, 1
                  )
                * spherical_harmonic_rmes[0][1] * spherical_harmonic_rmes[1][1]
                * am::RacahReductionFactor12Rose(
                    bra_L, bra_S, bra_J, ket_L, ket_S, ket_J, 1, 2, 1
                  )
                * spin_rmes[2];

          if (
              am::AllowedTriangle(bra_S, ket_S, 2)
              && am::AllowedTriangle(bra_lr, ket_lr, 1)
              && am::AllowedTriangle(bra_lc, ket_lc, 1)
              && am::AllowedTriangle(bra_L, ket_L, 2)
            )
            // "A_3"
            tensor_product_rmes[3] =
                std::sqrt(9. / 5.)
                * am::RacahReductionFactor12Rose(
                    bra_lr, bra_lc, bra_L, ket_lr, ket_lc, ket_L, 1, 1, 2
                  )
                * spherical_harmonic_rmes[0][1] * spherical_harmonic_rmes[1][1]
                * am::RacahReductionFactor12Rose(
                    bra_L, bra_S, bra_J, ket_L, ket_S, ket_J, 2, 2, 1
                  )
                * spin_rmes[2];

          if (
              am::AllowedTriangle(bra_S, ket_S, 2)
              && am::AllowedTriangle(bra_lr, ket_lr, 3)
              && am::AllowedTriangle(bra_lc, ket_lc, 1)
              && am::AllowedTriangle(bra_L, ket_L, 2)
            )
            // "A_4"
            tensor_product_rmes[4] =
                std::sqrt(14. / 5.)
                * am::RacahReductionFactor12Rose(
                    bra_lr, bra_lc, bra_L, ket_lr, ket_lc, ket_L, 3, 1, 2
                  )
                * spherical_harmonic_rmes[0][3] * spherical_harmonic_rmes[1][1]
                * am::RacahReductionFactor12Rose(
                    bra_L, bra_S, bra_J, ket_L, ket_S, ket_J, 2, 2, 1
                  )
                * spin_rmes[2];

          if (
              am::AllowedTriangle(bra_S, ket_S, 2)
              && am::AllowedTriangle(bra_lr, ket_lr, 3)
              && am::AllowedTriangle(bra_lc, ket_lc, 1)
              && am::AllowedTriangle(bra_L, ket_L, 3)
            )
            // "A_5"
            tensor_product_rmes[5] =
                std::sqrt(28. / 5.)
                * am::RacahReductionFactor12Rose(
                    bra_lr, bra_lc, bra_L, ket_lr, ket_lc, ket_L, 3, 1, 3
                  )
                * spherical_harmonic_rmes[0][3] * spherical_harmonic_rmes[1][1]
                * am::RacahReductionFactor12Rose(
                    bra_L, bra_S, bra_J, ket_L, ket_S, ket_J, 3, 2, 1
                  )
                * spin_rmes[2];
        }
        if (cm_preserving)
        {
          // "A_6(S_1)"
          if (
              am::AllowedTriangle(bra_S, ket_S, 1) && (bra_lc == ket_lc)
              && am::AllowedTriangle(bra_lr, ket_lr, 2)
              && am::AllowedTriangle(bra_L, ket_L, 2)
            )
            tensor_product_rmes[6] =
                std::sqrt(10.)
                * am::RacahReductionFactor1Rose(
                    bra_lr, bra_lc, bra_L, ket_lr, ket_lc, ket_L, 2
                  )
                * spherical_harmonic_rmes[0][2] * spherical_harmonic_rmes[1][0]
                * am::RacahReductionFactor12Rose(
                    bra_L, bra_S, bra_J, ket_L, ket_S, ket_J, 2, 1, 1
                  )
                * spin_rmes[1];
        }

        ////////////////////////////////////////////////////////////////
        // integrals
        //
        // f[n] is exp(-m_pi r)/r^n, where the r^n factor is absorbed into the
        // oscillator functions
        ////////////////////////////////////////////////////////////////
        double f0_integral = 0, f1_integral = 0, f2_integral = 0;
        Eigen::ArrayXd f0_integrand =
            oscillator_functions_grid[bra_lr][bra_nr] * exp_mpi_r_grid
            * oscillator_functions_grid[ket_lr][ket_nr] * regulator_grid
            * jacobian_grid;
        shell::FixEndpointSingularitiesToZero(f0_integrand);
        f0_integral = spline::CubicIntegrate(z_grid, f0_integrand);
        if (std::isnan(f0_integral))
          throw std::runtime_error("f0_integral is nan");

        Eigen::ArrayXd f1_integrand =
            oscillator_functions_rm1_grid[bra_lr][bra_nr] * exp_mpi_r_grid
            * oscillator_functions_rm1_grid[ket_lr][ket_nr] * regulator_grid
            * jacobian_grid;
        shell::FixEndpointSingularitiesToZero(f1_integrand);
        f1_integral = spline::CubicIntegrate(z_grid, f1_integrand);
        if (std::isnan(f1_integral))
          throw std::runtime_error("f1_integral is nan");
        if (cm_changing)
        {
          Eigen::ArrayXd f2_integrand =
              oscillator_functions_rm2_grid[bra_lr][bra_nr] * exp_mpi_r_grid
              * oscillator_functions_rm2_grid[ket_lr][ket_nr] * regulator_grid
              * jacobian_grid;
          shell::FixEndpointSingularitiesToZero(f2_integrand);
          f2_integral = spline::CubicIntegrate(z_grid, f2_integrand);
          if (std::isnan(f2_integral))
            throw std::runtime_error("f2_integral is nan");
        }
        ////////////////////////////////////////////////////////////////

        // pure relative part
        if (cm_preserving)
        {
          matrix(bra_index, ket_index) +=
              f0_integral * (tensor_product_rmes[6] + 2 * tensor_product_rmes[0])
              + (f1_integral / constants::pion_mass_fm)
                    * (tensor_product_rmes[6] - tensor_product_rmes[0]);
        }
        // center-of-mass changing part
        if (cm_changing)
        {
          double cm_radial_me =
              b_cm * constants::pion_mass_fm
              * analytic::CoordinateOscillatorMatrixElement(
                  bra_Nc, bra_lc, ket_Nc, ket_lc, +1
                );

          matrix(bra_index, ket_index) +=
              (f0_integral * cm_radial_me
               * (tensor_product_rmes[1] + tensor_product_rmes[2]
                  + tensor_product_rmes[3] + tensor_product_rmes[4]
                  + tensor_product_rmes[5]))
              + (3 * (f1_integral / constants::pion_mass_fm) * cm_radial_me
                 * (tensor_product_rmes[2] + tensor_product_rmes[3]
                    + tensor_product_rmes[4] + tensor_product_rmes[5]))
              + (3
                 * (f2_integral
                    / (constants::pion_mass_fm * constants::pion_mass_fm))
                 * cm_radial_me
                 * (tensor_product_rmes[2] + tensor_product_rmes[3]
                    + tensor_product_rmes[4] + tensor_product_rmes[5]));
        }
      }
    }
    matrix *= -constants::gA * constants::gA * constants::nucleon_mass_fm
              * constants::pion_mass_fm
              / (24 * constants::pi * constants::pion_decay_constant_fm
                 * constants::pion_decay_constant_fm)
              * isospin_rme;
  }
  std::cout << "done." << std::endl;
}


};  // namespace relative::lenpic
