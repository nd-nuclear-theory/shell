/****************************************************************

  lenpic_relative_me.cpp

  Patrick J. Fasano
  University of Notre Dame

****************************************************************/

#include "lenpic_relative_me.h"

#include <Eigen/Dense>
#include <stdexcept>
#include <vector>

#include "am/racah_reduction.h"
#include "am/rme.h"
#include "am/wigner_gsl.h"
#include "analytic/radial_oscillator_me.h"
#include "fmt/format.h"
#include "lenpic_constants.h"
#include "quadrature/mesh.h"
#include "quadrature/oscillator_wf.h"
#include "relative/relative_me.h"
#include "spline/spline.h"

namespace relative::lenpic
{

////////////////////////////////////////////////////////////////
// LO Gamow-Teller operator (as two-body operator)
////////////////////////////////////////////////////////////////

void ConstructLOGTOperator(
    const basis::OperatorLabelsJT& operator_labels,
    const basis::RelativeSpaceLSJT& relative_space,
    std::array<basis::RelativeSectorsLSJT, 3>& relative_component_sectors,
    std::array<basis::OperatorBlocks<double>, 3>& relative_component_matrices
  )
{
  ConstructSpinAMOperator(
      operator_labels,
      relative_space,
      relative_component_sectors,
      relative_component_matrices,
      /*T0=*/1
    );

  basis::ScalarMultiplyOperator(
      relative_component_sectors[1], relative_component_matrices[1], -constants::gA
    );
}

////////////////////////////////////////////////////////////////
// N2LO Gamow-Teller operator
////////////////////////////////////////////////////////////////

void ConstructN2LOGTOperator(
    const basis::OperatorLabelsJT& operator_labels,
    const basis::RelativeSpaceLSJT& relative_space,
    std::array<basis::RelativeSectorsLSJT, 3>& relative_component_sectors,
    std::array<basis::OperatorBlocks<double>, 3>& relative_component_matrices,
    double regulator_R,
    double single_particle_b,
    std::size_t num_points
  )
{
  // diagnostic output
  fmt::print("  Generating LENPIC SCS N2LO GT two-body operator");

  // overall constants
  //
  // opep_prefactor = g_A M_pi^3 / 12 pi F_pi^2
  // contact_prefactor = (-1/4) * D
  const double opep_prefactor =
      constants::gA * std::pow(constants::pion_mass_fm, 3)
      / (12. * constants::pi * std::pow(constants::pion_decay_constant_fm, 2));
  const double contact_prefactor = -constants::SCS_D_fm3(regulator_R) / 4.;

  // (iso)spin coefficients for [sigma_1 x sigma_2]_S0
  //
  // evaluated using an immediately-invoked lambda function so kSpinFactors can
  // be const
  static const auto kSpinFactors{
    []() {
      static constexpr int S0 = 1;
      std::array<std::array<double, 2>, 2> temp{};
      for (int bra_S : {0, 1})
        for (int ket_S : {0, 1})
            if (am::AllowedTriangle(bra_S, ket_S, S0))
              temp[bra_S][ket_S] = am::RacahReductionFactor12Rose(
                  0.5_hi, 0.5_hi, bra_S,
                  0.5_hi, 0.5_hi, ket_S,
                  1, 1, S0
                )
                * 2 * am::AngularMomentumJRME(0.5_hi, 0.5_hi)
                * 2 * am::AngularMomentumJRME(0.5_hi, 0.5_hi);
      return temp;
    }()
  };

  // coefficients for (tau_1 sigma_1 + tau_2 sigma_2)
  //
  // labeled by bra_S, ket_S, bra_T, ket_T
  //
  // evaluated using an immediately-invoked lambda function so kSpinFactors can
  // be const
  static const auto kSU4Factors{
    []() {
      std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2> temp{};
      for (int bra_S : {0, 1})
        for (int ket_S : {0, 1})
          for (int bra_T : {0, 1})
            for (int ket_T : {0, 1})
              if (am::AllowedTriangle(bra_S, ket_S, 1) && am::AllowedTriangle(bra_T, ket_T, 1))
                temp[bra_S][ket_S][bra_T][ket_T] =
                  2*am::jjJCoupledAngularMomentumJ1RME(0.5_hi, 0.5_hi, bra_S, 0.5_hi, 0.5_hi, ket_S)
                  * 2*am::jjJCoupledAngularMomentumJ1RME(0.5_hi, 0.5_hi, bra_T, 0.5_hi, 0.5_hi, ket_T)
                  + 2*am::jjJCoupledAngularMomentumJ2RME(0.5_hi, 0.5_hi, bra_S, 0.5_hi, 0.5_hi, ket_S)
                    * 2*am::jjJCoupledAngularMomentumJ2RME(0.5_hi, 0.5_hi, bra_T, 0.5_hi, 0.5_hi, ket_T);
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
  basis::ConstructZeroOperatorRelativeLSJT(
      operator_labels, relative_space, relative_component_sectors, relative_component_matrices
    );

  // only T0=1 component exists
  const int T0 = 1;
  const auto& sectors = relative_component_sectors[T0];
  auto& matrices = relative_component_matrices[T0];

  // common setup
  const double b_rel = single_particle_b * std::sqrt(2);
  // z \in [0,1]
  const Eigen::ArrayXd z_grid = shell::UniformMesh(num_points);
  // r \in [0,infty]
  const Eigen::ArrayXd r_grid = shell::TransformedCoordinate<+1>(b_rel, z_grid);
  const Eigen::ArrayXd jacobian_grid =
      shell::TransformedCoordinateJacobian<+1>(b_rel, z_grid);
  // exp(- M_pi * r)
  const Eigen::ArrayXd exp_mpi_r_grid = exp(-constants::pion_mass_fm * r_grid);
  // W(r) = 1+3*(1+M_pi*r)/(M_pi r)^2
  const Eigen::ArrayXd w_grid =
      (1.
       + 3. * (1. + constants::pion_mass_fm * r_grid)
             / (constants::pion_mass_fm * constants::pion_mass_fm * r_grid
                * r_grid));
  // [1-exp(-r^2/R^2)]^6
  const Eigen::ArrayXd local_regulator_grid =
      (regulator_R == 0)
          ? Eigen::ArrayXd::Ones(num_points).eval()
          : (pow(1. - exp(-r_grid * r_grid / (regulator_R * regulator_R)), 6)).eval();
  // exp(-r^2/R^2)/(pi R^2)^3/2
  const Eigen::ArrayXd nonlocal_regulator_grid =
      (regulator_R == 0)
          ? [&](){ auto tmp = Eigen::ArrayXd::Zero(num_points).eval(); tmp(0) = 1; return tmp; }()
          : (exp(-r_grid * r_grid / (regulator_R * regulator_R))
             / pow(constants::pi * regulator_R * regulator_R, 1.5))
                .eval();

  std::vector<std::vector<Eigen::ArrayXd>> oscillator_functions_grid,
      oscillator_functions_rm1_grid;
  // psi_{l,n}(r)
  oscillator_functions_grid.reserve(relative_space.Nmax() + 1);
  // psi_{l,n}(r)/r^{1/2} -- factor out 1/r symmetrically between bra and ket w.f.
  oscillator_functions_rm1_grid.reserve(relative_space.Nmax() + 1);
  for (int l = 0; l <= relative_space.Nmax(); ++l)
  {
    oscillator_functions_grid.push_back(shell::OscillatorWaveFunctions<+1>(
        (relative_space.Nmax() - l) / 2, l, b_rel, r_grid, /*power_shift=*/0.
      ));
    oscillator_functions_rm1_grid.push_back(shell::OscillatorWaveFunctions<+1>(
        (relative_space.Nmax() - l) / 2, l, b_rel, r_grid, /*power_shift=*/-0.5
      ));
  }

// iterate over sectors
#pragma omp parallel for default(none) schedule(dynamic) shared( \
    std::cout,                                                   \
    constants::pion_mass_fm,                                     \
    sectors,                                                     \
    matrices,                                                    \
    kSpinFactors,                                                \
    kSU4Factors,                                                 \
    regulator_R,                                                 \
    opep_prefactor,                                              \
    contact_prefactor,                                           \
    z_grid,                                                      \
    jacobian_grid,                                               \
    local_regulator_grid,                                        \
    nonlocal_regulator_grid,                                     \
    oscillator_functions_grid,                                   \
    oscillator_functions_rm1_grid,                               \
    exp_mpi_r_grid,                                              \
    w_grid                                                       \
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

    // short circuit on orbital angular momentum selection rules
    if (!((bra_L == ket_L) || am::AllowedTriangle(bra_L, ket_L, 2)))
      continue;
    if (!(am::AllowedTriangle(bra_S, ket_S, 1)
          && am::AllowedTriangle(bra_T, ket_T, 1)))
      continue;

    const double C2_rme = am::AllowedTriangle(bra_L, ket_L, 2)
                              ? am::SphericalHarmonicCRME(bra_L, ket_L, 2)
                              : 0;


    ////////////////////////////////////////////////////////////////
    // angular momentum/isospin factors in the LSJT basis
    //
    // * su4_generator_rme is (sigma_1 tau_1 + sigma_2 tau_2)
    // * spinisospin_rme is [tau_1 x tau_2]_1 * [sigma_1 x sigma_2]_1
    // * C2_su4_generator_rme is [C_2 x (sigma_1 tau_1 + sigma_2 tau_2)]_1
    // * C2_spinisospin_rme is [tau_1 x tau_2]_1 [C_2 x [sigma_1 x sigma_2]_1]_1
    ////////////////////////////////////////////////////////////////
    const double su4_generator_rme =
        am::RacahReductionFactor2Rose(bra_L, bra_S, bra_J, ket_L, ket_S, ket_J, 1)
        * kSU4Factors[bra_S][ket_S][bra_T][ket_T];
    const double spinisospin_rme =
        am::RacahReductionFactor2Rose(bra_L, bra_S, bra_J, ket_L, ket_S, ket_J, 1)
        * kSpinFactors[bra_S][ket_S] * kSpinFactors[bra_T][ket_T];
    const double C2_su4_generator_rme =
        am::RacahReductionFactor12Rose(
            bra_L, bra_S, bra_J, ket_L, ket_S, ket_J, 2, 1, 1
          )
        * C2_rme * kSU4Factors[bra_S][ket_S][bra_T][ket_T];
    const double C2_spinisospin_rme =
        am::RacahReductionFactor12Rose(
            bra_L, bra_S, bra_J, ket_L, ket_S, ket_J, 2, 1, 1
          )
        * C2_rme * kSpinFactors[bra_S][ket_S] * kSpinFactors[bra_T][ket_T];

    for (std::size_t bra_index = 0; bra_index < bra_subspace.size(); ++bra_index)
    {
      for (std::size_t ket_index = 0; ket_index < ket_subspace.size(); ++ket_index)
      {
        const auto bra_state = bra_subspace.GetState(bra_index);
        const auto ket_state = ket_subspace.GetState(ket_index);
        const int bra_Nr = bra_state.N();
        const int ket_Nr = ket_state.N();
        const int bra_nr = bra_state.n();
        const int ket_nr = ket_state.n();

        ////////////////////////////////////////////////////////////////
        // integrals
        //
        // * yukawa_integral is integral of e^(M_pi r)/(M_pi r)
        // * yukawa_w_integral is integral of e^(M_pi r)/(M_pi r) *
        // (1+3*(1+M_pi*r)/(M_pi r)^2)
        // * nonlocal_regulator_integral is integral
        // of e^(-(r'^2+r^2)/R^2)
        //
        // note: additional factor of 1/M_pi is here in yukawa_integral and
        // yukawa_w_integral because only 1/r is factored out into oscillator
        // functions
        ////////////////////////////////////////////////////////////////
        double yukawa_integral = 0, yukawa_w_integral = 0,
               nonlocal_regulator_integral = 0;
        Eigen::ArrayXd yukawa_integrand =
            oscillator_functions_rm1_grid[bra_L][bra_nr] * exp_mpi_r_grid
            * oscillator_functions_rm1_grid[ket_L][ket_nr]
            * local_regulator_grid * jacobian_grid;
        shell::FixEndpointSingularitiesToZero(yukawa_integrand);
        yukawa_integral = spline::CubicIntegrate(z_grid, yukawa_integrand)
                          / constants::pion_mass_fm;
        if (std::isnan(yukawa_integral))
          throw std::runtime_error("yukawa_integral is nan");

        Eigen::ArrayXd yukawa_w_integrand =
            oscillator_functions_rm1_grid[bra_L][bra_nr] * exp_mpi_r_grid
            * w_grid * oscillator_functions_rm1_grid[ket_L][ket_nr]
            * local_regulator_grid * jacobian_grid;
        shell::FixEndpointSingularitiesToZero(yukawa_w_integrand);
        yukawa_w_integral = spline::CubicIntegrate(z_grid, yukawa_w_integrand)
                            / constants::pion_mass_fm;
        if (std::isnan(yukawa_w_integral))
          throw std::runtime_error("yukawa_w_integral is nan");

        if (bra_L == ket_L)
        {
          Eigen::ArrayXd nonlocal_regulator_integrand_bra =
              oscillator_functions_grid[bra_L][bra_nr] * nonlocal_regulator_grid
              * jacobian_grid;
          shell::FixEndpointSingularitiesToZero(nonlocal_regulator_integrand_bra);
          Eigen::ArrayXd nonlocal_regulator_integrand_ket =
              oscillator_functions_grid[ket_L][ket_nr] * nonlocal_regulator_grid
              * jacobian_grid;
          shell::FixEndpointSingularitiesToZero(nonlocal_regulator_integrand_ket);
          nonlocal_regulator_integral =
              spline::CubicIntegrate(z_grid, nonlocal_regulator_integrand_bra)
              * spline::CubicIntegrate(z_grid, nonlocal_regulator_integrand_ket);
          if (std::isnan(nonlocal_regulator_integral))
            throw std::runtime_error("nonlocal_regulator_integral is nan");
        }
        ////////////////////////////////////////////////////////////////
        matrix(bra_index, ket_index) =
            opep_prefactor * yukawa_integral
            * (-constants::SCS_c3_fm(regulator_R) * su4_generator_rme
               + 2 * constants::SCS_c4_fm(regulator_R) * spinisospin_rme);
        matrix(bra_index, ket_index) +=
            opep_prefactor * std::sqrt(10) * yukawa_w_integral
            * (constants::SCS_c3_fm(regulator_R) * C2_su4_generator_rme
               + constants::SCS_c4_fm(regulator_R) * C2_spinisospin_rme);
        matrix(bra_index, ket_index) +=
            contact_prefactor * su4_generator_rme * nonlocal_regulator_integral;
      }
    }
    // matrix *= opep_prefactor;
  }
  std::cout << "done." << std::endl;
}


};  // namespace relative::lenpic
