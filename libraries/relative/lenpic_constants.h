/****************************************************************
  lenpic_relative_me.h

  Physical constants for use with LENPIC operator matrix element generation.

  Constants modified from chime/programs/constants.h (S. Pal, MIT license)

  Note: some alternative values of constants are included as comments below

  updated with data from:
  [1] E. Tiesinga, P. J. Mohr, D. B. Newell, and B. N. Taylor (2020), "The 2018
      CODATA Recommended Values of the Fundamental Physical Constants" (Web
      Version 8.1). Database developed by J. Baker, M. Douma, and S.
      Kotochigova. Available at http://physics.nist.gov/constants, National
      Institute of Standards and Technology, Gaithersburg, MD
      20899.
  [2] P.A. Zyla et al. (Particle Data Group), Prog. Theor. Exp. Phys. 2020,
      083C01 (2020).
  [3] E. Epelbaum, H. Krebs, U.-G. Meissner, Eur. Phys. J. A 51, 53 (2015).
  [4] E. Epelbaum et al., Phys. Rev. C 99, 024313 (2019).

  Language: C++17

  Patrick J. Fasano
  University of Notre Dame

  + 05/25/22 (pjf): Created, split from lenpic_relcm_me.h.
  + 06/02/22 (pjf): Added LECs for LENPIC SCS interactions.

****************************************************************/

#ifndef LENPIC_CONSTANTS_H_
#define LENPIC_CONSTANTS_H_

#include <cmath>
#include <stdexcept>
#if __has_include(<numbers>)
#  include <numbers>
#endif

namespace relative::lenpic::constants
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
constexpr double hbarc_GeV = hbarc/1000;  // (in GeV fm)

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

// LECs for LENPIC SCS

// chiral symmetry breaking scale \Lambda_\chi = 700MeV
constexpr double chiral_breaking_scale_MeV = 700.;  // [4]
constexpr double chiral_breaking_scale_GeV = chiral_breaking_scale_MeV/1000.;
constexpr double chiral_breaking_scale_fm = chiral_breaking_scale_MeV / hbarc;

// L_{\pi N}^{(2)} LECs, in GeV^-1, MeV^-1, or fm
constexpr double SCS_c1_GeV(double R) noexcept { return -0.81; }  // [3,4]
constexpr double SCS_c1_MeV(double R) noexcept { return SCS_c1_GeV(R)/1000.; }
constexpr double SCS_c1_fm(double R) noexcept { return SCS_c1_GeV(R)*hbarc_GeV; }
constexpr double SCS_c2_GeV(double R) noexcept { return 3.28; }  // [3,4]
constexpr double SCS_c2_MeV(double R) noexcept { return SCS_c2_GeV(R)/1000.; }
constexpr double SCS_c2_fm(double R) noexcept { return SCS_c2_GeV(R)*hbarc_GeV; }
constexpr double SCS_c3_GeV(double R) noexcept { return -4.69; }  // [3,4]
constexpr double SCS_c3_MeV(double R) noexcept { return SCS_c3_GeV(R)/1000.; }
constexpr double SCS_c3_fm(double R) noexcept { return SCS_c3_GeV(R)*hbarc_GeV; }
constexpr double SCS_c4_GeV(double R) noexcept { return 3.40; }  // [3,4]
constexpr double SCS_c4_MeV(double R) noexcept { return SCS_c4_GeV(R)/1000.; }
constexpr double SCS_c4_fm(double R) noexcept { return SCS_c4_GeV(R)*hbarc_GeV; }

// L_{\pi NN}^{(1)} LEC, unitless for c_D, in units of fm^3 or MeV^-3 for D
constexpr double SCS_cD(double R)  // [4]
{
  if (std::abs(R - 0.9) < 1e-4)
    return 2.1;
  else if (std::abs(R - 1.0) < 1e-4)
    return 7.2;
  throw std::invalid_argument("invalid regulator R");
};
constexpr double SCS_D_fm3(double R)  // D in units of fm^3
{
  return SCS_cD(R)
         / (pion_decay_constant_fm * pion_decay_constant_fm
            * chiral_breaking_scale_fm);
}
constexpr double SCS_D_MeV3(double R)  // D in units of MeV^-3
{
  return SCS_cD(R)
         / (pion_decay_constant_MeV * pion_decay_constant_MeV
            * chiral_breaking_scale_MeV);
}

// L_{NNN}^{(1)} LEC, unitless for c_E, in units of fm^5 or MeV^-5 for E
constexpr double SCS_cE(double R)  // [4]
{
  if (std::abs(R - 0.9) < 1e-4)
    return -0.329;
  else if (std::abs(R - 1.0) < 1e-4)
    return -0.652;
  throw std::invalid_argument("invalid regulator R");
};
constexpr double SCS_E_fm3(double R)  // E in units of fm^5
{
  return SCS_cE(R)
         / (pion_decay_constant_fm * pion_decay_constant_fm
            * pion_decay_constant_fm * pion_decay_constant_fm
            * chiral_breaking_scale_fm);
}
constexpr double SCS_E_MeV3(double R)  // E in units of MeV^-5
{
  return SCS_cE(R)
         / (pion_decay_constant_MeV * pion_decay_constant_MeV
            * pion_decay_constant_MeV * pion_decay_constant_MeV
            * chiral_breaking_scale_MeV);
}

}  // namespace relative::lenpic::constants

#endif  // LENPIC_CONSTANTS_H_
