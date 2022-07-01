/****************************************************************
  oscillator_wf.h

  Defines recurrence relations for harmonic oscillator wave functions.

  Language: C++17

  Patrick J. Fasano
  University of Notre Dame

  + 03/03/22 (pjf): Created.
  + 03/20/22 (pjf): Add OscillatorWaveFunction().

****************************************************************/

#ifndef OSCILLATOR_WF_H_
#define OSCILLATOR_WF_H_

#include <cmath>
#include <utility>
#include <vector>

namespace shell
{

template<int coord_sign, typename T>
inline auto UnitlessCoordSqr(double b, const T& coord)
{
  static_assert((coord_sign == +1) || (coord_sign == -1));
  if constexpr (coord_sign == +1)
    return ((coord * coord) / (b * b));
  else  // if (coord_sign == -1)
    return ((coord * coord) * (b * b));
}
template<int coord_sign, typename T>
inline auto OscillatorWaveFunctionRecurrence0(
    int l, double b, const T& x2, double power_shift = 0.0
  )
{
  static_assert((coord_sign == +1) || (coord_sign == -1));
  // n = 0 explicit form
  if constexpr (coord_sign == +1)
  {
    return (std::sqrt(2. / (b * std::tgamma(l + 1.5))) * std::pow(b, power_shift))
           * exp(-x2 / 2) * pow(x2, (l + 1 + power_shift) / 2);
  }
  else  // if (coord_sign == -1)
  {
    return (std::sqrt((2. * b) / std::tgamma(l + 1.5)) * std::pow(b, -power_shift))
           * exp(-x2 / 2) * pow(x2, (l + 1 + power_shift) / 2);
  }
}

template<int coord_sign, typename T>
inline auto OscillatorWaveFunctionRecurrence1(int l, const T& x2, const T& wf_m1)
{
  static_assert((coord_sign == +1) || (coord_sign == -1));
  // n = 1 in terms of n = 0; only first term of recurrence
  return (coord_sign * ((1.5 + l - x2) / std::sqrt(l + 1.5))) * wf_m1;
}

template<int coord_sign, typename T>
inline auto OscillatorWaveFunctionRecurrenceN(
    int n, int l, const T& x2, const T& wf_m1, const T& wf_m2
  )
{
  static_assert((coord_sign == +1) || (coord_sign == -1));
  return coord_sign
         * ((-x2 + (2 * n + l - 0.5)) * wf_m1
            - std::sqrt((n - 1) * (n + l - 0.5)) * wf_m2)
         / std::sqrt(n * (n + l + 0.5));
}

template<int coord_sign, typename T>
auto OscillatorWaveFunctions(
    int nmax, int l, double b, const T& coord, double power_shift = 0.
  )
// Generate harmonic oscillator wave functions up to nmax.
//
// This template function can evaluate the oscillator wave functions on any type
// for which coefficient-wise +, *, pow() and exp() are defined. For example, if
// T is double, this function will evaluate all the oscillator functions up to
// nmax at the value of coord. If T is, e.g. Eigen::ArrayXd, then the oscillator
// functions will be evaluated at all values contained in the 1-d array coord.
//
// See Mathematica notebook 220303-oscillator-recurrence.nb for derivation.
//
// Note: this function *returns* a std::vector instead of having a non-const
// reference parameter for output, as C++17 has made return value optimization
// mandatory, and thus there is no performance benefit to passing a non-const
// reference parameter for output.
//
// Template arguments:
//   coord_sign (int): sign selecting position- (+1) or momentum- (-1) space
//   T (typename): type of input coordinate values
//
// Arguments:
//   nmax (int): maximum radial quantum number
//   l (int): orbital angular momentum
//   b (double): oscillator length parameter
//   coord (T): value(s) of coordinate at which to evaluate wave function
//   power_shift (double, optional): additional power of r or k to include
//
// Returns:
//   (std::vector<T> of length nmax+1): wave functions for n=0 to n=nmax
{
  static_assert((coord_sign == +1) || (coord_sign == -1));
  if (nmax < 0)
    return std::vector<T>{};
  // output storage
  std::vector<T> wf(nmax + 1);

  // convenience variable
  //   x^2 = (r/b)^2 or x^2 = (k b)^2
  T x2 = UnitlessCoordSqr<coord_sign>(b, coord);

  // n = 0 base case
  wf[0] = OscillatorWaveFunctionRecurrence0<coord_sign>(l, b, x2, power_shift);
  if (nmax == 0)
    return wf;

  // n = 1 in terms of n = 0; only first term of recurrence
  wf[1] = OscillatorWaveFunctionRecurrence1<coord_sign>(l, x2, wf[0]);
  if (nmax == 1)
    return wf;

  // recurrence from n = 2 to n = nmax
  for (int n = 2; n <= nmax; ++n)
  {
    wf[n] = OscillatorWaveFunctionRecurrenceN<coord_sign>(
        n, l, x2, wf[n - 1], wf[n - 2]
      );
  }

  return wf;
}

template<int coord_sign, typename T>
auto OscillatorWaveFunction(
    int n, int l, double b, const T& coord, double power_shift = 0.
  )
// Generate harmonic oscillator wave function.
//
// This template function can evaluate the oscillator wave functions on any type
// for which coefficient-wise +, *, pow() and exp() are defined. For example, if
// T is double, this function will evaluate all the oscillator function at the
// value of coord. If T is, e.g. Eigen::ArrayXd, then the oscillator functions
// will be evaluated at all values contained in the 1-d array coord.
//
// See Mathematica notebook 220303-oscillator-recurrence.nb for derivation.
//
// Note: this function *returns* T instead of having a non-const reference
// parameter for output, as C++17 has made return value optimization mandatory,
// and thus there is no performance benefit to passing a non-const reference
// parameter for output.
//
// Template arguments:
//   coord_sign (int): sign selecting position- (+1) or momentum- (-1) space
//   T (typename): type of input coordinate values
//
// Arguments:
//   n (int): radial quantum number
//   l (int): orbital angular momentum
//   b (double): oscillator length parameter
//   coord (T): value(s) of coordinate at which to evaluate wave function
//   power_shift (double, optional): additional power of r or k to include
//
// Returns:
//   (std::vector<T> of length nmax+1): wave functions for n=0 to n=nmax
{
  static_assert((coord_sign == +1) || (coord_sign == -1));
  if (n < 0)
    return T{};
  // output storage
  T wf, wf_m1, wf_m2;

  // convenience variable
  //   x^2 = (r/b)^2 or x^2 = (k b)^2
  T x2 = UnitlessCoordSqr<coord_sign>(b, coord);

  // n = 0 base case
  wf = OscillatorWaveFunctionRecurrence0<coord_sign>(l, b, x2, power_shift);
  if (n == 0)
    return wf;

  // n = 1 in terms of n = 0; only first term of recurrence
  std::swap(wf, wf_m1);
  wf = OscillatorWaveFunctionRecurrence1<coord_sign>(l, x2, wf_m1);
  if (n == 1)
    return wf;

  // recurrence from n = 2 to n = nmax
  for (int i = 2; i <= n; ++i)
  {
    std::swap(wf_m2, wf_m1);
    std::swap(wf, wf_m1);
    wf = OscillatorWaveFunctionRecurrenceN<coord_sign>(i, l, x2, wf_m1, wf_m2);
  }

  return wf;
}
}  // namespace shell

#endif  // OSCILLATOR_WF_H_
