/****************************************************************
  oscillator_wf.h

  Defines recurrence relations for harmonic oscillator wave functions.

  Language: C++17

  Patrick J. Fasano
  University of Notre Dame

  + 03/03/22 (pjf): Created.

****************************************************************/

#ifndef OSCILLATOR_WF_H_
#define OSCILLATOR_WF_H_

#include <cmath>
#include <vector>

namespace shell
{
template<int coord_sign, typename T>
std::vector<T> OscillatorWaveFunctions(
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
    return {};
  // output storage
  std::vector<T> wf(nmax + 1);

  // convenience variable
  //   x^2 = (r/b)^2 or x^2 = (k b)^2
  T x2;
  if constexpr (coord_sign == +1)
    x2 = ((coord * coord) / (b * b));
  else  // if (coord_sign == -1)
    x2 = ((coord * coord) * (b * b));

  // n = 0 explicit form
  if constexpr (coord_sign == +1)
  {
    wf[0] =
        (std::sqrt(2. / (b * std::tgamma(l + 1.5))) * std::pow(b, power_shift))
        * exp(-x2 / 2) * pow(x2, (l + 1 + power_shift) / 2);
  }
  else  // if (coord_sign == -1)
  {
    wf[0] =
        (std::sqrt((2. * b) / std::tgamma(l + 1.5)) * std::pow(b, -power_shift))
        * exp(-x2 / 2) * pow(x2, (l + 1 + power_shift) / 2);
  }
  if (nmax == 0)
    return wf;

  // n = 1 in terms of n = 0; only first term of recurrence
  wf[1] = (coord_sign * (((1.5 + l) - x2) / std::sqrt(l + 1.5))) * wf[0];
  if (nmax == 1)
    return wf;

  // recurrence from n = 2 to n = nmax
  for (int n = 2; n <= nmax; ++n)
  {
    wf[n] = coord_sign
            * (((2 * n + l - 0.5) - x2) * wf[n - 1]
               - std::sqrt((n - 1) * (n + l - 0.5)) * wf[n - 2])
            / std::sqrt(n * (n + l + 0.5));
  }
  return wf;
}
}  // namespace shell

#endif  // OSCILLATOR_WF_H_
