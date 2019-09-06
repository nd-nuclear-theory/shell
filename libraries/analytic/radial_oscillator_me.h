/****************************************************************
  @file radial_oscillator_me.h

  Defines functions for generating matrix elements of radial operators
  in the harmonic oscillator basis.

  Language: C++11

  Patrick J. Fasano
  University of Notre Dame

  + 08/07/18 (pjf): Created.
  + 10/17/18 (pjf): Fix signs on k matrix elements.
  + 01/11/19 (pjf):
    - Define k_sqrt_2.
    - Move ladder oscillator matrix element functions to top and define
      coordinate matrix elements in terms of them.
  + 09/03/19 (pjf):
    - Remove unused variables from CoordinateOscillatorMatrixElement.
    - Fix operator signs on CoordinateSqrOscillatorMatrixElement.
****************************************************************/

#ifndef RADIAL_OSCILLATOR_ME_H_
#define RADIAL_OSCILLATOR_ME_H_

#include <cassert>
#include <cmath>

namespace analytic
{
constexpr double k_sqrt_2 = 1.414213562373095;

inline double CDaggerOscillatorMatrixElement(int bra_N, int bra_l, int ket_N, int ket_l)
// Calculate matrix elements of $C^\dagger$ operator in oscillator radial basis.
//
// Arguments:
//   bra_N, ket_N (input): bra and ket N(=2n+l) quantum numbers
//   bra_l, ket_l (input): bra and ket angular momentum quantum numbers
//
// Returns:
//   radial matrix element
{
  const int delta_l = bra_l - ket_l;
  const int delta_N = bra_N - ket_N;
  // cast to double to ensure floating-point arithmetic
  const double N = ket_N;
  const double l = ket_l;
  double matrix_element = 0.;
  if (delta_N == +1) {
    if (delta_l == +1) {
      matrix_element = std::sqrt(N + l + 3);
    } else if (delta_l == -1) {
      matrix_element = -std::sqrt(N - l + 2);
    }
  }

  return matrix_element;
}

inline double COscillatorMatrixElement(int bra_N, int bra_l, int ket_N, int ket_l)
// Calculate matrix elements of $c$ (lowering) operator in oscillator radial basis.
//
// Arguments:
//   bra_N, ket_N (input): bra and ket N(=2n+l) quantum numbers
//   bra_l, ket_l (input): bra and ket angular momentum quantum numbers
//
// Returns:
//   radial matrix element
{
  const int delta_l = bra_l - ket_l;
  const int delta_N = bra_N - ket_N;
  // cast to double to ensure floating-point arithmetic
  const double N = ket_N;
  const double l = ket_l;
  double matrix_element = 0.;
  if (delta_N == -1) {
    if (delta_l == +1) {
      matrix_element = -std::sqrt(N - l);
    } else if (delta_l == -1) {
      matrix_element = std::sqrt(N + l + 1);
    }
  }

  return matrix_element;
}

inline double CoordinateOscillatorMatrixElement(
    int bra_N, int bra_l, int ket_N, int ket_l, int operator_sign)
// Calculate matrix elements for r and ik operator in oscillator radial basis.
//
// Calculated from SU(1,1) algebraic expressions in (64) & (65) of D. J. Rowe,
// JPA 38, 10181 (2005), but with phase converted to "positive at origin"
// convention on the radial wave functions.  We use lambda=v+N/2 for
// oscillator functions, with N=3 and v=l in three dimensions.  Radial matrix
// elements must be complemented with angular matrix elements given by Rowe
// (97).
//
// While Rowe uses "positive at infinity" convention for the radial wave
// functions, we commonly use "positive at origin" convention.  Conversion
// introduces a factor (-)^(bra_n+ket_n), i.e., adding a (-) sign on the
// delta_N+delta_l=0 terms.
//
// Arguments:
//   bra_N, ket_N (input): bra and ket N(=2n+l) quantum numbers
//   bra_l, ket_l (input): bra and ket angular momentum quantum numbers
//   operator_sign (int): sign selecting coordinate or momentum
//     representation (+1 for "r", -1 for "k")
//
// Returns:
//   radial matrix element
{
  const int delta_l = bra_l - ket_l;
  const int delta_N = bra_N - ket_N;
  assert((delta_l == -1) || (delta_l == +1));

  double matrix_element = 0.;
  matrix_element =
    (operator_sign*CDaggerOscillatorMatrixElement(bra_N, bra_l, ket_N, ket_l)
    + COscillatorMatrixElement(bra_N, bra_l, ket_N, ket_l))/k_sqrt_2;

  return matrix_element;
}

inline double CoordinateSqrOscillatorMatrixElement(
    int bra_N, int bra_l, int ket_N, int ket_l, int operator_sign)
// Calculate matrix elements of r^2 and k^2 operators in oscillator radial basis.
//
// For r^2: Calculated from SU(1,1) algebraic expression in (39) of
// D. J. Rowe, JPA 38, 10181 (2005), but with phase converted to "positive at
// origin" convention on the radial wave functions.  We use lambda=v+N/2 for
// oscillator functions, with N=3 and v=l in three dimensions.
//
// - for (delta L)=0, use (delta l)=0 matrix elements of r^2
//   directly from (39)
//
// - for (delta L)=2, use resolution of identity over
//   intermediate L space, by double application of the (delta
//   l)=1 matrix elements of r from (64)&(65)
//
// For k^2: This is the "reduced" kinetic energy operatator, or k^2=-del^2
// operator, where p = hbar * k.  Calculated from SU(1,1) algebraic expression
// in (48) [making use of (41) and (47)] of D. J. Rowe, JPA 38, 10181 (2005),
// but with phase converted to "positive at origin" convention on the radial
// wave functions.  Note that the last term, or 1/r^2 term, of (48) drops out
// since lambda=v+N/2 for oscillator functions.  We use lambda=v+N/2 for
// oscillator functions, with N=3 and v=l in three dimensions.
//
// While Rowe uses "positive at infinity" convention for the radial wave
// functions, we commonly use "positive at origin" convention.  Conversion
// introduces a factor (-)^(bra_n+ket_n), i.e., adding a (-) sign on the
// |delta_N+delta_l|=2 terms.
//
// Arguments:
//   bra_N, ket_N (input): bra and ket N(=2n+l) quantum numbers
//   bra_l, ket_l (input): bra and ket angular momentum quantum numbers
//   operator_sign (int): sign selecting coordinate or momentum
//     representation (+1 for "r", -1 for "k")
//
// Returns:
//   radial matrix element
{
  const int delta_l = bra_l - ket_l;
  const int delta_N = bra_N - ket_N;

  // cast to double to ensure floating-point arithmetic
  const double N = ket_N;
  const double l = ket_l;
  double matrix_element = 0.;
  if (delta_N == +2) {
    if (delta_l == +2) {
      matrix_element =  0.5 * std::sqrt((N + l + 5) * (N + l + 3));
    } else if (delta_l == 0) {
      matrix_element = -0.5 * std::sqrt((N + l + 3) * (N - l + 2));
    } else if (delta_l == -2) {
      matrix_element =  0.5 * std::sqrt((N - l + 2) * (N - l + 4));
    }
  } else if (delta_N == 0) {
    if (delta_l == +2) {
      matrix_element = -operator_sign * std::sqrt((N - l) * (N + l + 3));
    } else if (delta_l == 0) {
      matrix_element =  operator_sign * (N + 1.5);
    } else if (delta_l == -2) {
      matrix_element = -operator_sign * std::sqrt((N - l + 2) * (N + l + 1));
    }
  } else if (delta_N == -2) {
    if (delta_l == +2) {
      matrix_element =  0.5 * std::sqrt((N - l - 2) * (N - l));
    } else if (delta_l == 0) {
      matrix_element = -0.5 *std::sqrt((N - l) * (N + l + 1));
    } else if (delta_l == -2) {
      matrix_element =  0.5 * std::sqrt((N + l + 1) * (N + l - 1));
    }
  }

  return matrix_element;
}

};      // namespace analytic
#endif  // RADIAL_OSCILLATOR_ME_H_
