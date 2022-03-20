/****************************************************************
  mesh.h

  Defines recurrence relations for harmonic oscillator wave functions.

  Language: C++17

  Patrick J. Fasano
  University of Notre Dame

  + 03/07/22 (pjf): Created.

****************************************************************/

#ifndef MESH_H_
#define MESH_H_

#include <Eigen/Dense>
#include <cmath>
#include <type_traits>

namespace shell
{
inline auto UniformMesh(std::size_t num_points)
// Generate uniform mesh on closed interval [0,1].
//
// Note: this function *returns* an Eigen::ArrayXd instead of having a non-const
// reference parameter for output, as C++17 has made return value optimization
// mandatory, and thus there is no performance benefit to passing a non-const
// reference parameter for output.
//
// Arguments:
//   num_points (std::size_t): number of points
//
// Returns:
//   (Eigen::ArrayXd of length num_points): uniform mesh on interval [0,1]
{
  return Eigen::ArrayXd::LinSpaced(num_points, 0.0, 1.0).eval();
}

template<int coord_sign, typename T>
inline auto TransformedCoordinate(double b, const T& z)
// Generate coordinate transformed from [0,1] to [0,infinity].
//
// Note: this function *returns* auto instead of having a non-const reference
// parameter for output, as C++17 has made return value optimization mandatory,
// and thus there is no performance benefit to passing a non-const reference
// parameter for output.
//
// Template arguments:
//   coord_sign (int): +1 for position, -1 for momentum
//   T (typename): type of input point(s)
//
// Arguments:
//   b (double): oscillator length scale
//   z (T): point(s) in the interval [0,1]
//
// Returns:
//   (auto): transformed point(s)
{
  if constexpr (coord_sign == +1)
    return ((z / (1. - z)) * b);
  else  // if (coord_sign == -1)
    return ((z / (1. - z)) / b);
}

template<int coord_sign, typename T>
inline auto TransformedCoordinateJacobian(double b, const T& z)
// Generate Jacobian associated with transforming coordinate from [0,1] to [0,infinity].
//
// Note: this function *returns* auto instead of having a non-const reference
// parameter for output, as C++17 has made return value optimization mandatory,
// and thus there is no performance benefit to passing a non-const reference
// parameter for output.
//
// Template arguments:
//   coord_sign (int): +1 for position, -1 for momentum
//   T (typename): type of input point(s)
//
// Arguments:
//   b (double): oscillator length scale
//   z (T): point(s) in the interval [0,1]
//
// Returns:
//   (auto): Jacobian for transformation, evaluated at point(s)
{
  if constexpr (coord_sign == +1)
    return ((1. / ((1. - z) * (1. - z))) * b);
  else  // if (coord_sign == -1)
    return ((1. / ((1. - z) * (1. - z))) / b);
}

template<typename T, std::enable_if_t<std::is_floating_point_v<typename T::value_type>>* = nullptr>
inline void FixEndpointSingularitiesToZero(T& values)
// Fix endpoint singularities to zero.
//
// This function assumes that a value of NaN at the endpoint of the array-like
// corresponds to a removable point singularity with limit zero.
//
// Arguments:
//   values (T): array-like of floating-point values to be modified
{
  if (std::isnan(values[0]))
    values[0] = 0.;
  if (std::isnan(values[values.size() - 1]))
    values[values.size() - 1] = 0.;
}

}  // namespace shell

#endif  // OSCILLATOR_WF_H_
