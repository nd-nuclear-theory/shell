/****************************************************************
  @file obme.h

  Defines functions for generating matrix elements of one-body operators.

  Language: C++11

  Patrick J. Fasano
  University of Notre Dame

  + 08/13/18 (pjf): Created.
  + 12/11/18 (pjf): Add additional operator generation functions:
    - LadderOneBodyOperator
    - AngularMomentumOneBodyOperator
    - AngularMomentumSquaredOneBodyOperator
    - IsospinOneBodyOperator
  + 05/09/19 (pjf): Use std::size_t for basis indices and sizes.
  + 04/19/22 (mac): Expand docstrings.
  + 04/19/22 (pjf): Add additional comments/documentation to docstrings.
****************************************************************/

#ifndef OBME_OBME_H_
#define OBME_OBME_H_

#include "am/rme.h"
#include "basis/nlj_orbital.h"
#include "obme/obme_operator.h"

namespace shell
{
void SolidHarmonicOneBodyOperator(
    shell::RadialBasisType basis_type,
    shell::RadialOperatorType operator_type,
    int order,
    const basis::OrbitalSpaceLJPN& space,
    const basis::OrbitalSectorsLJPN& sectors,
    basis::OperatorBlocks<double>& matrices
  );
// Generate reduced matrix elements of r^n C_J0 and (ik)^n C_J0.
//
// Spherical harmonic is in Racah's normalization [G. Racah, Phys. Rev. 62,
// 438 (1942), eqn. (46); Brink & Satchler (1993), eqn. (2.9)].
//
//   C_J0 = ((4*pi)/(2*J0+1))^(1/2) * Y_J0
//
// RMEs are in Rose convention.
//
// Examples:
//   - For order=0, J0=0, this is the identity operator.
//   - For order=1, J0=1, this is $\vec{r}$.
//   - For order=2, J0=0, this is $\vec{r}\cdot\vec{r}$.
//   - For order=2, J0=2, this is $\sqrt{3/2} [\vec{r}\otimes\vec{r}]_2$ [see
//       sec. 3.2, eqn. (23) of Varshalovich (1988)]; alternatively, this is
//       (4*pi/5)^(1/2) times the transition E2 operator Q=r^2Y_2, or 1/2 the
//       "quadrupole moment" operator Qmom [see footnote 10 of "intrinsic"
//       JPG 47, 122001 (2020)].
//
// Arguments:
//   basis_type (shell::RadialBasisType): radial basis type
//   operator_type (shell::RadialOperatorType): radial operator type
//   order (int): radial power
//   space (OrbitalSpaceLJPN): one-body space
//   sectors (OrbitalSectorsLJPN): one-body sectors
//   matrices (OperatorBlocks): operator matrices

void LadderOneBodyOperator(
    shell::RadialBasisType basis_type,
    shell::LadderOperatorType operator_type,
    const basis::OrbitalSpaceLJPN& space,
    const basis::OrbitalSectorsLJPN& sectors,
    basis::OperatorBlocks<double>& matrices
  );
// Generate reduced matrix elements of harmonic oscillator ladder operator
// $c^\dagger$ or $\tilde{c}$.
//
// TODO mac (3/24/23): Docstring still says c-tilde, but underlying code
// indicates c.  Confirm and document which it is.  These differ by a minus
// sign.  See Appendix F.1 of intrinsic.
//
// RMEs are in Rose convention.
//
// Arguments:
//   basis_type (shell::RadialBasisType): radial basis type
//   operator_type (shell::LadderOperatorType): radial operator type
//   space (OrbitalSpaceLJPN): one-body space
//   sectors (OrbitalSectorsLJPN): one-body sectors
//   matrices (OperatorBlocks): operator matrices

void AngularMomentumOneBodyOperator(
    am::AngularMomentumOperatorType operator_type,
    const basis::OrbitalSpaceLJPN& space,
    const basis::OrbitalSectorsLJPN& sectors,
    basis::OperatorBlocks<double>& matrices
  );
// Generate reduced matrix elements of angular momentum operator.
//
// RMEs are in Rose convention.
//
// Arguments:
//   operator_type (am::AngularMomentumOperatorType): angular momentum type
//   space (OrbitalSpaceLJPN): one-body space
//   sectors (OrbitalSectorsLJPN): one-body sectors
//   matrices (OperatorBlocks): operator matrices

void AngularMomentumSquaredOneBodyOperator(
    am::AngularMomentumOperatorType operator_type,
    const basis::OrbitalSpaceLJPN& space,
    const basis::OrbitalSectorsLJPN& sectors,
    basis::OperatorBlocks<double>& matrices
  );
// Generate matrix elements of angular momentum squared operator.
//
// RMEs are in Rose convention.
//
// Arguments:
//   operator_type (am::AngularMomentumOperatorType): angular momentum type
//   space (OrbitalSpaceLJPN): one-body space
//   sectors (OrbitalSectorsLJPN): one-body sectors
//   matrices (OperatorBlocks): operator matrices

void IsospinOneBodyOperator(
    const basis::OrbitalSpaceLJPN& space,
    const basis::OrbitalSectorsLJPN& sectors,
    basis::OperatorBlocks<double>& matrices
  );
// Generate reduced matrix elements of isospin operator.
//
// RMEs are in Rose convention.
//
// The spherical tensor component is deduced from Tz0 of input sectors, e.g.
// if Tz0=+1, matrix elements of T_+ are computed.
//
// Arguments:
//   space (OrbitalSpaceLJPN): one-body space
//   sectors (OrbitalSectorsLJPN): one-body sectors
//   matrices (OperatorBlocks): operator matrices

void OneBodyOperatorTensorProduct(
    const basis::OrbitalSpaceLJPN& space,
    const basis::OrbitalSectorsLJPN& sectors_a,
    const basis::OperatorBlocks<double>& matrices_a,
    const basis::OrbitalSectorsLJPN& sectors_b,
    const basis::OperatorBlocks<double>& matrices_b,
    const basis::OrbitalSectorsLJPN& sectors,
    basis::OperatorBlocks<double>& matrices);
// Compute reduced matrix elements of coupled tensor product [AB]_J of two
// one-body operators.
//
// RMEs are in Rose convention.
//
// Arguments:
//   space (basis::OrbitalSpaceLJPN): one-body space
//   sectors_a (basis::OrbitalSectors)

};      // namespace shell
#endif  // OBME_OBME_H_
