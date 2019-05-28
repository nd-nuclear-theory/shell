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
// Generate matrix elements of r^n C_j0 and (ik)^n C_j0 in Rose convention.
//
// For j0=1, order=1, this is \vec{r}; for j0=0, order=2, this
// coincides with r.r
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
// Generate matrix elements of $b^\dagger$ and $\tilde{b}$ in Rose convention.
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
// Generate matrix elements of angular momentum operator in Rose convention.
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
// Generate matrix elements of angular momentum squared operator in Rose convention.
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
// Generate matrix elements of isospin operator in Rose convention.
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
// Compute coupled tensor product [AB]_J of two one-body operators.
//
// Arguments:
//   space (basis::OrbitalSpaceLJPN): one-body space
//   sectors_a (basis::OrbitalSectors)

};      // namespace shell
#endif  // OBME_OBME_H_
