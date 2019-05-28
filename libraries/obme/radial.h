/****************************************************************
  @file radial.h

  Defines functions for generating matrix elements of radial operators.

  Language: C++11

  Patrick J. Fasano
  University of Notre Dame

  + 06/15/18 (pjf): Created.
  + 02/01/19 (pjf): Add documentation.
  + 05/09/19 (pjf): Use std::size_t for basis indices and sizes.
****************************************************************/

#ifndef OBME_RADIAL_H_
#define OBME_RADIAL_H_

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "basis/nlj_operator.h"
#include "basis/nlj_orbital.h"
#include "obme/obme_operator.h"

namespace shell
{
void GenerateRadialOperator(
    shell::RadialBasisType basis_type,
    shell::RadialOperatorType operator_type,
    int order,
    const basis::OrbitalSpaceLJPN& space,
    const basis::OrbitalSectorsLJPN& sectors,
    basis::OperatorBlocks<double>& matrices
  );
// Generate radial matrix elements of r^n and (ik)^n.
//
// Arguments:
//   basis_type (shell::RadialBasisType): radial basis type
//   operator_type (shell::RadialOperatorType): radial operator type
//   order (int): radial power
//   space (OrbitalSpaceLJPN): one-body space
//   sectors (OrbitalSectorsLJPN): one-body sectors
//   matrices (OperatorBlocks): operator matrices

void GenerateRadialOverlaps(
    shell::RadialBasisType bra_basis_type,
    shell::RadialBasisType ket_basis_type,
    double scale_factor,
    const basis::OrbitalSpaceLJPN& bra_space,
    const basis::OrbitalSpaceLJPN& ket_space,
    const basis::OrbitalSectorsLJPN& sectors,
    basis::OperatorBlocks<double>& matrices
  );
// Generate radial overlap matrix elements.
//
// Calculates overlaps using spline integration. Scale factor is
// as defined in mac thesis.
//
// Arguments:
//   bra_basis_type (shell::RadialBasisType): bra radial basis type
//   ket_basis_type (shell::RadialBasisType): ket radial basis type
//   scale_factor (double): scale factor between ket and bra radial functions
//   bra_space (OrbitalSpaceLJPN): bra one-body space
//   ket_space (OrbitalSpaceLJPN): ket one-body space
//   sectors (OrbitalSectorsLJPN): one-body sectors
//   matrices (OperatorBlocks): operator matrices

void ComposeRadialOperators(
    const basis::OrbitalSpaceLJPN& bra_space_a,
    const basis::OrbitalSpaceLJPN& ket_space_a,
    const basis::OrbitalSectorsLJPN& sectors_a,
    const basis::OperatorBlocks<double>& matrices_a,
    const basis::OrbitalSpaceLJPN& bra_space_b,
    const basis::OrbitalSpaceLJPN& ket_space_b,
    const basis::OrbitalSectorsLJPN& sectors_b,
    const basis::OperatorBlocks<double>& matrices_b,
    const basis::OrbitalSectorsLJPN& sectors,
    basis::OperatorBlocks<double>& matrices
  );
// Generate product of two radial operators.
//
// Calculates a.b product as effective matrix-matrix multiplication,
// with appropriate adaptation to sparse (by sector) storage.
//
// Arguments:
//   bra_space_a (OrbitalSpaceLJPN): bra one-body space for operator a
//   ket_space_a (OrbitalSpaceLJPN): ket one-body space for operator a
//   sectors_a (OrbitalSectorsLJPN): one-body sectors for operator a
//   matrices_a (OperatorBlocks): operator matrices for operator a
//   bra_space_b (OrbitalSpaceLJPN): bra one-body space for operator b
//   ket_space_b (OrbitalSpaceLJPN): ket one-body space for operator b
//   sectors_b (OrbitalSectorsLJPN): one-body sectors for operator b
//   matrices_b (OperatorBlocks): operator matrices for operator b
//   sectors (OrbitalSectorsLJPN): output one-body sectors
//   matrices (OperatorBlocks): output operator matrices

void SimilarityTransformOperator(
    const basis::OrbitalSpaceLJPN& source_space,
    const basis::OrbitalSpaceLJPN& target_space,
    const basis::OrbitalSectorsLJPN& xform_sectors,
    const basis::OperatorBlocks<double>& xform_matrices,
    const basis::OrbitalSectorsLJPN& operator_sectors,
    const basis::OperatorBlocks<double>& operator_matrices,
    const basis::OrbitalSectorsLJPN& output_sectors,
    basis::OperatorBlocks<double>& output_matrices
  );
// Generate similarity-transformed operator.
//
// Transforms operator via $U^T O U$ and constructs in new ouput matrices.
//
// Arguments:
//   source_space (OrbitalSpaceLJPN): space for source operator
//   target_space (OrbitalSpaceLJPN): space for operator after transformation
//   xform_sectors (OrbitalSectorsLJPN): sectors for transformation
//   xform_matrices (OperatorBlocks): matrices for transformation
//   operator_sectors (OrbitalSectorsLJPN): sectors for source operator
//   operator_matrices (OperatorBlocks): matrices for source operator
//   output_sectors (OrbitalSectorsLJPN): output sectors for transformed operator
//   output_matrices (OperatorBlocks): output matrices for transformed operator

};      // namespace shell
#endif  // OBME_RADIAL_H_
