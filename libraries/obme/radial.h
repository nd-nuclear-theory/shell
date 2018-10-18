/****************************************************************
  @file radial.h

  Defines functions for generating matrix elements of radial operators.

  Language: C++11

  Patrick J. Fasano
  University of Notre Dame

  + 06/15/18 (pjf): Created.
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
    basis::OperatorBlocks<double>& matrices);

void GenerateRadialOverlaps(
    shell::RadialBasisType bra_basis_type,
    shell::RadialBasisType ket_basis_type,
    double scale_factor,
    const basis::OrbitalSpaceLJPN& bra_space,
    const basis::OrbitalSpaceLJPN& ket_space,
    const basis::OrbitalSectorsLJPN& sectors,
    basis::OperatorBlocks<double>& matrices);

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
    basis::OperatorBlocks<double>& matrices);

void SimilarityTransformOperator(
    const basis::OrbitalSpaceLJPN& source_space,
    const basis::OrbitalSpaceLJPN& target_space,
    const basis::OrbitalSectorsLJPN& xform_sectors,
    const basis::OperatorBlocks<double>& xform_matrices,
    const basis::OrbitalSectorsLJPN& operator_sectors,
    const basis::OperatorBlocks<double>& operator_matrices,
    const basis::OrbitalSectorsLJPN& output_sectors,
    basis::OperatorBlocks<double>& output_matrices);

};      // namespace shell
#endif  // OBME_RADIAL_H_
