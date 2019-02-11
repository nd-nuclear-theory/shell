/****************************************************************
  @file obme_operator.h

  Defines storage classes for one-body operators.

  Language: C++11

  Patrick J. Fasano
  University of Notre Dame

  + 02/20/18 (pjf): Created.
  + 12/11/18 (pjf): Add LadderOperatorType enum.
****************************************************************/

#ifndef OBME_OPERATOR_H_
#define OBME_OPERATOR_H_

#include "basis/nlj_operator.h"
#include "basis/operator.h"

namespace shell {

/**
 * Radial IDs
 */
enum class RadialOperatorType : char {
  kR = 'r',       // radius
  kK = 'k',       // momentum
  kO = 'o',       // overlaps
  kGeneric = 'g'  // generic
};

/**
 * Ladder IDs
 */
enum class LadderOperatorType {
  kRaising,
  kLowering,
};

/**
 * Basis IDs
 */
enum class RadialBasisType : char {
  kGeneric = 'g',     // general
  kOscillator = 'o',  // harmonic oscillator
  kLaguerre = 'l'     // Laguerre
};

}
#endif  // OBME_OPERATOR_H_
