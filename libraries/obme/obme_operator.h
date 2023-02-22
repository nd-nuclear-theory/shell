/****************************************************************
  @file obme_operator.h

  Defines storage classes for one-body operators.

  Language: C++11

  Patrick J. Fasano
  University of Notre Dame

  + 02/20/18 (pjf): Created.
  + 12/11/18 (pjf): Add LadderOperatorType enum.
  + 08/26/19 (pjf): Add reverse definition for RadialOperatorType.
  + 09/08/20 (pjf): Add "ik" as alias for RadialOperatorType::kK.
****************************************************************/

#ifndef OBME_OPERATOR_H_
#define OBME_OPERATOR_H_

#include <string>
#include <unordered_map>

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

// notational reverse definitions for operator types
//
// Example:
//   std::string operator_type_code = "p";
//   ...
//   os << shell::kCharCodeRadialOperatorType[operator_type_code];
extern const std::unordered_map<std::string, RadialOperatorType> kCharCodeRadialOperatorType;

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
