/****************************************************************
  @file obme_operator.h

  Defines storage classes for one-body operators.

  Language: C++11

  Patrick J. Fasano
  University of Notre Dame

  + 02/20/18 (pjf): Created.
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
 * One-body operator storage
 */
struct OneBodyOperator : public basis::OneBodyOperatorLJPN {
  RadialOperatorType radial_operator_type;
  int radial_operator_power;
};

}
#endif  // OBME_OPERATOR_H_
