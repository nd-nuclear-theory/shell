/****************************************************************
  obme_operator.cpp

  Patrick J. Fasano
  University of Notre Dame

****************************************************************/

#include "obme/obme_operator.h"

namespace shell
{
  const std::unordered_map<std::string, RadialOperatorType> kCharCodeRadialOperatorType({
    {"r", RadialOperatorType::kR},
    {"k", RadialOperatorType::kK},
    {"o", RadialOperatorType::kO},
    {"g", RadialOperatorType::kGeneric}
    });
}  // namespace shell
