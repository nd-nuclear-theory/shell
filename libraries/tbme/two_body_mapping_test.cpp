/****************************************************************
  two_body_mapping_test.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "tbme/two_body_mapping.h"

#include "cppformat/format.h"

void TestMappingEqual()
{

  std::cout << "TestMappingEqual" << std::endl;

  // source space -- triangle
  const int source_Nmax = 2;
  const basis::OrbitalSpacePN source_orbital_space(source_Nmax);
  const basis::TwoBodySpaceJJJPN source_space(source_orbital_space,basis::WeightMax(basis::Rank::kTwoBody,source_Nmax));

  // target space -- triangle
  const int target_Nmax = 2;
  const basis::OrbitalSpacePN target_orbital_space(target_Nmax);
  const basis::TwoBodySpaceJJJPN target_space(source_orbital_space,basis::WeightMax(basis::Rank::kTwoBody,target_Nmax));

  // set up mapping
  const shell::TwoBodyMapping two_body_mapping(
      source_orbital_space,source_space,
      target_orbital_space,target_space
    );

  std::cout << two_body_mapping.DebugStr() << std::endl;

}

void TestMappingSubset()
{

  std::cout << "TestMappingSubset" << std::endl;

  // source space -- square
  const int source_Nmax = 4;
  const basis::OrbitalSpacePN source_orbital_space(source_Nmax);
  const basis::TwoBodySpaceJJJPN source_space(source_orbital_space,basis::WeightMax(basis::Rank::kTwoBody,source_Nmax));

  // target space -- subtriangle
  const int target_Nmax = 2;
  const basis::OrbitalSpacePN target_orbital_space(target_Nmax);
  const basis::TwoBodySpaceJJJPN target_space(source_orbital_space,basis::WeightMax(basis::Rank::kTwoBody,target_Nmax));

  // set up mapping
  const shell::TwoBodyMapping two_body_mapping(
      source_orbital_space,source_space,
      target_orbital_space,target_space
    );

  std::cout << two_body_mapping.DebugStr() << std::endl;

}

void TestMappingSuperset()
{

  std::cout << "TestMappingSuperset" << std::endl;

  // source space -- square
  const int source_Nmax = 2;
  const basis::OrbitalSpacePN source_orbital_space(source_Nmax);
  const basis::TwoBodySpaceJJJPN source_space(source_orbital_space,basis::WeightMax(basis::Rank::kTwoBody,source_Nmax));

  // target space -- subtriangle
  const int target_Nmax = 4;
  const basis::OrbitalSpacePN target_orbital_space(target_Nmax);
  const basis::TwoBodySpaceJJJPN target_space(source_orbital_space,basis::WeightMax(basis::Rank::kTwoBody,target_Nmax));

  // set up mapping
  const shell::TwoBodyMapping two_body_mapping(
      source_orbital_space,source_space,
      target_orbital_space,target_space
    );

  std::cout << two_body_mapping.DebugStr() << std::endl;

}

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  TestMappingEqual();
  TestMappingSubset();
  TestMappingSuperset();

  // termination
  return 0;
}
