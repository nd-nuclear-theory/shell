/****************************************************************
  tbme_mapping_test.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "tbme/tbme_mapping.h"

#include <algorithm>

#include "fmt/format.h"

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

void TestMappingRenumber()
{

  std::cout << "TestMappingRenumber" << std::endl;

  // source space -- square
  const int source_Nmax = 2;
  const basis::OrbitalSpacePN source_orbital_space(source_Nmax);
  const basis::TwoBodySpaceJJJPN source_space(
      source_orbital_space,basis::WeightMax(basis::Rank::kTwoBody,source_Nmax)
    );
  std::cout << source_orbital_space.DebugStr();
  std::cout << basis::OrbitalDefinitionStr(source_orbital_space.OrbitalInfo()) << std::endl;
  std::cout << source_space.DebugStr();

  // target space -- permuted
  //
  // reverse all orbitals (but it's okay that proton and neutron
  // orbitals are no longer in order in the orbital info lists, as
  // they will be partitioned back into their appropriate subspaces)

  const int target_Nmax = 2;
  std::vector<basis::OrbitalPNInfo> source_orbital_info = source_orbital_space.OrbitalInfo();
  std::vector<basis::OrbitalPNInfo> target_orbital_info;
  target_orbital_info.resize(source_orbital_info.size());
  std::reverse_copy(
      std::begin(source_orbital_info),std::end(source_orbital_info),
      std::begin(target_orbital_info)
    );
  const basis::OrbitalSpacePN target_orbital_space(target_orbital_info);
  const basis::TwoBodySpaceJJJPN target_space(
      target_orbital_space,basis::WeightMax(basis::Rank::kTwoBody,target_Nmax)
    );
  std::cout << target_orbital_space.DebugStr() << std::endl;
  std::cout << basis::OrbitalDefinitionStr(target_orbital_space.OrbitalInfo()) << std::endl;
  std::cout << target_space.DebugStr();

  // const basis::OrbitalSpacePN target_orbital_space(target_Nmax);
  // const basis::TwoBodySpaceJJJPN target_space(source_orbital_space,basis::WeightMax(basis::Rank::kTwoBody,target_Nmax));

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
  TestMappingRenumber();

  // termination
  return 0;
}
