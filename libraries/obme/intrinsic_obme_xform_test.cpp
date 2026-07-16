/****************************************************************
  intrinsic_obme_xform_test.cpp

  Mark A. Caprio and Victor Dumenil
  University of Notre Dame

****************************************************************/

#include <Eigen/Core>

#include "basis/operator.h"

#include "obme/intrinsic_obme_xform.h"

////////////////////////////////////////////////////////////////
// test code
////////////////////////////////////////////////////////////////

void TestOneBodyOperatorDeltaNSubspace()
{
  // subspace construction
  //
  // Normally subspaces are always constructed as part of a space.  We would not
  // construct a standalone subspace.  But, for testing purposes, it is
  // convenient to do so, since the code for subspace and state can be written
  // and tested independently of the code for space and/or sectors.

  std::cout << "Subspace construction" << std::endl;
  std::cout << std::endl;

  const shell::OneBodyOperatorDeltaNSubspace subspace(0, 0, 0, 2, 4);  // J0=0, g0=0, Delta_N=0, N1max=2, N2max=4
  std::cout << subspace.LabelStr() << std::endl;
  std::cout << subspace.DebugStr();

  // // state construction
  // //
  // // Let us give the state type a more thorough workout.  Note that the state
  // // type was already used in the implementation of the subspace's DebugStr()
  // // above.
  // 
  // std::cout << std::endl;
  // std::cout << "State construction" << std::endl;
  // std::cout << std::endl;
  // 
  // // illustrate construction by index vs. state labels
  // 
  // // construct by index -- invalid cases
  // 
  // if (false)
  //   {
  //     // Index -1 should be trapped and lead to assertion failure.
  //     //
  //     // Note: Negative index cast to size_t will be large positive index, by
  //     // two-complement.
  // 
  //     std::size_t test_index = (std::size_t)(-1);
  //     std::cout << "Constructing with index " << test_index << std::endl;
  //     const basis::OscillatorOrbitalState state_from_index(subspace, test_index);
  //   }
  // 
  // if (false)
  //   {
  //     // But a very large negative integer could wrap and become a valid index again!
  //     //
  //     // std::size_t test_index = (std::size_t)(-18446744073709551615);
  // 
  //     std::size_t max_size = std::numeric_limits<std::size_t>().max();
  //     std::size_t test_index = -max_size;  // should be equivalent to +1
  //     std::cout << "Constructing with index " << test_index << std::endl;
  //     const basis::OscillatorOrbitalState state_from_index(subspace, test_index);
  //   }
  // 
  // // construct by index -- valid case
  // std::size_t test_index = 1;
  // std::cout << "Constructing with index " << test_index << std::endl;
  // const basis::OscillatorOrbitalState state_from_index(subspace, test_index);
  // std::cout << "  label string " << state_from_index.LabelStr() << std::endl;
  // std::cout << "  quantum numbers"
  //           << " l " << state_from_index.l()
  //           << " g " << state_from_index.g()
  //           << " n " << state_from_index.n()
  //           << " N " << state_from_index.N()
  //           << std::endl;
  // 
  // // construct by state labels
  // int n = 0;
  // std::cout << "Constructing with n " << n << std::endl;
  // const basis::OscillatorOrbitalState state_from_labels(subspace, basis::OscillatorOrbitalState::StateLabelsType(n));
  // std::cout << "  label string " << state_from_labels.LabelStr() << std::endl;
  // std::cout << "  quantum numbers"
  //           << " l " << state_from_labels.l()
  //           << " g " << state_from_labels.g()
  //           << " n " << state_from_labels.n()
  //           << " N " << state_from_labels.N()
  //           << std::endl;
  // 
  // // try out accessors
  // std::cout << "Try out accessors" << std::endl;
  // std::cout << "subspace " << state_from_labels.subspace().LabelStr() << std::endl;
  // std::cout << "labels " << std::get<0>(state_from_labels.labels()) << std::endl;
  // std::cout << "index " << state_from_labels.index() << std::endl;
  // 
  // // text equality operator
  // std::cout << "  equality test " << (state_from_index == state_from_index) << std::endl;
  // std::cout << "  equality test " << (state_from_index == state_from_labels) << std::endl;
  // 
  // // iterate over states within subspace
  // //
  // // here we also try out the "state factory" member function of the subspace
  // std::cout << std::endl;
  // std::cout << "Iterate over states in subspace" << std::endl;
  // std::cout << std::endl;
  // 
  // // iterate by state index
  // std::cout << "by index" << std::endl;
  // for (std::size_t state_index=0; state_index<subspace.size(); ++state_index)
  //   {
  //     // const basis::OscillatorOrbitalState state(subspace,state_index);
  //     const basis::OscillatorOrbitalState state = subspace.GetState(state_index);
  //     std::cout << "index " << state.index() << " N " << state.N() << std::endl;
  //   };
  // 
  // // iterate by labels
  // std::cout << "by labels" << std::endl;
  // int n_max = (subspace.Nmax()-subspace.l())/2;
  // for (int n=0; n<=n_max; ++n)
  //   {
  //     const basis::OscillatorOrbitalState state = subspace.GetState(basis::OscillatorOrbitalState::StateLabelsType(n));
  //     std::cout << "index " << state.index() << " N " << state.N() << std::endl;
  //   };
  // 
  // std::cout << std::endl;

}


void TestOneBodyOperatorDeltaNSpace()
{

  std::cout << "Space construction" << std::endl;
  std::cout << std::endl;

  // construct space
  const shell::OneBodyOperatorDeltaNSpace space(0, 0, 2, 2, 4);  // J0=0, g0=0, Delta_N_max=2, N1max=2, N2max=4

  // print diagnostics
  std::cout << space.DebugStr();

  std::cout << std::endl;

  // try out accessors
  // std::cout << "Try out accessors" << std::endl;
  // std::cout << "Nmax " << space.Nmax() << std::endl;
  // std::cout << "ContainsSubspace " << space.ContainsSubspace(basis::OneBodyOperatorDeltaNSubspace::SubspaceLabelsType(5)) << std::endl;
  // std::cout << "LookUpSubspaceIndex " << space.LookUpSubspaceIndex(basis::OneBodyOperatorDeltaNSubspace::SubspaceLabelsType(2)) << std::endl;
  // std::cout << "LookUpSubspace " << space.LookUpSubspace(basis::OneBodyOperatorDeltaNSubspace::SubspaceLabelsType(2)).LabelStr() << std::endl;
  // std::cout << "size " << space.size() << std::endl;
  // std::cout << "dimension " << space.dimension() << std::endl;
  // 
  // const basis::OneBodyOperatorDeltaNSubspace& subspace = space.LookUpSubspace(basis::OneBodyOperatorDeltaNSubspace::SubspaceLabelsType(2));
  // std::cout << "LookUpSubspace with a reference variable " << subspace.LabelStr() << std::endl;
}

// void TestOneBodyOperatorDeltaNSectors()
// {
// 
//   std::cout << "Sectors construction" << std::endl;
//   std::cout << std::endl;
// 
//   // construct space
//   int Nmax = 4;
//   basis::OneBodyOperatorDeltaNSpace space(Nmax);
// 
//   // construct sectors -- L0=0, g0=0 (Hamiltonian-like)
//   basis::OneBodyOperatorDeltaNSectors hamiltonian_sectors(space, 0, 0);
//   std::cout << "hamiltonian_sectors" << std::endl
//             << hamiltonian_sectors.DebugStr()
//             << std::endl;
// 
//   // find indices for a matrix element from labels
//   std::cout << "Find indices for a matrix element from labels" << std::endl;
//   std::size_t bra_subspace_index=space.LookUpSubspaceIndex(basis::OneBodyOperatorDeltaNSubspace::SubspaceLabelsType(2));
//   std::size_t ket_subspace_index=space.LookUpSubspaceIndex(basis::OneBodyOperatorDeltaNSubspace::SubspaceLabelsType(2));
//   std::size_t sector_index=hamiltonian_sectors.LookUpSectorIndex(bra_subspace_index,ket_subspace_index);
//   std::cout << "sector index " << sector_index << std::endl;
//   basis::OneBodyOperatorDeltaNSectors::SectorType sector=hamiltonian_sectors.GetSector(sector_index);
//   // const basis::OneBodyOperatorDeltaNSectors::SectorType::BraSubspaceType& bra_subspace=sector.bra_subspace();
//   // const basis::OneBodyOperatorDeltaNSectors::SectorType::KetSubspaceType& ket_subspace=sector.ket_subspace();
//   const basis::OneBodyOperatorDeltaNSubspace& bra_subspace=sector.bra_subspace();
//   const basis::OneBodyOperatorDeltaNSubspace& ket_subspace=sector.ket_subspace();
//   // auto& bra_subspace=sector.bra_subspace();
//   // auto& ket_subspace=sector.ket_subspace();
//   std::cout << "bra_subspace " << bra_subspace.LabelStr() << std::endl;
//   std::cout << "ket_subspace " << ket_subspace.LabelStr() << std::endl;
//   std::cout << "bra_subspace " << std::endl << bra_subspace.DebugStr() << std::endl;
//   std::cout << "ket_subspace " << std::endl << ket_subspace.DebugStr() << std::endl;
//   std::size_t bra_state_index=bra_subspace.LookUpStateIndex(basis::OneBodyOperatorDeltaNState::StateLabelsType(0));
//   std::size_t ket_state_index=ket_subspace.LookUpStateIndex(basis::OneBodyOperatorDeltaNState::StateLabelsType(1));
//   std::cout << "bra state index " << bra_state_index << std::endl;
//   std::cout << "ket state index " << ket_state_index << std::endl;
// 
//   // construct sectors -- L0=1, g0=1 (E1-like)
//   basis::OneBodyOperatorDeltaNSectors e1_sectors(space, 1, 1);
//   std::cout << "e1_sectors" << std::endl
//             << e1_sectors.DebugStr()
//             << std::endl;
// 
//   // construct sectors -- L0=1, g0=1 (E1-like) -- but including noncanonical ("lower-triangle") sectors
//   basis::OneBodyOperatorDeltaNSectors e1_sectors_both_ways(
//       space, 1, 1,
//       basis::SectorDirection::kBoth
//     );
//   std::cout << "e1_sectors_both_ways" << std::endl
//             << e1_sectors_both_ways.DebugStr()
//             << std::endl;
// 
// 
//   // construct sectors -- L0=2, g0=0 (E2-like)
//   basis::OneBodyOperatorDeltaNSectors e2_sectors(space, 2, 0);
//   std::cout << "e2_sectors" << std::endl
//             << e2_sectors.DebugStr()
//             << std::endl;
// 
//   // construct sectors -- L0=2, g0=1 (M2-like)
//   basis::OneBodyOperatorDeltaNSectors m2_sectors(space, 2, 1);
//   std::cout << "n2_sectors" << std::endl
//             << m2_sectors.DebugStr()
//             << std::endl;
// 
// }


////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  TestOneBodyOperatorDeltaNSubspace();
  TestOneBodyOperatorDeltaNSpace();
    
  // termination
  return EXIT_SUCCESS;
}
