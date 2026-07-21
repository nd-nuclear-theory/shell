/****************************************************************
  intrinsic_obme_xform_test.cpp

  Mark A. Caprio and Victor Dumenil
  University of Notre Dame

****************************************************************/

#include <Eigen/Core>

#include "basis/nlj_orbital.h"
#include "basis/operator.h"
#include "mcutils/eigen.h"

#include "obme/intrinsic_obme_xform.h"
#include "moshinsky/moshinsky_bracket.h"

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

  std::cout << std::endl;

}


void TestOneBodyOperatorDeltaNSpace()
{

  std::cout << "Space construction" << std::endl;
  std::cout << std::endl;

  // construct space
  const shell::OneBodyOperatorDeltaNSpace space(0, 0, 2, 2, 4);  // J0=0, g0=0, Delta_N_max=2, N1max=2, N2max=4
  // const shell::OneBodyOperatorDeltaNSpace space(0, 0, 1, 2, 4);  // J0=0, g0=0, Delta_N_max=1, N1max=2, N2max=4  -- parity inconsistent (assertion fails)
  // const shell::OneBodyOperatorDeltaNSpace space(0, 1, 1, 2, 4);  // J0=0, g0=1, Delta_N_max=1, N1max=2, N2max=4
  // const shell::OneBodyOperatorDeltaNSpace space(1, 0, 2, 2, 4);  // J0=1, g0=0, Delta_N_max=2, N1max=2, N2max=4
  // const shell::OneBodyOperatorDeltaNSpace space(2, 0, 2, 2, 4);  // J0=2, g0=0, Delta_N_max=2, N1max=2, N2max=4
 
  // print diagnostics
  std::cout << space.DebugStr();

  std::cout << std::endl;

  // dump subspace contents
  for (std::size_t subspace_index=0; subspace_index<space.size(); ++subspace_index)
    {
      const auto& subspace = space.GetSubspace(subspace_index);
      std::cout << subspace.LabelStr() << std::endl
                << subspace.DebugStr()
                << std::endl;
    }
    
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

void TestOneBodyOperatorDeltaNSectors()
{

  std::cout << "Sectors construction" << std::endl;
  std::cout << std::endl;

  // construct space
  const shell::OneBodyOperatorDeltaNSpace space(0, 0, 2, 2, 4);  // J0=0, g0=0, Delta_N_max=2, N1max=2, N2max=4
  std::cout << space.DebugStr()
            << std::endl;

  // construct sectors
  const shell::OneBodyOperatorDeltaNSectors sectors(space);
  std::cout << sectors.DebugStr()
            << std::endl;

  // // find indices for a matrix element from labels
  // std::cout << "Find indices for a matrix element from labels" << std::endl;
  // std::size_t bra_subspace_index=space.LookUpSubspaceIndex(basis::OneBodyOperatorDeltaNSubspace::SubspaceLabelsType(2));
  // std::size_t ket_subspace_index=space.LookUpSubspaceIndex(basis::OneBodyOperatorDeltaNSubspace::SubspaceLabelsType(2));
  // std::size_t sector_index=hamiltonian_sectors.LookUpSectorIndex(bra_subspace_index,ket_subspace_index);
  // std::cout << "sector index " << sector_index << std::endl;
  // basis::OneBodyOperatorDeltaNSectors::SectorType sector=hamiltonian_sectors.GetSector(sector_index);
  // // const basis::OneBodyOperatorDeltaNSectors::SectorType::BraSubspaceType& bra_subspace=sector.bra_subspace();
  // // const basis::OneBodyOperatorDeltaNSectors::SectorType::KetSubspaceType& ket_subspace=sector.ket_subspace();
  // const basis::OneBodyOperatorDeltaNSubspace& bra_subspace=sector.bra_subspace();
  // const basis::OneBodyOperatorDeltaNSubspace& ket_subspace=sector.ket_subspace();
  // // auto& bra_subspace=sector.bra_subspace();
  // // auto& ket_subspace=sector.ket_subspace();
  // std::cout << "bra_subspace " << bra_subspace.LabelStr() << std::endl;
  // std::cout << "ket_subspace " << ket_subspace.LabelStr() << std::endl;
  // std::cout << "bra_subspace " << std::endl << bra_subspace.DebugStr() << std::endl;
  // std::cout << "ket_subspace " << std::endl << ket_subspace.DebugStr() << std::endl;
  // std::size_t bra_state_index=bra_subspace.LookUpStateIndex(basis::OneBodyOperatorDeltaNState::StateLabelsType(0));
  // std::size_t ket_state_index=ket_subspace.LookUpStateIndex(basis::OneBodyOperatorDeltaNState::StateLabelsType(1));
  // std::cout << "bra state index " << bra_state_index << std::endl;
  // std::cout << "ket state index " << ket_state_index << std::endl;
  // 
  // // construct sectors -- L0=1, g0=1 (E1-like)
  // basis::OneBodyOperatorDeltaNSectors e1_sectors(space, 1, 1);
  // std::cout << "e1_sectors" << std::endl
  //           << e1_sectors.DebugStr()
  //           << std::endl;
  // 
  // // construct sectors -- L0=1, g0=1 (E1-like) -- but including noncanonical ("lower-triangle") sectors
  // basis::OneBodyOperatorDeltaNSectors e1_sectors_both_ways(
  //     space, 1, 1,
  //     basis::SectorDirection::kBoth
  //   );
  // std::cout << "e1_sectors_both_ways" << std::endl
  //           << e1_sectors_both_ways.DebugStr()
  //           << std::endl;
  // 
  // 
  // // construct sectors -- L0=2, g0=0 (E2-like)
  // basis::OneBodyOperatorDeltaNSectors e2_sectors(space, 2, 0);
  // std::cout << "e2_sectors" << std::endl
  //           << e2_sectors.DebugStr()
  //           << std::endl;
  // 
  // // construct sectors -- L0=2, g0=1 (M2-like)
  // basis::OneBodyOperatorDeltaNSectors m2_sectors(space, 2, 1);
  // std::cout << "n2_sectors" << std::endl
  //           << m2_sectors.DebugStr()
  //           << std::endl;

}

void PopulateOperator()
// Populate M matrices with dummy values
{
  std::cout << "Populating operator" << std::endl;
  std::cout << std::endl;

  // set multipolarity
  int J0 = 0;
  int g0 = 0;

  // set truncation
  int Delta_N_max=2;
  int N1max=2;
  int N2max=4;

  // set up data structures
  const shell::OneBodyOperatorDeltaNSpace space(J0, g0, Delta_N_max, N1max, N2max);
  const shell::OneBodyOperatorDeltaNSectors sectors(space);
  basis::OperatorBlocks<double> matrices;
  basis::SetOperatorToZero(sectors, matrices);

  // populate blocks
  for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
  {
    // make aliases for sector and block
    const auto& sector = sectors.GetSector(sector_index);
    auto& sector_matrix = matrices[sector_index];

    // get subspaces
    const auto& bra_subspace = sector.bra_subspace();
    const auto& ket_subspace = sector.ket_subspace();

    // loop over matrix elements in block
    //
    // #pragma omp parallel for collapse(2)
    for (std::size_t bra_index = 0; bra_index < bra_subspace.size(); ++bra_index)
      {
        for (std::size_t ket_index = 0; ket_index < ket_subspace.size(); ++ket_index)
          {
            // get states
            const auto& bra_state = bra_subspace.GetState(bra_index);
            const auto& ket_state = ket_subspace.GetState(ket_index);

            // extract labels
            int bra_n1, bra_l1, bra_n2, bra_l2; HalfInt bra_j1, bra_j2;
            std::tie(bra_n1, bra_l1, bra_j1, bra_n2, bra_l2, bra_j2) = bra_state.labels();
            // TODO: and similarly for ket...
            const int bra_N1 = bra_state.N1();
            const int bra_N2 = bra_state.N2();
            const int ket_N1 = ket_state.N1();
            const int ket_N2 = ket_state.N2();
            
            // calculate matrix element
            float matrix_element = 42.;  // TODO implement matrix element calculations
              
            // save full matrix element
            sector_matrix(bra_index, ket_index) = matrix_element;
          }
      }
    
    // print diagnostic
    std::cout << mcutils::FormatMatrix(sector_matrix, "+.8e") << std::endl
              << std::endl;

  }

}

void PopulateOperator2()
// Populate M matrices with actual computed values
{

  std::cout << "Populating operator (M)" << std::endl;
  std::cout << std::endl;

  // set multipolarity
  int J0 = 2;
  int g0 = 0;

  // set truncation
  int Delta_N_max=2;
  int N1max=2;
  int N2max=4;

  // set number of nucleons (needed for Moshinsky bracket mass ratio 1/(A-1))
  int A = 6;

  // set up data structures
  const shell::OneBodyOperatorDeltaNSpace space(J0, g0, Delta_N_max, N1max, N2max);
  const shell::OneBodyOperatorDeltaNSectors sectors(space);
  basis::OperatorBlocks<double> matrices;

  // populate blocks -- calls ConstructOneBodyOperatorDeltaNMatrix directly,
  // which fills in matrices for all sectors at once
  shell::ConstructOneBodyOperatorDeltaNMatrix(space, sectors, A, matrices);

  // print diagnostic
  for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
    {
      const auto& sector_matrix = matrices[sector_index];
      std::cout << mcutils::FormatMatrix(sector_matrix, "+.8e") << std::endl
                << std::endl;
    }

}

void SerializeOneBodyOperator()
// Testbed code for serializing OBMEs (i.e., packing them into a vector) subject
// to the M matrix indexing scheme.
{

  // set multipolarity
  int J0 = 0;
  int g0 = 0;
  int Tz0 = 0;

  // set up OBO in native format
  int orbital_Nmax = 2;
  basis::OrbitalSpaceLJPN orbital_space(orbital_Nmax);
  basis::OrbitalSectorsLJPN obo_sectors(orbital_space, orbital_space, J0, g0, Tz0);
  basis::OperatorBlocks<double> obo_matrices;
  basis::SetOperatorToZero(obo_sectors, obo_matrices);
  std::cout << "Orbital space" << std::endl
            << orbital_space.DebugStr()
            << std::endl
            << "Native OBO sectors" << std::endl
            << obo_sectors.DebugStr()
            << std::endl;
    
  // set M truncation
  int Delta_N_max=2;
  int N1max=2;
  int N2max=4;
  
  // set up M matrix data structures
  const shell::OneBodyOperatorDeltaNSpace space(J0, g0, Delta_N_max, N1max, N2max);
  const shell::OneBodyOperatorDeltaNSectors sectors(space);
  basis::OperatorBlocks<double> matrices;
  std::cout << "M matrix space" << std::endl
            << space.DebugStr()
            << std::endl
            << "M matrix sectors" << std::endl
            << sectors.DebugStr()
            << std::endl;

// Orbital space
//  index   0 species   0 dim   2 
//  index   1 species   0 dim   1 
//  index   2 species   0 dim   1 
//  index   3 species   0 dim   1 
//  index   4 species   0 dim   1 
//  index   5 species   1 dim   2 
//  index   6 species   1 dim   1 
//  index   7 species   1 dim   1 
//  index   8 species   1 dim   1 
//  index   9 species   1 dim   1 
// 
// Native OBO sectors
//   0 bra   0 (0, 0, 1/2) ket   0 (0, 0, 1/2)
//   1 bra   1 (0, 1, 1/2) ket   1 (0, 1, 1/2)
//   2 bra   2 (0, 1, 3/2) ket   2 (0, 1, 3/2)
//   3 bra   3 (0, 2, 3/2) ket   3 (0, 2, 3/2)
//   4 bra   4 (0, 2, 5/2) ket   4 (0, 2, 5/2)
//   5 bra   5 (1, 0, 1/2) ket   5 (1, 0, 1/2)
//   6 bra   6 (1, 1, 1/2) ket   6 (1, 1, 1/2)
//   7 bra   7 (1, 1, 3/2) ket   7 (1, 1, 3/2)
//   8 bra   8 (1, 2, 3/2) ket   8 (1, 2, 3/2)
//   9 bra   9 (1, 2, 5/2) ket   9 (1, 2, 5/2)
// 
// M matrix space
// J0 0 g0 0 Delta_N_max 2 N1max 2 N2max 4
//   index   0  dim    1  labels [-2]
//   index   1  dim    6  labels [0]
//   index   2  dim    1  labels [2]
// 
// M matrix sectors
//   sector 0  bra index 0 labels [-2] size 1 dim 1  ket index 0 labels [-2] size 1 dim 1  multiplicity index 1  elements 1
//   sector 1  bra index 1 labels [0] size 6 dim 6  ket index 1 labels [0] size 6 dim 6  multiplicity index 1  elements 36
//   sector 2  bra index 2 labels [2] size 1 dim 1  ket index 2 labels [2] size 1 dim 1  multiplicity index 1  elements 1

  // What we need to do...
  //
  // Pick which species we are dealing with.
  //
  // Approach #1: Iterating over target, and retrieving source matrix elements on demand... 
  //
  // Approach #2: Iterating over source, and inserting target matrix elements on demand... 
  //
  // Lookups in the source operator seem more tedious, so maybe take Approach #2?
  //
  // Do we ever expect one of these structures to "overrrun" the other?  Yes,
  // perhaps the M matrix indexing, but we will probably want to truncate the M
  // matrix to match the operator before applying the transformation?  Or maybe
  // that is not necessary?  (The matvec will be cheap compared to the prior
  // matrix inversion.)  So iterate over the smaller structure, which will be
  // the native OBO.

  // select species for transformation (Tz conserving)
  basis::OrbitalSpeciesPN species = basis::OrbitalSpeciesPN::kP;
  
  for (std::size_t obo_sector_index=0; obo_sector_index < obo_sectors.size(); ++obo_sector_index)
    {
    // get sector
    const auto obo_sector = obo_sectors.GetSector(obo_sector_index);
    const auto& obo_bra_subspace = obo_sector.bra_subspace();
    const auto& obo_ket_subspace = obo_sector.ket_subspace();

    // short circuit select for species of interest
    if ((obo_bra_subspace.orbital_species() != species) || (obo_ket_subspace.orbital_species() != species))
      continue;
    
    for (std::size_t obo_bra_index = 0; obo_bra_index < obo_bra_subspace.size(); ++obo_bra_index)
      for (std::size_t obo_ket_index = 0; obo_ket_index < obo_ket_subspace.size(); ++obo_ket_index)
          {
            // unpack orbital info
            const auto& obo_bra_state = obo_bra_subspace.GetState(obo_bra_index);
            const auto& obo_ket_state = obo_ket_subspace.GetState(obo_ket_index);
            int n1 = obo_bra_state.n();
            int l1 = obo_bra_state.l();
            HalfInt j1 = obo_bra_state.j();
            int N1 = obo_bra_state.N();
            int n2 = obo_ket_state.n();
            int l2 = obo_ket_state.l();
            HalfInt j2 = obo_ket_state.j();
            int N2 = obo_ket_state.N();

            // look up target index

            // copy value
          }
      }
    

}


////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  // TestOneBodyOperatorDeltaNSubspace();
  // TestOneBodyOperatorDeltaNSpace();
  // TestOneBodyOperatorDeltaNSectors();
  // PopulateOperator();
  // PopulateOperator2();

  SerializeOneBodyOperator();
  
  // termination
  return EXIT_SUCCESS;
}
