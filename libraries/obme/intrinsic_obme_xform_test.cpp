/****************************************************************
  intrinsic_obme_xform_test.cpp

  Mark A. Caprio and Victor Dumenil
  University of Notre Dame

****************************************************************/

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <vector>

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
  int J0 = 1;
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
// tests of M^K, its inverse, and translationally invariant OBDMEs
////////////////////////////////////////////////////////////////

int DeltaNMaxForParity(int N1max, int g0)
// Largest Delta_N <= N1max compatible with the parity grade g0
// (Delta_N must have the same parity as g0).
{
  return ((N1max-g0)%2==0) ? N1max : N1max-1;
}

void TestOneBodyOperatorDeltaNMatrixInverse()
// Quick checks on M and M^{-1}, for a range of (J0,g0).
//
//  (1) Inverse:  max|M*Minv-1| and max|Minv*M-1| should be ~1e-12 or better
//      (checks inversion only).
//
//  (2) Structure of M (checks M itself): by the oscillator quanta conservation
//      in the brackets <n l 0 0|N1 L1 n1 l1> of (13), the Jacobi-type labels
//      (n,l,n',l') must carry at least as many quanta as the orbital labels:
//      N_orb <= N_Jac for each of the two orbitals.  We count nonzero
//      elements with N_orb<N_Jac ("expected") and with N_orb>N_Jac
//      ("reversed").  Expect reversed=0.
//
//  (3) Limit A->infinity (d=1/(A-1)->0): the c.m. becomes infinitely heavy, so
//      xi_{A-1} -> -r_A and M must tend to +-identity (a uniform sign per
//      block).  Checks the hat factors, 6-j symbols and phase in (13).
{

  std::cout << "\n Test of M and M^-1" << std::endl << std::endl;

  const int N1max = 8;
  const int N2max = 16;

  for (int J0=0; J0<=2; ++J0)
    for (int g0=0; g0<=1; ++g0)
      {
        int Delta_N_max = DeltaNMaxForParity(N1max,g0);
        const shell::OneBodyOperatorDeltaNSpace space(J0,g0,Delta_N_max,N1max,N2max);
        const shell::OneBodyOperatorDeltaNSectors sectors(space);

        // (1) inverse, at A=19 (finite-A brackets)
        int A = 19;
        basis::OperatorBlocks<double> matrices, inverse_matrices;
        shell::ConstructOneBodyOperatorDeltaNMatrix(space,sectors,A,matrices);
        double cond = 0.;
        shell::InvertOneBodyOperatorDeltaNMatrix(matrices,inverse_matrices,&cond);

        double err_right = 0., err_left = 0.;
        for (std::size_t i=0; i<matrices.size(); ++i)
          {
            if (matrices[i].size()==0) continue;
            Eigen::MatrixXd id = Eigen::MatrixXd::Identity(matrices[i].rows(),matrices[i].cols());
            err_right = std::max(err_right,(matrices[i]*inverse_matrices[i]-id).cwiseAbs().maxCoeff());
            err_left  = std::max(err_left, (inverse_matrices[i]*matrices[i]-id).cwiseAbs().maxCoeff());
          }

        // (2) structure of M
        long n_expected = 0, n_reversed = 0;
        for (std::size_t i=0; i<matrices.size(); ++i)
          {
            const auto& sector = sectors.GetSector(i);
            for (int row=0; row<matrices[i].rows(); ++row)
              for (int col=0; col<matrices[i].cols(); ++col)
                {
                  if (std::abs(matrices[i](row,col))<1e-12) continue;
                  const shell::OneBodyOperatorDeltaNState orb(sector.bra_subspace(),row);
                  const shell::OneBodyOperatorDeltaNState jac(sector.ket_subspace(),col);
                  if ((orb.N1()<jac.N1())||(orb.N2()<jac.N2())) ++n_expected;
                  if ((orb.N1()>jac.N1())||(orb.N2()>jac.N2())) ++n_reversed;
                }
          }

        // (3) A -> infinity limit
        basis::OperatorBlocks<double> matrices_inf;
        shell::ConstructOneBodyOperatorDeltaNMatrix(space,sectors,10000000,matrices_inf);
        double err_inf = 0.;
        bool uniform_sign = true;
        for (std::size_t i=0; i<matrices_inf.size(); ++i)
          {
            if (matrices_inf[i].size()==0) continue;
            double sign = (matrices_inf[i](0,0)>=0.) ? 1. : -1.;
            Eigen::MatrixXd id = Eigen::MatrixXd::Identity(matrices_inf[i].rows(),matrices_inf[i].cols());
            err_inf = std::max(err_inf,(matrices_inf[i]-sign*id).cwiseAbs().maxCoeff());
          }

        std::cout << "J0 " << J0 << " g0 " << g0
                  << " | max|M Minv-1| " << err_right
                  << "  max|Minv M-1| " << err_left
                  << "  cond " << cond
                  << " | nonzero N_orb<N_Jac " << n_expected
                  << "  N_orb>N_Jac " << n_reversed << " (expect 0)"
                  << " | A->inf: max|M-(+-1)| " << err_inf << " (expect ~1e-5 or less)"
                  << std::endl;
      }
  std::cout << std::endl;
}

////////////////////////////////////////////////////////////////
// reading density files (trdens_* and obd_trinv_*)
////////////////////////////////////////////////////////////////

struct DensityOrbital { int n; int l; HalfInt j; };
struct DensityEntry { int a; int b; double value[2]; };  // a,b: 1-based orbital indices as in file
struct DensityBlock { int J0; std::vector<DensityEntry> entries; };
struct DensityFile
{
  int A = 0;
  int N1max = -1;
  int N12max = -1;
  std::vector<DensityOrbital> orbitals;
  std::vector<DensityBlock> blocks;
};

DensityFile ReadDensityFile(const std::string& filename)
// Parse header (A, N1_max, N12_max, orbital list) and the "Jtrans=" blocks of
// lines "a b value_1 value_2".
{
  DensityFile f;
  std::ifstream stream(filename);
  if (!stream)
    {
      std::cerr << "ERROR: cannot open " << filename << std::endl;
      std::exit(EXIT_FAILURE);
    }

  std::string line;
  bool in_block = false;
  while (std::getline(stream,line))
    {
      int i1,i2,i3,i4,i5;
      double v1,v2;
      const char* s = line.c_str();

      if ((f.A==0) && (std::sscanf(s," A= %d",&i1)==1))
        {
          f.A = i1;
        }
      else if ((f.N1max<0) && (std::sscanf(s," N1_max= %d N12_max= %d",&i1,&i2)==2))
        {
          f.N1max = i1;
          f.N12max = i2;
        }
      else if (std::sscanf(s," # %d n= %d l= %d j= %d/%d",&i1,&i2,&i3,&i4,&i5)==5)
        {
          if (i1!=int(f.orbitals.size())+1)
            {
              std::cerr << "ERROR: orbital list out of order" << std::endl;
              std::exit(EXIT_FAILURE);
            }
          f.orbitals.push_back({i2,i3,HalfInt(i4,i5)});
        }
      else if (std::sscanf(s," Jtrans= %d",&i1)==1)
        {
          DensityBlock block;
          block.J0 = i1;
          f.blocks.push_back(block);
          in_block = true;
        }
      else if (in_block && (std::sscanf(s,"%d %d %lf %lf",&i1,&i2,&v1,&v2)==4))
        {
          DensityEntry e;
          e.a = i1; e.b = i2; e.value[0] = v1; e.value[1] = v2;
          f.blocks.back().entries.push_back(e);
        }
    }

  if ((f.A==0) || f.orbitals.empty() || f.blocks.empty())
    {
      std::cerr << "ERROR: could not parse " << filename << std::endl;
      std::exit(EXIT_FAILURE);
    }
  return f;
}

bool LocateOrbitalPair(
    const shell::OneBodyOperatorDeltaNSpace& space,
    const DensityOrbital& o1, const DensityOrbital& o2,
    std::size_t& subspace_index, std::size_t& state_index
  )
// Find the (subspace,state) of the M indexing for orbital pair (o1,o2), i.e.,
// labels (n1,l1,j1,n2,l2,j2).  Returns false if the pair is not in the space.
{
  int Delta_N = (2*o1.n+o1.l) - (2*o2.n+o2.l);
  for (std::size_t si=0; si<space.size(); ++si)
    {
      const auto& subspace = space.GetSubspace(si);
      if (subspace.Delta_N()!=Delta_N)
        continue;
      for (std::size_t k=0; k<subspace.size(); ++k)
        {
          const shell::OneBodyOperatorDeltaNState st(subspace,k);
          if ((st.n1()==o1.n)&&(st.l1()==o1.l)&&(st.j1()==o1.j)
              &&(st.n2()==o2.n)&&(st.l2()==o2.l)&&(st.j2()==o2.j))
            {
              subspace_index = si;
              state_index = k;
              return true;
            }
        }
    }
  return false;
}

typedef std::map<std::tuple<int,int,int>,std::array<double,2>> DensityMap;  // (J0,a,b) -> values

DensityMap TransformDensity(const DensityFile& f, bool swap_indices, int extra_N)
// Apply (14): rho_intrinsic = M^{-1} rho, block by block in (J0,g0,Delta_N).
//
// swap_indices: interpret file pair (a,b) as (n2,n1) rather than (n1,n2)
// extra_N: enlarge the M truncation beyond the density's (N1_max,N12_max)
{
  int Nmax_op = f.N1max;
  if (Nmax_op<0)
    for (const auto& o : f.orbitals) Nmax_op = std::max(Nmax_op,2*o.n+o.l);
  int Ntot_max = (f.N12max>=0) ? f.N12max : 2*Nmax_op;
  int N1max_mat = Nmax_op + extra_N;
  int N2max_mat = Ntot_max + 2*extra_N;

  DensityMap result;

  for (const auto& block : f.blocks)
    for (int g0=0; g0<=1; ++g0)
      {
        // entries of this parity
        std::vector<const DensityEntry*> sel;
        for (const auto& e : block.entries)
          {
            const auto& oa = f.orbitals.at(e.a-1);
            const auto& ob = f.orbitals.at(e.b-1);
            if ((oa.l+ob.l)%2==g0) sel.push_back(&e);
          }
        if (sel.empty()) continue;

        int Delta_N_max = DeltaNMaxForParity(N1max_mat,g0);
        const shell::OneBodyOperatorDeltaNSpace space(block.J0,g0,Delta_N_max,N1max_mat,N2max_mat);
        const shell::OneBodyOperatorDeltaNSectors sectors(space);
        basis::OperatorBlocks<double> matrices, inverse_matrices;
        shell::ConstructOneBodyOperatorDeltaNMatrix(space,sectors,f.A,matrices);
        shell::InvertOneBodyOperatorDeltaNMatrix(matrices,inverse_matrices);

        for (int col=0; col<2; ++col)
          {
            // source vectors, one per subspace (= sector)
            std::vector<Eigen::VectorXd> source(space.size());
            for (std::size_t si=0; si<space.size(); ++si)
              source[si] = Eigen::VectorXd::Zero(space.GetSubspace(si).size());

            std::vector<std::array<std::size_t,2>> location(sel.size());
            std::vector<bool> found(sel.size(),false);

            for (std::size_t k=0; k<sel.size(); ++k)
              {
                const auto& oa = f.orbitals.at(sel[k]->a-1);
                const auto& ob = f.orbitals.at(sel[k]->b-1);
                std::size_t si, st;
                bool ok = swap_indices ? LocateOrbitalPair(space,ob,oa,si,st)
                                       : LocateOrbitalPair(space,oa,ob,si,st);
                if (!ok)
                  {
                    if (col==0)
                      std::cerr << "warning: pair (" << sel[k]->a << "," << sel[k]->b
                                << ") J0=" << block.J0 << " not in M space; skipped" << std::endl;
                    continue;
                  }
                found[k] = true;
                location[k] = {si,st};
                source[si](st) = sel[k]->value[col];
              }

            // apply M^{-1} on each subspace
            std::vector<Eigen::VectorXd> target(space.size());
            for (std::size_t si=0; si<space.size(); ++si)
              target[si] = inverse_matrices[si]*source[si];

            for (std::size_t k=0; k<sel.size(); ++k)
              {
                if (!found[k]) continue;
                auto key = std::make_tuple(block.J0,sel[k]->a,sel[k]->b);
                result[key][col] = target[location[k][0]](location[k][1]);
              }
          }
      }

  return result;
}

bool TestTranslationallyInvariantOBDME(
    const std::string& trdens_filename,
    const std::string& obd_trinv_filename,
    double tolerance = 1e-6
  )
// Transform the OBDMEs of a "trdens" file (c.m. not removed) with M^{-1}, as in
// (14), and compare with the "obd_trinv" file (translationally invariant
// OBDMEs from the established code).
//
// Both possible readings of the pair indices (a,b) in the file are tried
// ((a,b)=(n1,n2) or (n2,n1)); the one with smaller discrepancy is shown.
{

  std::cout << "Translationally invariant OBDMEs: M^-1 * trdens vs. obd_trinv" << std::endl
            << "  trdens   " << trdens_filename << std::endl
            << "  obd_trinv " << obd_trinv_filename << std::endl << std::endl;

  const DensityFile raw = ReadDensityFile(trdens_filename);
  const DensityFile ref = ReadDensityFile(obd_trinv_filename);
  std::cout << "A=" << raw.A << "  N1_max=" << raw.N1max << "  N12_max=" << raw.N12max
            << "  orbitals " << raw.orbitals.size() << std::endl << std::endl;

  // reference and raw values by key
  DensityMap ref_map, raw_map;
  for (const auto& b : ref.blocks)
    for (const auto& e : b.entries)
      ref_map[std::make_tuple(b.J0,e.a,e.b)] = {e.value[0],e.value[1]};
  for (const auto& b : raw.blocks)
    for (const auto& e : b.entries)
      raw_map[std::make_tuple(b.J0,e.a,e.b)] = {e.value[0],e.value[1]};

  // baseline: no transformation
  double dev_raw = 0.;
  for (const auto& kv : raw_map)
    if (ref_map.count(kv.first))
      for (int c=0; c<2; ++c)
        dev_raw = std::max(dev_raw,std::abs(kv.second[c]-ref_map[kv.first][c]));
  std::cout << "max|trdens - obd_trinv| (no transformation)  = " << dev_raw << std::endl;

  // try both index conventions
  double best_dev = 1e300;
  DensityMap best;
  bool best_swap = false;
  for (int swap=0; swap<=1; ++swap)
    {
      // with this configuration (raw,swap==1,6), we succes to reproduce the same values as in 
      // obd_trinv at Nmax=6 (and below). 
      // While this configuration (raw,swap==1,0) works only for Nmax=0. 
      //DensityMap t = TransformDensity(raw,swap==1,0);
      DensityMap t = TransformDensity(raw,swap==1,6);
      double dev = 0.;
      for (const auto& kv : t)
        if (ref_map.count(kv.first))
          for (int c=0; c<2; ++c)
            dev = std::max(dev,std::abs(kv.second[c]-ref_map[kv.first][c]));
      std::cout << "max|M^-1 trdens - obd_trinv|, (a,b)=" << (swap ? "(n2,n1)" : "(n1,n2)")
                << "  = " << dev << std::endl;
      if (dev<best_dev) { best_dev = dev; best = t; best_swap = (swap==1); }
    }
  std::cout << "-> best convention: (a,b)=" << (best_swap ? "(n2,n1)" : "(n1,n2)") << std::endl << std::endl;

  // detailed table
  std::cout << std::setw(3) << "J0" << std::setw(4) << "a" << std::setw(4) << "b" << std::setw(4) << "col"
            << std::setw(16) << "trdens" << std::setw(16) << "M^-1 trdens"
            << std::setw(16) << "obd_trinv" << std::setw(12) << "diff" << std::endl;
  for (const auto& kv : best)
    {
      if (!ref_map.count(kv.first)) continue;
      for (int c=0; c<2; ++c)
        {
          double r = raw_map[kv.first][c], t = kv.second[c], x = ref_map[kv.first][c];
          if ((r==0.)&&(t==0.)&&(x==0.)) continue;  // skip all-zero rows (e.g. isoscalar column)
          std::cout << std::setw(3) << std::get<0>(kv.first)
                    << std::setw(4) << std::get<1>(kv.first)
                    << std::setw(4) << std::get<2>(kv.first)
                    << std::setw(4) << c
                    << std::scientific << std::setprecision(8)
                    << std::setw(16) << r << std::setw(16) << t << std::setw(16) << x
                    << std::setw(12) << std::setprecision(2) << (t-x)
                    << std::defaultfloat << std::endl;
        }
    }

  bool pass = (best_dev<tolerance);
  std::cout << std::endl << (pass ? "PASS" : "FAIL")
            << ": max deviation " << best_dev << " (tolerance " << tolerance << ")" << std::endl;
  return pass;
}


////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
  // Usage:
  //   intrinsic_obme_xform_test
  //       -> checks on M and M^-1 (no input files needed)
  //   intrinsic_obme_xform_test <trdens_file> <obd_trinv_file>
  //       -> additionally, M^-1 * trdens vs. obd_trinv comparison

{

  // TestOneBodyOperatorDeltaNSubspace();
  // TestOneBodyOperatorDeltaNSpace();
  // TestOneBodyOperatorDeltaNSectors();
  // PopulateOperator();
  //PopulateOperator2();

  //SerializeOneBodyOperator();
  
  TestOneBodyOperatorDeltaNMatrixInverse();

  if (argc>=3)
    {
      bool ok = TestTranslationallyInvariantOBDME(argv[1],argv[2]);
      return ok ? EXIT_SUCCESS : EXIT_FAILURE;
    }

  // termination
  return EXIT_SUCCESS;
}
