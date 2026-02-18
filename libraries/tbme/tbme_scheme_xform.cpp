/****************************************************************
  tbme_scheme_xform.cpp

  Zhou Zhou
  University of Notre Dame

****************************************************************/

#include "tbme/tbme_scheme_xform.h"

#include "fmt/format.h"

namespace shell {
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  void TransformOperatorTwoBodyJJJTToTwoBodyJJJTTz(
      const basis::TwoBodySpaceJJJT& two_body_jjjt_space,
      const std::array<basis::TwoBodySectorsJJJT,3>& two_body_jjjt_component_sectors,
      const std::array<basis::OperatorBlocks<double>,3>& two_body_jjjt_component_matrices,
      const basis::TwoBodySpaceJJJTTz& two_body_jjjttz_space,
      basis::TwoBodySectorsJJJTTz& two_body_jjjttz_sectors,
      basis::OperatorBlocks<double>& two_body_jjjttz_matrices
    )
  // assumptions: J, T, g, Tz are the same between bra and ket for a given matrix element
  {
    // enforce assumed operator labels
    for (int T0 = 0; T0 <= 2; T0++) {
      if (two_body_jjjt_component_sectors[T0].J0()!=0 || two_body_jjjt_component_sectors[T0].g0()!=0) {
        std::cerr << "ERROR: Provided operator has unsupported (J0,g0)!=(0,0)." << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    
    int J0=0;
    int g0=0;
    int Tz0=0;
    
    // enumerate target sectors
    two_body_jjjttz_sectors = basis::TwoBodySectorsJJJTTz(two_body_jjjttz_space,J0,g0,Tz0);

    // populate matrices
    basis::SetOperatorToZero(two_body_jjjttz_sectors,two_body_jjjttz_matrices);
    for (std::size_t two_body_jjjttz_sector_index=0; two_body_jjjttz_sector_index<two_body_jjjttz_sectors.size(); two_body_jjjttz_sector_index++)
      {
        // make reference to target sector
        const basis::TwoBodySectorsJJJTTz::SectorType& two_body_jjjttz_sector
          = two_body_jjjttz_sectors.GetSector(two_body_jjjttz_sector_index);
        
        // matrix elements between states with different T are assumed to vanish (and do not exist in me2j files)
        if (two_body_jjjttz_sector.bra_subspace().T()!=two_body_jjjttz_sector.ket_subspace().T()) {
          continue;
        }

        int J = two_body_jjjttz_sector.bra_subspace().J();
        int T = two_body_jjjttz_sector.bra_subspace().T();
        int g = two_body_jjjttz_sector.bra_subspace().g();
        int Tz = two_body_jjjttz_sector.bra_subspace().Tz();

        // find subspace indices for corresponding source sector
        std::size_t two_body_jjjt_bra_subspace_index = two_body_jjjt_space.LookUpSubspaceIndex(basis::TwoBodySubspaceJJJT::SubspaceLabelsType(J,T,g));
        std::size_t two_body_jjjt_ket_subspace_index = two_body_jjjt_space.LookUpSubspaceIndex(basis::TwoBodySubspaceJJJT::SubspaceLabelsType(J,T,g));
        // std::cout << "two_body_jjjt_bra_subspace_index" << two_body_jjjt_bra_subspace_index << std::endl;

        // populate matrix elements
        Eigen::MatrixXd& matrix = two_body_jjjttz_matrices[two_body_jjjttz_sector_index];
        for (std::size_t two_body_jjjttz_bra_state_index=0; two_body_jjjttz_bra_state_index<two_body_jjjttz_sector.bra_subspace().size(); two_body_jjjttz_bra_state_index++)
          {
            for (std::size_t two_body_jjjttz_ket_state_index=two_body_jjjttz_bra_state_index; two_body_jjjttz_ket_state_index<two_body_jjjttz_sector.ket_subspace().size(); two_body_jjjttz_ket_state_index++)
              {
                // find bra and ket state labels for jjjttz and pass to jjjt to find bra and ket states indices for jjjt
                // going from target indices to source indices
                const basis::TwoBodySubspaceJJJTTz::StateLabelsType& two_body_jjjttz_bra_state_labels = two_body_jjjttz_sector.bra_subspace().GetStateLabels(two_body_jjjttz_bra_state_index);
                const basis::TwoBodySubspaceJJJTTz::StateLabelsType& two_body_jjjttz_ket_state_labels = two_body_jjjttz_sector.ket_subspace().GetStateLabels(two_body_jjjttz_ket_state_index);
                std::array<double,3> matrix_element_by_T0;
                // find source sector indices from subspace indices
                for (int T0 = 0; T0 <= 2; T0++) {
                  std::size_t two_body_jjjt_sector_index = two_body_jjjt_component_sectors[T0].LookUpSectorIndex(two_body_jjjt_bra_subspace_index,two_body_jjjt_ket_subspace_index);
                  if (two_body_jjjt_sector_index == basis::kNone) {
                    matrix_element_by_T0[T0] = 0;
                  } else {
                    basis::TwoBodySectorsJJJT::SectorType two_body_jjjt_sector = two_body_jjjt_component_sectors[T0].GetSector(two_body_jjjt_sector_index);
                    std::size_t two_body_jjjt_bra_state_index = two_body_jjjt_sector.bra_subspace().LookUpStateIndex(two_body_jjjttz_bra_state_labels);
                    std::size_t two_body_jjjt_ket_state_index = two_body_jjjt_sector.ket_subspace().LookUpStateIndex(two_body_jjjttz_ket_state_labels);
                    // std::cout << "two_body_jjjt_bra/ket_state_index" << two_body_jjjt_bra_state_index << " " << two_body_jjjt_ket_state_index << std::endl;
                    matrix_element_by_T0[T0] = two_body_jjjt_component_matrices[T0][two_body_jjjt_sector_index](two_body_jjjt_bra_state_index,two_body_jjjt_ket_state_index);
                  }
                }
                if (T==0) {
                  matrix(two_body_jjjttz_bra_state_index,two_body_jjjttz_ket_state_index)
                    = matrix_element_by_T0[0];
                } else if (T==1) {
                  if (Tz==1) {
                    matrix(two_body_jjjttz_bra_state_index,two_body_jjjttz_ket_state_index)
                      = matrix_element_by_T0[0] + 1.0/std::sqrt(2.0)*matrix_element_by_T0[1] + 1.0/std::sqrt(10.0)*matrix_element_by_T0[2];
                  } else if (Tz==0) {
                    matrix(two_body_jjjttz_bra_state_index,two_body_jjjttz_ket_state_index)
                      = matrix_element_by_T0[0] - std::sqrt(2.0/5.0)*matrix_element_by_T0[2];
                  } else if (Tz==-1) {
                    matrix(two_body_jjjttz_bra_state_index,two_body_jjjttz_ket_state_index)
                      = matrix_element_by_T0[0] - 1.0/std::sqrt(2.0)*matrix_element_by_T0[1] + 1.0/std::sqrt(10.0)*matrix_element_by_T0[2];
                  }
                }
              }
          }
      }
  }

  void TransformOperatorTwoBodyJJJTTzToTwoBodyJJJT(
      const basis::TwoBodySpaceJJJTTz& two_body_jjjttz_space,
      const basis::TwoBodySectorsJJJTTz& two_body_jjjttz_sectors,
      const basis::OperatorBlocks<double>& two_body_jjjttz_matrices,
      const basis::TwoBodySpaceJJJT& two_body_jjjt_space,
      std::array<basis::TwoBodySectorsJJJT,3>& two_body_jjjt_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& two_body_jjjt_component_matrices
    )
  {

    // enforce assumed operator labels
    if (two_body_jjjttz_sectors.J0()!=0 || two_body_jjjttz_sectors.g0()!=0 || two_body_jjjttz_sectors.Tz0()!=0) {
      std::cerr << "ERROR: Provided operator has unsupported (J0,g0,Tz0)!=(0,0,0)." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    int J0=0;
    int g0=0;
    
    for (int T0 = 0; T0 <= 2; T0++) {
      // enumerate target sectors
      two_body_jjjt_component_sectors[T0] = basis::TwoBodySectorsJJJT(two_body_jjjt_space,J0,T0,g0);
      // populate matrices
      basis::SetOperatorToZero(two_body_jjjt_component_sectors[T0],two_body_jjjt_component_matrices[T0]);
      // two_body_jjjt_component_matrices[T0].resize(two_body_jjjt_component_sectors[T0].size());
      for (std::size_t two_body_jjjt_sector_index=0; two_body_jjjt_sector_index<two_body_jjjt_component_sectors[T0].size(); two_body_jjjt_sector_index++)
        {
          // make reference to target sector
          const basis::TwoBodySectorsJJJT::SectorType& two_body_jjjt_sector
            = two_body_jjjt_component_sectors[T0].GetSector(two_body_jjjt_sector_index);

          if (two_body_jjjt_sector.bra_subspace().T()!=two_body_jjjt_sector.ket_subspace().T()) {
            continue;
          }

          int J = two_body_jjjt_sector.bra_subspace().J();
          int T = two_body_jjjt_sector.bra_subspace().T();
          int g = two_body_jjjt_sector.bra_subspace().g();

          Eigen::MatrixXd& matrix = two_body_jjjt_component_matrices[T0][two_body_jjjt_sector_index];
          for (std::size_t two_body_jjjt_bra_state_index=0; two_body_jjjt_bra_state_index<two_body_jjjt_sector.bra_subspace().size(); two_body_jjjt_bra_state_index++)
            {
              for (std::size_t two_body_jjjt_ket_state_index=two_body_jjjt_bra_state_index; two_body_jjjt_ket_state_index<two_body_jjjt_sector.ket_subspace().size(); two_body_jjjt_ket_state_index++)
                {
                  const basis::TwoBodySubspaceJJJT::StateLabelsType& two_body_jjjt_bra_state_labels = two_body_jjjt_sector.bra_subspace().GetStateLabels(two_body_jjjt_bra_state_index);
                  const basis::TwoBodySubspaceJJJT::StateLabelsType& two_body_jjjt_ket_state_labels = two_body_jjjt_sector.ket_subspace().GetStateLabels(two_body_jjjt_ket_state_index);
                  std::array<double,3> matrix_element_by_Tz; // for Tz = -1, 0, 1
                  for (int Tz = -T; Tz <= T; Tz++) {
                    std::size_t two_body_jjjttz_bra_subspace_index_Tz = two_body_jjjttz_space.LookUpSubspaceIndex(basis::TwoBodySubspaceJJJTTz::SubspaceLabelsType(J,T,g,Tz));
                    std::size_t two_body_jjjttz_sector_index_Tz = two_body_jjjttz_sectors.LookUpSectorIndex(two_body_jjjttz_bra_subspace_index_Tz,two_body_jjjttz_bra_subspace_index_Tz);
                    if (two_body_jjjttz_sector_index_Tz == basis::kNone) {
                      matrix_element_by_Tz[Tz+1] = 0;
                    } else {
                      basis::TwoBodySectorsJJJTTz::SectorType two_body_jjjttz_sector_Tz = two_body_jjjttz_sectors.GetSector(two_body_jjjttz_sector_index_Tz);
                      std::size_t two_body_jjjttz_bra_state_index_Tz = two_body_jjjttz_sector_Tz.bra_subspace().LookUpStateIndex(two_body_jjjt_bra_state_labels);
                      std::size_t two_body_jjjttz_ket_state_index_Tz = two_body_jjjttz_sector_Tz.ket_subspace().LookUpStateIndex(two_body_jjjt_ket_state_labels);
                      matrix_element_by_Tz[Tz+1] = two_body_jjjttz_matrices[two_body_jjjttz_sector_index_Tz](two_body_jjjttz_bra_state_index_Tz,two_body_jjjttz_ket_state_index_Tz);
                    }
                  }
                  if (T==0) {
                    matrix(two_body_jjjt_bra_state_index,two_body_jjjt_ket_state_index)
                      = matrix_element_by_Tz[1];
                  } else if (T==1) {
                    if (T0==0) {
                      matrix(two_body_jjjt_bra_state_index,two_body_jjjt_ket_state_index)
                        = 1/3.0*(matrix_element_by_Tz[0]+matrix_element_by_Tz[1]+matrix_element_by_Tz[2]);
                    } else if (T0==1) {
                      matrix(two_body_jjjt_bra_state_index,two_body_jjjt_ket_state_index)
                        = 1/3.0*(-std::sqrt(9.0/2.0)*matrix_element_by_Tz[0]+std::sqrt(9.0/2.0)*matrix_element_by_Tz[2]);
                    } else if (T0==2) {
                      matrix(two_body_jjjt_bra_state_index,two_body_jjjt_ket_state_index)
                        = 1/3.0*(std::sqrt(5.0/2.0)*matrix_element_by_Tz[0]-std::sqrt(10.0)*matrix_element_by_Tz[1]+std::sqrt(5.0/2.0)*matrix_element_by_Tz[2]);
                    }
                  }
                }
            }
        }
    }

  }

  basis::TwoBodySpeciesPN TwoBodySpeciesPNByTz(int Tz)
  // Look up two body species code
  {
    assert((-1<=Tz) && (Tz<=1));
    switch (Tz)
      {
      case +1: return basis::TwoBodySpeciesPN::kPP;
      case -1: return basis::TwoBodySpeciesPN::kNN;
      case 0: return basis::TwoBodySpeciesPN::kPN;
      }
    std::exit(EXIT_FAILURE);  // to suppress warning that control reaches end of non-void function
  };
        

  std::size_t OscillatorOrbitalIndexFromNj(int N, HalfInt j)
  // Index (0-based) for oscillator orbital (N,j) in standard lexicographic
  // ordering.
  {
    return N*(N+1)/2 + (j.TwiceValue()-1)/2;
  }
  
  void TransformOperatorTwoBodyJJJPNToTwoBodyJJJTTz(
      const basis::TwoBodySpaceJJJPN& two_body_jjjpn_space,
      const basis::TwoBodySectorsJJJPN& two_body_jjjpn_sectors,
      const basis::OperatorBlocks<double>& two_body_jjjpn_matrices,
      const basis::TwoBodySpaceJJJTTz& two_body_jjjttz_space,
      basis::TwoBodySectorsJJJTTz& two_body_jjjttz_sectors,
      basis::OperatorBlocks<double>& two_body_jjjttz_matrices
    )
  {
    // check space truncation
    //
    // TODO debug:
    //
    // terminate called after throwing an instance of 'std::bad_alloc'
    //  what():  std::bad_alloc
    //
    //  if (!(two_body_jjjpn_space.orbital_space().is_oscillator_like()))
    //    {
    //    std::cout << "ERROR: Two-body space not constructed from oscillator-like orbital space."
    //              << std::endl;
    //    std::exit(EXIT_FAILURE);
    //    }
    
    // extract operator labels
    int J0 = two_body_jjjpn_sectors.J0();
    int g0 = two_body_jjjpn_sectors.g0();
    int Tz0  = two_body_jjjpn_sectors.Tz0();

    // initialize target operator
    two_body_jjjttz_sectors = basis::TwoBodySectorsJJJTTz(two_body_jjjttz_space, J0, g0, Tz0);
    basis::SetOperatorToZero(two_body_jjjttz_sectors, two_body_jjjttz_matrices);

    // populate matrices
    for (std::size_t two_body_jjjttz_sector_index=0; two_body_jjjttz_sector_index<two_body_jjjttz_sectors.size(); two_body_jjjttz_sector_index++)
      {
        // make references to target sector
        const basis::TwoBodySectorsJJJTTz::SectorType& two_body_jjjttz_sector
          = two_body_jjjttz_sectors.GetSector(two_body_jjjttz_sector_index);
        //std::cout << two_body_jjjttz_sector.DebugStr() << std::endl;
        Eigen::MatrixXd& two_body_jjjttz_matrix = two_body_jjjttz_matrices[two_body_jjjttz_sector_index];

        // extract target sector labels
        int bra_J, bra_T, bra_g, bra_Tz;
        std::tie(bra_J, bra_T, bra_g, bra_Tz) = two_body_jjjttz_sector.bra_subspace().labels();
        int ket_J, ket_T, ket_g, ket_Tz;
        std::tie(ket_J, ket_T, ket_g, ket_Tz) = two_body_jjjttz_sector.ket_subspace().labels();

        // look up corresponding source subspaces
        basis::TwoBodySpeciesPN bra_two_body_species = basis::kTzCodeTwoBodySpeciesPN.at(bra_Tz);
        basis::TwoBodySpeciesPN ket_two_body_species = basis::kTzCodeTwoBodySpeciesPN.at(ket_Tz);
        std::size_t two_body_jjjpn_bra_subspace_index
          = two_body_jjjpn_space.LookUpSubspaceIndex(basis::TwoBodySubspaceJJJPN::SubspaceLabelsType(bra_two_body_species, bra_J, bra_g));
        const basis::TwoBodySubspaceJJJPN& two_body_jjjpn_bra_subspace = two_body_jjjpn_space.GetSubspace(two_body_jjjpn_bra_subspace_index);
        std::size_t two_body_jjjpn_ket_subspace_index
          = two_body_jjjpn_space.LookUpSubspaceIndex(basis::TwoBodySubspaceJJJPN::SubspaceLabelsType(ket_two_body_species, ket_J, ket_g));
        const basis::TwoBodySubspaceJJJPN& two_body_jjjpn_ket_subspace = two_body_jjjpn_space.GetSubspace(two_body_jjjpn_ket_subspace_index);

        // validate source sector orbitals
        assert(
            two_body_jjjpn_bra_subspace.orbital_subspace1().is_oscillator_like()
            &&
            two_body_jjjpn_bra_subspace.orbital_subspace2().is_oscillator_like()
            &&
            two_body_jjjpn_ket_subspace.orbital_subspace1().is_oscillator_like()
            &&
            two_body_jjjpn_ket_subspace.orbital_subspace2().is_oscillator_like()
          );
        
        // Relations of canonicalized NAS states
        //
        // (Omit delta factors to obtain relations of canonicalized AS states.)
        //
        // In the following expressions, orbitals are ordered a<=b.
        //
        // jjJTTz <-> jjJpn
        //
        // Tz=+1
        //
        //   |ab;J;T,Tz=+1> = |ab;J;pp>
        //
        // Tz=-1
        //
        //   |ab;J;T,Tz=-1> = |ab;J;nn>
        //
        // Tz=0
        //
        //   |ab;J;T={1,0},Tz=0> = (1+delta_ab)^(-1/2) * (1/sqrt(2))
        //      * [ |ab;J;pn> {+,-} (-)(-)^(J-ja-jb) |ba;J;pn>]
        //
        //   Inverting the Tz=0 relation, to give the pn states in terms of TTz
        //   states, we find
        //
        //     {|ab;J;pn> , (-)(-)^(J-ja-jb) |ba;J;pn>} = (1+delta_ab)^(+1/2) * (1/sqrt(2))
        //       * [ |ab;J;T=1,Tz=0> {+,-} |ab;J;T=-,Tz=0>]
        //
        //   or, splitting out these cases and moving the overall sign to the rhs,
        //
        //     |ab;J;pn> = (1+delta_ab)^(+1/2) * (1/sqrt(2)) * [ |ab;J;T=1,Tz=0> + |ab;J;T=-,Tz=0>]
        //
        //     |ba;J;pn> = (1+delta_ab)^(+1/2) * (1/sqrt(2)) * (-)(-)^(J-ja-jb)
        //       * [ |ab;J;T=1,Tz=0> - |ab;J;T=-,Tz=0>]
        //
        //   In the case a=b, this relation reduces to
        //
        //     |aa;J;T=1,Tz=0> = |aa;J;pn>
        //
        //   in conjunction with the constraint for like-orbital states that J
        //   must be even for the pn state to be nonvanishing, and that J+T must
        //   be odd for the T-coupled state to be nonvanishing.
        
        // populate matrix elements
        for (std::size_t two_body_jjjttz_bra_state_index=0; two_body_jjjttz_bra_state_index<two_body_jjjttz_sector.bra_subspace().size(); two_body_jjjttz_bra_state_index++)
          for (std::size_t two_body_jjjttz_ket_state_index=two_body_jjjttz_bra_state_index; two_body_jjjttz_ket_state_index<two_body_jjjttz_sector.ket_subspace().size(); two_body_jjjttz_ket_state_index++)
            {

              // diagonal sector: restrict to upper triangle
              if (two_body_jjjttz_sector.IsDiagonal())
                if (!(two_body_jjjttz_bra_state_index<=two_body_jjjttz_ket_state_index))
                  continue;
              
              // extract target matrix element orbital labels
              int bra_N1; HalfInt bra_j1; int bra_N2; HalfInt bra_j2;
              int ket_N1; HalfInt ket_j1; int ket_N2; HalfInt ket_j2;
              std::tie(bra_N1, bra_j1, bra_N2, bra_j2) = two_body_jjjttz_sector.bra_subspace().GetStateLabels(two_body_jjjttz_bra_state_index);
              std::tie(ket_N1, ket_j1, ket_N2, ket_j2) = two_body_jjjttz_sector.ket_subspace().GetStateLabels(two_body_jjjttz_ket_state_index);
              // std::cout << fmt::format("target rme indices {:2d} {:2d}", two_body_jjjttz_bra_state_index, two_body_jjjttz_ket_state_index) << std::endl;

              // accumulate contributions from source matrix elements
              bool bra_must_swap_orbitals = bra_Tz==0;
              bool ket_must_swap_orbitals = ket_Tz==0;
              double target_rme = 0.;
              for (int bra_swap_orbitals=0; bra_swap_orbitals <= int(bra_must_swap_orbitals); ++bra_swap_orbitals)
                for (int ket_swap_orbitals=0; ket_swap_orbitals <= int(ket_must_swap_orbitals); ++ket_swap_orbitals)
                  {
                    
                    // std::cout << fmt::format("source swap iteration {:2d} {:2d}", bra_swap_orbitals, ket_swap_orbitals) << std::endl;

                    // deduce orbital indices from (N,j) labels
                    //
                    // index = N*(N+1)/2 + (j-1/2)
                    std::size_t bra_index_1 = OscillatorOrbitalIndexFromNj(bra_N1, bra_j1);
                    std::size_t bra_index_2 = OscillatorOrbitalIndexFromNj(bra_N2, bra_j2);
                    std::size_t ket_index_1 = OscillatorOrbitalIndexFromNj(ket_N1, ket_j1);
                    std::size_t ket_index_2 = OscillatorOrbitalIndexFromNj(ket_N2, ket_j2);
                    
                    // look up two-body state indices
                    if (bra_swap_orbitals)
                      std::swap(bra_index_1, bra_index_2);
                    std::size_t two_body_jjjpn_bra_state_index = two_body_jjjpn_bra_subspace.LookUpStateIndex(basis::TwoBodySubspaceJJJPN::StateLabelsType(bra_index_1, bra_index_2));
                    if (ket_swap_orbitals)
                      std::swap(ket_index_1, ket_index_2);
                    std::size_t two_body_jjjpn_ket_state_index = two_body_jjjpn_ket_subspace.LookUpStateIndex(basis::TwoBodySubspaceJJJPN::StateLabelsType(ket_index_1, ket_index_2));

                    // look up canonicalized source matrix element
                    std::size_t two_body_jjjpn_bra_subspace_index_canonical, two_body_jjjpn_ket_subspace_index_canonical,
                      two_body_jjjpn_bra_state_index_canonical, two_body_jjjpn_ket_state_index_canonical;
                    double canonicalization_factor;
                    std::tie(
                        two_body_jjjpn_bra_subspace_index_canonical, two_body_jjjpn_ket_subspace_index_canonical,
                        two_body_jjjpn_bra_state_index_canonical, two_body_jjjpn_ket_state_index_canonical,
                        canonicalization_factor
                      ) =
                      CanonicalizeIndicesJJJPN(
                          two_body_jjjpn_space, J0, g0,
                          two_body_jjjpn_bra_subspace_index, two_body_jjjpn_ket_subspace_index,
                          two_body_jjjpn_bra_state_index, two_body_jjjpn_ket_state_index
                        );
                    std::size_t two_body_jjjpn_sector_index
                      = two_body_jjjpn_sectors.LookUpSectorIndex(two_body_jjjpn_bra_subspace_index_canonical, two_body_jjjpn_ket_subspace_index_canonical);
                    double source_rme = canonicalization_factor * two_body_jjjpn_matrices[two_body_jjjpn_sector_index](two_body_jjjpn_bra_state_index_canonical, two_body_jjjpn_ket_state_index_canonical);

                    // accumulate to target matrix element
                    double bra_prefactor = 1.;
                    if (bra_Tz == 0)
                      {
                        bra_prefactor *= 1/std::sqrt(2); // 1/sqrt(2)
                        if (bra_swap_orbitals && bra_T==0)
                          bra_prefactor *= -1;  // {+,-} 
                        if (bra_swap_orbitals)
                          bra_prefactor *= (-1)*ParitySign(bra_J-bra_j1-bra_j2);  // (-)(-)^(J-ja-jb)
                      }
                    double ket_prefactor = 1.;
                    if (ket_Tz == 0)
                      {
                        ket_prefactor *= 1/std::sqrt(2); // 1/sqrt(2)
                        if (ket_swap_orbitals && ket_T==0)
                          ket_prefactor *= -1;  // {+,-} 
                        if (ket_swap_orbitals)
                          ket_prefactor *= (-1)*ParitySign(ket_J-ket_j1-ket_j2);  // (-)(-)^(J-ja-jb)
                      }
                    
                    target_rme += bra_prefactor * ket_prefactor * source_rme;
                  }
              
              two_body_jjjttz_matrix(two_body_jjjttz_bra_state_index, two_body_jjjttz_ket_state_index) = target_rme;

            }
      }
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
