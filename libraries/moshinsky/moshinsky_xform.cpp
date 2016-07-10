/****************************************************************
  moshinsky_xform.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "moshinsky/moshinsky_xform.h"

#include "moshinsky/moshinsky_bracket.h"

namespace moshinsky {

  ////////////////////////////////////////////////////////////////
  // Obtaining relative-cm matrix elements from relative
  ////////////////////////////////////////////////////////////////

  double RacahReductionFactorFirstSystemGT(
      const HalfInt& J1p, const HalfInt& J2p, const HalfInt& Jp, 
      const HalfInt& J1, const HalfInt& J2, const HalfInt& J, 
      const HalfInt& J0
    )
  {

    assert(J2p==J2);

    //std::cout << "Wigner6J "
    //          << J1p << Jp << J2 << J << J1 << J0
    //          << " " << am::Wigner6J(J1p,Jp,J2,J,J1,J0)
    //          << std::endl;

    double value = ParitySign(J1p+J2+J+J0)
      *Hat(J1p)*Hat(J)
      *am::Wigner6J(J1p,Jp,J2,J,J1,J0);
    return value;
  }

  Eigen::MatrixXd
  RelativeCMMatrixNLSJT(
      const basis::RelativeSpaceLSJT& relative_space,
      const basis::RelativeSectorsLSJT& relative_sectors,
      const basis::MatrixVector& relative_matrices,
      const typename basis::RelativeCMSectorsNLSJT::SectorType& relative_cm_sector,
      int J0, int T0, int g0,
      basis::SymmetryPhaseMode symmetry_phase_mode
    )
  {
  
    // initialize matrix
    Eigen::MatrixXd relative_cm_matrix = Eigen::MatrixXd::Zero(
        relative_cm_sector.bra_subspace().size(),
        relative_cm_sector.ket_subspace().size()
      );

    // extract relative-cm target subspace labels
    // int Np = relative_cm_sector.ket_subspace().N();
    // int Lp = relative_cm_sector.ket_subspace().L();
    // int Sp = relative_cm_sector.ket_subspace().S();
    // int Jp = relative_cm_sector.ket_subspace().J();
    // int Tp = relative_cm_sector.ket_subspace().T();
    // int gp = relative_cm_sector.ket_subspace().g();
    // int N = relative_cm_sector.bra_subspace().N();
    // int L = relative_cm_sector.bra_subspace().L();
    // int S = relative_cm_sector.bra_subspace().S();
    // int J = relative_cm_sector.bra_subspace().J();
    // int T = relative_cm_sector.bra_subspace().T();
    // int g = relative_cm_sector.bra_subspace().g();

    // iterate over spurious relative (target) matrix elements
    for (int relative_cm_bra_index=0; relative_cm_bra_index<relative_cm_sector.bra_subspace().size(); ++relative_cm_bra_index)
      for (int relative_cm_ket_index=0; relative_cm_ket_index<relative_cm_sector.ket_subspace().size(); ++relative_cm_ket_index)
        {

          // DISABLED: short circuit if target matrix element is non-canonical
          //
          // Actually, it is helpful to populate the entire diagonal
          // sector, not just the upper triangle.  We will need to carry
          // out a similarity transformation, and it is easiest to do
          // that as a simple matrix multiplication "sandwich", for
          // which we need the full matrix to be populated.

          // if (relative_cm_sector.IsDiagonal())
          //   if (!(relative_cm_bra_index<=relative_cm_ket_index))
          //     continue;

          // extract needed relative-cm state labels
          basis::RelativeCMStateNLSJT relative_cm_bra(relative_cm_sector.bra_subspace(),relative_cm_bra_index);
          int Nrp = relative_cm_bra.Nr();
          int lrp = relative_cm_bra.lr();
          int Ncp = relative_cm_bra.Nc();
          int lcp = relative_cm_bra.lc();
          int Lp = relative_cm_bra.L();
          int Sp = relative_cm_bra.S();
          int Jp = relative_cm_bra.J();
          int Tp = relative_cm_bra.T();

          basis::RelativeCMStateNLSJT relative_cm_ket(relative_cm_sector.ket_subspace(),relative_cm_ket_index);
          int Nr = relative_cm_ket.Nr();
          int lr = relative_cm_ket.lr();
          int Nc = relative_cm_ket.Nc();
          int lc = relative_cm_ket.lc();
          int L = relative_cm_ket.L();
          int S = relative_cm_ket.S();
          int J = relative_cm_ket.J();
          int T = relative_cm_ket.T();

          // short circuit if cm motion not conserved
          if (!((Ncp==Nc)&&(lcp==lc)))
            continue;

          // check for adequate relative N truncation on source interaction
          assert(std::max(Nrp,Nr)<=relative_space.Nmax());

          // determine relative (source) angular momenta
          //   by triangle selection between relative and cm motion
          HalfInt::pair relative_bra_Jr_range = am::AngularMomentumRangeIntersection(
              am::ProductAngularMomentumRange(relative_cm_bra.lc(),relative_cm_bra.J()),
              am::ProductAngularMomentumRange(relative_cm_bra.lr(),relative_cm_bra.S())
            );
          HalfInt::pair relative_ket_Jr_range = am::AngularMomentumRangeIntersection(
              am::ProductAngularMomentumRange(relative_cm_ket.lc(),relative_cm_ket.J()),
              am::ProductAngularMomentumRange(relative_cm_ket.lr(),relative_cm_ket.S())
            );

          // sum relative angular momentum contributions
          double relative_cm_matrix_element = 0;
          for (int Jrp=int(relative_bra_Jr_range.first); Jrp<=int(relative_bra_Jr_range.second); ++Jrp)
            for (int Jr=int(relative_ket_Jr_range.first); Jr<=int(relative_ket_Jr_range.second); ++Jr)
              {

                // short circuit if relative J truncation on source
                // interaction eliminates this source sector
                if (!(std::max(Jrp,Jr)<=relative_space.Jmax()))
                  continue;

                // short circuit if relative J values are not triangular
                //
                // Note that coservation of CM motion guarantees that
                // the relative parity condition will be satisfied
                // (grp+g0+gr~0), given that the relative_cm parity
                // condition is presumed satisfied as a precondition.
                if (!(am::AllowedTriangle(Jrp,J0,Jr)))
                  continue;

                // look up relative subspace indices
                int grp = lrp%2;
                int relative_bra_subspace_index = relative_space.LookUpSubspaceIndex(
                    basis::RelativeSubspaceLSJT::SubspaceLabelsType(lrp,Sp,Jrp,Tp,grp)
                  );
                int gr = lr%2;
                int relative_ket_subspace_index = relative_space.LookUpSubspaceIndex(
                    basis::RelativeSubspaceLSJT::SubspaceLabelsType(lr,S,Jr,T,gr)
                  );

                // look up relative matrix element indices
                const basis::RelativeSubspaceLSJT& relative_bra_subspace = relative_space.GetSubspace(
                    relative_bra_subspace_index
                  );
                const basis::RelativeSubspaceLSJT& relative_ket_subspace = relative_space.GetSubspace(
                    relative_ket_subspace_index
                  );
                int relative_bra_index = relative_bra_subspace.LookUpStateIndex(
                    basis::RelativeStateLSJT::StateLabelsType(Nrp)
                  );
                int relative_ket_index = relative_ket_subspace.LookUpStateIndex(
                    basis::RelativeStateLSJT::StateLabelsType(Nr)
                  );
                // std::cout << relative_bra_subspace.LabelStr() << relative_ket_subspace.LabelStr() << std::endl;


                // canonicalize indices for matrix element lookup
                double canonicalization_factor;
                basis::CanonicalizeIndicesRelativeLSJT(
                    relative_space,
                    relative_bra_subspace_index, relative_ket_subspace_index,
                    relative_bra_index, relative_ket_index,
                    canonicalization_factor,
                    J0, T0, g0,
                    symmetry_phase_mode
                  );

                // look up matrix element
                int relative_sector_index = relative_sectors.LookUpSectorIndex(
                    relative_bra_subspace_index,
                    relative_ket_subspace_index
                  );
                double relative_matrix_element = canonicalization_factor
                  * relative_matrices[relative_sector_index](relative_bra_index,relative_ket_index);
              
                // accumulate contribution to relative-cm matrix element
                double contribution 
                  = am::Unitary6JZ(lrp,lc,Lp,Jp,Sp,Jrp) * am::Unitary6JZ(lr,lc,L,J,S,Jr)
                  * RacahReductionFactorFirstSystemGT(Jrp,lc,Jp,Jr,lc,J,J0)
                  * relative_matrix_element;
                // std::cout << " Jrp " << Jrp << " Jr " << Jr
                //           << " " << relative_bra_subspace.LabelStr()
                //           << " " << relative_ket_subspace.LabelStr()
                //           << " " << am::Unitary6JZ(lrp,lc,Lp,Jp,Sp,Jrp)
                //           << " " <<  am::Unitary6JZ(lr,lc,L,J,S,Jr)
                //           << " " << RacahReductionFactorFirstSystemGT(Jrp,lc,Jp,Jr,lc,J,J0)
                //           << " " << relative_matrix_element
                //           << " contribution " << contribution
                //           << std::endl;
                relative_cm_matrix_element += contribution;

              }	

          // save matrix element
          relative_cm_matrix(relative_cm_bra_index,relative_cm_ket_index) = relative_cm_matrix_element;
        }

    // return matrix
    return relative_cm_matrix;

  }

  void TransformOperatorRelativeLSJTToRelativeCMNLSJT(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      const std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      const std::array<basis::MatrixVector,3>& relative_component_matrices,
      const basis::RelativeCMSpaceNLSJT& relative_cm_nlsjt_space,
      std::array<basis::RelativeCMSectorsNLSJT,3>& relative_cm_nlsjt_component_sectors,
      std::array<basis::MatrixVector,3>& relative_cm_nlsjt_component_matrices
    )
  {
    for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
      // for each isospin component
      {

        // enumerate sectors
        relative_cm_nlsjt_component_sectors[T0]
          = basis::RelativeCMSectorsNLSJT(relative_cm_nlsjt_space,operator_labels.J0,T0,operator_labels.g0);

        // populate matrices
        relative_cm_nlsjt_component_matrices[T0].resize(relative_cm_nlsjt_component_sectors[T0].size());
        for (int sector_index=0; sector_index<relative_cm_nlsjt_component_sectors[T0].size(); ++sector_index)
          {
            const basis::RelativeCMSectorsNLSJT::SectorType& relative_cm_nlsjt_sector
              = relative_cm_nlsjt_component_sectors[T0].GetSector(sector_index);
            relative_cm_nlsjt_component_matrices[T0][sector_index] = moshinsky::RelativeCMMatrixNLSJT(
                relative_space,relative_component_sectors[T0],relative_component_matrices[T0],
                relative_cm_nlsjt_sector,
                operator_labels.J0, T0, operator_labels.g0,
                operator_labels.symmetry_phase_mode
              );
          }
      }
  }

  ////////////////////////////////////////////////////////////////
  // Moshinsky transformation
  ////////////////////////////////////////////////////////////////

  Eigen::MatrixXd
  TransformationMatrixRelativeCMTwoBodyNLSJT(
      const basis::RelativeCMSubspaceNLSJT& relative_cm_subspace,
      const basis::TwoBodySubspaceNLSJT& two_body_subspace
    )
  {

    // set up matrix to hold results
    Eigen::MatrixXd matrix(relative_cm_subspace.size(),two_body_subspace.size());

    // guard against call with mismatched subspaces
    assert(relative_cm_subspace.N()==two_body_subspace.N());
    assert(relative_cm_subspace.L()==two_body_subspace.L());
    assert(relative_cm_subspace.S()==two_body_subspace.S());
    assert(relative_cm_subspace.J()==two_body_subspace.J());
    assert(relative_cm_subspace.T()==two_body_subspace.T());

    // retrieve L for sector
    int L = relative_cm_subspace.L();

    // populate matrix
    for (int index_relative_cm=0; index_relative_cm<relative_cm_subspace.size(); ++index_relative_cm)
      for (int index_two_body=0; index_two_body<two_body_subspace.size(); ++index_two_body)
        {

          // retrieve relative-cm state
          basis::RelativeCMStateNLSJT relative_cm_state(relative_cm_subspace,index_relative_cm);
          int Nr = relative_cm_state.Nr();
          int lr = relative_cm_state.lr();
          int nr = (Nr - lr)/2;
          int Nc = relative_cm_state.Nc();
          int lc = relative_cm_state.lc();
          int nc = (Nc - lc)/2;
          int L = relative_cm_state.L();

          // retrieve two-body state
          basis::TwoBodyStateNLSJT two_body_state(two_body_subspace,index_two_body);
          int N1 = two_body_state.N1();
          int l1 = two_body_state.l1();
          int n1 = (N1 - l1)/2;
          int N2 = two_body_state.N2();
          int l2 = two_body_state.l2();
          int n2 = (N2 - l2)/2;

          // evaluate bracket
          matrix(index_relative_cm,index_two_body) = sqrt(2.)*moshinsky::MoshinskyBracket(nr,lr,nc,lc,n1,l1,n2,l2,L);
				  
        }

    // return matrix
    //  std::cout << "Moshinsky " << matrix << std::endl;
    return matrix;
  }

  Eigen::MatrixXd
  TwoBodyMatrixNLSJT(
      const basis::RelativeCMSectorsNLSJT::SectorType& relative_cm_sector,
      const basis::TwoBodySectorsNLSJT::SectorType& two_body_sector,
      const Eigen::MatrixXd& relative_cm_matrix
    )
  {

    // obtain transformation matrices
    Eigen::MatrixXd bra_transformation_matrix = TransformationMatrixRelativeCMTwoBodyNLSJT(
        relative_cm_sector.bra_subspace(),
        two_body_sector.bra_subspace()
      );
    Eigen::MatrixXd ket_transformation_matrix = TransformationMatrixRelativeCMTwoBodyNLSJT(
        relative_cm_sector.ket_subspace(),
        two_body_sector.ket_subspace()
      );

    // obtain target matrix
    Eigen::MatrixXd two_body_matrix
      = bra_transformation_matrix.transpose()
      * relative_cm_matrix
      * ket_transformation_matrix;

    return two_body_matrix;

  }

  void TransformOperatorRelativeCMNLSJTToTwoBodyNLSJT(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeCMSpaceNLSJT& relative_cm_nlsjt_space,
      const std::array<basis::RelativeCMSectorsNLSJT,3>& relative_cm_nlsjt_component_sectors,
      const std::array<basis::MatrixVector,3>& relative_cm_nlsjt_component_matrices,
      const basis::TwoBodySpaceNLSJT& two_body_nlsjt_space,
      std::array<basis::TwoBodySectorsNLSJT,3>& two_body_nlsjt_component_sectors,
      std::array<basis::MatrixVector,3>& two_body_nlsjt_component_matrices
    )
  {
    for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
      // for each isospin component
      {

        // enumerate sectors
        two_body_nlsjt_component_sectors[T0]
          = basis::TwoBodySectorsNLSJT(two_body_nlsjt_space,operator_labels.J0,T0,operator_labels.g0);
         
        // std::cout << " T0 " << T0
        //           << " sectors " << two_body_nlsjt_component_sectors[T0].size()
        //           << std::endl;

        // populate matrices
        two_body_nlsjt_component_matrices[T0].resize(two_body_nlsjt_component_sectors[T0].size());
        for (int sector_index=0; sector_index<two_body_nlsjt_component_sectors[T0].size(); ++sector_index)
          {
            // target sector
            const basis::TwoBodySectorsNLSJT::SectorType& two_body_nlsjt_sector
              = two_body_nlsjt_component_sectors[T0].GetSector(sector_index);

            // look up source sector
            int relative_cm_nlsjt_bra_subspace_index = relative_cm_nlsjt_space.LookUpSubspaceIndex(
                basis::RelativeCMSubspaceNLSJTLabels(
                    two_body_nlsjt_sector.bra_subspace().L(),
                    two_body_nlsjt_sector.bra_subspace().S(),
                    two_body_nlsjt_sector.bra_subspace().J(),
                    two_body_nlsjt_sector.bra_subspace().T(),
                    two_body_nlsjt_sector.bra_subspace().g(),
                    two_body_nlsjt_sector.bra_subspace().N()
                  )
              );
            int relative_cm_nlsjt_ket_subspace_index = relative_cm_nlsjt_space.LookUpSubspaceIndex(
                basis::RelativeCMSubspaceNLSJTLabels(
                    two_body_nlsjt_sector.ket_subspace().L(),
                    two_body_nlsjt_sector.ket_subspace().S(),
                    two_body_nlsjt_sector.ket_subspace().J(),
                    two_body_nlsjt_sector.ket_subspace().T(),
                    two_body_nlsjt_sector.ket_subspace().g(),
                    two_body_nlsjt_sector.ket_subspace().N()
                  )
              );
            int relative_cm_nlsjt_sector_index = relative_cm_nlsjt_component_sectors[T0].LookUpSectorIndex(
                relative_cm_nlsjt_bra_subspace_index,
                relative_cm_nlsjt_ket_subspace_index
              );
            const basis::RelativeCMSectorsNLSJT::SectorType& relative_cm_nlsjt_sector
              = relative_cm_nlsjt_component_sectors[T0].GetSector(relative_cm_nlsjt_sector_index);
            const Eigen::MatrixXd& relative_cm_nlsjt_matrix
              = relative_cm_nlsjt_component_matrices[T0][relative_cm_nlsjt_sector_index];

            // transform
            two_body_nlsjt_component_matrices[T0][sector_index] = moshinsky::TwoBodyMatrixNLSJT(
                relative_cm_nlsjt_sector,
                two_body_nlsjt_sector,
                relative_cm_nlsjt_matrix
              );
          }
      }

  }


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
