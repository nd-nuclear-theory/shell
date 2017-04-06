/****************************************************************
  moshinsky_xform.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

// #include "cppformat/format.h"  // debugging
#include "moshinsky/moshinsky_bracket.h"
#include "moshinsky/moshinsky_xform.h"

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
  RelativeCMMatrixLSJTN(
      const basis::RelativeSpaceLSJT& relative_space,
      const basis::RelativeSectorsLSJT& relative_sectors,
      const basis::MatrixVector& relative_matrices,
      const typename basis::RelativeCMSectorsLSJTN::SectorType& relative_cm_sector,
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
          basis::RelativeCMStateLSJTN relative_cm_bra(relative_cm_sector.bra_subspace(),relative_cm_bra_index);
          int Nrp = relative_cm_bra.Nr();
          int lrp = relative_cm_bra.lr();
          int Ncp = relative_cm_bra.Nc();
          int lcp = relative_cm_bra.lc();
          int Lp = relative_cm_bra.L();
          int Sp = relative_cm_bra.S();
          int Jp = relative_cm_bra.J();
          int Tp = relative_cm_bra.T();

          basis::RelativeCMStateLSJTN relative_cm_ket(relative_cm_sector.ket_subspace(),relative_cm_ket_index);
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
                int relative_subspace_index_bra = relative_space.LookUpSubspaceIndex(
                    basis::RelativeSubspaceLSJT::SubspaceLabelsType(lrp,Sp,Jrp,Tp,grp)
                  );
                int gr = lr%2;
                int relative_subspace_index_ket = relative_space.LookUpSubspaceIndex(
                    basis::RelativeSubspaceLSJT::SubspaceLabelsType(lr,S,Jr,T,gr)
                  );

                // look up relative matrix element indices
                const basis::RelativeSubspaceLSJT& relative_bra_subspace = relative_space.GetSubspace(
                    relative_subspace_index_bra
                  );
                const basis::RelativeSubspaceLSJT& relative_ket_subspace = relative_space.GetSubspace(
                    relative_subspace_index_ket
                  );
                int relative_state_index_bra = relative_bra_subspace.LookUpStateIndex(
                    basis::RelativeStateLSJT::StateLabelsType(Nrp)
                  );
                int relative_state_index_ket = relative_ket_subspace.LookUpStateIndex(
                    basis::RelativeStateLSJT::StateLabelsType(Nr)
                  );
                // std::cout << relative_bra_subspace.LabelStr() << relative_ket_subspace.LabelStr() << std::endl;


                // canonicalize indices for matrix element lookup
                int canonical_relative_subspace_index_bra, canonical_relative_subspace_index_ket;
                int canonical_relative_state_index_bra, canonical_relative_state_index_ket;
                double canonicalization_factor;
                std::tie(
                    canonical_relative_subspace_index_bra, canonical_relative_subspace_index_ket,
                    canonical_relative_state_index_bra, canonical_relative_state_index_ket,
                    std::ignore,
                    canonicalization_factor
                  )
                  = basis::CanonicalizeIndicesJT(
                      relative_space,
                      J0, T0, g0,
                      symmetry_phase_mode,
                      relative_subspace_index_bra, relative_subspace_index_ket,
                      relative_state_index_bra, relative_state_index_ket
                    );

                // look up matrix element
                int relative_sector_index = relative_sectors.LookUpSectorIndex(
                    canonical_relative_subspace_index_bra,
                    canonical_relative_subspace_index_ket
                  );
                double relative_matrix_element = canonicalization_factor
                  * relative_matrices[relative_sector_index](
                      canonical_relative_state_index_bra, canonical_relative_state_index_ket
                    );
              
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

  void TransformOperatorRelativeLSJTToRelativeCMLSJTN(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      const std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      const std::array<basis::MatrixVector,3>& relative_component_matrices,
      const basis::RelativeCMSpaceLSJTN& relative_cm_lsjtn_space,
      std::array<basis::RelativeCMSectorsLSJTN,3>& relative_cm_lsjtn_component_sectors,
      std::array<basis::MatrixVector,3>& relative_cm_lsjtn_component_matrices
    )
  {
    for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
      // for each isospin component
      {

        // enumerate sectors
        relative_cm_lsjtn_component_sectors[T0]
          = basis::RelativeCMSectorsLSJTN(relative_cm_lsjtn_space,operator_labels.J0,T0,operator_labels.g0);

        // populate matrices
        relative_cm_lsjtn_component_matrices[T0].resize(relative_cm_lsjtn_component_sectors[T0].size());
        for (int sector_index=0; sector_index<relative_cm_lsjtn_component_sectors[T0].size(); ++sector_index)
          {
            const basis::RelativeCMSectorsLSJTN::SectorType& relative_cm_lsjtn_sector
              = relative_cm_lsjtn_component_sectors[T0].GetSector(sector_index);
            relative_cm_lsjtn_component_matrices[T0][sector_index] = moshinsky::RelativeCMMatrixLSJTN(
                relative_space,relative_component_sectors[T0],relative_component_matrices[T0],
                relative_cm_lsjtn_sector,
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
  TransformationMatrixRelativeCMTwoBodyLSJTN(
      const basis::RelativeCMSubspaceLSJTN& relative_cm_subspace,
      const basis::TwoBodySubspaceLSJTN& two_body_subspace
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
    assert(relative_cm_subspace.g()==two_body_subspace.g());

    // retrieve L for sector
    int L = relative_cm_subspace.L();

    // populate matrix
    for (int index_relative_cm=0; index_relative_cm<relative_cm_subspace.size(); ++index_relative_cm)
      for (int index_two_body=0; index_two_body<two_body_subspace.size(); ++index_two_body)
        {

          // retrieve relative-cm state
          basis::RelativeCMStateLSJTN relative_cm_state(relative_cm_subspace,index_relative_cm);
          int Nr = relative_cm_state.Nr();
          int lr = relative_cm_state.lr();
          int nr = (Nr - lr)/2;
          int Nc = relative_cm_state.Nc();
          int lc = relative_cm_state.lc();
          int nc = (Nc - lc)/2;
          int L = relative_cm_state.L();

          // retrieve two-body state
          basis::TwoBodyStateLSJTN two_body_state(two_body_subspace,index_two_body);
          int N1 = two_body_state.N1();
          int l1 = two_body_state.l1();
          int n1 = (N1 - l1)/2;
          int N2 = two_body_state.N2();
          int l2 = two_body_state.l2();
          int n2 = (N2 - l2)/2;

          // evaluate bracket
          matrix(index_relative_cm,index_two_body)
            = sqrt(2.)*moshinsky::MoshinskyBracket(nr,lr,nc,lc,n1,l1,n2,l2,L);
				  
        }

    // return matrix
    //  std::cout << "Moshinsky " << matrix << std::endl;
    return matrix;
  }

  Eigen::MatrixXd
  TwoBodyMatrixLSJTN(
      const basis::RelativeCMSectorsLSJTN::SectorType& relative_cm_sector,
      const basis::TwoBodySectorsLSJTN::SectorType& two_body_sector,
      const Eigen::MatrixXd& relative_cm_matrix
    )
  {

    // obtain transformation matrices
    Eigen::MatrixXd bra_transformation_matrix = TransformationMatrixRelativeCMTwoBodyLSJTN(
        relative_cm_sector.bra_subspace(),
        two_body_sector.bra_subspace()
      );
    Eigen::MatrixXd ket_transformation_matrix = TransformationMatrixRelativeCMTwoBodyLSJTN(
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

  void TransformOperatorRelativeCMLSJTNToTwoBodyLSJTN(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeCMSpaceLSJTN& relative_cm_lsjtn_space,
      const std::array<basis::RelativeCMSectorsLSJTN,3>& relative_cm_lsjtn_component_sectors,
      const std::array<basis::MatrixVector,3>& relative_cm_lsjtn_component_matrices,
      const basis::TwoBodySpaceLSJTN& two_body_lsjtn_space,
      std::array<basis::TwoBodySectorsLSJTN,3>& two_body_lsjtn_component_sectors,
      std::array<basis::MatrixVector,3>& two_body_lsjtn_component_matrices
    )
  {
    for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
      // for each isospin component
      {

        // enumerate target sectors
        two_body_lsjtn_component_sectors[T0]
          = basis::TwoBodySectorsLSJTN(two_body_lsjtn_space,operator_labels.J0,T0,operator_labels.g0);
         
        // std::cout << " T0 " << T0
        //           << " sectors " << two_body_lsjtn_component_sectors[T0].size()
        //           << std::endl;

        // populate matrices
        two_body_lsjtn_component_matrices[T0].resize(two_body_lsjtn_component_sectors[T0].size());
        for (int sector_index=0; sector_index<two_body_lsjtn_component_sectors[T0].size(); ++sector_index)
          {
            // make reference to target sector
            const basis::TwoBodySectorsLSJTN::SectorType& two_body_lsjtn_sector
              = two_body_lsjtn_component_sectors[T0].GetSector(sector_index);

            // look up source sector
            int relative_cm_lsjtn_bra_subspace_index = relative_cm_lsjtn_space.LookUpSubspaceIndex(
                basis::RelativeCMSubspaceLSJTNLabels(
                    two_body_lsjtn_sector.bra_subspace().L(),
                    two_body_lsjtn_sector.bra_subspace().S(),
                    two_body_lsjtn_sector.bra_subspace().J(),
                    two_body_lsjtn_sector.bra_subspace().T(),
                    two_body_lsjtn_sector.bra_subspace().g(),
                    two_body_lsjtn_sector.bra_subspace().N()
                  )
              );
            int relative_cm_lsjtn_ket_subspace_index = relative_cm_lsjtn_space.LookUpSubspaceIndex(
                basis::RelativeCMSubspaceLSJTNLabels(
                    two_body_lsjtn_sector.ket_subspace().L(),
                    two_body_lsjtn_sector.ket_subspace().S(),
                    two_body_lsjtn_sector.ket_subspace().J(),
                    two_body_lsjtn_sector.ket_subspace().T(),
                    two_body_lsjtn_sector.ket_subspace().g(),
                    two_body_lsjtn_sector.ket_subspace().N()
                  )
              );
            int relative_cm_lsjtn_sector_index = relative_cm_lsjtn_component_sectors[T0].LookUpSectorIndex(
                relative_cm_lsjtn_bra_subspace_index,
                relative_cm_lsjtn_ket_subspace_index
              );
            const basis::RelativeCMSectorsLSJTN::SectorType& relative_cm_lsjtn_sector
              = relative_cm_lsjtn_component_sectors[T0].GetSector(relative_cm_lsjtn_sector_index);
            const Eigen::MatrixXd& relative_cm_lsjtn_matrix
              = relative_cm_lsjtn_component_matrices[T0][relative_cm_lsjtn_sector_index];

            // transform
            two_body_lsjtn_component_matrices[T0][sector_index] = moshinsky::TwoBodyMatrixLSJTN(
                relative_cm_lsjtn_sector,
                two_body_lsjtn_sector,
                relative_cm_lsjtn_matrix
              );
          }
      }

  }


  ////////////////////////////////////////////////////////////////
  // recoupling to jjJT scheme
  ////////////////////////////////////////////////////////////////

  Eigen::SparseMatrix<double>
    TransformationMatrixTwoBodyLSJTNToTwoBodyJJJTN(
        const basis::TwoBodySubspaceLSJTN& two_body_lsjtn_subspace,
        const basis::TwoBodySubspaceJJJTN& two_body_jjjtn_subspace
      )
  {

    // Since it is possible to know an upper bound on the number of
    // entries per column in advance and perform insertions in
    //  row order within each column, we can follow the pattern of a
    // "sorted insertion" as described in the Eigen 3.2.8
    // documentation:
    //
    //   eigen-doc-3.2.8/group__TutorialSparse.html
    //
    //     SparseMatrix<double> mat(rows,cols);         // default is column major
    //     mat.reserve(VectorXi::Constant(cols,6));
    //     for each i,j such that v_ij != 0
    //        mat.insert(i,j) = v_ij;                   // alternative: mat.coeffRef(i,j) += v_ij;
    //     mat.makeCompressed();                        // optional
    //
    // Is there a preferred row/column major order for multiplication
    // against a dense matrix?
    //
    // dv2 = sm1 * dv1;
    // dm2 = dm1 * sm1.adjoint();
    // dm2 = 2. * sm1 * dm1;

    // set up matrix to hold results
    Eigen::SparseMatrix<double> matrix(
        two_body_lsjtn_subspace.size(),
        two_body_jjjtn_subspace.size()
      );
    matrix.reserve(1);  // we expect maximum one entry per column

    // guard against call with mismatched subspaces
    assert(two_body_lsjtn_subspace.N()==two_body_jjjtn_subspace.N());
    assert(two_body_lsjtn_subspace.J()==two_body_jjjtn_subspace.J());
    assert(two_body_lsjtn_subspace.T()==two_body_jjjtn_subspace.T());
    assert(two_body_lsjtn_subspace.g()==two_body_jjjtn_subspace.g());

    // retrieve subspace angular momentum labels
    int L = two_body_lsjtn_subspace.L();
    int S = two_body_lsjtn_subspace.S();
    int J = two_body_lsjtn_subspace.J();

    // populate matrix -- scan columns
    for (int index_two_body_jjjtn=0; index_two_body_jjjtn<two_body_jjjtn_subspace.size(); ++index_two_body_jjjtn)
      // for each target JJJTN state
      {
        // retrieve target JJJTN state
        basis::TwoBodyStateJJJTN two_body_jjjtn_state(two_body_jjjtn_subspace,index_two_body_jjjtn);
        int N1 = two_body_jjjtn_state.N1();
        int l1 = two_body_jjjtn_state.l1();
        HalfInt j1 = two_body_jjjtn_state.j1();
        int N2 = two_body_jjjtn_state.N2();
        int l2 = two_body_jjjtn_state.l2();
        HalfInt j2 = two_body_jjjtn_state.j2();

        // check for corresponding source LSJTN state
        //
        // Note that canonicalization of orbital labels should not be
        // a concern.  We claim that canonical order (N1,j1)<=(N2,j2)
        // within (N,N1,j1,N2,j2) implies canonical order
        // (N1,l1)<=(N2,l2) within (N,N1,l1,N2,l2).  Argument: If
        // N1<N2, canonical ordering is ordering by N, which is
        // trivially preserved.  If N1==N2, then j1<=j2 has two cases.
        // If j1<j2, then l1<=l2 is trivially obtained, since j and l
        // differ by at most 1/2.  And, if j1==j2, then l1==l2, by
        // this same constraint combined with the realization that l
        // goes by steps of two within a major shell, so j values map
        // uniquely onto l values.
        int index_two_body_lsjtn = two_body_lsjtn_subspace.LookUpStateIndex(
            basis::TwoBodyStateLSJTNLabels(N1,l1,N2,l2)
          );
        if (index_two_body_lsjtn==basis::kNone)
          continue;


        // Alternately:
        //
        // basis::TwoBodyStateLSJTN two_body_lsjtn_state(
        //     two_body_lsjtn_subspace,
        //     basis::TwoBodyStateLSJTNLabels(N1,l1,N2,l2)
        //   );
        // int index_two_body_lsjtnt = two_body_lsjtn_state.index();

        // evaluate bracket
        const HalfInt s(1,2);
        double matrix_element = am::Unitary9J(
            l1,s,j1,
            l2,s,j2,
            L,S,J
          );

        matrix.insert(index_two_body_lsjtn,index_two_body_jjjtn) = matrix_element;
        }

    // compress matrix
    //   (does this matter?)
    matrix.makeCompressed();

    // return matrix
    return matrix;
  }

  bool SubspacesConnectedLSJTNToJJJTN(
        const basis::TwoBodySubspaceLSJTN& two_body_lsjtn_subspace,
        const basis::TwoBodySubspaceJJJTN& two_body_jjjtn_subspace
    )
  // Check if given LSJTN subspace contributes to given JJJTN subspace.
  //
  // Note that triangularity of (L,S,J) is already enforced by
  // construction in the source sector, so it need not be checked.
  //
  // Arguments:
  //     two_body_lsjtn_subspace (input): LSJTN subspace
  //     two_body_jjjtn_subspace (input): JJJTN subspace
  //
  // Returns:
  //     whether or not subspaces are connected
  {
    bool connected = true;

    // check equality of spectator labels
    connected &= (two_body_lsjtn_subspace.N()==two_body_jjjtn_subspace.N());
    connected &= (two_body_lsjtn_subspace.J()==two_body_jjjtn_subspace.J());
    connected &= (two_body_lsjtn_subspace.T()==two_body_jjjtn_subspace.T());
    connected &= (two_body_lsjtn_subspace.g()==two_body_jjjtn_subspace.g());

    return connected;
  }

  Eigen::MatrixXd 
    TwoBodyMatrixJJJTN(
        const basis::TwoBodySpaceLSJTN& two_body_lsjtn_space,
        const basis::TwoBodySectorsLSJTN& two_body_lsjtn_sectors,
        const basis::MatrixVector& two_body_lsjtn_matrices,
        const basis::TwoBodySectorsJJJTN::SectorType& two_body_jjjtn_sector,
        int J0, int T0, int g0,
        basis::SymmetryPhaseMode symmetry_phase_mode
      )
  {
    // We have two options for how we complete the quadruple sum over
    // source (L',S')(L,S) labels:
    //
    // (1) Source subspace lookup: (a) For the given target NJTg ket
    // subspace, we can scan over all S(=0,1,2) and L which are
    // triangular with J, and look up the corresponding source ket
    // subspace.  (b) Similarly, the given target NJTg bra subspace we
    // can scan over all S'(=0,1,2) and L' which are triangular with
    // J', and look up the corresponding source bra subspace.  These
    // lookups require access to the TwoBodyLSJTN space (unless lookup
    // by labels is implemented as part of BaseSectors).
    //
    // (2) Source sector scan: We can just scan through all source
    // sectors and pick off all sectors where the source bra/ket NJTg
    // labels match the target bra/ket NJTg labels.  As the source
    // sector list grows long, this brute force approach might become
    // slightly inefficient, but it is the simplest to code and does
    // avoid lookups.
    //
    // Note however that the separate canonicalization of (bra,ket)
    // subspace order in the LSJTN and JJJPN schemes means that the
    // comparison of bra/ket labels in the source LSJTN sector with
    // those in the target JJJTN sector may require bra and ket
    // subspaces to be interchanged.  This somewhat compromises the
    // simplicity of the source sector scan scheme.
    //
    // So we abandoned this option for...
    //
    // (3) Source *subspace* scan: We can just scan through all source
    // subspaces and pick off all subspaces where the source bra/ket
    // NJTg labels match the target bra/ket NJTg labels.  We then look
    // up the corresponding *canonicalized* source sector.

    // set up matrix to hold results
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(
        two_body_jjjtn_sector.bra_subspace().size(),
        two_body_jjjtn_sector.ket_subspace().size()
      );
    //std::cout << fmt::format("target {}x{}",matrix.rows(),matrix.cols()) << std::endl; 

    // identify source subspaces which contribute (and collect their
    // transformation matrices)
    //
    // We take the "easy" approach of scanning once through all source
    // subspaces.  But one could also explicitly attempt lookups of the small set
    // of possibly subspaces given by (1) S=0,1 and (2) L in triangle(J,S).

    std::vector<int> bra_two_body_lsjtn_subspace_indices, ket_two_body_lsjtn_subspace_indices;
    std::vector<Eigen::SparseMatrix<double>> bra_transformation_matrices, ket_transformation_matrices;

    // scan for bra subspace
    for (
        int bra_two_body_lsjtn_subspace_index=0;
        bra_two_body_lsjtn_subspace_index<two_body_lsjtn_space.size();
        ++bra_two_body_lsjtn_subspace_index
      )
      {
        // check if source subspace contributes
        const basis::TwoBodySubspaceLSJTN& bra_two_body_lsjtn_subspace = two_body_lsjtn_space.GetSubspace(bra_two_body_lsjtn_subspace_index);
        if (!SubspacesConnectedLSJTNToJJJTN(bra_two_body_lsjtn_subspace,two_body_jjjtn_sector.bra_subspace()))
          continue;

        // save source subspace index and transformation matrix
        bra_two_body_lsjtn_subspace_indices.push_back(bra_two_body_lsjtn_subspace_index);
        Eigen::SparseMatrix<double> bra_transformation_matrix = TransformationMatrixTwoBodyLSJTNToTwoBodyJJJTN(
            bra_two_body_lsjtn_subspace,
            two_body_jjjtn_sector.bra_subspace()
          );
        bra_transformation_matrices.push_back(bra_transformation_matrix);
      }

    // scan for ket subspace
    for (
        int ket_two_body_lsjtn_subspace_index=0;
        ket_two_body_lsjtn_subspace_index<two_body_lsjtn_space.size();
        ++ket_two_body_lsjtn_subspace_index
      )
      {
        // check if source subspace contributes
        const basis::TwoBodySubspaceLSJTN& ket_two_body_lsjtn_subspace = two_body_lsjtn_space.GetSubspace(ket_two_body_lsjtn_subspace_index);
        if (!SubspacesConnectedLSJTNToJJJTN(ket_two_body_lsjtn_subspace,two_body_jjjtn_sector.ket_subspace()))
          continue;

        // save source subspace index and transformation matrix
        ket_two_body_lsjtn_subspace_indices.push_back(ket_two_body_lsjtn_subspace_index);
        Eigen::SparseMatrix<double> ket_transformation_matrix = TransformationMatrixTwoBodyLSJTNToTwoBodyJJJTN(
            ket_two_body_lsjtn_subspace,
            two_body_jjjtn_sector.ket_subspace()
          );
        ket_transformation_matrices.push_back(ket_transformation_matrix);
      }

    for (int bra_subspace_metaindex=0; bra_subspace_metaindex<bra_two_body_lsjtn_subspace_indices.size(); ++bra_subspace_metaindex)
      for (int ket_subspace_metaindex=0; ket_subspace_metaindex<ket_two_body_lsjtn_subspace_indices.size(); ++ket_subspace_metaindex)
        {
          // alias source subspace information -- bra
          int bra_two_body_lsjtn_subspace_index = bra_two_body_lsjtn_subspace_indices[bra_subspace_metaindex];
          const basis::TwoBodySubspaceLSJTN& bra_two_body_lsjtn_subspace = two_body_lsjtn_space.GetSubspace(bra_two_body_lsjtn_subspace_index);
          Eigen::SparseMatrix<double>& bra_transformation_matrix = bra_transformation_matrices[bra_subspace_metaindex];

          // alias source subspace information -- ket
          int ket_two_body_lsjtn_subspace_index = ket_two_body_lsjtn_subspace_indices[ket_subspace_metaindex];
          const basis::TwoBodySubspaceLSJTN& ket_two_body_lsjtn_subspace = two_body_lsjtn_space.GetSubspace(ket_two_body_lsjtn_subspace_index);
          Eigen::SparseMatrix<double>& ket_transformation_matrix = ket_transformation_matrices[ket_subspace_metaindex];

          //std::cout << fmt::format("  source {}x{}",bra_two_body_lsjtn_subspace.size(),ket_two_body_lsjtn_subspace.size()) << std::endl; 

          // identify source sector

          int bra_canonical_two_body_lsjtn_subspace_index,ket_canonical_two_body_lsjtn_subspace_index;
          bool swapped_subspaces;
          double canonicalization_factor;
          std::tie(
              bra_canonical_two_body_lsjtn_subspace_index,ket_canonical_two_body_lsjtn_subspace_index,
              swapped_subspaces,
              canonicalization_factor
            )
            = CanonicalizeIndicesJT(
                two_body_lsjtn_space,
                J0,T0,g0,symmetry_phase_mode,
                bra_two_body_lsjtn_subspace_index,ket_two_body_lsjtn_subspace_index
              );
          int two_body_lsjtn_sector_index = two_body_lsjtn_sectors.LookUpSectorIndex(
              bra_canonical_two_body_lsjtn_subspace_index,
              ket_canonical_two_body_lsjtn_subspace_index
            );

          // accumulate contribution from source matrix
          //
          // We are assuming that diagonal sector matrices are fully-populated square matrices.

          const Eigen::MatrixXd& two_body_lsjtn_matrix = two_body_lsjtn_matrices[two_body_lsjtn_sector_index];

          if (swapped_subspaces)
            {
              matrix
                += bra_transformation_matrix.transpose()
                * canonicalization_factor * two_body_lsjtn_matrix.transpose()
                * ket_transformation_matrix;
            }
          else
            {
              matrix
                += bra_transformation_matrix.transpose()
                * two_body_lsjtn_matrix
                * ket_transformation_matrix;
            }
        }

    // return target matrix
    return matrix;


    
  }

  void TransformOperatorTwoBodyLSJTNToTwoBodyJJJTN(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::TwoBodySpaceLSJTN& two_body_lsjtn_space,
      const std::array<basis::TwoBodySectorsLSJTN,3>& two_body_lsjtn_component_sectors,
      const std::array<basis::MatrixVector,3>& two_body_lsjtn_component_matrices,
      const basis::TwoBodySpaceJJJTN& two_body_jjjtn_space,
      std::array<basis::TwoBodySectorsJJJTN,3>& two_body_jjjtn_component_sectors,
      std::array<basis::MatrixVector,3>& two_body_jjjtn_component_matrices
    )
  {
    for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
      // for each isospin component
      {

        // enumerate target sectors
        two_body_jjjtn_component_sectors[T0]
          = basis::TwoBodySectorsJJJTN(two_body_jjjtn_space,operator_labels.J0,T0,operator_labels.g0);

        // populate matrices
        two_body_jjjtn_component_matrices[T0].resize(two_body_jjjtn_component_sectors[T0].size());
        for (int sector_index=0; sector_index<two_body_jjjtn_component_sectors[T0].size(); ++sector_index)
          // for each target sector
          {
            // make reference to target sector
            const basis::TwoBodySectorsJJJTN::SectorType& two_body_jjjtn_sector
              = two_body_jjjtn_component_sectors[T0].GetSector(sector_index);

            // make references to isospin component of source operator
            const basis::TwoBodySectorsLSJTN& two_body_lsjtn_sectors = two_body_lsjtn_component_sectors[T0];
            const basis::MatrixVector& two_body_lsjtn_matrices = two_body_lsjtn_component_matrices[T0];
         
            // transform
            Eigen::MatrixXd& matrix = two_body_jjjtn_component_matrices[T0][sector_index];
            matrix = TwoBodyMatrixJJJTN(
                two_body_lsjtn_space,
                two_body_lsjtn_sectors,
                two_body_lsjtn_matrices,
                two_body_jjjtn_sector,
                operator_labels.J0,T0,operator_labels.g0,operator_labels.symmetry_phase_mode
              );
          }
      }

  }

  ////////////////////////////////////////////////////////////////
  // branching to jjJpn scheme
  ////////////////////////////////////////////////////////////////

  Eigen::MatrixXd 
    TwoBodyMatrixJJJPN(
        const basis::OperatorLabelsJT& operator_labels,
        const basis::TwoBodySpaceJJJT& two_body_jjjt_space,
        const std::array<basis::TwoBodySectorsJJJT,3>& two_body_jjjt_component_sectors,
        const std::array<basis::MatrixVector,3>& two_body_jjjt_component_matrices,
        const basis::TwoBodySectorsJJJPN::SectorType& two_body_jjjpn_sector
      )
  {

    // define isospin Clebsch-Gordan coefficients
    //
    // as static arrays over T0
    static const double kPPCoefficients[] = {+1,+sqrt(1/2.),+sqrt(1/10.)};
    static const double kNNCoefficients[] = {+1,-sqrt(1/2.),+sqrt(1/10.)};
    static const double kPNCoefficients11[] = {+1,0,-sqrt(2./5.)};
    static const double kPNCoefficient10 = 1.;
    static const double kPNCoefficient01 = -sqrt(1/3.);
    static const double kPNCoefficient00 = 1.;

    // set up matrix to hold results
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(
        two_body_jjjpn_sector.bra_subspace().size(),
        two_body_jjjpn_sector.ket_subspace().size()
      );

    // extract sector labels
    basis::TwoBodySpeciesPN two_body_species = two_body_jjjpn_sector.bra_subspace().two_body_species();
    int Jp = two_body_jjjpn_sector.bra_subspace().J();
    int J = two_body_jjjpn_sector.ket_subspace().J();

    // scan source sectors
    for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
      // for each isospin component
      {
        for (int Tp=0; Tp<=1; ++Tp)
          // for bra isospin
          for (int T=0; T<=1; ++T)
            // for ket isospin
            {
              // impose isospin triangularity
              if (!am::AllowedTriangle(Tp,T0,T))
                continue;

              // impose Tz sufficiency for current sector two-body species
              if (!(
                      (two_body_species==basis::TwoBodySpeciesPN::kPN)
                      || ((Tp==1)&&(T==1))
                    ))
                continue;

              // look up source subspaces
              int two_body_jjjt_subspace_index_bra = two_body_jjjt_space.LookUpSubspaceIndex(
                  basis::TwoBodySubspaceJJJTLabels(
                      two_body_jjjpn_sector.bra_subspace().J(),
                      Tp,
                      two_body_jjjpn_sector.bra_subspace().g()
                    )
                );
              const basis::TwoBodySubspaceJJJT& two_body_jjjt_subspace_bra
                = two_body_jjjt_space.GetSubspace(two_body_jjjt_subspace_index_bra);
              int two_body_jjjt_subspace_index_ket = two_body_jjjt_space.LookUpSubspaceIndex(
                  basis::TwoBodySubspaceJJJTLabels(
                      two_body_jjjpn_sector.ket_subspace().J(),
                      T,
                      two_body_jjjpn_sector.ket_subspace().g()
                    )
                );
              const basis::TwoBodySubspaceJJJT& two_body_jjjt_subspace_ket
                = two_body_jjjt_space.GetSubspace(two_body_jjjt_subspace_index_ket);

              // std::cout
              //   << " JJJT bra subspace " << two_body_jjjt_subspace_bra.LabelStr()
              //   << " JJJT ket subspace " << two_body_jjjt_subspace_ket.LabelStr()
              //   << std::endl;
              // std::cout
              //   << " JJJT bra subspace contents " << std::endl
              //   << two_body_jjjt_subspace_bra.DebugStr()
              //   << " JJJT ket subspace contents " << std::endl
              //   << two_body_jjjt_subspace_ket.DebugStr();

              // look up *canonicalized* source sector
              //
              // Note on canonicalization:
              //
              // With an isoscalar operator, there are no
              // <T=0|...|T=1> cross sectors.  When bra and ket have
              // the same T, (J,g,s) canonicality of subspaces (fixed
              // s between bra and ket) ensures (J,T,g) canonicality
              // of subspaces (fixed T between bra and ket).  We can
              // therefore assume canonicality of source bra and ket
              // subspaces.
              //
              // However, when <T=0|...|T=1> sectors are present, we
              // must canonicalize the sector bra and ket subspaces
              // before looking up the sector.
              //
              // We can do this "early", outside the loop over matrix
              // elements, for efficiency, or "late", inside the usual
              // full canonicalization of the matrix element.
              //
              // In either case, we will have to go back and
              // canonicalize the matrix element indices as well,
              // since two-body state indices in the jjJpn subspaces
              // do not map trivially to indices in the jjJT
              // subspaces.  In fact, subspace dimensions can be quite
              // different, and two-body state labels are not in 1-1
              // correspondence between jjJT and jjJpn subspaces.
              //
              // Notably, orbitals (i1,i2) can appear in the order
              // i1>i2 in the pn subspaces of the jjJpn scheme but
              // only in the order i1<i2 in the jjJT subspaces.
              // (Slight differences in the sets of state labels also
              // arise due to the antisymmetry constraint on J and T
              // for like-orbital states.)  This ordering constraint
              // (and swap) must be taken into account in the state
              // index lookup in the pn sectors.

              int two_body_jjjt_sector_index = two_body_jjjt_component_sectors[T0].LookUpSectorIndex(
                  std::min(two_body_jjjt_subspace_index_bra,two_body_jjjt_subspace_index_ket),
                  std::max(two_body_jjjt_subspace_index_bra,two_body_jjjt_subspace_index_ket)
                );
              const basis::TwoBodySectorsJJJT::SectorType& two_body_jjjt_sector
                = two_body_jjjt_component_sectors[T0].GetSector(two_body_jjjt_sector_index);
              const Eigen::MatrixXd& two_body_jjjt_matrix
                = two_body_jjjt_component_matrices[T0][two_body_jjjt_sector_index];

              // identify isospin Clebsch-Gordan coefficient to apply
              // to this source sector
              double isospin_coefficient;
              if (two_body_species==basis::TwoBodySpeciesPN::kPP)
                {
                  isospin_coefficient = kPPCoefficients[T0];
                }
              else if (two_body_species==basis::TwoBodySpeciesPN::kPP)
                {
                  isospin_coefficient = kNNCoefficients[T0];
                }
              else if (two_body_species==basis::TwoBodySpeciesPN::kPN)
                {
                  // Note: There are just the isospin Clebsch-Gordan
                  // coefficients for the Wigner-Eckhart branching of
                  // the isospin reduced matrix elements.  These do not
                  // include the Clebsch/normalization factors in the
                  // expansion of the pn states in terms of T=0/1
                  // states, which are treated below at the level of
                  // individual matrix elements.
                  if ((Tp==1)&&(T==1))
                    isospin_coefficient = kPNCoefficients11[T0];
                  else if ((Tp==0)&&(T==1))
                    isospin_coefficient = kPNCoefficient01;
                  else if ((Tp==1)&&(T==0))
                    isospin_coefficient = kPNCoefficient10;
                  else if ((Tp==0)&&(T==0))
                    isospin_coefficient = kPNCoefficient00;
                }


              // accumulate contributions to target sector matrix elements
              for (int bra_index = 0; bra_index < two_body_jjjpn_sector.bra_subspace().size(); ++bra_index)
                for (int ket_index = 0; ket_index < two_body_jjjpn_sector.ket_subspace().size(); ++ket_index)
                  // for each target matrix element
                  {
                
                    // ensure target matrix element is canonical, if
                    // diagonal sector
                    if (two_body_jjjpn_sector.IsDiagonal())
                      if (!(bra_index<=ket_index))
                        continue;

                    // retrieve target states
                    basis::TwoBodyStateJJJPN two_body_jjjpn_bra(two_body_jjjpn_sector.bra_subspace(),bra_index);
                    basis::TwoBodyStateJJJPN two_body_jjjpn_ket(two_body_jjjpn_sector.ket_subspace(),ket_index);

                    // std::cout << " JJJPN target states "
                    //           << two_body_jjjpn_bra.LabelStr() 
                    //           << two_body_jjjpn_ket.LabelStr() << std::endl;
                
                    // verify source states allowed by antisymmetry
                    //
                    // Recall jjJT antisymmetry constraint:
                    //   J+T~1 if (N1,j1)==(N2,j2)   
                    if (two_body_jjjpn_bra.index1()==two_body_jjjpn_bra.index2())
                      if (!((Jp+Tp)%2==1))
                        continue;
                    if (two_body_jjjpn_ket.index1()==two_body_jjjpn_ket.index2())
                      if (!((J+T)%2==1))
                        continue;

                    // retrieve source state indices
                    int index1_bra = two_body_jjjpn_bra.index1();
                    int index2_bra = two_body_jjjpn_bra.index2();
                    int index1_ket = two_body_jjjpn_ket.index1();
                    int index2_ket = two_body_jjjpn_ket.index2();
                    int N1_bra = two_body_jjjpn_bra.GetOrbital1().N();
                    int N2_bra = two_body_jjjpn_bra.GetOrbital2().N();
                    int N1_ket = two_body_jjjpn_ket.GetOrbital1().N();
                    int N2_ket = two_body_jjjpn_ket.GetOrbital2().N();
                    HalfInt j1_bra = two_body_jjjpn_bra.GetOrbital1().j();
                    HalfInt j2_bra = two_body_jjjpn_bra.GetOrbital2().j();
                    HalfInt j1_ket = two_body_jjjpn_ket.GetOrbital1().j();
                    HalfInt j2_ket = two_body_jjjpn_ket.GetOrbital2().j();

                    // calculate pn state expansion coefficients
                    double pn_normalization_factor = 1.0;
                    if (two_body_species==basis::TwoBodySpeciesPN::kPN)
                      {
                        pn_normalization_factor = 0.5;
                        if (index1_bra==index2_bra)
                          pn_normalization_factor *= sqrt(2.);
                        if (index1_ket==index2_ket)
                          pn_normalization_factor *= sqrt(2.);
                      }

                    // canonicalize indices of orbitals within
                    // two-body states
                    //
                    // Orbitals may need to be swapped going from
                    // |ab;J>_pn to |ab;JT>, if a>b, inducing a phase
                    // factor ~(ja+jb+J+T).
                    double canonicalization_factor_bra;
                    basis::TwoBodyStateJJJTLabels two_body_jjjt_state_labels_bra;
                    if (index1_bra <= index2_bra)
                      // OR: (std::pair<int,HalfInt>(N1_bra,j1_bra) <= std::pair<int,HalfInt>(N2_bra,j2_bra))
                      // state is canonical
                      {
                        canonicalization_factor_bra = 1.;
                        two_body_jjjt_state_labels_bra = basis::TwoBodyStateJJJTLabels(
                            N1_bra,j1_bra,N2_bra,j2_bra
                          );
                      }
                    else
                      // state is non-canonical
                      {
                        canonicalization_factor_bra = ParitySign(int(j1_bra+j2_bra)+Jp+Tp);
                        two_body_jjjt_state_labels_bra = basis::TwoBodyStateJJJTLabels(
                            N2_bra,j2_bra,N1_bra,j1_bra
                          );
                      }

                    // std::cout << " JJJT bra lookup "
                    //           << " " << two_body_jjjpn_bra.GetOrbital1().N()
                    //           << " " << two_body_jjjpn_bra.GetOrbital1().j()
                    //           << " " << two_body_jjjpn_bra.GetOrbital2().N()
                    //           << " " << two_body_jjjpn_bra.GetOrbital2().j()
                    //           << std::endl;

                    int two_body_jjjt_state_index_bra
                      = two_body_jjjt_subspace_bra.LookUpStateIndex(
                          two_body_jjjt_state_labels_bra
                        );

                    double canonicalization_factor_ket;
                    basis::TwoBodyStateJJJTLabels two_body_jjjt_state_labels_ket;
                    if (index1_ket <= index2_ket)
                      // OR: (std::pair<int,HalfInt>(N1_ket,j1_ket) <= std::pair<int,HalfInt>(N2_ket,j2_ket))
                      // state is canonical
                      {
                        canonicalization_factor_ket = 1.;
                        two_body_jjjt_state_labels_ket = basis::TwoBodyStateJJJTLabels(
                            N1_ket,j1_ket,N2_ket,j2_ket
                          );
                      }
                    else
                      // state is non-canonical
                      {
                        canonicalization_factor_ket = ParitySign(int(j1_ket+j2_ket)+Jp+Tp);
                        two_body_jjjt_state_labels_ket = basis::TwoBodyStateJJJTLabels(
                            N2_ket,j2_ket,N1_ket,j1_ket
                          );
                      }

                    // std::cout << " JJJT ket lookup "
                    //           << " " << two_body_jjjpn_ket.GetOrbital1().N()
                    //           << " " << two_body_jjjpn_ket.GetOrbital1().j()
                    //           << " " << two_body_jjjpn_ket.GetOrbital2().N()
                    //           << " " << two_body_jjjpn_ket.GetOrbital2().j()
                    //           << std::endl;

                    int two_body_jjjt_state_index_ket
                      = two_body_jjjt_subspace_ket.LookUpStateIndex(
                          two_body_jjjt_state_labels_ket
                        );

                    // canonicalize indices for matrix element lookup
                    //
                    // Recall that we have already looked up the
                    // canonicalized sector index, so we do not care
                    // any more about the subspace indices returned by
                    // this call.
                    int canonical_two_body_jjjt_subspace_index_bra, canonical_two_body_jjjt_subspace_index_ket;
                    int canonical_two_body_jjjt_state_index_bra, canonical_two_body_jjjt_state_index_ket;
                    double canonicalization_factor;
                    std::tie(
                        canonical_two_body_jjjt_subspace_index_bra, canonical_two_body_jjjt_subspace_index_ket,
                        canonical_two_body_jjjt_state_index_bra, canonical_two_body_jjjt_state_index_ket,
                        std::ignore,
                        canonicalization_factor
                      )
                      = basis::CanonicalizeIndicesJT(
                          two_body_jjjt_space,
                          operator_labels.J0, T0, operator_labels.g0,
                          operator_labels.symmetry_phase_mode,
                          two_body_jjjt_subspace_index_bra, two_body_jjjt_subspace_index_ket,
                          two_body_jjjt_state_index_bra, two_body_jjjt_state_index_ket
                        );


                    // look up source matrix element
                    double two_body_jjjt_matrix_element
                      = pn_normalization_factor
                      * canonicalization_factor_bra * canonicalization_factor_ket * canonicalization_factor
                      * two_body_jjjt_matrix(
                        canonical_two_body_jjjt_state_index_bra,canonical_two_body_jjjt_state_index_ket
                      );

                    // incorporate contribution
                    matrix(bra_index,ket_index) += isospin_coefficient * two_body_jjjt_matrix_element;

                  }
            }
      }

    // return target matrix
    return matrix;

  }

  void TransformOperatorTwoBodyJJJTToTwoBodyJJJPN(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::TwoBodySpaceJJJT& two_body_jjjt_space,
      const std::array<basis::TwoBodySectorsJJJT,3>& two_body_jjjt_component_sectors,
      const std::array<basis::MatrixVector,3>& two_body_jjjt_component_matrices,
      const basis::TwoBodySpaceJJJPN& two_body_jjjpn_space,
      basis::TwoBodySectorsJJJPN& two_body_jjjpn_sectors,
      basis::MatrixVector& two_body_jjjpn_matrices
    )
  {
    // enumerate target sectors
    int Tz0 = 0;
    two_body_jjjpn_sectors
      = basis::TwoBodySectorsJJJPN(two_body_jjjpn_space,operator_labels.J0,operator_labels.g0,Tz0);

    // populate matrices
    two_body_jjjpn_matrices.resize(two_body_jjjpn_sectors.size());
    for (int sector_index=0; sector_index<two_body_jjjpn_sectors.size(); ++sector_index)
      {
        // make reference to target sector
        const basis::TwoBodySectorsJJJPN::SectorType& two_body_jjjpn_sector
          = two_body_jjjpn_sectors.GetSector(sector_index);

        // transform
        Eigen::MatrixXd& matrix = two_body_jjjpn_matrices[sector_index];
        matrix = TwoBodyMatrixJJJPN(
            operator_labels,
            two_body_jjjt_space,
            two_body_jjjt_component_sectors,
            two_body_jjjt_component_matrices,
            two_body_jjjpn_sector
          );
      }

  }



  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
