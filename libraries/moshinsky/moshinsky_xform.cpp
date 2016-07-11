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
    assert(relative_cm_subspace.g()==two_body_subspace.g());

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
  // recoupling to jjJT scheme
  ////////////////////////////////////////////////////////////////

  Eigen::SparseMatrix<double>
    TransformationMatrixTwoBodyNLSJTToTwoBodyNJJJT(
        const basis::TwoBodySubspaceNLSJT& two_body_nlsjt_subspace,
        const basis::TwoBodySubspaceNJJJT& two_body_njjjt_subspace
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
        two_body_nlsjt_subspace.size(),
        two_body_njjjt_subspace.size()
      );
    matrix.reserve(1);  // we expect maximum one entry per column

    // guard against call with mismatched subspaces
    assert(two_body_nlsjt_subspace.N()==two_body_njjjt_subspace.N());
    assert(two_body_nlsjt_subspace.J()==two_body_njjjt_subspace.J());
    assert(two_body_nlsjt_subspace.T()==two_body_njjjt_subspace.T());
    assert(two_body_nlsjt_subspace.g()==two_body_njjjt_subspace.g());

    // retrieve subspace angular momentum labels
    int L = two_body_nlsjt_subspace.L();
    int S = two_body_nlsjt_subspace.S();
    int J = two_body_nlsjt_subspace.J();

    // populate matrix -- scan columns
    for (int index_two_body_njjjt=0; index_two_body_njjjt<two_body_njjjt_subspace.size(); ++index_two_body_njjjt)
      {
        // retrieve target NJJJT state
        basis::TwoBodyStateNJJJT two_body_njjjt_state(two_body_njjjt_subspace,index_two_body_njjjt);
        int N1 = two_body_njjjt_state.N1();
        int l1 = two_body_njjjt_state.l1();
        HalfInt j1 = two_body_njjjt_state.j1();
        int N2 = two_body_njjjt_state.N2();
        int l2 = two_body_njjjt_state.l2();
        HalfInt j2 = two_body_njjjt_state.j2();

        // look up corresponding source NLSJT state

        // basis::TwoBodyStateNLSJT two_body_nlsjt_state(
        //     two_body_nlsjt_subspace,
        //     basis::TwoBodyStateNLSJTLabels(N1,l1,N2,l2)
        //   );
        // int index_two_body_nlsjtt = two_body_nlsjt_state.index();

        int index_two_body_nlsjt = two_body_nlsjt_subspace.LookUpStateIndex(
            basis::TwoBodyStateNLSJTLabels(N1,l1,N2,l2)
          );

        // evaluate bracket
        const HalfInt s(1,2);
        double matrix_element = am::Unitary9J(
            l1,s,j1,
            l2,s,j2,
            L,S,J
          );

        matrix.insert(index_two_body_nlsjt,index_two_body_njjjt) = matrix_element;
        }

    // compress matrix
    //   (does this matter?)
    matrix.makeCompressed();

    // return matrix
    return matrix;
  }

  Eigen::MatrixXd 
    TwoBodyMatrixNJJJT(
        const basis::TwoBodySectorsNLSJT& two_body_nlsjt_sectors,
        const basis::MatrixVector& two_body_nlsjt_matrices,
        const basis::TwoBodySectorsNJJJT::SectorType& two_body_njjjt_sector
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
    // lookups require access to the TwoBodyNLSJT space (unless lookup
    // by labels is implemented as part of BaseSectors).
    //
    // (2) Source sector scan: We can just scan through all source
    // sectors and pick off all sectors where the source bra/ket NJTg
    // labels match the target bra/ket NJTg labels.  As the source
    // sector list grows long, this brute force approach might become
    // slightly inefficient, but it is the simplest to code and does
    // avoid lookups.

    // set up matrix to hold results
    Eigen::MatrixXd matrix(
        two_body_njjjt_sector.bra_subspace().size(),
        two_body_njjjt_sector.ket_subspace().size()
      );

    // scan source sectors
    for (int two_body_nlsjt_sector_index=0; two_body_nlsjt_sector_index<two_body_nlsjt_sectors.size(); ++two_body_nlsjt_sector_index)
      {
        // retrieve source sector
        const basis::TwoBodySectorsNLSJT::SectorType& two_body_nlsjt_sector
          = two_body_nlsjt_sectors.GetSector(two_body_nlsjt_sector_index);

        // check for sector match
        if (!(
                true
                // bra subspace match
                && (two_body_nlsjt_sector.bra_subspace().N()==two_body_njjjt_sector.bra_subspace().N())
                && (two_body_nlsjt_sector.bra_subspace().J()==two_body_njjjt_sector.bra_subspace().J())
                && (two_body_nlsjt_sector.bra_subspace().T()==two_body_njjjt_sector.bra_subspace().T())
                && (two_body_nlsjt_sector.bra_subspace().g()==two_body_njjjt_sector.bra_subspace().g())
                // ket subspace match
                && (two_body_nlsjt_sector.ket_subspace().N()==two_body_njjjt_sector.ket_subspace().N())
                && (two_body_nlsjt_sector.ket_subspace().J()==two_body_njjjt_sector.ket_subspace().J())
                && (two_body_nlsjt_sector.ket_subspace().T()==two_body_njjjt_sector.ket_subspace().T())
                && (two_body_nlsjt_sector.ket_subspace().g()==two_body_njjjt_sector.ket_subspace().g())
              ))
          continue;

        // obtain transformation matrices
        Eigen::SparseMatrix<double> bra_transformation_matrix = TransformationMatrixTwoBodyNLSJTToTwoBodyNJJJT(
            two_body_nlsjt_sector.bra_subspace(),
            two_body_njjjt_sector.bra_subspace()
          );
        Eigen::SparseMatrix<double> ket_transformation_matrix = TransformationMatrixTwoBodyNLSJTToTwoBodyNJJJT(
            two_body_nlsjt_sector.ket_subspace(),
            two_body_njjjt_sector.ket_subspace()
          );

        // obtain target matrix
        const Eigen::MatrixXd& two_body_nlsjt_cm_matrix = two_body_nlsjt_matrices[two_body_nlsjt_sector_index];
        matrix
          += bra_transformation_matrix.transpose()
          * two_body_nlsjt_cm_matrix
          * ket_transformation_matrix;
      }

    // return target matrix
    return matrix;
  }


  void TransformOperatorTwoBodyNLSJTToTwoBodyNJJJT(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::TwoBodySpaceNLSJT& two_body_nlsjt_space,
      const std::array<basis::TwoBodySectorsNLSJT,3>& two_body_nlsjt_component_sectors,
      const std::array<basis::MatrixVector,3>& two_body_nlsjt_component_matrices,
      const basis::TwoBodySpaceNJJJT& two_body_njjjt_space,
      std::array<basis::TwoBodySectorsNJJJT,3>& two_body_njjjt_component_sectors,
      std::array<basis::MatrixVector,3>& two_body_njjjt_component_matrices
    )
  {
    // TODO
  }



  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
