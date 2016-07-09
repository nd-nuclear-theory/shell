/****************************************************************
  moshinsky_xform.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "moshinsky/moshinsky_xform.h"

#include "moshinsky/moshinsky_bracket.h"

namespace moshinsky {

double RacahReductionFactorFirstSystemGT(
    const HalfInt& J1p, const HalfInt& Jp, 
    const HalfInt& J1, const HalfInt& J, 
    const HalfInt& J2, const HalfInt& J0
  )
{
  double value = ParitySign(J1p+J2+J+J0)
    *Hat(J1p)*Hat(J)
    *am::Wigner6J(J1p,Jp,J2,J,J1,J0);
  return value;
}

////////////////////////////////////////////////////////////////
// Moshinsky transformation
////////////////////////////////////////////////////////////////

Eigen::MatrixXd
MoshinskyMatrixNLSJT(
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
        matrix(index_relative_cm,index_two_body) = moshinsky::MoshinskyBracket(nr,lr,nc,lc,n1,l1,n2,l2,L);
				  
      }

  // return matrix
  //  std::cout << "Moshinsky " << matrix << std::endl;
  return matrix;
}


Eigen::MatrixXd
RelativeCMMatrixNLSJT(
    const basis::RelativeSpaceLSJT& relative_space,
    const basis::RelativeSectorsLSJT& relative_sectors,
    const basis::MatrixVector& relative_matrices,
    const typename basis::RelativeCMSectorsNLSJT::SectorType& relative_cm_sector,
    int J0
  )
{
  
  // initialize matrix
  Eigen::MatrixXd relative_cm_matrix = Eigen::MatrixXd::Zero(
      relative_cm_sector.ket_subspace().size(),
      relative_cm_sector.bra_subspace().size()
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

              // look up relative sector
              int grp = lrp%2;
              int relative_bra_subspace_index = relative_space.LookUpSubspaceIndex(
                  basis::RelativeSubspaceLSJT::SubspaceLabelsType(lrp,Sp,Jrp,Tp,grp)
                );
              int gr = lr%2;
              int relative_ket_subspace_index = relative_space.LookUpSubspaceIndex(
                  basis::RelativeSubspaceLSJT::SubspaceLabelsType(lr,S,Jr,T,gr)
                );
              int relative_sector_index = relative_sectors.LookUpSectorIndex(
                  relative_bra_subspace_index,
                  relative_ket_subspace_index
                );

              // look up relative matrix element
              const basis::RelativeSubspaceLSJT& relative_bra_subspace = relative_space.GetSubspace(
                  relative_bra_subspace_index
                );
              int relative_bra_index = relative_bra_subspace.LookUpStateIndex(
                  basis::RelativeStateLSJT::StateLabelsType(Nrp)
                );
              const basis::RelativeSubspaceLSJT& relative_ket_subspace = relative_space.GetSubspace(
                  relative_ket_subspace_index
                );
              int relative_ket_index = relative_ket_subspace.LookUpStateIndex(
                  basis::RelativeStateLSJT::StateLabelsType(Nr)
                );
              double relative_matrix_element = relative_matrices[relative_sector_index](relative_bra_index,relative_ket_index);
              
              // accumulate contribution to relative-cm matrix element
              relative_cm_matrix_element
                += am::Unitary6JZ(lrp,lc,Lp,Jp,Sp,Jrp) * am::Unitary6JZ(lr,lc,L,J,S,Jr)
                * RacahReductionFactorFirstSystemGT(Jrp,Jp,Jr,J,lr,J0)
                * relative_matrix_element;
            }	

        // save matrix element
        relative_cm_matrix(relative_cm_bra_index,relative_cm_ket_index) = relative_cm_matrix_element;
      }

  // return matrix
  return relative_cm_matrix;

}

Eigen::MatrixXd 
TransformedSector(
    const basis::TwoBodySubspaceNLSJT& two_body_subspace2,
    const basis::TwoBodySubspaceNLSJT& two_body_subspace1,
    const basis::RelativeSpaceLSJT& relative_space,
    const basis::RelativeSectorsLSJT& relative_sectors,
    const basis::MatrixVector& relative_matrices,
    int J0
  )
{
  // define target matrix
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(two_body_subspace2.size(),two_body_subspace1.size());

  // extract subspace labels
  int L2 = two_body_subspace2.L();
  int S2 = two_body_subspace2.S();
  int J2 = two_body_subspace2.J();
  int T2 = two_body_subspace2.T();
  int g2 = two_body_subspace2.g();

  int L1 = two_body_subspace1.L();
  int S1 = two_body_subspace1.S();
  int J1 = two_body_subspace1.J();
  int T1 = two_body_subspace1.T();
  int g1 = two_body_subspace1.g();


  // obtain N truncations
  //
  // spurious spaces contribute for Nmax only up to lesser of the bra
  // and ket two-body space Nmax values
//  int Nmax2 = two_body_subspace2.Nmax();
//  int Nmax1 = two_body_subspace1.Nmax();
//  int Nc_max = std::min(Nmax2,Nmax1); 

//  // accumulate cm contributions
//  for (int Nc=0; Nc <= Nc_max; ++Nc)
//    for (int lc = Nc%2; lc <= Nc; lc+=2)
//      {
//        // define intermediate spurious relative subspaces
//        const basis::RelativeCMSubspaceLSJT relative_cm_subspace2(Nc,lc,L2,S2,J2,T2,g2,Nmax2);
//        const basis::RelativeCMSubspaceLSJT relative_cm_subspace1(Nc,lc,L1,S1,J1,T1,g1,Nmax1);
//
//        // check that spurious relative subspace is nonempty
//        if (relative_cm_subspace2.size()==0)
//          continue;
//        if (relative_cm_subspace1.size()==0)
//          continue;
//
//        // construct Ar_spurious
//        // PLACEHOLDER: Eigen::MatrixXd As = Eigen::MatrixXd::Zero(relative_cm_subspace2.size(),relative_cm_subspace1.size());
//        Eigen::MatrixXd As = RelativeCMMatrix(
//            relative_cm_subspace2,relative_cm_subspace1,
//            relative_space,relative_sectors,relative_matrices,J0
//          );
//        // std::cout << "As " << As << std::endl;
//
//        // construct M2
//        // PLACEHOLDER: Eigen::MatrixXd M2 = Eigen::MatrixXd::Zero(relative_cm_subspace2.size(),two_body_subspace2.size());
//        // TO UPDATE: Eigen::MatrixXd M2 = MoshinskyMatrix(relative_cm_subspace2,two_body_subspace2);
//        // std::cout << "M2 " << M2 << std::endl;
//
//        // construct M1
//        // PLACEHOLDER: Eigen::MatrixXd M1 = Eigen::MatrixXd::Zero(relative_cm_subspace1.size(),two_body_subspace1.size());
//        // TO UPDATE: Eigen::MatrixXd M1 = MoshinskyMatrix(relative_cm_subspace1,two_body_subspace1);
//        // std::cout << "M1 " << M1 << std::endl;
//      
//        // generate contribution to A
//        Eigen::MatrixXd A_contribution = 2 * M2.transpose() * As * M1;
//
//        A += A_contribution;
//      }


  return A;

}


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
