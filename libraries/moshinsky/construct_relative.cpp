/****************************************************************

  construct_relative.cpp

  Mark A. Caprio, University of Notre Dame.

****************************************************************/

#include "construct_relative.h"

namespace relative {

  void ConstructDiagonalConstantOperator(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::MatrixVector,3>& relative_component_matrices,
      double c
    )
  {

    // validate operator labels
    assert(operator_labels.J0==0);
    assert(operator_labels.g0==0);
    assert(operator_labels.symmetry_phase_mode==basis::SymmetryPhaseMode::kHermitian);
    assert(operator_labels.T0_min==0);

    for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
      // for each isospin component
      {

        // enumerate sectors
        relative_component_sectors[T0]
          = basis::RelativeSectorsLSJT(relative_space,operator_labels.J0,T0,operator_labels.g0);
         
        // populate matrices
        relative_component_matrices[T0].resize(relative_component_sectors[T0].size());
        if (T0==0)
          // identity matrices in
          basis::SetOperatorToDiagonalConstant(relative_component_sectors[T0],relative_component_matrices[T0],c);
        else
          basis::SetOperatorToZero(relative_component_sectors[T0],relative_component_matrices[T0]);
      }
  }


  void ConstructKinematicOperator(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::MatrixVector,3>& relative_component_matrices,
      relative::KinematicOperator kinematic_operator
    )
  {

    // zero initialize operator
    ConstructDiagonalConstantOperator(
        operator_labels,relative_space,relative_component_sectors,relative_component_matrices,0.
      );

    // select T0=0 component
    const basis::RelativeSectorsLSJT& sectors = relative_component_sectors[0];
    basis::MatrixVector& matrices = relative_component_matrices[0];

    // iterate over sectors
    for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
      {

        // extract sector
	const basis::RelativeSectorsLSJT::SectorType& sector = sectors.GetSector(sector_index);

        // short-circuit select diagonal sectors
	if (!sector.IsDiagonal())
          continue;

        // extract subspace and labels
        const basis::RelativeSubspaceLSJT& subspace = sector.ket_subspace();
        int Nmax = subspace.Nmax();
        int L = subspace.L();
        int nmax = (Nmax-L)/2;   // max radial quantum number
        std::cout << " Nmax " << Nmax << " L " << L << " nmax " << nmax << " size " << subspace.size() << std::endl;
        std::cout << subspace.LabelStr() << std::endl;
        std::cout << subspace.DebugStr();

        // Although actually we could get nmax just from the
        // subspace dimension...
        assert(nmax==subspace.size()-1);

        // select which kinematic operator
        //
        // Note: Signs on off-diagonal terms are for "positive at
        // origin" convention, would be reversed for "positive at
        // infinity" convention.
        int operator_sign;
        if (kinematic_operator == relative::KinematicOperator::kRSqr)
          operator_sign = -1;
        else if (kinematic_operator == relative::KinematicOperator::kKSqr)
          operator_sign = +1;

        // populate nonzero entries
        //
        // We make use of the known indexing scheme for a
        // RelativeLSJT basis, that the radial quantum number n is
        // just the 0-based state index.

        Eigen::MatrixXd& matrix = matrices[sector_index];
        for (int n=0; n<=nmax; ++n)
          {
            if (n>0)
              matrix(n-1,n) = operator_sign*sqrt(n*(n+L+0.5));
            matrix(n,n) = (2*n+L+1.5);
            if (n<nmax)
              matrix(n+1,n) = operator_sign*sqrt((n+1)*(n+L+1.5));
	  }
      }
  }

  // std::vector<Eigen::MatrixXd>
  // ImportInteraction_Kinetic(const basis::RelativeSpaceLSJT& space,const basis::RelativeSectorsLSJT& sectors)
  // {
  //   int Nmax=40;
  //   std::vector<Eigen::MatrixXd> sector_vector(sectors.size());
  // 
  //   for(int i=0; i<sectors.size(); ++i)
  //   {
  //     
  //     const basis::RelativeSubspaceLSJT& bra_subspace = sectors.GetSector(i).bra_subspace();
  //     const basis::RelativeSubspaceLSJT& ket_subspace = sectors.GetSector(i).ket_subspace();
  //     // Extract ket labels 
  //     int L=ket_subspace.L();
  //     int S=ket_subspace.S();
  //     int J=ket_subspace.J();
  //     int T=ket_subspace.T();
  //     int g=ket_subspace.g();
  //     // Extract bra labels
  //     int Lp=bra_subspace.L();
  //     int Sp=bra_subspace.S();
  //     int Jp=bra_subspace.J();
  //     int Tp=bra_subspace.T();
  //     int gp=bra_subspace.g();
  // 
  //     int npmax=(Nmax-Lp)/2;
  //     int nmax=(Nmax-L)/2;
  //     sector_vector[i]=Eigen::MatrixXd::Constant(nmax+1, nmax+1, 0);
  // 
  //     if((Lp==L)&&(Sp==S)&&(Jp==J)&&(Tp==T)&&(gp==g))
  //       {
  //         for(int n=0; n<=nmax; ++n)
  //           {
  //             if((0<=(n-1))&&((n-1)<=npmax))
  //               sector_vector[i](n-1,n)=-sqrt(n*(n+L+.5))*10;
  //             if(n<=npmax)
  //               sector_vector[i](n,n)=(2*n+L+1.5)*10;
  //             if((0<=(n+1))&&((n+1)<=npmax))
  //               sector_vector[i](n+1,n)=-sqrt((n+1)*(n+L+1.5))*10;
  //           }
  //       }
  //   }
  //   return sector_vector;
  // }
  // 
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
