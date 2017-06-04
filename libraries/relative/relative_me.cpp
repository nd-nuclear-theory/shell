/****************************************************************

  relative_me.cpp

  Mark A. Caprio, University of Notre Dame.

****************************************************************/

#include "relative_me.h"

#include "spline/wavefunction_class.h"

namespace relative {

  void ConstructKinematicOperator(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::MatrixVector,3>& relative_component_matrices,
      relative::KinematicOperator kinematic_operator
    )
  {

    // zero initialize operators
    basis::ConstructZeroOperatorRelativeLSJT(
        operator_labels,relative_space,relative_component_sectors,relative_component_matrices
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

  void ConstructCoulombOperator(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::MatrixVector,3>& relative_component_matrices,
      int num_steps
    )
  {

    // zero initialize operator
    basis::ConstructZeroOperatorRelativeLSJT(
        operator_labels,relative_space,relative_component_sectors,relative_component_matrices
      );

    for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
      // for each isospin component
      {
        // select T0 component
        const basis::RelativeSectorsLSJT& sectors = relative_component_sectors[T0];
        basis::MatrixVector& matrices = relative_component_matrices[T0];

        // iterate over sectors
        for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
          {

            // extract sector
            const basis::RelativeSectorsLSJT::SectorType& sector = sectors.GetSector(sector_index);


            // short-circuit select diagonal sectors
            //
            // The a.m. scalar nature imposes diagonal on (L,S,J),
            // positive parity imposes diagonal on (g), and then there
            // is that apparent parity selection imposing diagonal in
            // T.
            if (!sector.IsDiagonal())
              continue;

            // extract subspace and labels
            const basis::RelativeSubspaceLSJT& subspace = sector.ket_subspace();
            int Nmax = subspace.Nmax();
            int L = subspace.L();
            int S = subspace.S();
            int J = subspace.J();
            int T = subspace.T();
            int nmax = (Nmax-L)/2;   // max radial quantum number

            // short circuit-select T=1 sectors (for proton-only operator)


            // Although actually we could get nmax just from the
            // subspace dimension...
            assert(nmax==subspace.size()-1);

            // populate nonzero entries
            //
            // We make use of the known indexing scheme for a
            // RelativeLSJT basis, that the radial quantum number n is
            // just the 0-based state index.

            Eigen::MatrixXd& matrix = matrices[sector_index];

            double isospin_factor;
            if (T0==0)
              isospin_factor = 1/3.;
            else if (T0==1)
              isospin_factor = 1/std::sqrt(2.);
            else if (T0==2)
              isospin_factor = std::sqrt(10.)/6.;

            for (int bra_n=0; bra_n<=nmax; ++bra_n)
              for (int ket_n=0; ket_n<=nmax; ++ket_n)
                {

                  // get bra and ket states
                  spline::WaveFunction bra_wavefunction(bra_n,L,1,spline::Basis::HC);
                  spline::WaveFunction ket_wavefunction(ket_n,L,1,spline::Basis::HC);
            
                  // evaluate radial integral
                  const int num_size = num_steps+1;
                  double radial_integral = bra_wavefunction.MatrixElement(num_size,ket_wavefunction,-1);

                  // impose isospin factors
                  matrix(bra_n,ket_n) = isospin_factor * radial_integral;
                }


          }
      }

  }


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
