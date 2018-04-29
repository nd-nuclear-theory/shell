/****************************************************************

  relative_me.cpp

  Mark A. Caprio, University of Notre Dame.

****************************************************************/


#include "relative_me.h"

#include "am/racah_reduction.h"
#include "am/wigner_gsl.h"
#include "cppformat/format.h"
#include "spline/wavefunction_class.h"

namespace relative {

  // isospin coefficients for operator which can see only protons,
  // only neutrons
  //
  // outer index: operator type (OperatorTypePN) = kPP=0,kNN=1,kTotal=2 (short circuited)
  // inner index: isospin component T0 = 0,1,2
  std::array<std::array<double,3>,2>
  kSingleSpeciesOperatorIsospinFactors({
      std::array<double,3>({
        1/3.,1/std::sqrt(2.),std::sqrt(10.)/6.
          }),
        std::array<double,3>({
        1/3.,-1/std::sqrt(2.),std::sqrt(10.)/6.
          })
        });

  void ConstructKinematicOperator(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
      relative::KinematicOperator kinematic_operator
    )
  {
    // validate operator parameters
    assert(operator_labels.J0==0);
    assert(operator_labels.g0==0);

    // select phase for coordinate or momentum space
    int operator_sign;
    if (kinematic_operator == relative::KinematicOperator::kRSqr)
      operator_sign = +1;
    else if (kinematic_operator == relative::KinematicOperator::kKSqr)
      operator_sign = -1;

    // zero initialize operators
    basis::ConstructZeroOperatorRelativeLSJT(
        operator_labels,relative_space,relative_component_sectors,relative_component_matrices
      );

    // select T0=0 component
    const basis::RelativeSectorsLSJT& sectors = relative_component_sectors[0];
    basis::OperatorBlocks<double>& matrices = relative_component_matrices[0];

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

        // relative coordinate dilation factor
        //
        // See "Note on oscillator length" at start of header file.

        double relative_oscillator_scale_factor;
        if (kinematic_operator == relative::KinematicOperator::kRSqr)
          relative_oscillator_scale_factor = 2.;
        else if (kinematic_operator == relative::KinematicOperator::kKSqr)
          relative_oscillator_scale_factor = 1/2.;

        // populate nonzero entries
        //
        // We make use of the known indexing scheme for a
        // RelativeLSJT basis, that the radial quantum number n is
        // just the 0-based state index.
        //
        // While Rowe uses "positive at infinity" convention
        // for the radial wave functions, we use "positive
        // at origin" convention.  Conversion introduces a
        // factor (-)^(bra_n+ket_n), i.e., adding a (-) sign
        // on the bra.n()=n+1 or n-1 terms.

        basis::OperatorBlock<double>& matrix = matrices[sector_index];
        for (int n=0; n<=nmax; ++n)
          {
            if (n>0)
              matrix(n-1,n) = -operator_sign*relative_oscillator_scale_factor*sqrt(n*(n+L+0.5));
            matrix(n,n) = relative_oscillator_scale_factor*(2*n+L+1.5);
            if (n<nmax)
              matrix(n+1,n) = -operator_sign*relative_oscillator_scale_factor*sqrt((n+1)*(n+L+1.5));
	  }
      }
  }

  void ConstructQuadrupoleOperator(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
      relative::KinematicOperator kinematic_operator,
      basis::OperatorTypePN operator_type
    )
  {

    // validate operator parameters
    assert(operator_labels.J0==2);
    assert(operator_labels.g0==0);

    assert(operator_type == basis::OperatorTypePN::kTotal);


    // select phase for coordinate or momentum space
    // TODO: check phase for momentum space version and insert into appropriate terms below
    assert(kinematic_operator == relative::KinematicOperator::kRSqr);
    int operator_sign;
    if (kinematic_operator == relative::KinematicOperator::kRSqr)
      operator_sign = +1;
    else if (kinematic_operator == relative::KinematicOperator::kKSqr)
      operator_sign = -1;

    // zero initialize operator
    basis::ConstructZeroOperatorRelativeLSJT(
        operator_labels,relative_space,relative_component_sectors,relative_component_matrices
      );

 
    for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
      // for each isospin component
      {

        // select T0 component
        const basis::RelativeSectorsLSJT& sectors = relative_component_sectors[T0];
        basis::OperatorBlocks<double>& matrices = relative_component_matrices[T0];

        // only implemented for total operator so far
        // double isospin_factor = kSingleSpeciesOperatorIsospinFactors[int(operator_type)][T0];
        if (T0!=0)
          continue;

        // iterate over sectors
        for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
          {

            // set us aliases -- for sector and subspaces
            const basis::RelativeSectorsLSJT::SectorType& sector = sectors.GetSector(sector_index);
            const basis::RelativeSubspaceLSJT& bra_subspace = sector.bra_subspace();
            const basis::RelativeSubspaceLSJT& ket_subspace = sector.ket_subspace();

            // short-circuit select allowed sectors
            //
            // Quadrupole (J0=2) and positive parity (g0=0) selection rules
            // should already have been enforced in sector construction, but
            // here we also impose orbital quadrupole (L0=2) and spin scalar
            // selection (S0=0).
            if (!(
                    am::AllowedTriangle(bra_subspace.L(),2,ket_subspace.L())
                    && am::AllowedTriangle(bra_subspace.S(),0,ket_subspace.S())
                  )
              )
              continue;

            // extract subspaces
            // int Nmax = subspace.Nmax();
            // int L = subspace.L();
            // int S = subspace.S();
            // int J = subspace.J();
            // int T = subspace.T();
            // int nmax = (Nmax-L)/2;   // max radial quantum number

            // determine angular factor
            //
            // We evaluate <L'S'L'||Y_2||LSJ> (S'=S) using the known
            // <L'||Y_2||L> and Racah's two-system reduction formula for the
            // case where one operator is the identity.  This is given for the
            // special case of j orbitals (s=1/2) in Suhonen (2.57), but here we
            // have must generalize to S=0,1.
            //
            // We then convert the RME normalization from Edmonds convention to
            // Rose convention.
            //
            // TODO: recode more cleanly all in Rose convention and using
            // am::RacahReductionFactor1Rose

            int bra_L = bra_subspace.L();
            int ket_L = ket_subspace.L();
            int S = ket_subspace.S();
            int bra_J = bra_subspace.J();
            int ket_J = ket_subspace.J();
            double angular_factor
              = std::sqrt(5./(4*M_PI))
              *ParitySign(bra_L+ket_L+S+ket_J)
              *Hat(bra_J)*Hat(ket_J)*Hat(bra_L)*Hat(ket_L)
              *am::Wigner3J(ket_L,2,bra_L,0,0,0)
              *am::Wigner6J(bra_L,bra_J,S,ket_J,ket_L,2);
            angular_factor /= Hat(bra_J);  // convert to Rose convention

            // alias to matrix
            basis::OperatorBlock<double>& matrix = matrices[sector_index];

            // populate nonzero entries
            //
            // Restrict to a tri-diagonal loop over radial labels.
            //
            // We make use of the known indexing scheme for a RelativeLSJT
            // basis, that the radial quantum number n is just the 0-based state
            // index.

            int bra_subspace_size = bra_subspace.size();
            int ket_subspace_size = ket_subspace.size();
            for (int ket_n=0; ket_n<ket_subspace_size; ++ket_n)
              {
                // determine bounds on bra indices by tridiagonal nature of r^2 matrices
                int bra_n_min, bra_n_max;
                if (bra_subspace.L()==ket_subspace.L()-2)
                  {
                    bra_n_min=ket_n;
                    bra_n_max=ket_n+2;
                  }
                else if (bra_subspace.L()==ket_subspace.L())
                  {
                    bra_n_min=ket_n-1;
                    bra_n_max=ket_n+1;
                  }
                else if (bra_subspace.L()==ket_subspace.L()+2)
                  {
                    bra_n_min=ket_n-1;
                    bra_n_max=ket_n;
                  }
                bra_n_min=std::max(bra_n_min,0);
                bra_n_max=std::min(bra_n_max,bra_subspace_size-1);

                // restricted loop over bra indices
                for (int bra_n=bra_n_min; bra_n<=bra_n_max; ++bra_n)
                  {

                  const basis::RelativeStateLSJT bra(bra_subspace,bra_n);
                  const basis::RelativeStateLSJT ket(ket_subspace,ket_n);

                  // alias ket quantum numbers
                  //   since formulas refer to ket quantum numbers
                  int L = ket.L();
                  int n = ket.n();
                  assert(n==ket_n);

                  // radial matrix element
                  //
                  // Apply SU(1,1) radial matrix element formulas from Rowe JPA
                  // 38, 10181 (2015):
                  //
                  // - for (delta L)=0, use (delta l)=0 matrix elements of r^2
                  // - from (39)
                  //
                  // - for (delta L)=2, use resolution of identity over
                  //   intermediate L space, by double application of the (delta
                  //   l)=1 matrix elements of r from (64)&(65)
                  //
                  // While Rowe uses "positive at infinity" convention for the
                  // radial wave functions, we use "positive at origin"
                  // convention.  Conversion introduces a factor
                  // (-)^(bra_n+ket_n), i.e., adding a (-) sign on the
                  // bra.n()=n+1 or n-1 terms.

                  // validation of radial formulas: cross check with numerical
                  // evaluation
                  //
                  // Carried out with unrestricted loop over bra_n to
                  // numerically check zero entries as well.
                  //
                  //   for (int bra_n=0; bra_n<bra_subspace_size; ++bra_n)
                  //
                  // (delta L)=-2
                  // ...
                  //  bra L 0 n 1 ; ket L 2 n 9 : analytic      0.00000000 numerical     -0.00000000
                  //  bra L 0 n 2 ; ket L 2 n 0 : analytic      1.41421356 numerical      1.41421356
                  //  bra L 0 n 2 ; ket L 2 n 1 : analytic     -5.29150262 numerical     -5.29150262
                  //  bra L 0 n 2 ; ket L 2 n 2 : analytic      3.96862697 numerical      3.96862697
                  //  bra L 0 n 2 ; ket L 2 n 3 : analytic      0.00000000 numerical      0.00000000
                  // ...
                  // (delta L)=0
                  // ...
                  //  bra L 1 n 0 ; ket L 1 n 9 : analytic      0.00000000 numerical      0.00000000
                  //  bra L 1 n 1 ; ket L 1 n 0 : analytic     -1.58113883 numerical     -1.58113883
                  //  bra L 1 n 1 ; ket L 1 n 1 : analytic      4.50000000 numerical      4.50000000
                  //  bra L 1 n 1 ; ket L 1 n 2 : analytic     -2.64575131 numerical     -2.64575131
                  //  bra L 1 n 1 ; ket L 1 n 3 : analytic      0.00000000 numerical      0.00000000
                  // ...
                  // (delta L)=+2
                  //
                  // This case does not appear naturally for canonical ordering
                  // of sectors, but it wast forced for testing purposes:
                  //
                  //   // debugging: override canonical sector ordering
                  //   relative_component_sectors[T0]
                  //     = basis::RelativeSectorsLSJT(relative_space,operator_labels.J0,T0,operator_labels.g0,basis::SectorDirection::kBoth);
                  //
                  // ...
                  // bra L 2 n 1 ; ket L 0 n 0 : analytic      0.00000000 numerical     -0.00000000
                  // bra L 2 n 1 ; ket L 0 n 1 : analytic      2.95803989 numerical      2.95803989
                  // bra L 2 n 1 ; ket L 0 n 2 : analytic     -5.29150262 numerical     -5.29150262
                  // bra L 2 n 1 ; ket L 0 n 3 : analytic      2.44948974 numerical      2.44948974
                  // bra L 2 n 1 ; ket L 0 n 4 : analytic      0.00000000 numerical      0.00000000
                  // ...

                  double radial_integral = 0.;
                  if (bra.L()==L-2)
                    // (delta L) = -2
                    {
                      if (bra.n()==n)
                        radial_integral = std::sqrt((L+n-0.5)*(L+n+0.5));
                      else if (bra.n()==n+1)
                        radial_integral = -2*std::sqrt((n+1)*(L+n+0.5));
                      else if (bra.n()==n+2)
                        radial_integral = std::sqrt((n+2)*(n+1));
                    }
                  else if (bra.L()==L)
                    // (delta L) = 0
                    {
                      if (bra.n()==n-1)
                        radial_integral = -std::sqrt((L+n+0.5)*n);
                      else if (bra.n()==n)
                        radial_integral = (L+2*n+1.5);
                      else if (bra.n()==n+1)
                        radial_integral = -std::sqrt((L+n+1.5)*(n+1));
                    }
                  else if (bra.L()==L+2)
                    // (delta L) = +2
                    {
                      if (bra.n()==n-2)
                        radial_integral = std::sqrt((n-1)*n);
                      else if (bra.n()==n-1)
                        radial_integral = -2*std::sqrt(n*(L+n+1.5));
                      else if (bra.n()==n)
                        radial_integral = std::sqrt((L+n+2.5)*(L+n+1.5));
                    }

                  // // numerical validation code
                  // spline::WaveFunction bra_wavefunction(bra_n,bra.L(),1,spline::Basis::HC);
                  // spline::WaveFunction ket_wavefunction(ket_n,ket.L(),1,spline::Basis::HC);
                  // const int num_steps = 500;
                  // const int num_size = num_steps+1;
                  // double radial_integral_numerical = bra_wavefunction.MatrixElement(num_size,ket_wavefunction,2);
                  // std::cout
                  //   << fmt::format(
                  //       " bra L {:d} n {:d} ; ket L {:d} n {:d} : analytic {:15.8f} numerical {:15.8f}",
                  //       bra.L(),bra.n(),ket.L(),ket.n(),radial_integral,radial_integral_numerical
                  //     )
                  //   << std::endl;

                  // relative coordinate dilation factor
                  //
                  // See "Note on oscillator length" at start of header file.

                  const double relative_oscillator_scale_factor = 2.;

                  // combine factors
                  matrix(bra_n,ket_n) = angular_factor * relative_oscillator_scale_factor * radial_integral;
                }

              }
          }
      }
  }

  void ConstructOrbitalAMOperator(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
      basis::OperatorTypePN operator_type
    )
  {

    // validate operator parameters
    assert(operator_labels.J0==1);
    assert(operator_labels.g0==0);

    assert(operator_type == basis::OperatorTypePN::kTotal);

    // zero initialize operator
    basis::ConstructZeroOperatorRelativeLSJT(
        operator_labels,relative_space,relative_component_sectors,relative_component_matrices
      );

    for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
      // for each isospin component
      {

        // select T0 component
        const basis::RelativeSectorsLSJT& sectors = relative_component_sectors[T0];
        basis::OperatorBlocks<double>& matrices = relative_component_matrices[T0];

        // only implemented for total operator so far
        // double isospin_factor = kSingleSpeciesOperatorIsospinFactors[int(operator_type)][T0];
        if (T0!=0)
          continue;

        // iterate over sectors
        for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
          {

            // set us aliases -- for sector and subspaces
            const basis::RelativeSectorsLSJT::SectorType& sector = sectors.GetSector(sector_index);
            const basis::RelativeSubspaceLSJT& bra_subspace = sector.bra_subspace();
            const basis::RelativeSubspaceLSJT& ket_subspace = sector.ket_subspace();

            // short-circuit select allowed sectors
            //
            // Dipole (J0=1) and positive parity (g0=0) selection rules should
            // already have been enforced in sector construction, but here we
            // also impose orbital angular momentum conservation by the
            // generator L (L'=L) and spin scalar selection (S0=0).
            if (!(
                    (bra_subspace.L()==ket_subspace.L())
                    && am::AllowedTriangle(bra_subspace.S(),0,ket_subspace.S())
                  )
              )
              continue;

            // determine angular factor
            //
            // Obtained from Racah's two-system reduction formula and the basic
            // RME for the angular momentum generator.  See "intrinsic dipole
            // operators" notes (p. 2).

            int L = ket_subspace.L();
            int S = ket_subspace.S();
            int bra_J = bra_subspace.J();
            int ket_J = ket_subspace.J();

            //TODO
            double angular_factor
              = am::RacahReductionFactor1Rose(L,S,bra_J,L,S,ket_J,1)
              *std::sqrt(L*(L+1));

            // alias to matrix
            basis::OperatorBlock<double>& matrix = matrices[sector_index];

            // populate nonzero entries
            //
            // We make use of the known indexing scheme for a
            // RelativeLSJT basis, that the radial quantum number n is
            // just the 0-based state index.  
            //
            // The (LSJ') and (LSJ) subspaces are the same size, so this is a
            // square block, and the block is simply proportional to the
            // identity matrix.

            int bra_subspace_size = bra_subspace.size();
            int ket_subspace_size = ket_subspace.size();
            assert(bra_subspace_size==ket_subspace_size);  // (LSJ') and (LSJ) subspaces same size
            // matrix = angular_factor*basis::OperatorBlock<double>:Identity(bra_subspace.size(),ket_subspace.size());
            for (int ket_n=0; ket_n<ket_subspace_size; ++ket_n)
              {
                int bra_n = ket_n;
                matrix(bra_n,ket_n) = angular_factor;
              }

          }
      }
  }

  void ConstructSpinAMOperator(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
      basis::OperatorTypePN operator_type
    )
  {

    // validate operator parameters
    assert(operator_labels.J0==1);
    assert(operator_labels.g0==0);

    assert(operator_type == basis::OperatorTypePN::kTotal);

    // zero initialize operator
    basis::ConstructZeroOperatorRelativeLSJT(
        operator_labels,relative_space,relative_component_sectors,relative_component_matrices
      );

    for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
      // for each isospin component
      {

        // select T0 component
        const basis::RelativeSectorsLSJT& sectors = relative_component_sectors[T0];
        basis::OperatorBlocks<double>& matrices = relative_component_matrices[T0];

        // only implemented for total operator so far
        // double isospin_factor = kSingleSpeciesOperatorIsospinFactors[int(operator_type)][T0];
        if (T0!=0)
          continue;

        // iterate over sectors
        for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
          {

            // set us aliases -- for sector and subspaces
            const basis::RelativeSectorsLSJT::SectorType& sector = sectors.GetSector(sector_index);
            const basis::RelativeSubspaceLSJT& bra_subspace = sector.bra_subspace();
            const basis::RelativeSubspaceLSJT& ket_subspace = sector.ket_subspace();

            // short-circuit select allowed sectors
            //
            // Dipole (J0=1) and positive parity (g0=0) selection rules should
            // already have been enforced in sector construction, but here we
            // also impose spin angular momentum conservation by the generator S
            // (S'=S) and orbital angular momentum scalar selection (L0=0).
            if (!(
                    (bra_subspace.S()==ket_subspace.S())
                    && am::AllowedTriangle(bra_subspace.L(),0,ket_subspace.L())
                  )
              )
              continue;

            // determine angular factor
            //
            // Obtained from Racah's two-system reduction formula and the basic
            // RME for the angular momentum generator.  See "intrinsic dipole
            // operators" notes (p. 2).

            int L = ket_subspace.L();
            int S = ket_subspace.S();
            int bra_J = bra_subspace.J();
            int ket_J = ket_subspace.J();

            //TODO
            double angular_factor
              = am::RacahReductionFactor2Rose(L,S,bra_J,L,S,ket_J,1)
              *std::sqrt(S*(S+1));

            // alias to matrix
            basis::OperatorBlock<double>& matrix = matrices[sector_index];

            // populate nonzero entries
            //
            // We make use of the known indexing scheme for a
            // RelativeLSJT basis, that the radial quantum number n is
            // just the 0-based state index.  
            //
            // The (LSJ') and (LSJ) subspaces are the same size, so this is a
            // square block, and the block is simply proportional to the
            // identity matrix.

            int bra_subspace_size = bra_subspace.size();
            int ket_subspace_size = ket_subspace.size();
            assert(bra_subspace_size==ket_subspace_size);  // (LSJ') and (LSJ) subspaces same size
            // matrix = angular_factor*basis::OperatorBlock<double>:Identity(bra_subspace.size(),ket_subspace.size());
            for (int ket_n=0; ket_n<ket_subspace_size; ++ket_n)
              {
                int bra_n = ket_n;
                matrix(bra_n,ket_n) = angular_factor;
              }

          }
      }
  }
  

  void ConstructCoulombOperator(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
      basis::OperatorTypePN operator_type,
      int num_steps
    )
  {

    // validate operator parameters
    assert(operator_labels.J0==0);
    assert(operator_labels.g0==0);

    // zero initialize operator
    basis::ConstructZeroOperatorRelativeLSJT(
        operator_labels,relative_space,relative_component_sectors,relative_component_matrices
      );

    for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
      // for each isospin component
      {
        // select T0 component
        const basis::RelativeSectorsLSJT& sectors = relative_component_sectors[T0];
        basis::OperatorBlocks<double>& matrices = relative_component_matrices[T0];

        // iterate over sectors
        for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
          {

            // extract sector
            const basis::RelativeSectorsLSJT::SectorType& sector = sectors.GetSector(sector_index);


            // short-circuit select diagonal sectors
            //
            // The a.m. scalar nature imposes diagonal on (L,S,J),
            // positive parity imposes diagonal on (g), and then
            // consequently T is also fixed to be the same on both
            // sides.
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

            // although actually we could get nmax just from the
            // subspace dimension...
            assert(nmax==subspace.size()-1);

            // populate nonzero entries
            //
            // We make use of the known indexing scheme for a
            // RelativeLSJT basis, that the radial quantum number n is
            // just the 0-based state index.

            basis::OperatorBlock<double>& matrix = matrices[sector_index];

            // calculate isospin factors
            //   and short-circuit evaluation of non-contributing sectors
            double isospin_factor;
            if (operator_type==basis::OperatorTypePN::kTotal)
              {
                if (T0==0)
                  isospin_factor = 1.;
                else
                  continue;
              }
            else
              {
                if (T==1)
                  isospin_factor = kSingleSpeciesOperatorIsospinFactors[int(operator_type)][T0];
                else
                  continue;
              }
             
            for (int bra_n=0; bra_n<=nmax; ++bra_n)
              for (int ket_n=0; ket_n<=nmax; ++ket_n)
                {

                  // get bra and ket states
                  spline::WaveFunction bra_wavefunction(bra_n,L,1,spline::Basis::HC);
                  spline::WaveFunction ket_wavefunction(ket_n,L,1,spline::Basis::HC);
            
                  // evaluate radial integral
                  const int num_size = num_steps+1;
                  double radial_integral = bra_wavefunction.MatrixElement(num_size,ket_wavefunction,-1);

                  // relative coordinate dilation factor
                  //
                  // See "Note on oscillator length" at start of header file.

                  const double relative_oscillator_scale_factor = 1/std::sqrt(2.);

                  // impose isospin factors
                  matrix(bra_n,ket_n) = isospin_factor * relative_oscillator_scale_factor * radial_integral;
                }

          }
      }

  }


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
