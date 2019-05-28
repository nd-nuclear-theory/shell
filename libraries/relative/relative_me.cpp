/****************************************************************

  relative_me.cpp

  Mark A. Caprio, University of Notre Dame.

****************************************************************/


#include "relative_me.h"

#include "am/racah_reduction.h"
#include "am/wigner_gsl.h"
#include "fmt/format.h"
#include "analytic/radial_oscillator_me.h"
#include "spline/spline_me.h"

namespace relative {

  ////////////////////////////////////////////////////////////////
  // oscillator radial matrix
  ////////////////////////////////////////////////////////////////


  basis::OperatorBlock<double> RadialCoordinateSqrMatrix(
      int bra_dim, int ket_dim,
      int L, int delta_L,
      int operator_sign=+1
    )
  // Construct radial matrices for r^2 and k^2 operators in oscillator radial
  // basis.
  //
  // Arguments:
  //   bra_dim, ket_dim (input): dimensions of radial subspaces
  //   L (input): ket radial orbital angular momentum quantum number
  //   delta_L (input): bra radial orbital angular momentum quantum
  //     number relative to ket
  //   operator_sign (int): sign selecting coordinate or momentum
  //     representation (+1 for "r", -1 for "k")
  //
  // Returns:
  //   radial matrix
  {

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
    // of sectors, but it was forced for testing purposes:
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

    // TODO: need to check signs for L-changing k^2 matrix elements
    assert((delta_L==0)||(operator_sign==+1));
    const int ket_L = L, bra_L = L+delta_L;

    // set up matrix
    basis::OperatorBlock<double> matrix = basis::OperatorBlock<double>::Zero(bra_dim,ket_dim);

    // populate nonzero entries
    for (int ket_n=0; ket_n<ket_dim; ++ket_n)
      {
        const int ket_N = 2*ket_n + ket_L;

        // restricted loop over bra indices
        int bra_n_min, bra_n_max;
        if (delta_L==-2)
          {
            bra_n_min=ket_n;
            bra_n_max=ket_n+2;
          }
        else if (delta_L==0)
          {
            bra_n_min=ket_n-1;
            bra_n_max=ket_n+1;
          }
        else if (delta_L==+2)
          {
            bra_n_min=ket_n-1;
            bra_n_max=ket_n;
          }
        bra_n_min=std::max(bra_n_min,0);
        bra_n_max=std::min(bra_n_max,bra_dim-1);
        for (int bra_n=bra_n_min; bra_n<=bra_n_max; ++bra_n)
          {
            const int bra_N = 2*bra_n + bra_L;
            double matrix_element =
              analytic::CoordinateSqrOscillatorMatrixElement(
                bra_N, bra_L, ket_N, ket_L, operator_sign
              );
            // // numerical validation code
            // spline::WaveFunction bra_wavefunction(bra_n,bra.L(),1,spline::Basis::HC);
            // spline::WaveFunction ket_wavefunction(ket_n,ket.L(),1,spline::Basis::HC);
            // const int num_steps = 500;
            // const int num_size = num_steps+1;
            // double matrix_element_numerical = bra_wavefunction.MatrixElement(num_size,ket_wavefunction,2);
            // std::cout
            //   << fmt::format(
            //       " bra L {:d} n {:d} ; ket L {:d} n {:d} : analytic {:15.8f} numerical {:15.8f}",
            //       bra.L(),bra.n(),ket.L(),ket.n(),matrix_element,matrix_element_numerical
            //     )
            //   << std::endl;

            matrix(bra_n,ket_n) = matrix_element;
          }
      }

    return matrix;
  }

  basis::OperatorBlock<double> RadialCoordinateMatrix(
      int bra_dim, int ket_dim,
      int L, int delta_L,
      int operator_sign=+1
    )
  // Construct radial matrix for r operator in oscillator radial basis.
  //
  // Calculated from SU(1,1) algebraic expressions in (64) & (65) of D. J. Rowe,
  // JPA 38, 10181 (2005), but with phase converted to "positive at origin"
  // convention on the radial wave functions.  We use lambda=v+N/2 for
  // oscillator functions, with N=3 and v=l in three dimensions.  Radial matrix
  // elements must be complemented with angular matrix elements given by Rowe
  // (97).
  //
  // While Rowe uses "positive at infinity" convention for the radial wave
  // functions, we commonly use "positive at origin" convention.  Conversion
  // introduces a factor (-)^(bra_n+ket_n), i.e., adding a (-) sign on the
  // bra_n=ket_n+1 or ket_n-1 terms.
  //
  // Arguments:
  //   bra_dim, ket_dim (input): dimensions of radial subspaces
  //   L (input): ket radial orbital angular momentum quantum number
  //   delta_L (input): bra radial orbital angular momentum quantum
  //     number relative to ket
  //
  // Returns:
  //   radial matrix
  {

    // validate arguments
    assert((delta_L==-1)||(delta_L==+1));
    const int ket_L = L, bra_L = L+delta_L;

    // set up matrix
    basis::OperatorBlock<double> matrix = basis::OperatorBlock<double>::Zero(bra_dim,ket_dim);

    // populate nonzero entries
    for (int ket_n=0; ket_n<ket_dim; ++ket_n)
      {
        const int ket_N = 2*ket_n + ket_L;

        // restricted loop over bra indices
        int bra_n_min, bra_n_max;
        if (delta_L==-1)
          {
            bra_n_min=ket_n;
            bra_n_max=ket_n+1;
          }
        else if (delta_L==+1)
          {
            bra_n_min=ket_n-1;
            bra_n_max=ket_n;
          }
        bra_n_min=std::max(bra_n_min,0);
        bra_n_max=std::min(bra_n_max,bra_dim-1);
        for (int bra_n=bra_n_min; bra_n<=bra_n_max; ++bra_n)
          {
            const int bra_N = 2*bra_n + bra_L;
            double matrix_element =
              analytic::CoordinateOscillatorMatrixElement(
                bra_N, bra_L, ket_N, ket_L, operator_sign
              );

            matrix(bra_n,ket_n) = matrix_element;
          }
      }

    return matrix;
  }

  ////////////////////////////////////////////////////////////////
  // spatial kinematic/transition operators
  ////////////////////////////////////////////////////////////////

  void ConstructCoordinateSqr(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
      relative::CoordinateType coordinate_type,
      int T0
    )
  {
    // validate operator parameters
    assert(operator_labels.J0==0);
    assert(operator_labels.g0==0);
    assert((T0==0)||(T0==1));
    assert((operator_labels.T0_min<=T0)&&(T0<=operator_labels.T0_max));

    // select phase for coordinate or momentum space
    int operator_sign;
    if (coordinate_type == relative::CoordinateType::kR)
      operator_sign = +1;
    else if (coordinate_type == relative::CoordinateType::kK)
      operator_sign = -1;

    // zero initialize operators
    basis::ConstructZeroOperatorRelativeLSJT(
        operator_labels,relative_space,relative_component_sectors,relative_component_matrices
      );

    // select T0 component
    const basis::RelativeSectorsLSJT& sectors = relative_component_sectors[T0];
    basis::OperatorBlocks<double>& matrices = relative_component_matrices[T0];

    // iterate over sectors
    for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
      {

        // set us aliases -- for sector and subspaces
        const basis::RelativeSectorsLSJT::SectorType& sector = sectors.GetSector(sector_index);
        const basis::RelativeSubspaceLSJT& bra_subspace = sector.bra_subspace();
        const basis::RelativeSubspaceLSJT& ket_subspace = sector.ket_subspace();

        // short-circuit select diagonal sectors -- only valid in isoscalar case
	// if (!sector.IsDiagonal())
        //   continue;

        // short-circuit select allowed sectors
        //
        // Scalar (J0=0) and positive parity (g0=0) selection rules
        // should already have been enforced in sector construction.
        //
        // But here we also impose (L0,S0)=(0,0) and T0 triangularity.
        //
        // Furthermore, even in isovector case, we actually have T conservation,
        // by generator action of T.
        if (!(
                am::AllowedTriangle(bra_subspace.L(),0,ket_subspace.L())
                && am::AllowedTriangle(bra_subspace.S(),0,ket_subspace.S())
                // && am::AllowedTriangle(bra_subspace.T(),T0,ket_subspace.T())
                && (bra_subspace.T()==ket_subspace.T())
              )
          )
          continue;

        // extract subspace labels

        int bra_L = bra_subspace.L();
        int ket_L = ket_subspace.L();
        int bra_S = bra_subspace.S();
        int ket_S = ket_subspace.S();
        int bra_J = bra_subspace.J();
        int ket_J = ket_subspace.J();
        int bra_T = bra_subspace.T();
        int ket_T = ket_subspace.T();

        // determine angular/isospin factor
        //
        // Since operators are scalar, only angular factor comes from isospin.

        const int T = ket_T;
        double angular_factor = 1.;
        if (T0==1)
          angular_factor *= 2*std::sqrt(T*(T+1));

        // relative coordinate dilation factor
        //
        // See "Note on oscillator length" at start of header file.

        double relative_oscillator_scale_factor;
        if (coordinate_type == relative::CoordinateType::kR)
          relative_oscillator_scale_factor = 2.;
        else if (coordinate_type == relative::CoordinateType::kK)
          relative_oscillator_scale_factor = 1/2.;

        // alias to matrix
        basis::OperatorBlock<double>& matrix = matrices[sector_index];

        // construct matrix
        //
        // We make use of the known indexing scheme for a
        // RelativeLSJT basis, that the radial quantum number n is
        // just the 0-based state index.
        std::size_t dim = ket_subspace.size();
        int L = ket_L;
        int delta_L = 0;
        matrix = angular_factor * relative_oscillator_scale_factor
          *relative::RadialCoordinateSqrMatrix(dim,dim,L,delta_L,operator_sign);
      }
  }

  void ConstructDipoleOperator(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
      relative::CoordinateType coordinate_type,
      int T0
    )
  {
    // validate operator parameters
    assert(operator_labels.J0==1);
    assert(operator_labels.g0==1);
    assert((T0==1));
    assert((operator_labels.T0_min<=T0)&&(T0<=operator_labels.T0_max));

    // zero initialize operators
    basis::ConstructZeroOperatorRelativeLSJT(
        operator_labels,relative_space,relative_component_sectors,relative_component_matrices
      );

    // select T0 component
    const basis::RelativeSectorsLSJT& sectors = relative_component_sectors[T0];
    basis::OperatorBlocks<double>& matrices = relative_component_matrices[T0];

    // iterate over sectors
    for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
      {

        // set us aliases -- for sector and subspaces
        const basis::RelativeSectorsLSJT::SectorType& sector = sectors.GetSector(sector_index);
        const basis::RelativeSubspaceLSJT& bra_subspace = sector.bra_subspace();
        const basis::RelativeSubspaceLSJT& ket_subspace = sector.ket_subspace();

        // short-circuit select allowed sectors
        //
        // Dipole (J0=1) and negative parity (g0=1) selection rules
        // should already have been enforced in sector construction.
        //
        // But here we also impose (L0,S0)=(1,0) and T0 triangularity.
        if (!(
                am::AllowedTriangle(bra_subspace.L(),1,ket_subspace.L())
                && am::AllowedTriangle(bra_subspace.S(),0,ket_subspace.S())
                && am::AllowedTriangle(bra_subspace.T(),T0,ket_subspace.T())
              )
          )
          continue;

        // extract subspace labels

        int bra_L = bra_subspace.L();
        int ket_L = ket_subspace.L();
        int bra_S = bra_subspace.S();
        int ket_S = ket_subspace.S();
        int bra_J = bra_subspace.J();
        int ket_J = ket_subspace.J();
        int bra_T = bra_subspace.T();
        int ket_T = ket_subspace.T();

        // determine angular/isospin factor
        //
        // Since operators are scalar, only angular factor comes from isospin.

        const HalfInt half = HalfInt(1,2);
        double angular_factor
          = std::sqrt(3.)
          * am::RacahReductionFactor1Rose(bra_L,bra_S,bra_J,ket_L,ket_S,ket_J,1)
          *(
              am::RacahReductionFactor1Rose(half,half,bra_T,half,half,ket_T,1)
              -
              am::RacahReductionFactor2Rose(half,half,bra_T,half,half,ket_T,1)
            );

        // orbital angular RME
        int L = ket_L;
        int delta_L = bra_L-ket_L;
        if (delta_L==-1)
          angular_factor *= std::sqrt((L)/(2.*L+1));
        else if (delta_L==+1)
          angular_factor *= std::sqrt((L+1)/(2.*L+3));

        // relative coordinate dilation factor
        //
        // See "Note on oscillator length" at start of header file.

        double relative_oscillator_scale_factor;
        if (coordinate_type == relative::CoordinateType::kR)
          relative_oscillator_scale_factor = std::sqrt(2.);
        else if (coordinate_type == relative::CoordinateType::kK)
          relative_oscillator_scale_factor = std::sqrt(1/2.);

        // alias to matrix
        basis::OperatorBlock<double>& matrix = matrices[sector_index];

        // populate nonzero entries
        //
        // We make use of the known indexing scheme for a
        // RelativeLSJT basis, that the radial quantum number n is
        // just the 0-based state index.

        std::size_t bra_subspace_size = bra_subspace.size();
        std::size_t ket_subspace_size = ket_subspace.size();

        matrix = angular_factor * relative_oscillator_scale_factor
          *relative::RadialCoordinateMatrix(bra_subspace_size,ket_subspace_size,ket_L,delta_L);
      }

  };

  void ConstructQuadrupoleOperator(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
      relative::CoordinateType coordinate_type,
      int T0
    )
  {

    // validate operator parameters
    assert(operator_labels.J0==2);
    assert(operator_labels.g0==0);
    assert((T0==0)||(T0==1));
    assert((operator_labels.T0_min<=T0)&&(T0<=operator_labels.T0_max));

    // select phase for coordinate or momentum space
    // TODO: check phase for momentum space version and insert into appropriate terms below
    assert(coordinate_type == relative::CoordinateType::kR);
    int operator_sign;
    if (coordinate_type == relative::CoordinateType::kR)
      operator_sign = +1;
    // else if (coordinate_type == relative::CoordinateType::kK)
    //   operator_sign = -1;

    // zero initialize operator
    basis::ConstructZeroOperatorRelativeLSJT(
        operator_labels,relative_space,relative_component_sectors,relative_component_matrices
      );

    // select T0 component
    const basis::RelativeSectorsLSJT& sectors = relative_component_sectors[T0];
    basis::OperatorBlocks<double>& matrices = relative_component_matrices[T0];

    // iterate over sectors
    for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
      {

        // set us aliases -- for sector and subspaces
        const basis::RelativeSectorsLSJT::SectorType& sector = sectors.GetSector(sector_index);
        const basis::RelativeSubspaceLSJT& bra_subspace = sector.bra_subspace();
        const basis::RelativeSubspaceLSJT& ket_subspace = sector.ket_subspace();

        // short-circuit select allowed sectors
        //
        // Quadrupole (J0=2) and positive parity (g0=0) selection rules
        // should already have been enforced in sector construction.
        //
        // But here we also impose (L0,S0)=(2,0) and T0 triangularity.
        //
        // Furthermore, even in isovector case, we actually have T conservation,
        // by generator action of T.
        if (!(
                am::AllowedTriangle(bra_subspace.L(),2,ket_subspace.L())
                && am::AllowedTriangle(bra_subspace.S(),0,ket_subspace.S())
                // && am::AllowedTriangle(bra_subspace.T(),T0,ket_subspace.T())
                && (bra_subspace.T()==ket_subspace.T())
              )
          )
          continue;

        // extract subspace labels

        int bra_L = bra_subspace.L();
        int ket_L = ket_subspace.L();
        int bra_S = bra_subspace.S();
        int ket_S = ket_subspace.S();
        int bra_J = bra_subspace.J();
        int ket_J = ket_subspace.J();
        int bra_T = bra_subspace.T();
        int ket_T = ket_subspace.T();

        // determine angular/isospin factor
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

        const int S = ket_S;
        const int T = ket_T;
        double angular_factor
          = std::sqrt(5./(4*M_PI))
          *ParitySign(bra_L+ket_L+S+ket_J)
          *Hat(bra_J)*Hat(ket_J)*Hat(bra_L)*Hat(ket_L)
          *am::Wigner3J(ket_L,2,bra_L,0,0,0)
          *am::Wigner6J(bra_L,bra_J,S,ket_J,ket_L,2);
        angular_factor /= Hat(bra_J);  // convert to Rose convention
        if (T0==1)
          angular_factor *= 2*std::sqrt(T*(T+1));

        // relative coordinate dilation factor
        //
        // See "Note on oscillator length" at start of header file.

        double relative_oscillator_scale_factor;
        if (coordinate_type == relative::CoordinateType::kR)
          relative_oscillator_scale_factor = 2.;
        else if (coordinate_type == relative::CoordinateType::kK)
          relative_oscillator_scale_factor = 1/2.;

        // alias to matrix
        basis::OperatorBlock<double>& matrix = matrices[sector_index];

        // populate nonzero entries
        //
        // Restrict to a tri-diagonal loop over radial labels.
        //
        // We make use of the known indexing scheme for a RelativeLSJT
        // basis, that the radial quantum number n is just the 0-based state
        // index.

        std::size_t bra_subspace_size = bra_subspace.size();
        std::size_t ket_subspace_size = ket_subspace.size();
        int delta_L = bra_L-ket_L;
        matrix = angular_factor * relative_oscillator_scale_factor
          *relative::RadialCoordinateSqrMatrix(bra_subspace_size,ket_subspace_size,ket_L,delta_L,operator_sign);

      }
  }

  void ConstructOrbitalAMOperator(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
      int T0
    )
  {

    // validate operator parameters
    assert(operator_labels.J0==1);
    assert(operator_labels.g0==0);
    assert((T0==0)||(T0==1));
    assert((operator_labels.T0_min<=T0)&&(T0<=operator_labels.T0_max));

    // zero initialize operator
    basis::ConstructZeroOperatorRelativeLSJT(
        operator_labels,relative_space,relative_component_sectors,relative_component_matrices
      );

    // select T0 component
    const basis::RelativeSectorsLSJT& sectors = relative_component_sectors[T0];
    basis::OperatorBlocks<double>& matrices = relative_component_matrices[T0];

    // iterate over sectors
    for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
      {

        // set us aliases -- for sector and subspaces
        const basis::RelativeSectorsLSJT::SectorType& sector = sectors.GetSector(sector_index);
        const basis::RelativeSubspaceLSJT& bra_subspace = sector.bra_subspace();
        const basis::RelativeSubspaceLSJT& ket_subspace = sector.ket_subspace();

        // short-circuit select allowed sectors
        //
        // Dipole (J0=1) and positive parity (g0=0) selection rules should
        // already have been enforced in sector construction.
        //
        // But here we also impose (L0,S0)=(1,0) and T0 triangularity.
        //
        // For both isoscalar and isovector: We furthermore impose orbital
        // angular momentum conservation by the orbital angular momentum
        // generator action (L'=L).
        //
        // Furthermore, even in isovector case, we actually have T conservation,
        // by generator action of T.
        if (
            !(
                // am::AllowedTriangle(bra_subspace.L(),1,ket_subspace.L())
                am::AllowedTriangle(bra_subspace.S(),0,ket_subspace.S())
                // && am::AllowedTriangle(bra_subspace.T(),T0,ket_subspace.T())
                && (bra_subspace.L()==ket_subspace.L())
                && (bra_subspace.T()==ket_subspace.T())
              )
          )
          continue;

        // extract subspace labels

        int bra_L = bra_subspace.L();
        int ket_L = ket_subspace.L();
        int bra_S = bra_subspace.S();
        int ket_S = ket_subspace.S();
        int bra_J = bra_subspace.J();
        int ket_J = ket_subspace.J();
        int bra_T = bra_subspace.T();
        int ket_T = ket_subspace.T();

        // determine angular/isospin factor

        const int L = ket_L;
        const int T = ket_T;
        double angular_factor
          = am::RacahReductionFactor1Rose(bra_L,bra_S,bra_J,ket_L,ket_S,ket_J,1)
          *std::sqrt(L*(L+1));
        if (T0==1)
          angular_factor *= 2*std::sqrt(T*(T+1));

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

        std::size_t bra_subspace_size = bra_subspace.size();
        std::size_t ket_subspace_size = ket_subspace.size();
        assert(bra_subspace_size==ket_subspace_size);  // (LSJ') and (LSJ) subspaces same size
        // matrix = angular_factor*basis::OperatorBlock<double>:Identity(bra_subspace.size(),ket_subspace.size());
        for (std::size_t ket_n=0; ket_n<ket_subspace_size; ++ket_n)
          {
            std::size_t bra_n = ket_n;
            matrix(bra_n,ket_n) = angular_factor;
          }

      }
  }

  ////////////////////////////////////////////////////////////////
  // spin transition operators
  ////////////////////////////////////////////////////////////////

  void ConstructSpinAMOperator(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
      int T0
    )
  {

    // validate operator parameters
    assert(operator_labels.J0==1);
    assert(operator_labels.g0==0);
    assert((T0==0)||(T0==1));
    assert((operator_labels.T0_min<=T0)&&(T0<=operator_labels.T0_max));

    // zero initialize operator
    basis::ConstructZeroOperatorRelativeLSJT(
        operator_labels,relative_space,relative_component_sectors,relative_component_matrices
      );

    // select T0 component
    const basis::RelativeSectorsLSJT& sectors = relative_component_sectors[T0];
    basis::OperatorBlocks<double>& matrices = relative_component_matrices[T0];

    // iterate over sectors
    for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
      {
        // set us aliases -- for sector and subspaces
        const basis::RelativeSectorsLSJT::SectorType& sector = sectors.GetSector(sector_index);
        const basis::RelativeSubspaceLSJT& bra_subspace = sector.bra_subspace();
        const basis::RelativeSubspaceLSJT& ket_subspace = sector.ket_subspace();

        // short-circuit select allowed sectors
        //
        // Dipole (J0=1) and positive parity (g0=0) selection rules should
        // already have been enforced in sector construction.
        //
        // But here we also impose (L0,S0)=(0,1) and T0 triangularity.
        //
        // For isoscalar: We furthermore impose spin angular momentum conservation by
        // the spin generator action (S'=S).
        if (
            !(
                am::AllowedTriangle(bra_subspace.L(),0,ket_subspace.L())
                && am::AllowedTriangle(bra_subspace.S(),1,ket_subspace.S())
                && am::AllowedTriangle(bra_subspace.T(),T0,ket_subspace.T())
              )
            ||
            ((T0==0)&&!(bra_subspace.S()==ket_subspace.S()))
          )
          continue;

        // extract subspace labels

        int bra_L = bra_subspace.L();
        int ket_L = ket_subspace.L();
        int bra_S = bra_subspace.S();
        int ket_S = ket_subspace.S();
        int bra_J = bra_subspace.J();
        int ket_J = ket_subspace.J();
        int bra_T = bra_subspace.T();
        int ket_T = ket_subspace.T();

        // determine angular/isospin factor

        double angular_factor;
        if (T0==0)
          {
            const int S = ket_S;
            angular_factor = am::RacahReductionFactor2Rose(bra_L,bra_S,bra_J,ket_L,ket_S,ket_J,1)
              *std::sqrt(S*(S+1));
          }
        else  // (T0==1)
          {
            const HalfInt half = HalfInt(1,2);
            angular_factor = 3./2.*am::RacahReductionFactor2Rose(bra_L,bra_S,bra_J,ket_L,ket_S,ket_J,1)
              *(
                  am::RacahReductionFactor1Rose(half,half,bra_S,half,half,ket_S,1)
                  *am::RacahReductionFactor1Rose(half,half,bra_T,half,half,ket_T,1)
                  +
                  am::RacahReductionFactor2Rose(half,half,bra_S,half,half,ket_S,1)
                  *am::RacahReductionFactor2Rose(half,half,bra_T,half,half,ket_T,1)
                );
          }

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

        std::size_t bra_subspace_size = bra_subspace.size();
        std::size_t ket_subspace_size = ket_subspace.size();
        assert(bra_subspace_size==ket_subspace_size);  // (LSJ') and (LSJ) subspaces same size
        // matrix = angular_factor*basis::OperatorBlock<double>:Identity(bra_subspace.size(),ket_subspace.size());
        for (std::size_t ket_n=0; ket_n<ket_subspace_size; ++ket_n)
          {
            std::size_t bra_n = ket_n;
            matrix(bra_n,ket_n) = angular_factor;
          }

      }
  }

  ////////////////////////////////////////////////////////////////
  // Coulomb
  ////////////////////////////////////////////////////////////////

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
        for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
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

            // alias to matrix
            basis::OperatorBlock<double>& matrix = matrices[sector_index];

            // populate nonzero entries
            //
            // We make use of the known indexing scheme for a
            // RelativeLSJT basis, that the radial quantum number n is
            // just the 0-based state index.

            for (int bra_n=0; bra_n<=nmax; ++bra_n)
              for (int ket_n=0; ket_n<=nmax; ++ket_n)
                {

                  // evaluate radial integral
                  const int num_size = num_steps+1;
                  double radial_integral =
                    spline::RadialMatrixElement(
                      bra_n, L, 1, spline::BasisType::kOscillator,
                      ket_n, L, 1, spline::BasisType::kOscillator,
                      spline::OperatorType::kR, -1,
                      num_size
                    );

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
