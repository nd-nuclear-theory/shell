/****************************************************************
  coulomb_su11.cpp

  RETAINED test code for disappointing SU11 evaluation of Coulomb
  integrals.

  Mark A. Caprio
  University of Notre Dame

  3/26/17 (mac): Extracted from coulomb -> contruct_radial.

****************************************************************/

#include "mcutils/gsl.h"
#include "spline/wavefunction_class.h"

#include "relative/construct_relative.h"
#include "relative/coulomb.h"

namespace relative {

  void ConstructCoulombOperatorSU11(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
      int nmax_extended
    );
  // Construct coulomb operator in relative LSJT basis.
  //
  // Makes use of analytic expression for radial integral of 1/r^2
  // from SU(1,1) factorization, from (48) of Rowe [JPA 38, 10181
  // (2005)], followed by matrix square root.  The operator 1/r is
  // simply sqrt(1/r^2) on coordinate space, but this square root
  // relation is only approximately true as an operator relation on
  // the Hilbert space spanned by a truncated oscillator basis.  This
  // introduces errors which are mitigated by adding a "buffer"
  // (n_extension) to the basis until after the square root is taken,
  // i.e., enlarging the basis for the resolution of the identity in
  // the relation
  //
  //    (1/r) (1/r) ~ (1/r^2)
  //
  // Note: Unfortunately, only poor (~1%) accuracy is obtained for
  // nmax_extended~50, and the code hangs for nmax_extended~100.
  //
  // See notes on "internal representation of an operator in JT
  // scheme" in lsjt_operator.h for the general principles of how the
  // operators are represented.
  //
  // Arguments:
  //   operator_labels (basis::OperatorLabelsJT): tensorial properties of operator
  //   relative_space (...): target space
  //   relative_component_sectors (..., output): target sectors
  //   relative_component_matrices (..., output): target matrices
  //   nmax_extended (int, input): expanded space in which to calculate 1/r^2 matrix before taking sqrt


  void ConstructCoulombOperatorSU11(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSpaceLSJT& relative_space,
      std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
      int nmax_extended
    )
  {

    // zero initialize operator
    ConstructDiagonalConstantOperator(
        operator_labels,relative_space,relative_component_sectors,relative_component_matrices,0.
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
            //
            // Make matrix of 1/r^2, somewhat larger than desired
            // target matrix of 1/r, take matrix square root, and
            // truncate down to size needed for target matrix.

            // set up 1/r^2 matrix
            Eigen::MatrixXd matrix_sqr(nmax_extended+1,nmax_extended+1);
            for (int bra_n=0; bra_n<=nmax_extended; ++bra_n)
              for (int ket_n=0; ket_n<=nmax_extended; ++ket_n)
                {
                  double radial_integral;
                  double lambda = L + 1.5;  // Rowe (17): lambda = v + N/2, N = dimensionality of space
                  if (bra_n <= ket_n)
                    {
                      radial_integral
                        = ParitySign(bra_n-ket_n) / (lambda-1)
                        * std::sqrt(
                            (mcutils::Factorial(ket_n)*mcutils::Factorial(lambda+bra_n-1))
                            / (mcutils::Factorial(bra_n)*mcutils::Factorial(lambda+ket_n-1))
                          );
                    }
                  else
                    {
                      radial_integral
                        = ParitySign(bra_n-ket_n) / (lambda-1)
                        * std::sqrt(
                            (mcutils::Factorial(bra_n)*mcutils::Factorial(lambda+ket_n-1))
                            / (mcutils::Factorial(ket_n)*mcutils::Factorial(lambda+bra_n-1))
                          );
                    }
                  matrix_sqr(bra_n,ket_n) = radial_integral;
                }

            // extract 1/r matrix as square root
            Eigen::MatrixXd& matrix = matrices[sector_index];
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> matrix_sqr_eigensystem(matrix_sqr);
            matrix = matrix_sqr_eigensystem.operatorSqrt().block(0,0,nmax+1,nmax+1);

            double isospin_factor;
            if (T0==0)
              isospin_factor = 1/3.;
            else if (T0==1)
              isospin_factor = 1/std::sqrt(2.);
            else if (T0==2)
              isospin_factor = std::sqrt(10.)/6.;
            matrix *= isospin_factor;


          }
      }

  }

void TestCoulombSU11()
{

  std::cout << "TestCoulombSU11" << std::endl;

  // set tensorial labels
  basis::OperatorLabelsJT operator_labels;
  operator_labels.J0=0;
  operator_labels.g0=0;
  operator_labels.symmetry_phase_mode=basis::SymmetryPhaseMode::kHermitian;
  operator_labels.T0_min=0;
  operator_labels.T0_max=2;

  // set basis parameters
  int Nmax=20;
  int Jmax = Nmax+1;

  // set up relative space
  basis::RelativeSpaceLSJT relative_space(Nmax,Jmax);

  // populate operator containers
  std::array<basis::RelativeSectorsLSJT,3> relative_component_sectors;
  std::array<basis::OperatorBlocks<double>,3> relative_component_matrices;
  basis::ConstructIdentityOperatorRelativeLSJT(
      operator_labels,
      relative_space,
      relative_component_sectors,
      relative_component_matrices
    );

  std::vector<int> nmax_extended_list({10,20,50,100});  // 50 runs easily in a blink, 100 hangs
  for (int nmax_extended : nmax_extended_list)
    {
      std::cout << fmt::format("Nmax {} nmax_extended {}",Nmax,nmax_extended) << std::endl;
      relative::ConstructCoulombOperatorSU11(
          operator_labels,
          relative_space,
          relative_component_sectors,
          relative_component_matrices,
          nmax_extended
        );
      
      std::string filename = fmt::format("coulomb_test_Nmax{}_nmaxe{}.dat",Nmax,nmax_extended);
      basis::WriteRelativeOperatorLSJT(
          filename,
          relative_space,
          operator_labels,relative_component_sectors,relative_component_matrices,
          true  // verbose
        );
    }

}

}
