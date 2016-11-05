/****************************************************************
  tbme_xform.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

// #include "cppformat/format.h"  // for debugging
#include "tbme/tbme_xform.h"

namespace shell {
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

#if 0
  void TwoBodyMatrixSectorTransform (const legacy::TwoBodyMatrixNljTzJP& source_matrix, legacy::TwoBodyMatrixNljTzJP& destination_matrix, 
                                     const XformRadialMatrices& xform_radial_matrices, 
                                     int N1b_cut, int N2b_cut,
                                     const legacy::SectorNljTzJP& sector
    )
  {
    // recover sector properties
    const legacy::TwoSpeciesStateType state_type = sector.GetStateType();
    const int J = sector.GetJ();
    const int g = sector.GetGrade();

    // Notation: unprimed variables are for untransformed
    // (oscillator) states, primed variables are for transformed
    // states

    // extract dimension of subspace
    const int dimension = destination_matrix.GetTwoBodyBasis().GetDimension(state_type,J,g);

    // for canonical pairs of states (in lexicographical order)

    // OMP: 0uter loop parallelization is "embarassingly parallel"
    // parallelization of target matrix elements involving
    // differen bra-ket pairs.  This avoids thread overhead of
    // parallelizing the reduction operation in the inner loops.
    // Optimal approach in tests based on a simple computational
    // load (i.e., no memory access) in inner loop.

    // #pragma omp parallel for collapse(2) 
    // DEBUGGING: under craype CC compiler with gcc 6, gives error:
    //   error: initializer expression refers to iteration variable 'k1p' 
    for (int k1p = 0; k1p < dimension; ++k1p)
      for (int k2p = k1p; k2p < dimension; ++k2p)
        {
          // identify target matrix element
          legacy::TwoBodyStateNlj s1p = destination_matrix.GetTwoBodyBasis().GetState(state_type,J,g,k1p);
          legacy::TwoBodyStateNlj s2p = destination_matrix.GetTwoBodyBasis().GetState(state_type,J,g,k2p);

          // DBG: std::cout << "  Evaluating:" 
          // DBG:      << " type " << state_type << " J " << J
          // DBG:      << " " << s1p.a1 << s1p.a2 << s2p.a1 << s2p.a2 
          // DBG:      << std::endl;

          // define shorthands for n and l quantum numbers (for more readable coding below)
          int n11p = s1p.a1.Getn();
          int n12p = s1p.a2.Getn();
          int n21p = s2p.a1.Getn();
          int n22p = s2p.a2.Getn();

          int l11p = s1p.a1.Getl();
          int l12p = s1p.a2.Getl();
          int l21p = s2p.a1.Getl();
          int l22p = s2p.a2.Getl();

          // determine summation limits under N1b truncation for source matrix elements
          //    n = (N-l)/2 and rely on integer division to give floor of (N1b - l) /2
          //    but if N1b_cut < l this cannot be relied upon, since integer division truncation with 
          //    negative argument is initially implementation defined with late standard towards zero
          int n11_max = (N1b_cut >= l11p) ? (N1b_cut - l11p) / 2 : -1;
          int n12_max = (N1b_cut >= l12p) ? (N1b_cut - l12p) / 2 : -1;
          int n21_max = (N1b_cut >= l21p) ? (N1b_cut - l21p) / 2 : -1;
          int n22_max = (N1b_cut >= l22p) ? (N1b_cut - l22p) / 2 : -1;

          // carry out sum

          // Note: Sum is explicitly over radial quantum numbers, rather than over all two-body 
          // states in unprimed basis, since this allows much greater a priori restriction 
          // of the sum, by imposition of the Kroenecker delta on l and j quantum numbers.
          // (Whether or not these reduced loops are actually warranted for efficiency is not
          // thoroughly explored and may be expected to depend on basis size.)  The for loops
          // directly implement the one-body cutoff on the sum over the unprimed basis (i.e., Ncut 
          // as defined in csbasis).  But the two-body cutoff is then implemented with a short-circuit
          // text inside the loop.

          double matrix_element = 0.;
          for (int n11 = 0; n11 <= n11_max; ++n11)
            for (int n12 = 0; n12 <= n12_max; ++n12)
              for (int n21 = 0; n21 <= n21_max; ++n21)
                for (int n22 = 0; n22 <= n22_max; ++n22)
                  {
                    // look up orbitals for source matrix element
                    int N11 = 2*n11 + l11p;
                    int N12	= 2*n12 + l12p;
                    int N21	= 2*n21 + l21p;
                    int N22	= 2*n22 + l22p;
                    legacy::SPOrbitalNlj a11(N11,s1p.a1.Getj());
                    legacy::SPOrbitalNlj a12(N12,s1p.a2.Getj());
                    legacy::SPOrbitalNlj a21(N21,s2p.a1.Getj());
                    legacy::SPOrbitalNlj a22(N22,s2p.a2.Getj());
							
                    // build states for bracket
                    legacy::TwoBodyStateNlj s1(a11,a12,J);
                    legacy::TwoBodyStateNlj s2(a21,a22,J);
                    // DBG: std::cout << "    Adding: " << s1.a1 << s1.a2 << s2.a1 << s2.a2 << std::endl;

                    // impose 2-body cutoff on unprimed states
                    if ( 
                        ( (s1.a1.GetN() + s1.a2.GetN()) > N2b_cut)
                        ||
                        ( (s2.a1.GetN() + s2.a2.GetN()) > N2b_cut)
                      )
                      continue;

                    // skip symmetry forbidden states
                    if ( !SymmetryAllowedState(state_type,s1) ||  !SymmetryAllowedState(state_type,s2) )
                      continue;

                    // revise bracket to canonical order
                    int phase = CanonicalizeTwoBodyBracket(state_type, s1, s2, 0);
                    // DBG: std::cout << "      ====> " << s1.a1 << s1.a2 << s2.a1 << s2.a2 << std::endl;

                    // accumulate source matrix element to target matrix element

                    // Note (2013): Transformation brackets in radial xform file are stored  
                    // with unprimed (HO) basis as bra (row index) and primed (new) basis as 
                    // ket (column).
							

                    // The patern for pn bases is...
                    // 	double coefficient 
                    // 		= transformation1.GetMatrixElement(l11p,l11p,n11,n11p)
                    // 		* transformation2.GetMatrixElement(l12p,l12p,n12,n12p)
                    // 		* transformation1.GetMatrixElement(l21p,l21p,n21,n21p)
                    // 		* transformation2.GetMatrixElement(l22p,l22p,n22,n22p);

                    double coefficient;
                    if (state_type == legacy::kPP)
                      {
                        coefficient
                          = xform_radial_matrices.xform_p.GetMatrixElement(l11p,l11p,n11,n11p)
                          * xform_radial_matrices.xform_p.GetMatrixElement(l12p,l12p,n12,n12p)
                          * xform_radial_matrices.xform_p.GetMatrixElement(l21p,l21p,n21,n21p)
                          * xform_radial_matrices.xform_p.GetMatrixElement(l22p,l22p,n22,n22p);
                      }
                    else if (state_type == legacy::kNN)
                      {
                        coefficient
                          = xform_radial_matrices.xform_n.GetMatrixElement(l11p,l11p,n11,n11p)
                          * xform_radial_matrices.xform_n.GetMatrixElement(l12p,l12p,n12,n12p)
                          * xform_radial_matrices.xform_n.GetMatrixElement(l21p,l21p,n21,n21p)
                          * xform_radial_matrices.xform_n.GetMatrixElement(l22p,l22p,n22,n22p);
                      }
                    else if (state_type == legacy::kPN)
                      {
                        coefficient
                          = xform_radial_matrices.xform_p.GetMatrixElement(l11p,l11p,n11,n11p)
                          * xform_radial_matrices.xform_n.GetMatrixElement(l12p,l12p,n12,n12p)
                          * xform_radial_matrices.xform_p.GetMatrixElement(l21p,l21p,n21,n21p)
                          * xform_radial_matrices.xform_n.GetMatrixElement(l22p,l22p,n22,n22p);
                      }

                    // DBG: std::cout << "      Contribution:" 
                    // DBG:      << " phase " << phase 
                    // DBG:      << " coefficient " << coefficient 
                    // DBG:      << " NAS " << source_matrix.GetMatrixElement(state_type, s1, s2) 
                    // DBG:      << " UNAS " << source_matrix.GetMatrixElementUNAS(state_type, s1, s2)
                    // DBG:      << " total " << phase * coefficient * source_matrix.GetMatrixElementUNAS(state_type, s1, s2)
                    // DBG:      << std::endl;

                    matrix_element += phase * coefficient * source_matrix.GetMatrixElementUNAS(state_type, s1, s2); 
                  }
			
          // save resultant matrix element

          // OMP: The "critical" designation for saving the matrix element is out of caution since
          // STL containers are supposedly not thread-safe.
          // #pragma omp critical
          destination_matrix.SetMatrixElementUNAS(state_type, s1p, s2p, matrix_element); 
			
          // DBG: std::cout << "    Resultant:" 
          // DBG:      << " calculated UNAS " << matrix_element 
          // DBG:      << " readback NAS " << destination_matrix.GetMatrixElement(state_type, s1p, s2p) 
          // DBG:      << std::endl;
			
        }
  }
#endif    

  Eigen::MatrixXd TwoBodyTransformedMatrix(
      // radial overlap data
      const basis::OrbitalSpaceLJPN& radial_source_orbital_space,
      const basis::OrbitalSpaceLJPN& radial_target_orbital_space,
      const basis::OrbitalSectorsLJPN& radial_sectors,
      const basis::MatrixVector& radial_matrices,
      // two-body indexing
      const typename basis::TwoBodySectorsJJJPN::SectorType& source_sector,
      const typename basis::TwoBodySectorsJJJPN::SectorType& target_sector,
      const shell::TwoBodyMapping& two_body_mapping,
      // matrix data
      const Eigen::MatrixXd& source_matrix
    )
  {

    // allocate xform matrices
    //
    // Note: For a diagonal sector, the bra and ket xform matrices are
    // identical (to within transposition).  For a scalar operator, only
    // diagonal sectors are involved.  So this is always the case when
    // considering Hamiltonian-like operators.  We assume diagonal
    // sectors for now, in the interest of efficiency, but we retain
    // some of the framework for later full generality.
    assert(source_sector.IsDiagonal());
    assert(target_sector.IsDiagonal());
    Eigen::MatrixXd ket_xform_matrix(source_sector.ket_subspace().size(),target_sector.ket_subspace().size());
    const Eigen::MatrixXd& bra_xform_matrix = ket_xform_matrix;

    // populate xform matrices
    //
    // Matrices contain <a1,a2;J|a1',a2';J> for normalized states.

    // carry out xform
    Eigen::MatrixXd target_matrix = bra_xform_matrix.transpose() * source_matrix * ket_xform_matrix;

    return target_matrix;
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
