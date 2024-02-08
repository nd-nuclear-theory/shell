/****************************************************************
  tbme_radial_xform.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "basis/nlj_operator.h"
#include "fmt/format.h"  // for debugging
#include "tbme/tbme_radial_xform.h"

namespace shell {
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  // transformation overlap matrix
  //
  // Note closely parallel structure to KinematicScalarTBME in
  // tbme_separable.
  ////////////////////////////////////////////////////////////////

  double OneBodyOverlap(
      // radial overlap data
      const basis::OrbitalSpaceLJPN& radial_source_orbital_space,
      const basis::OrbitalSpaceLJPN& radial_target_orbital_space,
      const basis::OrbitalSectorsLJPN& radial_sectors,
      const basis::OperatorBlocks<double>& radial_matrices,
      // orbital labels
      const basis::OrbitalStatePN& a, const basis::OrbitalStatePN& ap
    )
  // Evaluate one-body overlap <a|a'>.
  //
  // <a|a'> = <R_a|R_a'>[(la,ja)==(la',ja')]
  {

    int la = a.l();
    int lap = ap.l();
    HalfInt ja = a.j();
    HalfInt jap = ap.j();

    double matrix_element = 0.;
    if ( (la == lap) && (ja == jap) )
      {
        matrix_element += basis::MatrixElementLJPN(
            radial_source_orbital_space,radial_target_orbital_space,radial_sectors,radial_matrices,
            a,ap
          );
      }

    return matrix_element;
  }

  double TwoBodyOverlapProduct(
      // radial overlap data
      const basis::OrbitalSpaceLJPN& radial_source_orbital_space,
      const basis::OrbitalSpaceLJPN& radial_target_orbital_space,
      const basis::OrbitalSectorsLJPN& radial_sectors,
      const basis::OperatorBlocks<double>& radial_matrices,
      // two-body labels
      const basis::OrbitalStatePN& a1, const basis::OrbitalStatePN& a2,
      const basis::OrbitalStatePN& a1p, const basis::OrbitalStatePN& a2p
    )
  // Evaluate the unsymmetrized overlap (a1,a2;J|a1',a2';J).
  //
  // This is simply the product of one body overlaps, independent of
  // angular momentum:
  //
  //   (a1,a2;J|a1',a2';J) = <a1|a1'><a2|a2'>
  {

    // short circuit check equality of (l,j) on each one-body factor
    //
    // Note: This is redundant to (but preempts) the (l,j) equality
    // check in OneBodyOverlap.
    bool triangle_allowed = (
        ((a1.l()==a1p.l())&&(a1.j()==a1p.j()))
        && ((a2.l()==a2p.l())&&(a2.j()==a2p.j()))
      );
    if (!triangle_allowed )
      return 0.;

    // evaluate matrix element
    double matrix_element = OneBodyOverlap(
          radial_source_orbital_space,radial_target_orbital_space,radial_sectors,radial_matrices,
          a1,a1p
        )
      * OneBodyOverlap(
          radial_source_orbital_space,radial_target_orbital_space,radial_sectors,radial_matrices,
          a2,a2p
        );

    return matrix_element;
  }

  double TwoBodyOverlap(
      // radial overlap data
      const basis::OrbitalSpaceLJPN& radial_source_orbital_space,
      const basis::OrbitalSpaceLJPN& radial_target_orbital_space,
      const basis::OrbitalSectorsLJPN& radial_sectors,
      const basis::OperatorBlocks<double>& radial_matrices,
      // two-body labels
      const basis::TwoBodyStateJJJPN& bra, const basis::TwoBodyStateJJJPN& ket
    )
  // Evaluate the full (possibly antisymmetrized) overlap <cd;J|ab;J>.
  //
  // Computes <cd;J|ab;J>_AS for pp/nn states, or <cd;J|ab;J>_pn
  // for pn states.
  //
  // See csbasis (52) and (54).
  {
    // extract orbitals
    const basis::OrbitalStatePN& a1p = ket.GetOrbital1();
    const basis::OrbitalStatePN& a2p = ket.GetOrbital2();
    const basis::OrbitalStatePN& a1 = bra.GetOrbital1();
    const basis::OrbitalStatePN& a2 = bra.GetOrbital2();

    // extract sector parameters
    assert(bra.two_body_species()==ket.two_body_species());
    basis::TwoBodySpeciesPN two_body_species = ket.two_body_species();
    assert(bra.J()==ket.J());
    int J = ket.J();

    // evaluate matrix element
    double matrix_element = 0.;
    matrix_element += TwoBodyOverlapProduct(
          radial_source_orbital_space,radial_target_orbital_space,radial_sectors,radial_matrices,
          a1,a2,a1p,a2p
      );
    if (two_body_species != basis::TwoBodySpeciesPN::kPN)
      {
        int phase = - ParitySign(J-a1.j()-a2.j());
        matrix_element += phase * TwoBodyOverlapProduct(
            radial_source_orbital_space,radial_target_orbital_space,radial_sectors,radial_matrices,
            a1,a2,a2p,a1p
          );
      }

    return matrix_element;
  }


  Eigen::MatrixXd TwoBodyTransformationMatrix(
      // radial overlap data
      const basis::OrbitalSpaceLJPN& radial_source_orbital_space,
      const basis::OrbitalSpaceLJPN& radial_target_orbital_space,
      const basis::OrbitalSectorsLJPN& radial_sectors,
      const basis::OperatorBlocks<double>& radial_matrices,
      // two-body indexing
      const basis::TwoBodySubspaceJJJPN& bra_subspace,
      const basis::TwoBodySubspaceJJJPN& ket_subspace
    )
  // Generate transformation matrix from source subspace (bras) to
  // target subspace (kets).
  //
  // This generates a transformation matrix to act "on the right" of a
  // source sector matrix.  So don't be confused by the fact that the
  // kets of the operator matrix are the bras of this transformation
  // matrix.
  {
    // generate matrix for sector
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(bra_subspace.size(),ket_subspace.size());

    // recover sector properties
    assert(bra_subspace.two_body_species()==ket_subspace.two_body_species());
    const basis::TwoBodySpeciesPN two_body_species = bra_subspace.two_body_species();

    // for all pairs of states in sector
    const std::size_t bra_subspace_size = bra_subspace.size();
    const std::size_t ket_subspace_size = ket_subspace.size();
    // #pragma omp parallel for collapse(2)
    for (std::size_t bra_index = 0; bra_index < bra_subspace_size; ++bra_index)
      for (std::size_t ket_index = 0; ket_index < ket_subspace_size; ++ket_index)
        {

          // construct states
          basis::TwoBodyStateJJJPN bra(bra_subspace,bra_index);
          basis::TwoBodyStateJJJPN ket(ket_subspace,ket_index);

          // calculate matrix element (pn or AS)
          double matrix_element = TwoBodyOverlap(
              radial_source_orbital_space,radial_target_orbital_space,radial_sectors,radial_matrices,
              bra,ket
            );

          // convert to NAS if needed
          if (two_body_species!=basis::TwoBodySpeciesPN::kPN)
            {
              if (bra.index1()==bra.index2())
                matrix_element *= 1/(sqrt(2.));
              if (ket.index1()==ket.index2())
                matrix_element *= 1/(sqrt(2.));
            }

          // store matrix element
          matrix(bra_index,ket_index) = matrix_element;
        }

    return matrix;

  };

  ////////////////////////////////////////////////////////////////
  // similarity transformation
  ////////////////////////////////////////////////////////////////

  Eigen::MatrixXd TwoBodyTransformedMatrix(
      // radial overlap data
      const basis::OrbitalSpaceLJPN& radial_source_orbital_space,
      const basis::OrbitalSpaceLJPN& radial_target_orbital_space,
      const basis::OrbitalSectorsLJPN& radial_sectors,
      const basis::OperatorBlocks<double>& radial_matrices,
      // two-body indexing
      const typename basis::TwoBodySectorsJJJPN::SectorType& source_sector,
      const typename basis::TwoBodySectorsJJJPN::SectorType& target_sector,
      // matrix data
      const Eigen::MatrixXd& source_matrix
    )
  {

    // populate xform matrices
    //
    // Matrices contain <a1,a2;J|a1',a2';J>, for normalized (pn or
    // NAS) states.
    //
    // Note: For a diagonal sector, the bra and ket xform matrices are
    // identical (to within transposition).  For a scalar operator, only
    // diagonal sectors are involved.  So this is always the case when
    // considering Hamiltonian-like operators.  We exploit this symmetry
    // in the interest of efficiency, but we retain full generality for
    // nonscalar operators.
    Eigen::MatrixXd ket_xform_matrix = TwoBodyTransformationMatrix(
      radial_source_orbital_space,radial_target_orbital_space,radial_sectors,radial_matrices,
      source_sector.ket_subspace(),target_sector.ket_subspace()
      );
    Eigen::MatrixXd target_matrix;

    if (source_sector.IsDiagonal())
    {
      const Eigen::MatrixXd& bra_xform_matrix = ket_xform_matrix;

      // carry out xform
      target_matrix = bra_xform_matrix.transpose() * source_matrix * ket_xform_matrix;
    }
    else
    {
      Eigen::MatrixXd bra_xform_matrix = TwoBodyTransformationMatrix(
        radial_source_orbital_space,radial_target_orbital_space,radial_sectors,radial_matrices,
        source_sector.bra_subspace(),target_sector.bra_subspace()
        );

      // carry out xform
      target_matrix = bra_xform_matrix.transpose() * source_matrix * ket_xform_matrix;
    }

    // diagonstics for inspecting unitary transformations
    //
    // std::cout
    //   << fmt::format("Sector: {} {}",target_sector.bra_subspace().LabelStr(),target_sector.ket_subspace().LabelStr()) << std::endl
    //   << "Source" << std::endl
    //   << source_matrix  << std::endl
    //   << "Transformation" << std::endl
    //   << ket_xform_matrix  << std::endl
    //   << "Target" << std::endl
    //   << target_matrix  << std::endl;

    return target_matrix;
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
