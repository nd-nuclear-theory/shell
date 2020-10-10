/****************************************************************
  tbme_separable.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include <memory>

#include "am/racah_reduction.h"
#include "am/wigner_gsl.h"
#include "basis/nlj_operator.h"
#include "fmt/format.h"  // for debugging
#include "obme/obme.h"
#include "tbme/tbme_separable.h"

namespace shell {
  ////////////////////////////////////////////////////////////////
  // utility functions
  ////////////////////////////////////////////////////////////////

  inline bool one_body_sector_allowed(
      const basis::OrbitalSectorsLJPN& sectors,
      const basis::OrbitalStatePN& bra,
      const basis::OrbitalStatePN& ket
    )
  {
    bool allowed = true;
    allowed &= am::AllowedTriangle(ket.j(), sectors.j0(), bra.j());
    allowed &= ((ket.g()+sectors.g0()+bra.g())%2 == 0);
    allowed &= (bra.Tz() == ket.Tz() + sectors.Tz0());
    return allowed;
  }

  ////////////////////////////////////////////////////////////////
  // upgraded one-body operator
  ////////////////////////////////////////////////////////////////

  basis::OperatorBlock<double>
  UpgradeOneBodyOperatorJJJPN(
      const basis::OrbitalSpaceLJPN& ob_orbital_space,
      const basis::OrbitalSectorsLJPN& ob_sectors,
      const basis::OperatorBlocks<double>& ob_matrices,
      const basis::TwoBodySectorsJJJPN::SectorType& sector,
      const int A
    )
  {
    // set up aliases
    const int j0 = ob_sectors.j0();
    const int g0 = ob_sectors.g0();
    const int Tz0 = ob_sectors.Tz0();
    const basis::TwoBodySubspaceJJJPN& bra_subspace = sector.bra_subspace();
    const std::size_t bra_subspace_size = bra_subspace.size();
    const basis::TwoBodySubspaceJJJPN& ket_subspace = sector.ket_subspace();
    const std::size_t ket_subspace_size = ket_subspace.size();

    basis::OperatorBlock<double> matrix =
      basis::OperatorBlock<double>::Zero(bra_subspace_size, ket_subspace_size);

    // short circuit on triangle constraint
    if (!am::AllowedTriangle(bra_subspace.J(), ket_subspace.J(), j0))
      return matrix;

    // build sector matrix
    #pragma omp parallel for \
      collapse(2)            \
      default(none),         \
      shared(                \
        matrix, j0, g0, Tz0, \
        bra_subspace, bra_subspace_size, ket_subspace, ket_subspace_size, \
        ob_orbital_space, ob_sectors, ob_matrices, sector, A \
        )
    for (std::size_t bra_index = 0; bra_index < bra_subspace_size; ++bra_index)
      for (std::size_t ket_index = 0; ket_index < ket_subspace_size; ++ket_index)
      {
        // construct states
        const basis::TwoBodyStateJJJPN bra(bra_subspace, bra_index);
        const basis::TwoBodyStateJJJPN ket(ket_subspace, ket_index);
        const basis::OrbitalStatePN& a = ket.GetOrbital1();
        const basis::OrbitalStatePN& b = ket.GetOrbital2();
        const basis::OrbitalStatePN& c = bra.GetOrbital1();
        const basis::OrbitalStatePN& d = bra.GetOrbital2();

        // fetch one-body matrix elements
        double matel_ca = 0, matel_cb = 0, matel_da = 0, matel_db = 0;
        if (one_body_sector_allowed(ob_sectors, c, a))
          matel_ca =
            basis::MatrixElementLJPN(
                ob_orbital_space, ob_orbital_space, ob_sectors, ob_matrices, c, a
              );
        if (one_body_sector_allowed(ob_sectors, c, b))
          matel_cb =
            basis::MatrixElementLJPN(
                ob_orbital_space, ob_orbital_space, ob_sectors, ob_matrices, c, b
              );
        if (one_body_sector_allowed(ob_sectors, d, a))
          matel_da =
            basis::MatrixElementLJPN(
                ob_orbital_space, ob_orbital_space, ob_sectors, ob_matrices, d, a
              );
        if (one_body_sector_allowed(ob_sectors, d, b))
          matel_db =
            basis::MatrixElementLJPN(
                ob_orbital_space, ob_orbital_space, ob_sectors, ob_matrices, d, b
              );

        // calculate AS two-body element
        double matrix_element = 0.;
        int phase = 1;

        // (cd;J'||u1||ab;J)+(cd;J'||u2||ab;J)
        if (d == b)
        {
          matrix_element += phase
            * am::RacahReductionFactor1Rose(
                  c.j(), d.j(), bra.J(), a.j(), b.j(), ket.J(), j0
                )
            * matel_ca;
        }
        if (c == a)
        {
          matrix_element += phase
            * am::RacahReductionFactor2Rose(
                  c.j(), d.j(), bra.J(), a.j(), b.j(), ket.J(), j0
                )
            * matel_db;
        }

        // - (-)**(J-ja-jb) * [(cd;J'||u1||ba;J)+(cd;J'||u2||ba;J)]
        phase = - ParitySign(ket.J()-a.j()-b.j());
        if (d == a)
        {
          matrix_element += phase
            * am::RacahReductionFactor1Rose(
                  c.j(), d.j(), bra.J(), b.j(), a.j(), ket.J(), j0
                )
            * matel_cb;
        }
        if (c == b)
        {
          matrix_element += phase
            * am::RacahReductionFactor2Rose(
                  c.j(), d.j(), bra.J(), b.j(), a.j(), ket.J(), j0
                )
            * matel_da;
        }

        // - (-)**(J'-jc-jd) * [(dc;J'||u1||ab;J)+(dc;J'||u2||ab;J)]
        phase = - ParitySign(bra.J()-c.j()-d.j());
        if (c == b)
        {
          matrix_element += phase
            * am::RacahReductionFactor1Rose(
                  d.j(), c.j(), bra.J(), a.j(), b.j(), ket.J(), j0
                )
            * matel_da;
        }
        if (d == a)
        {
          matrix_element += phase
            * am::RacahReductionFactor2Rose(
                  d.j(), c.j(), bra.J(), a.j(), b.j(), ket.J(), j0
                )
            * matel_cb;
        }

        // (-)**(J'-jc-jd+J-ja-jb) * [(dc;J'||u1||ba;J)+(dc;J'||u2||ba;J)]
        phase = ParitySign(bra.J()-c.j()-d.j()+ket.J()-a.j()-b.j());
        if (c == a)
        {
          matrix_element += phase
            * am::RacahReductionFactor1Rose(
                  d.j(), c.j(), bra.J(), b.j(), a.j(), ket.J(), j0
                )
            * matel_db;
        }
        if (d == b)
        {
          matrix_element += phase
            * am::RacahReductionFactor2Rose(
                  d.j(), c.j(), bra.J(), b.j(), a.j(), ket.J(), j0
                )
            * matel_ca;
        }

        // normalize (from sqrt(2) in definition of AS states)
        matrix_element *= 0.5;

        // convert to NAS if needed
        if (a == b)
          matrix_element *= 1/(sqrt(2.));
        if (c == d)
          matrix_element *= 1/(sqrt(2.));

        // NOTE(pjf): uncertain if this write needs to be in a critical region,
        // but the synchronization cost is high
        // #pragma omp critical
        matrix(bra_index, ket_index) = matrix_element;
      }
    return 1./(A-1)*matrix;
  }


  ////////////////////////////////////////////////////////////////
  // Racah reduction formula
  ////////////////////////////////////////////////////////////////

  basis::OperatorBlock<double>
  RacahReduceTensorProductJJJPN(
      const basis::OrbitalSpaceLJPN& ob_orbital_space,
      const basis::OrbitalSectorsLJPN& ob_sectors1,
      const basis::OperatorBlocks<double>& ob_matrices1,
      const basis::OrbitalSectorsLJPN& ob_sectors2,
      const basis::OperatorBlocks<double>& ob_matrices2,
      const basis::TwoBodySectorsJJJPN::SectorType& sector,
      const int J0
    )
  {
    // sanity check on angular momenta
    assert(am::AllowedTriangle(ob_sectors1.j0(), ob_sectors2.j0(), J0));

    const basis::TwoBodySubspaceJJJPN& bra_subspace = sector.bra_subspace();
    const std::size_t bra_subspace_size = bra_subspace.size();
    const basis::TwoBodySubspaceJJJPN& ket_subspace = sector.ket_subspace();
    const std::size_t ket_subspace_size = ket_subspace.size();

    basis::OperatorBlock<double> matrix =
      basis::OperatorBlock<double>::Zero(bra_subspace_size, ket_subspace_size);

    // short circuit on triangle constraint
    if (!am::AllowedTriangle(bra_subspace.J(), ket_subspace.J(), J0))
      return matrix;

    // build sector matrix
    #pragma omp parallel for                                                    \
      collapse(2)                                                               \
      default(none),                                                            \
      shared(                                                                   \
        matrix, bra_subspace, bra_subspace_size, ket_subspace, ket_subspace_size,\
        ob_orbital_space, ob_sectors1, ob_matrices1, ob_sectors2, ob_matrices2, \
        sector, J0                                                              \
        )
    for (std::size_t bra_index = 0; bra_index < bra_subspace_size; ++bra_index)
      for (std::size_t ket_index = 0; ket_index < ket_subspace_size; ++ket_index)
      {
        // construct states
        const basis::TwoBodyStateJJJPN bra(bra_subspace, bra_index);
        const basis::TwoBodyStateJJJPN ket(ket_subspace, ket_index);
        const basis::OrbitalStatePN& a = ket.GetOrbital1();
        const basis::OrbitalStatePN& b = ket.GetOrbital2();
        const basis::OrbitalStatePN& c = bra.GetOrbital1();
        const basis::OrbitalStatePN& d = bra.GetOrbital2();

        // fetch one-body matrix elements
        double matel_c1a = 0., matel_c1b = 0., matel_d1a = 0., matel_d1b = 0.;
        double matel_c2a = 0., matel_c2b = 0., matel_d2a = 0., matel_d2b = 0.;
        if (one_body_sector_allowed(ob_sectors1, c, a))
          matel_c1a =
            basis::MatrixElementLJPN(
                ob_orbital_space, ob_orbital_space, ob_sectors1, ob_matrices1, c, a
              );
        if (one_body_sector_allowed(ob_sectors1, c, b))
          matel_c1b =
            basis::MatrixElementLJPN(
                ob_orbital_space, ob_orbital_space, ob_sectors1, ob_matrices1, c, b
              );
        if (one_body_sector_allowed(ob_sectors1, d, a))
          matel_d1a =
            basis::MatrixElementLJPN(
                ob_orbital_space, ob_orbital_space, ob_sectors1, ob_matrices1, d, a
              );
        if (one_body_sector_allowed(ob_sectors1, d, b))
          matel_d1b =
            basis::MatrixElementLJPN(
                ob_orbital_space, ob_orbital_space, ob_sectors1, ob_matrices1, d, b
              );
        if ((std::addressof(ob_sectors1) == std::addressof(ob_sectors2))
            && (std::addressof(ob_matrices1) == std::addressof(ob_matrices2)))
        {
          matel_c2a = matel_c1a;
          matel_c2b = matel_c1b;
          matel_d2a = matel_d1a;
          matel_d2b = matel_d1b;
        } else {
          if (one_body_sector_allowed(ob_sectors2, c, a))
            matel_c2a =
              basis::MatrixElementLJPN(
                  ob_orbital_space, ob_orbital_space, ob_sectors2, ob_matrices2, c, a
                );
          if (one_body_sector_allowed(ob_sectors2, c, b))
            matel_c2b =
              basis::MatrixElementLJPN(
                  ob_orbital_space, ob_orbital_space, ob_sectors2, ob_matrices2, c, b
                );
          if (one_body_sector_allowed(ob_sectors2, d, a))
            matel_d2a =
              basis::MatrixElementLJPN(
                  ob_orbital_space, ob_orbital_space, ob_sectors2, ob_matrices2, d, a
                );
          if (one_body_sector_allowed(ob_sectors2, d, b))
            matel_d2b =
              basis::MatrixElementLJPN(
                  ob_orbital_space, ob_orbital_space, ob_sectors2, ob_matrices2, d, b
                );
        }

        double matrix_element = 0.;
        int phase = 1;
        double factor_c1a_d2b, factor_d2b_c1a, factor_d1b_c2a, factor_c2a_d1b,
          factor_c1b_d2a, factor_d2a_c1b, factor_d1a_c2b, factor_c2b_d1a;

        factor_c1a_d2b =
          am::RacahReductionFactor12Rose(
              c.j(), d.j(), bra.J(),
              a.j(), b.j(), ket.J(),
              ob_sectors1.j0(), ob_sectors2.j0(), J0
            );
        // factor_d2b_c1a =
        //   ParitySign(bra.J()+c.j()+d.j()+ket.J()+a.j()+b.j())
        //   * factor_c1a_d2b;

        // phase = 1;
        // matrix_element += phase * factor_c1a_d2b * matel_c1a * matel_d2b;
        // phase = ParitySign(bra.J()-c.j()-d.j()+ket.J()-a.j()-b.j());  // AS phase
        // matrix_element += phase * factor_d2b_c1a * matel_d2b * matel_c1a;
        matrix_element += 2 * factor_c1a_d2b * matel_d2b * matel_c1a;


        if (ob_sectors1.j0() == ob_sectors2.j0())
          factor_d1b_c2a =
            ParitySign(bra.J()+c.j()+d.j()+ket.J()+a.j()+b.j()+ob_sectors1.j0()+ob_sectors2.j0()+J0)
            * factor_c1a_d2b;
        else
          factor_d1b_c2a =
            am::RacahReductionFactor12Rose(
                d.j(), c.j(), bra.J(),
                b.j(), a.j(), ket.J(),
                ob_sectors1.j0(), ob_sectors2.j0(), J0
              );
        // factor_c2a_d1b =
        //   ParitySign(bra.J()+c.j()+d.j()+ket.J()+a.j()+b.j())
        //   * factor_d1b_c2a;

        // phase = ParitySign(bra.J()-c.j()-d.j()+ket.J()-a.j()-b.j());  // AS phase
        // matrix_element += phase * factor_d1b_c2a * matel_d1b * matel_c2a;
        // phase = 1;
        // matrix_element += phase * factor_c2a_d1b * matel_c2a * matel_d1b;
        // matrix_element += (ParitySign(bra.J()+c.j()+d.j()+ket.J()+a.j()+b.j()) + ParitySign(bra.J()-c.j()-d.j()+ket.J()-a.j()-b.j())) * factor_d1b_c2a * matel_d1b * matel_c2a;
        matrix_element += 2 * ParitySign(bra.J()-c.j()-d.j()+ket.J()-a.j()-b.j()) * factor_d1b_c2a * matel_d1b * matel_c2a;


        factor_c1b_d2a =
          am::RacahReductionFactor12Rose(
              c.j(), d.j(), bra.J(),
              b.j(), a.j(), ket.J(),
              ob_sectors1.j0(), ob_sectors2.j0(), J0
            );
        // factor_d2a_c1b =
        //   ParitySign(bra.J()+c.j()+d.j()+ket.J()+a.j()+b.j())
        //   * factor_c1b_d2a;

        // phase = - ParitySign(ket.J()-a.j()-b.j());  // AS phase
        // matrix_element += phase * factor_c1b_d2a * matel_c1b * matel_d2a;
        // phase = - ParitySign(bra.J()-c.j()-d.j());  // AS phase
        // matrix_element += phase * factor_d2a_c1b * matel_d2a * matel_c1b;
        // matrix_element += -(ParitySign(ket.J()-a.j()-b.j()) + ParitySign(ket.J()+a.j()+b.j())) * factor_c1b_d2a * matel_d2a * matel_c1b;
        matrix_element += -2 * ParitySign(ket.J()-a.j()-b.j()) * factor_c1b_d2a * matel_d2a * matel_c1b;


        if (ob_sectors1.j0() == ob_sectors2.j0())
          factor_d1a_c2b =
            ParitySign(bra.J()+c.j()+d.j()+ket.J()+a.j()+b.j()+ob_sectors1.j0()+ob_sectors2.j0()+J0)
            * factor_c1b_d2a;
        else
          factor_d1a_c2b =
            am::RacahReductionFactor12Rose(
                d.j(), c.j(), bra.J(),
                a.j(), b.j(), ket.J(),
                ob_sectors1.j0(), ob_sectors2.j0(), J0
              );
        // factor_c2b_d1a =
        //   ParitySign(bra.J()+c.j()+d.j()+ket.J()+a.j()+b.j())
        //   * factor_d1a_c2b;

        // phase = - ParitySign(bra.J()-c.j()-d.j());  // AS phase
        // matrix_element += phase * factor_d1a_c2b * matel_d1a * matel_c2b;
        // phase = - ParitySign(ket.J()-a.j()-b.j());  // AS phase
        // matrix_element += phase * factor_c2b_d1a * matel_c2b * matel_d1a;
        matrix_element += -2 * ParitySign(bra.J()-c.j()-d.j()) * factor_d1a_c2b * matel_d1a * matel_c2b;


        // normalize (from symmetrization of v and sqrt(2) in AS state)
        matrix_element *= 0.25;

        // convert to NAS if needed
        if (a == b)
          matrix_element *= 1/(sqrt(2.));
        if (c == d)
          matrix_element *= 1/(sqrt(2.));

        // NOTE(pjf): uncertain if this write needs to be in a critical region,
        // but the synchronization cost is high
        // #pragma omp critical
        matrix(bra_index, ket_index) = matrix_element;
      }
    return matrix;
  }

  ////////////////////////////////////////////////////////////////
  // two-body identity operator
  ////////////////////////////////////////////////////////////////

  Eigen::MatrixXd
  IdentityOperatorMatrixJJJPN(
      const basis::TwoBodySectorsJJJPN::SectorType& sector,
      int A
    )
  {

    // set up aliases
    assert(sector.IsDiagonal());
    const basis::TwoBodySubspaceJJJPN& subspace = sector.ket_subspace();

    // generate matrix for sector
    Eigen::MatrixXd matrix;
    const double normalization_factor = 2./(A*(A-1));
    matrix = normalization_factor*Eigen::MatrixXd::Identity(subspace.size(),subspace.size());

    return matrix;
  }

  ////////////////////////////////////////////////////////////////
  // loop timing test
  ////////////////////////////////////////////////////////////////

  Eigen::MatrixXd
  TimingTestMatrixJJJPN(
      const basis::TwoBodySectorsJJJPN::SectorType& sector,
      bool loop,
      bool store
    )
  {

    // set up aliases
    assert(sector.IsDiagonal());
    const basis::TwoBodySubspaceJJJPN& subspace = sector.ket_subspace();

    // generate matrix for sector
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(subspace.size(),subspace.size());

    // recover sector properties
    const basis::TwoBodySpeciesPN two_body_species = subspace.two_body_species();
    const int J = subspace.J();
    const int g = subspace.g();

    if (loop)
      {
        // for upper-triangular pairs of states in sector
        const std::size_t subspace_size = subspace.size();
#pragma omp parallel for collapse(2) if (0)  // disabled until have chance to profile
        for (std::size_t bra_index = 0; bra_index < subspace_size; ++bra_index)
          for (std::size_t ket_index = 0; ket_index < subspace_size; ++ket_index)
            {

              // diagonal sector: restrict to upper triangle
              // if (sector.IsDiagonal())
              if (!(bra_index<=ket_index))
                continue;

              // construct states
              basis::TwoBodyStateJJJPN bra(subspace,bra_index);
              basis::TwoBodyStateJJJPN ket(subspace,ket_index);

              // store matrix element
              if (store)
                matrix(bra_index,ket_index) = -1.;
            }

      }
    return matrix;
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
}  // namespace shell
