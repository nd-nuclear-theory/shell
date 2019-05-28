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
#include "tbme/tbme_separable.h"

namespace shell {
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  bool one_body_sector_allowed(const basis::OrbitalSectorsLJPN& sectors,
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
  UpgradeOneBodyOperatorJJJPN(const basis::OrbitalSpaceLJPN& ob_orbital_space,
                              const basis::OrbitalSectorsLJPN& ob_sectors,
                              const basis::OperatorBlocks<double>& ob_matrices,
                              const basis::TwoBodySectorsJJJPN::SectorType& sector,
                              int A
                             )
  {
    const basis::TwoBodySubspaceJJJPN& bra_subspace = sector.bra_subspace();
    const std::size_t bra_subspace_size = bra_subspace.size();
    const basis::TwoBodySubspaceJJJPN& ket_subspace = sector.ket_subspace();
    const std::size_t ket_subspace_size = ket_subspace.size();

    basis::OperatorBlock<double> matrix =
      basis::OperatorBlock<double>::Zero(bra_subspace_size, ket_subspace_size);

    // short circuit on triangle constraint
    int j0 = ob_sectors.j0();
    int g0 = ob_sectors.g0();
    if (!am::AllowedTriangle(bra_subspace.J(), ket_subspace.J(), j0))
      return matrix;

    // build sector matrix
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
      int J0
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

        // (cd;J'||v12||ab;J)+(cd;J'||v21||ab;J)
        matrix_element += phase
          * am::RacahReductionFactor12Rose(
                c.j(), d.j(), bra.J(),
                a.j(), b.j(), ket.J(),
                ob_sectors1.j0(), ob_sectors2.j0(), J0
            )
          * matel_c1a * matel_d2b;
        matrix_element += phase
          * am::RacahReductionFactor21Rose(
                c.j(), d.j(), bra.J(),
                a.j(), b.j(), ket.J(),
                ob_sectors1.j0(), ob_sectors2.j0(), J0
            )
          * matel_c2a * matel_d1b;

        // - (-)**(J-ja-jb) * [(cd;J'||v12||ba;J)+(cd;J'||v21||ba;J)]
        phase = - ParitySign(ket.J()-a.j()-b.j());  // AS phase
        matrix_element += phase
          * am::RacahReductionFactor12Rose(
                c.j(), d.j(), bra.J(),
                b.j(), a.j(), ket.J(),
                ob_sectors1.j0(), ob_sectors2.j0(), J0
            )
          * matel_c1b * matel_d2a;
        matrix_element += phase
          * am::RacahReductionFactor21Rose(
                c.j(), d.j(), bra.J(),
                b.j(), a.j(), ket.J(),
                ob_sectors1.j0(), ob_sectors2.j0(), J0
            )
          * matel_c2b * matel_d1a;

        // - (-)**(J'-jc-jd) * [(dc;J'||v12||ab;J)+(dc;J'||v21||ab;J)]
        phase = - ParitySign(bra.J()-c.j()-d.j());  // AS phase
        matrix_element += phase
          * am::RacahReductionFactor12Rose(
                d.j(), c.j(), bra.J(),
                a.j(), b.j(), ket.J(),
                ob_sectors1.j0(), ob_sectors2.j0(), J0
            )
          * matel_d1a * matel_c2b;
        matrix_element += phase
          * am::RacahReductionFactor21Rose(
                d.j(), c.j(), bra.J(),
                a.j(), b.j(), ket.J(),
                ob_sectors1.j0(), ob_sectors2.j0(), J0
            )
          * matel_d2a * matel_c1b;

        // (-)**(J'-jc-jd+J-ja-jb) * [(dc;J'||v12||ba;J)+(dc;J'||v21||ba;J)]
        phase = ParitySign(bra.J()-c.j()-d.j()+ket.J()-a.j()-b.j());  // AS phase
        matrix_element += phase
          * am::RacahReductionFactor12Rose(
                d.j(), c.j(), bra.J(),
                b.j(), a.j(), ket.J(),
                ob_sectors1.j0(), ob_sectors2.j0(), J0
            )
          * matel_d1b * matel_c2a;
        matrix_element += phase
          * am::RacahReductionFactor21Rose(
                d.j(), c.j(), bra.J(),
                b.j(), a.j(), ket.J(),
                ob_sectors1.j0(), ob_sectors2.j0(), J0
            )
          * matel_d2b * matel_c1a;

        // normalize (from symmetrization of v and sqrt(2) in AS state)
        matrix_element *= 0.25;

        // convert to NAS if needed
        if (a == b)
          matrix_element *= 1/(sqrt(2.));
        if (c == d)
          matrix_element *= 1/(sqrt(2.));

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
  // kinematic operators
  ////////////////////////////////////////////////////////////////

  Eigen::MatrixXd
  CoordinateSqrMatrixJJJPN(
      const basis::OrbitalSpaceLJPN& radial_orbital_space,
      const basis::OrbitalSectorsLJPN& radial_sectors,
      const basis::OperatorBlocks<double>& radial_matrices,
      const RadialOperatorType& radial_operator_type,
      CoordinateSqrOperatorType coordinate_sqr_operator_type,
      const basis::TwoBodySectorsJJJPN::SectorType& sector,
      int A, int J0
    )
  {
    int momentum_space_phase = ParitySign(radial_operator_type == shell::RadialOperatorType::kK);

    basis::OperatorBlock<double> matrix;
    if (coordinate_sqr_operator_type == CoordinateSqrOperatorType::kUTSqr)
    {
      assert(J0==radial_sectors.j0());
      matrix = UpgradeOneBodyOperatorJJJPN(
        radial_orbital_space, radial_sectors, radial_matrices, sector, A
      );
    }
    else // if (coordinate_sqr_operator_type==CoordinateSqrOperatorType::kVT1T2)
    {
      assert(am::AllowedTriangle(radial_sectors.j0(), radial_sectors.j0(), J0));
      matrix = momentum_space_phase*RacahReduceTensorProductJJJPN(
        radial_orbital_space,
        radial_sectors, radial_matrices,
        radial_sectors, radial_matrices,
        sector, J0
      );
      // use the dot product instead of scalar product
      if (J0==0)
        matrix *= ParitySign(radial_sectors.j0()) * Hat(radial_sectors.j0());
    }

    return matrix;
  };

  ////////////////////////////////////////////////////////////////
  // angular momentum operators
  ////////////////////////////////////////////////////////////////

  double AngularMomentumScalarOBME(
      shell::AngularMomentumOperatorFamily operator_family,
      const basis::OrbitalStatePN& b, const basis::OrbitalStatePN& a
    )
  // Give <b|T^2|a> for a squared angular momentum operator, i.e. T^2
  // = l^2, s^2, or j^2.
  //
  // See mac "spin operator" notes page 3.
  {

    int na = a.n();
    int nb = b.n();
    int la = a.l();
    int lb = b.l();
    HalfInt ja = a.j();
    HalfInt jb = b.j();

    double matrix_element = 0.;

    if ( b == a )
      {
        if (operator_family == shell::AngularMomentumOperatorFamily::kOrbital)
          {
            matrix_element += la*(la+1);
          }
        else if (operator_family == shell::AngularMomentumOperatorFamily::kSpin)
          {
            matrix_element += 3./4.;
          }
        else if (operator_family == shell::AngularMomentumOperatorFamily::kTotal)
          {
            matrix_element += double(ja)*(double(ja)+1);
          }
      }

    return matrix_element;
  }


  double AngularMomentumVectorOBRME(
      shell::AngularMomentumOperatorFamily operator_family,
      const basis::OrbitalStatePN& b, const basis::OrbitalStatePN& a
    )
  // Evaluate <b||T||a> for an angular momentum operator, i.e., T = l,
  // s, or j.
  //
  // Based on Suhonen "From nucleons to nucleus" (2.56) and (2.58).  See
  // mac "spin operator" notes page 3.
  {

    int na = a.n();
    int nb = b.n();
    int la = a.l();
    int lb = b.l();
    HalfInt ja = a.j();
    HalfInt jb = b.j();

    double matrix_element = 0.;

    if ( (nb == na) && (lb == la) )
      {
          if (
            (operator_family == shell::AngularMomentumOperatorFamily::kOrbital)
            || (operator_family == shell::AngularMomentumOperatorFamily::kTotal)
          )
            {
              int phase = ParitySign(la+ja+HalfInt(3,2));
              matrix_element += Hat(jb)*Hat(ja)*sqrt(la*(la+1)*(2*la+1))*phase*am::Wigner6J(la,la,1,ja,jb,HalfInt(1,2));
            }
          if (
            (operator_family == shell::AngularMomentumOperatorFamily::kSpin)
            || (operator_family == shell::AngularMomentumOperatorFamily::kTotal)
          )
            {
              int phase = ParitySign(la+jb+HalfInt(3,2));
              matrix_element += sqrt(3./2.)*Hat(jb)*Hat(ja)*phase*am::Wigner6J(HalfInt(1,2),HalfInt(1,2),1,ja,jb,la);
            }
      }

    return matrix_element;
  }

  double AngularMomentumScalarTBME (
      shell::AngularMomentumOperatorFamily operator_family,
      shell::AngularMomentumOperatorSpecies operator_species,
      basis::TwoBodySpeciesPN two_body_species,
      int J,
      const basis::TwoBodyStateJJJPN& s2, const basis::TwoBodyStateJJJPN& s1
    )
  // Evaluate <cd|V_(T^2)|ab> for the "upgraded" two-body operator
  // obtained from a scalar one-body squared angular momentum
  // operator, using <b|T^2|a>.
  //
  // Computes <cd|V|ab>_AS for pp/nn states, or (cd|V|ab)_pn for pn
  //    states.
  //
  // After csbasis (52) and (54).  See mac "spin operator" notes page 3.
  //
  // TODO for elegance and parallelism to vector case, refactor in
  // terms of AngularMomentumScalarTBMEProduct.  See
  // KinematicScalarTBME.

  {
    // extract orbitals
    const basis::OrbitalStatePN& a = s1.GetOrbital1();
    const basis::OrbitalStatePN& b = s1.GetOrbital2();
    const basis::OrbitalStatePN& c = s2.GetOrbital1();
    const basis::OrbitalStatePN& d = s2.GetOrbital2();

    double matrix_element = 0.;

    if (
        ((two_body_species == basis::TwoBodySpeciesPN::kPP) && (operator_species == shell::AngularMomentumOperatorSpecies::kP))
        || ((two_body_species == basis::TwoBodySpeciesPN::kNN) && (operator_species == shell::AngularMomentumOperatorSpecies::kN))
        || ((two_body_species != basis::TwoBodySpeciesPN::kPN) && (operator_species == shell::AngularMomentumOperatorSpecies::kTotal))
        )
      {
        // like-nucleon case
        //
        // short circuited to only apply if operator is for same species or is total operator
        if (d == b)
          matrix_element += AngularMomentumScalarOBME(operator_family,c,a);
        if (c == a)
          matrix_element += AngularMomentumScalarOBME(operator_family,d,b);
        int phase = - ParitySign(J - a.j() - b.j());
        if (d == a)
          matrix_element += phase * AngularMomentumScalarOBME(operator_family,c,b);
        if (c == b)
          matrix_element += phase * AngularMomentumScalarOBME(operator_family,d,a);
      }
    else if (two_body_species == basis::TwoBodySpeciesPN::kPN)
      {
        // proton-neutron case
        if ((d == b) && (operator_species != shell::AngularMomentumOperatorSpecies::kN))
          // term only contributes to proton or total operators, not neutron operator
          matrix_element += AngularMomentumScalarOBME(operator_family,c,a);
        if ((c == a) && (operator_species != shell::AngularMomentumOperatorSpecies::kP))
          // term only contributes to neutron or total operators, not proton operator
          matrix_element += AngularMomentumScalarOBME(operator_family,d,b);
      }

    return matrix_element;
  }

  double AngularMomentumVectorDotTBMEProduct(
      shell::AngularMomentumOperatorFamily operator_family,
      int J,
      const basis::OrbitalStatePN& c, const basis::OrbitalStatePN& d,
      const basis::OrbitalStatePN& a, const basis::OrbitalStatePN& b
    )
  // Evaluate the unsymmetrized matrix element (cd|T1.T2|ab) of a dot
  // product of one-body vector angular momentum operators (i.e., T =
  // l, s, or j), using <b||T||a>.
  //
  // After csbasis (57).  See mac "spin operator" notes page 2.
  {

    // short circuit check on triangularity of each one-body factor
    bool triangle_allowed = ( am::AllowedTriangle(a.j(),c.j(),1) && am::AllowedTriangle(b.j(),d.j(),1));
    if (!triangle_allowed )
      return 0.;

    // evaluate matrix element
    int phase = ParitySign(d.j() + a.j() + J);
    double matrix_element = phase * am::Wigner6J(c.j(),d.j(),J,b.j(),a.j(),1)
      * AngularMomentumVectorOBRME(operator_family,c,a) * AngularMomentumVectorOBRME(operator_family,d,b);

    return matrix_element;
  }

  double AngularMomentumVectorDotTBME(
      shell::AngularMomentumOperatorFamily operator_family,
      shell::AngularMomentumOperatorSpecies operator_species,
      basis::TwoBodySpeciesPN two_body_species,
      int J,
      const basis::TwoBodyStateJJJPN& s2, const basis::TwoBodyStateJJJPN& s1
    )
  // Evaluate the matrix element <cd|T1.T2|ab> of a dot product of
  // one-body vector angular momentum operators (i.e., T = l, s, or
  // j), using the unsymmetrized matrix element (cd|T1.T2|ab).
  //
  // Computes <cd|T1.T2|ab>_AS for pp/nn states, or (cd|T1.T2|ab)_pn
  // for pn states.
  //
  // Note: For the V_(T1.T2) two-body term, a proton operator only has
  // nonvanishing pp sectors, and a neutron operator only has
  // nonvanishing nn sectors.
  //
  // After csbasis (55) and (56).  See mac "spin operator" notes page 2.
  {
    // extract orbitals
    const basis::OrbitalStatePN& a = s1.GetOrbital1();
    const basis::OrbitalStatePN& b = s1.GetOrbital2();
    const basis::OrbitalStatePN& c = s2.GetOrbital1();
    const basis::OrbitalStatePN& d = s2.GetOrbital2();

    double matrix_element = 0.;

    if (
        ((two_body_species == basis::TwoBodySpeciesPN::kPP) && (operator_species == shell::AngularMomentumOperatorSpecies::kP))
        || ((two_body_species == basis::TwoBodySpeciesPN::kNN) && (operator_species == shell::AngularMomentumOperatorSpecies::kN))
        || ((two_body_species != basis::TwoBodySpeciesPN::kPN) && (operator_species == shell::AngularMomentumOperatorSpecies::kTotal))
        )
      {
        // like nucleon case
        matrix_element += AngularMomentumVectorDotTBMEProduct(operator_family,J,c,d,a,b);
        int phase = - ParitySign(J - a.j() - b.j());
        matrix_element += phase * AngularMomentumVectorDotTBMEProduct(operator_family,J,c,d,b,a);
      }
    else if (
        (two_body_species == basis::TwoBodySpeciesPN::kPN)
        && (operator_species == shell::AngularMomentumOperatorSpecies::kTotal)
      )
      {
        // proton-neutron case
        matrix_element += AngularMomentumVectorDotTBMEProduct(operator_family,J,c,d,a,b);
      }

    return matrix_element;
  }

  Eigen::MatrixXd
  AngularMomentumMatrixJJJPN(
      am::AngularMomentumOperatorType operator_family,
      basis::OperatorTypePN operator_species,
      const basis::TwoBodySectorsJJJPN::SectorType& sector,
      int A
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

          // calculate matrix element (pn or AS)
          double matrix_element_t2
            = AngularMomentumScalarTBME(operator_family,operator_species,two_body_species,J,bra,ket);
          double matrix_element_t1t2
            = AngularMomentumVectorDotTBME(operator_family,operator_species,two_body_species,J,bra,ket);
          double matrix_element = 1./(A-1)*matrix_element_t2 + 2*matrix_element_t1t2;

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
          // matrix(ket_index,bra_index) = matrix_element;  // lower triangle
        }

    return matrix;
  }


  ////////////////////////////////////////////////////////////////
  // isospin operators
  ////////////////////////////////////////////////////////////////

  double IsospinOBME(
      IsospinOperatorType operator_type,
      const basis::OrbitalSpaceLJPN& overlap_orbital_space,
      const basis::OrbitalSectorsLJPN& overlap_sectors,
      const basis::OperatorBlocks<double>& overlap_matrices,
      const basis::OrbitalStatePN& b, const basis::OrbitalStatePN& a
    )
  // Give <b|T|a> for a scalar isospin operator, i.e., T = tz, t+, or t-.
  //
  // See pjf "isospin" notes.
  {
    int na = a.n();
    int nb = b.n();
    int la = a.l();
    int lb = b.l();
    HalfInt ja = a.j();
    HalfInt jb = b.j();
    HalfInt Tza = a.Tz();
    HalfInt Tzb = b.Tz();

    double matrix_element = 0.;

    if (operator_type == IsospinOperatorType::kSquared) {
      if ( b == a ) {
        matrix_element += 3./4.;
      }
    } else if (operator_type == IsospinOperatorType::kTz) {
      if ( b == a ) {
        matrix_element += double(a.Tz());
      }
    } else if (operator_type == IsospinOperatorType::kRaising) {
      if ((jb == ja) && (lb == la) && (Tzb == Tza + 1)) {
        matrix_element += basis::MatrixElementLJPN(
            overlap_orbital_space, overlap_orbital_space,
            overlap_sectors, overlap_matrices,
            b, a
          );
      }
    } else if (operator_type == IsospinOperatorType::kLowering) {
      if ((jb == ja) && (lb == la) && (Tzb == Tza - 1)) {
        matrix_element += basis::MatrixElementLJPN(
            overlap_orbital_space, overlap_orbital_space,
            overlap_sectors, overlap_matrices,
            b, a
          );
      }
    }
    return matrix_element;
  }

  double IsospinScalarTBMESum(
      IsospinOperatorType operator_type,
      // radial matrix element data
      const basis::OrbitalSpaceLJPN& overlap_orbital_space,
      const basis::OrbitalSectorsLJPN& overlap_sectors,
      const basis::OperatorBlocks<double>& overlap_matrices,
      // two-body labels
      const basis::OrbitalStatePN& c, const basis::OrbitalStatePN& d,
      const basis::OrbitalStatePN& a, const basis::OrbitalStatePN& b
    )
  // Evaluate the unsymmetrized matrix element (cd;J|T^2|ab;J) of
  // the "upgraded" two-body operator obtained from a scalar one-body
  // operator T = t^2, tz, t+, or t-.
  //
  // See csbasis (52), which is written for an unsymmetrized pn state but
  // applies to any unsymmetrized product state.
  {
    // short circuit check equality of (l,j) on each one-body factor
    //
    // Note: This is redundant to (but preempts) the (l,j) equality
    // check in IsospinOBME.
    bool triangle_allowed = (
        ((c.l() == a.l()) && (c.j() == a.j()))
        && ((d.l() == b.l()) && (d.j() == b.j()))
      );
    if (!triangle_allowed)
      return 0.;

    // evaluate matrix element
    double matrix_element = 0.;
    if (d == b)
      matrix_element += IsospinOBME(
          operator_type,
          overlap_orbital_space, overlap_sectors, overlap_matrices,
          c, a
        );
    if (c == a)
      matrix_element += IsospinOBME(
          operator_type,
          overlap_orbital_space, overlap_sectors, overlap_matrices,
          d, b
        );

    return matrix_element;
  }

  double IsospinScalarTBME(
      IsospinOperatorType operator_type,
      // overlap matrix element data
      const basis::OrbitalSpaceLJPN& overlap_orbital_space,
      const basis::OrbitalSectorsLJPN& overlap_sectors,
      const basis::OperatorBlocks<double>& overlap_matrices,
      // two-body labels
      const basis::TwoBodyStateJJJPN& bra, const basis::TwoBodyStateJJJPN& ket
    )
  // Evaluate <cd|V_(T^2)|ab>_AS for the "upgraded" two-body operator obtained
  // from a scalar one-body isospin operator, using <b|T^2|a>.
  //
  // Computes <cd|V|ab>_AS for pp/nn/pn states.
  //
  // After csbasis (52) and (54).  See pjf "isospin" notes.
  //
  {
    // extract orbitals
    const basis::OrbitalStatePN& a = ket.GetOrbital1();
    const basis::OrbitalStatePN& b = ket.GetOrbital2();
    const basis::OrbitalStatePN& c = bra.GetOrbital1();
    const basis::OrbitalStatePN& d = bra.GetOrbital2();

    // extract sector parameters
    assert(bra.J() == ket.J());
    int J = ket.J();

    // evaluate matrix element
    double matrix_element = 0.;
    matrix_element += IsospinScalarTBMESum(
        operator_type,
        overlap_orbital_space, overlap_sectors, overlap_matrices,
        c, d, a, b
      );
    int phase = - ParitySign(J-a.j()-b.j());
    matrix_element += phase * IsospinScalarTBMESum(
        operator_type,
        overlap_orbital_space, overlap_sectors, overlap_matrices,
        c, d, b, a
      );

    return matrix_element;
  }

  double IsospinDotTBMEProduct(
      IsospinOperatorType operator_type_1, IsospinOperatorType operator_type_2,
      // overlap matrix element data
      const basis::OrbitalSpaceLJPN& overlap_orbital_space,
      const basis::OrbitalSectorsLJPN& overlap_sectors,
      const basis::OperatorBlocks<double>& overlap_matrices,
      // two-body labels
      const basis::OrbitalStatePN& c, const basis::OrbitalStatePN& d,
      const basis::OrbitalStatePN& a, const basis::OrbitalStatePN& b,
      int J
    )
    // Evaluate the unsymmetrized matrix element (cd|T1*T2|ab) of a
    // product of one-body isospin operators using <b||T||a>.
    //
    // After csbasis (57).  See pjf "isospin" notes.
  {

    // short circuit check on triangularity of each one-body factor
    //
    // Note: This is redundant to (but preempts) the triangularity
    // checks in KinematicVectorOBRME.

    bool triangle_allowed = (
        ((c.l() == a.l()) && (c.j() == a.j()))
        && ((d.l() == b.l()) && (d.j() == b.j()))
      );
    if (!triangle_allowed)
      return 0.;

    // evaluate matrix element
    int racah_phase = ParitySign(d.j()+a.j()+J);
    double matrix_element = racah_phase * am::Wigner6J(c.j(),d.j(),J,b.j(),a.j(),0)
      * Hat(a.j()) * Hat(b.j())
      * IsospinOBME(operator_type_1, overlap_orbital_space, overlap_sectors, overlap_matrices, c, a)
      * IsospinOBME(operator_type_2, overlap_orbital_space, overlap_sectors, overlap_matrices, d, b);

    return matrix_element;
  }

  double IsospinDotTBME(
      IsospinOperatorType operator_type_1, IsospinOperatorType operator_type_2,
      // overlap matrix element data
      const basis::OrbitalSpaceLJPN& overlap_orbital_space,
      const basis::OrbitalSectorsLJPN& overlap_sectors,
      const basis::OperatorBlocks<double>& overlap_matrices,
      // two-body labels
      const basis::TwoBodyStateJJJPN& bra, const basis::TwoBodyStateJJJPN& ket
    )
  // Evaluate the matrix element <cd|T1.T2|ab> of a dot product of
  // one-body vector angular momentum operators (i.e., T = l, s, or
  // j), using the unsymmetrized matrix element (cd|T1.T2|ab).
  //
  // Computes <cd|T1.T2|ab>_AS for pp/nn states, or (cd|T1.T2|ab)_pn
  // for pn states.
  //
  // Note: For the V_(T1.T2) two-body term, a proton operator only has
  // nonvanishing pp sectors, and a neutron operator only has
  // nonvanishing nn sectors.
  //
  // After csbasis (55) and (56).  See mac "spin operator" notes page 2.
  {
    // extract orbitals
    const basis::OrbitalStatePN& a = ket.GetOrbital1();
    const basis::OrbitalStatePN& b = ket.GetOrbital2();
    const basis::OrbitalStatePN& c = bra.GetOrbital1();
    const basis::OrbitalStatePN& d = bra.GetOrbital2();

    // extract sector parameters
    assert(bra.J() == ket.J());
    int J = ket.J();

    double matrix_element = 0.;

    matrix_element
      += IsospinDotTBMEProduct(
          operator_type_1, operator_type_2,
          overlap_orbital_space, overlap_sectors, overlap_matrices,
          c, d, a, b,
          J
        );
    int phase = - ParitySign(J - a.j() - b.j());
    matrix_element += phase
      * IsospinDotTBMEProduct(
          operator_type_1, operator_type_2,
          overlap_orbital_space, overlap_sectors, overlap_matrices,
          c, d, b, a,
          J
        );

    return matrix_element;
  }

  Eigen::MatrixXd
  IsospinMatrixJJJPN(
      // overlap matrix element data
      const basis::OrbitalSpaceLJPN& overlap_orbital_space,
      const basis::OrbitalSectorsLJPN& overlap_sectors,
      const basis::OperatorBlocks<double>& overlap_matrices,
      IsospinOperatorType operator_type,
      const basis::TwoBodySectorsJJJPN::SectorType& sector,
      int A
    )
  {
    // set up aliases
    const basis::TwoBodySubspaceJJJPN& ket_subspace = sector.ket_subspace();
    const basis::TwoBodySubspaceJJJPN& bra_subspace = sector.bra_subspace();

    // generate matrix for sector
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(bra_subspace.size(), ket_subspace.size());

    // recover sector properties
    assert(bra_subspace.J() == ket_subspace.J());
    const int J = ket_subspace.J();
    assert(bra_subspace.g() == ket_subspace.g());
    const int g = ket_subspace.g();

    // for upper-triangular pairs of states in sector
    const std::size_t bra_subspace_size = bra_subspace.size();
    const int ket_subspace_size = ket_subspace.size();
    #pragma omp parallel for collapse(2) if (0)  // disabled until have chance to profile
    for (int bra_index = 0; bra_index < bra_subspace_size; ++bra_index) {
      for (int ket_index = 0; ket_index < ket_subspace_size; ++ket_index) {
        // diagonal sector: restrict to upper triangle
        if (sector.IsDiagonal() && !(bra_index <= ket_index))
          continue;

        // construct states
        basis::TwoBodyStateJJJPN bra(bra_subspace,bra_index);
        basis::TwoBodyStateJJJPN ket(ket_subspace,ket_index);

        // initalize one-body and two-body pieces
        double matrix_element_ob = 0;
        double matrix_element_tb = 0;

        // calculate matrix element
        matrix_element_ob
          += IsospinScalarTBME(
              operator_type,
              overlap_orbital_space, overlap_sectors, overlap_matrices,
              bra, ket
            );
        if (operator_type == IsospinOperatorType::kSquared) {
          matrix_element_tb
            += 2*IsospinDotTBME(
              IsospinOperatorType::kTz, IsospinOperatorType::kTz,
              overlap_orbital_space, overlap_sectors, overlap_matrices,
              bra, ket
            );
          matrix_element_tb
            += IsospinDotTBME(
              IsospinOperatorType::kRaising, IsospinOperatorType::kLowering,
              overlap_orbital_space, overlap_sectors, overlap_matrices,
              bra, ket
            );
          matrix_element_tb
            += IsospinDotTBME(
                IsospinOperatorType::kLowering, IsospinOperatorType::kRaising,
                overlap_orbital_space, overlap_sectors, overlap_matrices,
                bra, ket
              );
        }
        double matrix_element = 1./(A-1)*matrix_element_ob + matrix_element_tb;

        // convert to NAS if needed
        if (bra.two_body_species()!=basis::TwoBodySpeciesPN::kPN)
          if (bra.index1()==bra.index2())
            matrix_element *= 1/(sqrt(2.));
        if (ket.two_body_species()!=basis::TwoBodySpeciesPN::kPN)
          if (ket.index1()==ket.index2())
            matrix_element *= 1/(sqrt(2.));

        // store matrix element
        matrix(bra_index,ket_index) = matrix_element;
        // matrix(ket_index,bra_index) = matrix_element;  // lower triangle
      }
    }

    return matrix;
  }


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
