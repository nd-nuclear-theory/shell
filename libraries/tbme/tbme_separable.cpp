/****************************************************************
  tbme_separable.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "am/wigner_gsl.h"
#include "basis/nlj_operator.h"
#include "fmt/format.h"  // for debugging
#include "tbme/tbme_separable.h"

namespace shell {
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

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
        const int subspace_size = subspace.size();
#pragma omp parallel for collapse(2) if (0)  // disabled until have chance to profile
        for (int bra_index = 0; bra_index < subspace_size; ++bra_index)
          for (int ket_index = 0; ket_index < subspace_size; ++ket_index)
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
  // kinematic operators -- scalar (T^2)
  ////////////////////////////////////////////////////////////////

  double KinematicScalarOBME(
      // radial matrix element data
      const basis::OrbitalSpaceLJPN& radial_orbital_space,
      const basis::OrbitalSectorsLJPN& radial_sectors,
      const basis::OperatorBlocks<double>& radial_matrices,
      // one-body labels
      const basis::OrbitalStatePN& b, const basis::OrbitalStatePN& a
    )
  // Evaluate <b|T^2|a> for a scalar one-body operator T^2, i.e., r^2
  // or k^2.
  //
  // <b|T^2|a> = <R_b|T^2|R_a>[(la,ja)==(lb,jb)]
  //
  // See, e.g., below csbasis (39).
  {

    // int na = a.n();
    // int nb = b.n();
    int la = a.l();
    int lb = b.l();
    HalfInt ja = a.j();
    HalfInt jb = b.j();

    double matrix_element = 0.;
    if ( (lb == la) && (jb == ja) )
      {
        matrix_element += basis::MatrixElementLJPN(
            radial_orbital_space,radial_orbital_space,radial_sectors,radial_matrices,
            b,a
          );
      }

    return matrix_element;
  }

  double KinematicScalarTBMEProduct(
      // radial matrix element data
      const basis::OrbitalSpaceLJPN& radial_orbital_space,
      const basis::OrbitalSectorsLJPN& radial_sectors,
      const basis::OperatorBlocks<double>& radial_matrices,
      // two-body labels
      const basis::OrbitalStatePN& c, const basis::OrbitalStatePN& d,
      const basis::OrbitalStatePN& a, const basis::OrbitalStatePN& b
    )
  // Evaluate the unsymmetrized matrix element (cd;J|V_(T^2)|ab;J) of
  // the "upgraded" two-body operator obtained from a scalar one-body
  // operator T^2, i.e., r^2 or k^2.
  //
  // See csbasis (52), which is written for a pn state but applies to
  // any unsymmetrized product state.
  {
    // short circuit check equality of (l,j) on each one-body factor
    //
    // Note: This is redundant to (but preempts) the (l,j) equality
    // check in KinematicScalarOBME.
    bool triangle_allowed = (
        ((c.l()==a.l())&&(c.j()==a.j()))
        && ((d.l()==b.l())&&(d.j()==b.j()))
      );
    if (!triangle_allowed )
      return 0.;

    // evaluate matrix element
    double matrix_element = 0.;
    if (d == b)
      matrix_element += KinematicScalarOBME(
          radial_orbital_space,radial_sectors,radial_matrices,
          c,a
        );
    if (c == a)
      matrix_element += KinematicScalarOBME(
          radial_orbital_space,radial_sectors,radial_matrices,
          d,b
        );

    return matrix_element;
  }

  double KinematicScalarTBME(
      // radial matrix element data
      const basis::OrbitalSpaceLJPN& radial_orbital_space,
      const basis::OrbitalSectorsLJPN& radial_sectors,
      const basis::OperatorBlocks<double>& radial_matrices,
      // two-body labels
      const basis::TwoBodyStateJJJPN& bra, const basis::TwoBodyStateJJJPN& ket
    )
  // Evaluate the full (possibly antisymmetrized) matrix element
  // <cd;J|V_(T^2)|ab;J> for the "upgraded" two-body operator obtained
  // from a scalar one-body operator T^2, i.e., r^2 or k^2.
  //
  // Computes <cd;J|V|ab;J>_AS for pp/nn states, or <cd;J|V|ab;J>_pn
  // for pn states.
  //
  // See csbasis (52) and (54).
  {
    // extract orbitals
    const basis::OrbitalStatePN& a = ket.GetOrbital1();
    const basis::OrbitalStatePN& b = ket.GetOrbital2();
    const basis::OrbitalStatePN& c = bra.GetOrbital1();
    const basis::OrbitalStatePN& d = bra.GetOrbital2();

    // extract sector parameters
    assert(bra.two_body_species()==ket.two_body_species());
    basis::TwoBodySpeciesPN two_body_species = ket.two_body_species();
    assert(bra.J()==ket.J());
    int J = ket.J();

    // evaluate matrix element
    double matrix_element = 0.;
    matrix_element += KinematicScalarTBMEProduct(
          radial_orbital_space,radial_sectors,radial_matrices,
          c,d,a,b
      );
    if (two_body_species != basis::TwoBodySpeciesPN::kPN)
      {
        int phase = - ParitySign(J-a.j()-b.j());
        matrix_element += phase * KinematicScalarTBMEProduct(
            radial_orbital_space,radial_sectors,radial_matrices,
            c,d,b,a
          );
      }

    return matrix_element;
  }

  ////////////////////////////////////////////////////////////////
  // kinematic operators -- vector (T1.T2)
  ////////////////////////////////////////////////////////////////

  double KinematicVectorOBRME(
      // radial matrix element data
      const basis::OrbitalSpaceLJPN& radial_orbital_space,
      const basis::OrbitalSectorsLJPN& radial_sectors,
      const basis::OperatorBlocks<double>& radial_matrices,
      // one-body labels
      const basis::OrbitalStatePN& b, const basis::OrbitalStatePN& a
    )
  // Evaluate <b||T||a> for a vector one-body operator T, i.e., r or
  // k.
  //
  // See csbasis (58)-(60).  Omits the "extra" phase factor on the
  // momentum operator in csbasis (59), to avoid carrying around
  // imaginary factors, as well as for uniformity between the r and k
  // cases.  This phase factor must instead be taken into account
  // later by the calling function KinematicVectorDotTBMEProduct
  // when calculating the TBME.
  //
  // QUERY (mac): Why the extra factor Hat(ja)???  Without that
  // factor, this RME would be in Edmonds convention.  But I think
  // this factor is the wrong way to go to group theory convention.
  {

    // int na = a.n();
    // int nb = b.n();
    int la = a.l();
    int lb = b.l();
    HalfInt ja = a.j();
    HalfInt jb = b.j();

    double matrix_element = 0.;

    if (
        am::AllowedTriangle(la,1,lb)  // might be redundant to j and parity checks?
        && am::AllowedTriangle(ja,1,jb)
        && ((la+lb+1)%2==0)
      )
      {
        double radial_matrix_element = basis::MatrixElementLJPN(
              radial_orbital_space,radial_orbital_space,radial_sectors,radial_matrices,
              b,a
            );
        matrix_element += ParitySign(jb-ja+1) * Hat(ja)
          * am::ClebschGordan(ja,HalfInt(1,2),1,0,jb,HalfInt(1,2))
          * radial_matrix_element;
      }

    return matrix_element;
  }

  double KinematicVectorDotTBMEProduct(
      // radial matrix element data
      const basis::OrbitalSpaceLJPN& radial_orbital_space,
      const basis::OrbitalSectorsLJPN& radial_sectors,
      const basis::OperatorBlocks<double>& radial_matrices,
      // mode
      bool momentum_space,
      // two-body labels
      const basis::OrbitalStatePN& c, const basis::OrbitalStatePN& d,
      const basis::OrbitalStatePN& a, const basis::OrbitalStatePN& b,
      int J
    )
  // Evaluate the unsymmetrized matrix element (cd|T1.T2|ab) of a dot
  // product of one-body vector operators (i.e., T = r or k), using
  // <b||T||a>.
  //
  // See csbasis (57).  Incorporates momentum space phase factor from
  // csbasis (59), which was omitted from KinematicVectorOBRME to
  // avoid carrying around imaginary factors.
  {

    // short circuit check on triangularity of each one-body factor
    //
    // Note: This is redundant to (but preempts) the triangularity
    // checks in KinematicVectorOBRME.

    bool triangle_allowed = (
        (am::AllowedTriangle(c.l(),1,a.l()) && am::AllowedTriangle(c.j(),1,a.j()) && ((c.l()+1+a.l())%2==0))
        &&
        (am::AllowedTriangle(d.l(),1,b.l()) && am::AllowedTriangle(d.j(),1,b.j()) && ((d.l()+1+b.l())%2==0))
      );
    if (!triangle_allowed )
      return 0.;

    // evaluate matrix element
    int racah_phase = ParitySign(d.j()+a.j()+J);
    double matrix_element = racah_phase * am::Wigner6J(c.j(),d.j(),J,b.j(),a.j(),1)
      * KinematicVectorOBRME(radial_orbital_space,radial_sectors,radial_matrices,c,a)
      * KinematicVectorOBRME(radial_orbital_space,radial_sectors,radial_matrices,d,b);

    // momentum space phase factor
    if (momentum_space)
      matrix_element *= ParitySign((c.l()+d.l()-a.l()-b.l())/2);

    // std::cout
    //   << std::endl
    //   << fmt::format(
    //       "KinematicVectorDotTBMEProduct {} {} {} {} J {} => {}",
    //       c.LabelStr(),d.LabelStr(),a.LabelStr(),b.LabelStr(),J,
    //       matrix_element
    //     )
    //   << std::endl;

    return matrix_element;
  }

  double KinematicVectorDotTBME(
      // two-body labels
      const basis::OrbitalSpaceLJPN& radial_orbital_space,
      const basis::OrbitalSectorsLJPN& radial_sectors,
      const basis::OperatorBlocks<double>& radial_matrices,
      // mode
      bool momentum_space,
      // two-body labels
      const basis::TwoBodyStateJJJPN& bra, const basis::TwoBodyStateJJJPN& ket
    )
  // Evaluate the full (possibly antisymmetrized) matrix element
  // <cd|T1.T2|ab> of a dot product of one-body vector operators
  // (i.e., T = r or k), using the unsymmetrized matrix element
  // (cd|T1.T2|ab).
  //
  // Computes <cd;J|T1.T2|ab;J>_AS for pp/nn states, or
  // <cd;J|T1.T2|ab;J>_pn for pn states.
  //
  // See csbasis (55) and (56).
  {
    // extract orbitals
    const basis::OrbitalStatePN& a = ket.GetOrbital1();
    const basis::OrbitalStatePN& b = ket.GetOrbital2();
    const basis::OrbitalStatePN& c = bra.GetOrbital1();
    const basis::OrbitalStatePN& d = bra.GetOrbital2();

    // extract sector parameters
    assert(bra.two_body_species()==ket.two_body_species());
    basis::TwoBodySpeciesPN two_body_species = ket.two_body_species();
    assert(bra.J()==ket.J());
    int J = ket.J();

    // std::cout
    //   << std::endl
    //   << fmt::format(
    //       "KinematicVectorDotTBME {} {} {} {} J {} : ({},{})",
    //       c.LabelStr(),d.LabelStr(),a.LabelStr(),b.LabelStr(),J,
    //       bra.index(),ket.index()
    //     )
    //   << std::endl;

    // evaluate matrix element
    double matrix_element = 0.;
    matrix_element += KinematicVectorDotTBMEProduct(
        radial_orbital_space,radial_sectors,radial_matrices,
        momentum_space,
        c,d,a,b,J
      );
    if (two_body_species!=basis::TwoBodySpeciesPN::kPN)
      {
        int phase = -ParitySign(J-a.j()-b.j());
        matrix_element += phase * KinematicVectorDotTBMEProduct(
            radial_orbital_space,radial_sectors,radial_matrices,
            momentum_space,
            c,d,b,a,J
          );
      }

    //std::cout
    //  << std::endl
    //  << fmt::format(
    //      "   => {}",
    //      matrix_element
    //    )
    //  << std::endl;

    return matrix_element;
  }

  Eigen::MatrixXd
  KinematicMatrixJJJPN(
      const basis::OrbitalSpaceLJPN& radial_orbital_space,
      const basis::OrbitalSectorsLJPN& radial_sectors,
      const basis::OperatorBlocks<double>& radial_matrices,
      KinematicOperatorType kinematic_operator_type,
      shell::RadialOperatorType radial_operator_type,
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
    bool momentum_space = (radial_operator_type == shell::RadialOperatorType::kK);

    // for upper-triangular pairs of states in sector
    const int subspace_size = subspace.size();
#pragma omp parallel for collapse(2) if (0)  // disabled until have chance to profile
    for (int bra_index = 0; bra_index < subspace_size; ++bra_index)
      for (int ket_index = 0; ket_index < subspace_size; ++ket_index)
        {

          // diagonal sector: restrict to upper triangle
          // if (sector.IsDiagonal())
          if (!(bra_index<=ket_index))
            continue;

          // construct states
          basis::TwoBodyStateJJJPN bra(subspace,bra_index);
          basis::TwoBodyStateJJJPN ket(subspace,ket_index);

          // calculate matrix element (pn or AS)
          double matrix_element;
          if (kinematic_operator_type==KinematicOperatorType::kUTSqr)
            matrix_element = 1./(A-1)*KinematicScalarTBME(
                radial_orbital_space,radial_sectors,radial_matrices,
                bra,ket
              );
          else // if (kinematic_operator_type==KinematicOperatorType::kVT1T2)
            matrix_element = KinematicVectorDotTBME(
                radial_orbital_space,radial_sectors,radial_matrices,
                momentum_space,
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
      shell::AngularMomentumOperatorFamily operator_family,
      shell::AngularMomentumOperatorSpecies operator_species,
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
    const int subspace_size = subspace.size();
#pragma omp parallel for collapse(2) if (0)  // disabled until have chance to profile
    for (int bra_index = 0; bra_index < subspace_size; ++bra_index)
      for (int ket_index = 0; ket_index < subspace_size; ++ket_index)
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
    const int bra_subspace_size = bra_subspace.size();
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
