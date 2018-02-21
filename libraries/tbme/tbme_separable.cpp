/****************************************************************
  tbme_separable.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "am/wigner_gsl.h"
#include "basis/nlj_operator.h"
#include "cppformat/format.h"  // for debugging
#include "tbme/tbme_separable.h"

namespace shell {
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  // upgraded one-body operator
  ////////////////////////////////////////////////////////////////

  basis::OperatorBlock<double>
  UpgradeOneBodyOperatorJJJPN(basis::OneBodyOperatorLJPN ob_operator,
                              const basis::TwoBodySectorsJJJPN::SectorType& sector,
                              int A
                             )
  {
    const basis::TwoBodySubspaceJJJPN& bra_subspace = sector.bra_subspace();
    const int bra_subspace_size = bra_subspace.size();
    const basis::TwoBodySubspaceJJJPN& ket_subspace = sector.ket_subspace();
    const int ket_subspace_size = ket_subspace.size();

    basis::OperatorBlock<double> matrix =
      basis::OperatorBlock<double>::Zero(bra_subspace_size, ket_subspace_size);

    // short circuit on triangle constraint
    if (!am::AllowedTriangle(bra_subspace.J(), ket_subspace.J(), ob_operator.sectors.j0()))
      return matrix;

    // build sector matrix
    for (int bra_index = 0; bra_index < bra_subspace_size; ++bra_index)
      for (int ket_index = 0; ket_index < ket_subspace_size; ++ket_index)
      {
        // construct states
        const basis::TwoBodyStateJJJPN bra(bra_subspace, bra_index);
        const basis::TwoBodyStateJJJPN ket(ket_subspace, ket_index);
        const basis::OrbitalStatePN& a = ket.GetOrbital1();
        const basis::OrbitalStatePN& b = ket.GetOrbital2();
        const basis::OrbitalStatePN& c = bra.GetOrbital1();
        const basis::OrbitalStatePN& d = bra.GetOrbital2();

        // calculate AS two-body element
        double matrix_element = 0.;
        if (d == b)
          matrix_element += ob_operator.get_matrix_element(c, a);
        if (c == a)
          matrix_element += ob_operator.get_matrix_element(d, b);

        int phase = - ParitySign(ket.J()-a.j()-b.j());
        if (d == a)
          matrix_element += phase*ob_operator.get_matrix_element(c, b);
        if (c == b)
          matrix_element += phase*ob_operator.get_matrix_element(d, a);

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
  double RacahDotProductFormula(const basis::OneBodyOperatorLJPN& ob_operator1,
                                const basis::OneBodyOperatorLJPN& ob_operator2,
                                int bra_J, int ket_J,
                                const basis::OrbitalStatePN& c, const basis::OrbitalStatePN& d,
                                const basis::OrbitalStatePN& a, const basis::OrbitalStatePN& b
                               )
  {
    int operator_j0 = ob_operator1.sectors.j0();
    // evaluate Racah reduction formula
    double matrix_element;
    matrix_element = ob_operator1.get_matrix_element(c, a) * ob_operator2.get_matrix_element(d, b);
    // TODO(pjf): enable below Hat when converting to Edmonds convention inside h2mixer
    // matrix_element *= Hat(bra.J())
    matrix_element *= ParitySign(d.j()+bra_J+a.j());
    matrix_element *= am::Wigner6J(c.j(), d.j(), bra_J, b.j(), a.j(), operator_j0);

    return matrix_element;
  }

  basis::OperatorBlock<double>
  RacahReduceDotProductJJJPN(const basis::OneBodyOperatorLJPN& ob_operator1,
                             const basis::OneBodyOperatorLJPN& ob_operator2,
                             const basis::TwoBodySectorsJJJPN::SectorType& sector)
  {
    // sanity check on
    assert(ob_operator1.sectors.j0() == ob_operator2.sectors.j0());
    int operator_j0 = ob_operator1.sectors.j0();

    const basis::TwoBodySubspaceJJJPN& bra_subspace = sector.bra_subspace();
    const int bra_subspace_size = bra_subspace.size();
    const basis::TwoBodySubspaceJJJPN& ket_subspace = sector.ket_subspace();
    const int ket_subspace_size = ket_subspace.size();

    basis::OperatorBlock<double> matrix =
      basis::OperatorBlock<double>::Zero(bra_subspace_size, ket_subspace_size);

    // short circuit on triangle constraint
    if (bra_subspace.J() != ket_subspace.J())
      return matrix;

    // build sector matrix
    for (int bra_index = 0; bra_index < bra_subspace_size; ++bra_index)
      for (int ket_index = 0; ket_index < ket_subspace_size; ++ket_index)
      {
        // construct states
        const basis::TwoBodyStateJJJPN bra(bra_subspace, bra_index);
        const basis::TwoBodyStateJJJPN ket(ket_subspace, ket_index);
        const basis::OrbitalStatePN& a = ket.GetOrbital1();
        const basis::OrbitalStatePN& b = ket.GetOrbital2();
        const basis::OrbitalStatePN& c = bra.GetOrbital1();
        const basis::OrbitalStatePN& d = bra.GetOrbital2();

        int phase = - ParitySign(ket.J()-a.j()-b.j());  // AS phase
        double matrix_element =
          RacahDotProductFormula(ob_operator1, ob_operator2, bra.J(), ket.J(), c, d, a, b)
          + phase*RacahDotProductFormula(ob_operator1, ob_operator2, bra.J(), ket.J(), c, d, b, a);

        // convert to NAS if needed
        if (a == b)
          matrix_element *= 1/(sqrt(2.));
        if (c == d)
          matrix_element *= 1/(sqrt(2.));

        matrix(bra_index, ket_index) = matrix_element;
      }
    return matrix;
  }

  double RacahTensorProductFormula(const basis::OneBodyOperatorLJPN& ob_operator1,
                                   const basis::OneBodyOperatorLJPN& ob_operator2,
                                   int bra_J, int ket_J, int J0,
                                   const basis::OrbitalStatePN& c, const basis::OrbitalStatePN& d,
                                   const basis::OrbitalStatePN& a, const basis::OrbitalStatePN& b
                                  )
  {
    int ob_operator1_j0 = ob_operator1.sectors.j0();
    int ob_operator2_j0 = ob_operator2.sectors.j0();
    double matrix_element;
    matrix_element = ob_operator1.get_matrix_element(c, a) * ob_operator2.get_matrix_element(d, b);
    // TODO(pjf): enable below Hat when converting to Edmonds convention
    // matrix_element *= Hat(bra_J)
    matrix_element *= Hat(J0) * Hat(ket_J);
    matrix_element *= am::Wigner9J(c.j(), d.j(), bra_J, a.j(), b.j(), ket_J, ob_operator1_j0, ob_operator2_j0, J0);

    return matrix_element;
  }

  basis::OperatorBlock<double>
  RacahReduceTensorProductJJJPN(const basis::OneBodyOperatorLJPN ob_operator1,
                                const basis::OneBodyOperatorLJPN ob_operator2,
                                const basis::TwoBodySectorsJJJPN::SectorType& sector,
                                int J0
                               )
  {
    // sanity check on
    assert(am::AllowedTriangle(ob_operator1.sectors.j0(), ob_operator2.sectors.j0(), J0));

    const basis::TwoBodySubspaceJJJPN& bra_subspace = sector.bra_subspace();
    const int bra_subspace_size = bra_subspace.size();
    const basis::TwoBodySubspaceJJJPN& ket_subspace = sector.ket_subspace();
    const int ket_subspace_size = ket_subspace.size();

    basis::OperatorBlock<double> matrix =
      basis::OperatorBlock<double>::Zero(bra_subspace_size, ket_subspace_size);

    // short circuit on triangle constraint
    if (!am::AllowedTriangle(bra_subspace.J(), ket_subspace.J(), J0))
      return matrix;

    // build sector matrix
    for (int bra_index = 0; bra_index < bra_subspace_size; ++bra_index)
      for (int ket_index = 0; ket_index < ket_subspace_size; ++ket_index)
      {
        // construct states
        const basis::TwoBodyStateJJJPN bra(bra_subspace, bra_index);
        const basis::TwoBodyStateJJJPN ket(ket_subspace, ket_index);
        const basis::OrbitalStatePN& a = ket.GetOrbital1();
        const basis::OrbitalStatePN& b = ket.GetOrbital2();
        const basis::OrbitalStatePN& c = bra.GetOrbital1();
        const basis::OrbitalStatePN& d = bra.GetOrbital2();

        int phase = - ParitySign(ket.J()-a.j()-b.j());  // AS phase
        // evaluate Racah reduction formula on AS states
        double matrix_element =
          RacahTensorProductFormula(ob_operator1, ob_operator2, bra.J(), ket.J(), J0, c, d, a, b)
          + phase*RacahTensorProductFormula(ob_operator1, ob_operator2, bra.J(), ket.J(), J0, c, d, b, a);

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
  // kinematic operators
  ////////////////////////////////////////////////////////////////

  Eigen::MatrixXd
  KinematicMatrixJJJPN(
      const shell::OneBodyOperator& radial_operator,
      KinematicOperatorType kinematic_operator_type,
      const basis::TwoBodySectorsJJJPN::SectorType& sector,
      int A
    )
  {
    // set up aliases
    assert(sector.IsDiagonal());
    const basis::TwoBodySubspaceJJJPN& subspace = sector.ket_subspace();

    // recover sector properties
    const basis::TwoBodySpeciesPN two_body_species = subspace.two_body_species();
    const int J = subspace.J();
    const int g = subspace.g();
    int momentum_space_phase = ParitySign(radial_operator.radial_operator_type == shell::RadialOperatorType::kK);

    basis::OperatorBlock<double> matrix;
    if (kinematic_operator_type == KinematicOperatorType::kUTSqr)
    {
      matrix = UpgradeOneBodyOperatorJJJPN(radial_operator, sector, A);
    }
    else // if (kinematic_operator_type==KinematicOperatorType::kVT1T2)
    {
      matrix = momentum_space_phase*RacahReduceDotProductJJJPN(radial_operator, radial_operator, sector);
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
