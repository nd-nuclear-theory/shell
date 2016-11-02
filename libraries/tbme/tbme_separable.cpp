/****************************************************************
  tbme_separable.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "am/wigner_gsl.h"
// #include "cppformat/format.h"  // for debugging
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
  // kinematic operators
  ////////////////////////////////////////////////////////////////

  double ShellKinematicScalarOBME(
      const basis::OrbitalSpaceLJPN& radial_orbital_space,
      const basis::OrbitalSectorsLJPN& radial_sectors,
      const basis::MatrixVector& radial_matrices,
      const basis::OrbitalStatePN& b, const basis::OrbitalStatePN& a
    )
  // Evaluate <b|T^2|a> for a scalar radial operator T^2, i.e., r^2 or
  // k^2.
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
            radial_orbital_space,radial_orbital_space,
            radial_sectors,
            radial_matrices,
            b,a
          );
      }

    return matrix_element;
  }

  ////////////////////////////////////////////////////////////////
  // angular momentum operators
  ////////////////////////////////////////////////////////////////

  double ShellAngularMomentumScalarOBME(
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


  // ShellAngularMomentumVectorOBRME gives <b||T||a> for an angular
  // momentum operator, i.e., T = l, s, or j.
  //
  // Based on Suhonen "From nucleons to nucleus" (2.56) and (2.58).  See
  // mac "spin operator" notes page 3.
  
  double ShellAngularMomentumVectorOBRME(
      shell::AngularMomentumOperatorFamily operator_family,
      const basis::OrbitalStatePN& b, const basis::OrbitalStatePN& a
    )
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

  // ShellAngularMomentumScalarTBME uses <b|T^2|a> to compute the TBMEs
  // <cd|V_(T^2)|ab> of the "upgraded one-body part" of a squared angular
  // momentum operator.
  //
  // Computes <cd|V|ab>_AS for pp/nn states, or (cd|V|ab)_pn for pn
  //    states.
  //
  // After csbasis (52) and (54).  See mac "spin operator" notes page 3.

  double ShellAngularMomentumScalarTBME (
      shell::AngularMomentumOperatorFamily operator_family,
      shell::AngularMomentumOperatorSpecies operator_species,
      basis::TwoBodySpeciesPN two_body_species,
      int J, 
      const basis::TwoBodyStateJJJPN& s2, const basis::TwoBodyStateJJJPN& s1
    )
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
	  matrix_element += ShellAngularMomentumScalarOBME(operator_family,c,a);
        if (c == a)
	  matrix_element += ShellAngularMomentumScalarOBME(operator_family,d,b);
        int phase = - ParitySign(J - a.j() - b.j());
        if (d == a)
	  matrix_element += phase * ShellAngularMomentumScalarOBME(operator_family,c,b);
        if (c == b)
	  matrix_element += phase * ShellAngularMomentumScalarOBME(operator_family,d,a);
      }
    else if (two_body_species == basis::TwoBodySpeciesPN::kPN)
      {
	// proton-neutron case
        if ((d == b) && (operator_species != shell::AngularMomentumOperatorSpecies::kN))
	  // term only contributes to proton or total operators, not neutron operator
	  matrix_element += ShellAngularMomentumScalarOBME(operator_family,c,a);
        if ((c == a) && (operator_species != shell::AngularMomentumOperatorSpecies::kP))
	  // term only contributes to neutron or total operators, not proton operator
	  matrix_element += ShellAngularMomentumScalarOBME(operator_family,d,b);
      }

    return matrix_element;
  	
  }

  // ShellAngularMomentumVectorDotTBMEProduct uses <b||T||a> to
  // compute the unsymmetrized matrix element (cd|T.T|ab), i.e., for a
  // distinguishable-particle direct product state, for an angular
  // momentum operator, i.e., T = l, s, or j.
  //
  // After csbasis (57).  See mac "spin operator" notes page 2.

  double ShellAngularMomentumVectorDotTBMEProduct(
      shell::AngularMomentumOperatorFamily operator_family, 
      int J, 
      const basis::OrbitalStatePN& c, const basis::OrbitalStatePN& d,
      const basis::OrbitalStatePN& a, const basis::OrbitalStatePN& b
    )
  {

    // short circuit check on triangularity of each one-body factor
    bool triangle_allowed = ( am::AllowedTriangle(a.j(),c.j(),1) && am::AllowedTriangle(b.j(),d.j(),1));
    if (!triangle_allowed )
      return 0.;

    // evaluate matrix element
    int phase = ParitySign(d.j() + a.j() + J);
    double matrix_element = phase * am::Wigner6J(c.j(),d.j(),J,b.j(),a.j(),1)
      * ShellAngularMomentumVectorOBRME(operator_family,c,a) * ShellAngularMomentumVectorOBRME(operator_family,d,b);
  	
    return matrix_element;
  }

  // ShellAngularMomentumVectorDotTBME uses <b||T||a> [via
  // (cd|T.T|ab)] to compute the TBMEs <cd|T1.T2|ab> of the "two-body
  // part" of a squared angular momentum operator.
  //
  // Computes <cd|V|ab>_AS for pp/nn states, or (cd|V|ab)_pn for pn
  //    states.
  //
  // For V_(T.T) two-body term, proton operator only has pp sector
  // and neutron operator only has nn sector.
  //
  // After csbasis (55) and (56).  See mac "spin operator" notes page 2.

  double ShellAngularMomentumVectorDotTBME(
      shell::AngularMomentumOperatorFamily operator_family,
      shell::AngularMomentumOperatorSpecies operator_species, 
      basis::TwoBodySpeciesPN two_body_species,
      int J, 
      const basis::TwoBodyStateJJJPN& s2, const basis::TwoBodyStateJJJPN& s1
    )
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
	matrix_element += ShellAngularMomentumVectorDotTBMEProduct(operator_family,J,c,d,a,b);
	int phase = - ParitySign(J - a.j() - b.j());
	matrix_element += phase * ShellAngularMomentumVectorDotTBMEProduct(operator_family,J,c,d,b,a);
      }
    else if (
        (two_body_species == basis::TwoBodySpeciesPN::kPN)
        && (operator_species == shell::AngularMomentumOperatorSpecies::kTotal)
      )
      {
	// proton-neutron case
	matrix_element += ShellAngularMomentumVectorDotTBMEProduct(operator_family,J,c,d,a,b);
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
    for (int bra_index = 0; bra_index < subspace.size(); ++bra_index)
      for (int ket_index = 0; ket_index < subspace.size(); ++ket_index)
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
	    = ShellAngularMomentumScalarTBME(operator_family,operator_species,two_body_species,J,bra,ket);
	  double matrix_element_t1t2 
	    = ShellAngularMomentumVectorDotTBME(operator_family,operator_species,two_body_species,J,bra,ket);
	  double matrix_element = 1./(A-1)*matrix_element_t2 + 2*matrix_element_t1t2;

          // convert to NAS if needed
          if (subspace.two_body_species()!=basis::TwoBodySpeciesPN::kPN)
            {
              double conversion_factor = 1.;
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
  ////////////////////////////////////////////////////////////////
} // namespace
