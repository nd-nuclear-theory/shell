/****************************************************************
  shell_separable.cpp

  Created by Mark A. Caprio, University of Notre Dame.
  Last modified 5/25/15.

****************************************************************/

#include <shell/shell_separable.h>
#include <cmath>

#include <am/wigner_gsl.h>

namespace shell {

  ////////////////////////////////////////////////////////////////
  // angular momentum calculation
  ////////////////////////////////////////////////////////////////


  // ShellAngularMomentumScalarOBME gives <b|T^2|a> for a squared
  // angular momentum operator, i.e. T^2 = l^2, s^2, or j^2.
  //
  // See mac "spin operator" notes page 3.

  double ShellAngularMomentumScalarOBME (const AngularMomentumType angular_momentum_operator,
					 const SPOrbitalNlj& b, const SPOrbitalNlj& a)
  {

    int na = a.Getn();
    int nb = b.Getn();
    int la = a.Getl();
    int lb = b.Getl();
    HalfInt ja = a.Getj();
    HalfInt jb = b.Getj();

    double matrix_element = 0.;

    if ( b == a )
      {
	if (angular_momentum_operator == kOrbital)
	  {
	    matrix_element += la*(la+1);
	  }
	else if (angular_momentum_operator == kSpin)
	  {
	    matrix_element += 3./4.;
	  }
	else if (angular_momentum_operator == kTotal)
	  {
 	    matrix_element += ja.DValue()*(ja.DValue()+1);
	  }
      }

    return matrix_element;
  }


  // ShellAngularMomentumVectorOBRME gives <b||T||a> for an angular
  // momentum operator, i.e., T = l, s, or j.
  //
  // Based on Suhonen "From nucleons to nucleus" (2.56) and (2.58).  See
  // mac "spin operator" notes page 3.
  
  double ShellAngularMomentumVectorOBRME (const AngularMomentumType angular_momentum_operator, 
  					  const SPOrbitalNlj& b, const SPOrbitalNlj& a)
  {
    
    int na = a.Getn();
    int nb = b.Getn();
    int la = a.Getl();
    int lb = b.Getl();
    HalfInt ja = a.Getj();
    HalfInt jb = b.Getj();
    
    double matrix_element = 0.;
    
    if ( (nb == na) && (lb == la) )
      {
  	if ((angular_momentum_operator == kOrbital) || (angular_momentum_operator == kTotal))
  	  {
  	    int phase = ParitySign(la+ja+HalfInt(3,2));
  	    matrix_element += Hat(jb)*Hat(ja)*sqrt(la*(la+1)*(2*la+1))*phase*Wigner6J(la,la,1,ja,jb,HalfInt(1,2));
  	  }
  	if ((angular_momentum_operator == kSpin) || (angular_momentum_operator == kTotal))
  	  {
  	    int phase = ParitySign(la+jb+HalfInt(3,2));
  	    matrix_element += sqrt(3./2.)*Hat(jb)*Hat(ja)*phase*Wigner6J(HalfInt(1,2),HalfInt(1,2),1,ja,jb,la);
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

  double ShellAngularMomentumScalarTBME (AngularMomentumType angular_momentum_operator,
					 TwoSpeciesStateType operator_species, 
					 TwoSpeciesStateType state_type, int J, 
					 const TwoBodyStateNlj& s2, const TwoBodyStateNlj& s1)
  {
    // extract orbitals
    SPOrbitalNlj a = s1.a1;
    SPOrbitalNlj b = s1.a2;
    SPOrbitalNlj c = s2.a1;
    SPOrbitalNlj d = s2.a2;

    double matrix_element = 0.;

    if (
	((state_type == kPP) || (state_type == kNN))
	&&
	((operator_species == state_type) || (operator_species == kPN))
	)
      {
	// like-nucleon case
	//
	// short circuited to only apply if operator is for same species or is total operator
        if (d == b)
	  matrix_element += ShellAngularMomentumScalarOBME(angular_momentum_operator, c, a);
        if (c == a)
	  matrix_element += ShellAngularMomentumScalarOBME(angular_momentum_operator, d, b);
        int phase = - ParitySign(J - a.Getj() - b.Getj());
        if (d == a)
	  matrix_element += phase * ShellAngularMomentumScalarOBME(angular_momentum_operator, c, b);
        if (c == b)
	  matrix_element += phase * ShellAngularMomentumScalarOBME(angular_momentum_operator, d, a);
      }
    else if (state_type == kPN)
      {
	// proton-neutron case
        if ((d == b) && (operator_species != kNN))
	  // term only contributes to proton or total operators, not neutron operator
	  matrix_element += ShellAngularMomentumScalarOBME(angular_momentum_operator, c, a);
        if ((c == a) && (operator_species != kPP))
	  // term only contributes to neutron or total operators, not proton operator
	  matrix_element += ShellAngularMomentumScalarOBME(angular_momentum_operator, d, b);
      }

    return matrix_element;
  	
  }

  // ShellAngularMomentumVectorDotTBMEProduct uses <b||T||a> to
  // compute the the matrix element to compute the unsymmetrized
  // matrix element (cd|T.T|ab), i.e., for a distinguishable-particle
  // direct product state, for an angular momentum operator, i.e., T =
  // l, s, or j.
  //
  // After csbasis (57).  See mac "spin operator" notes page 2.

  double ShellAngularMomentumVectorDotTBMEProduct (AngularMomentumType angular_momentum_operator, 
						   int J, 
						   const SPOrbitalNlj& c, const SPOrbitalNlj& d, 
						   const SPOrbitalNlj& a, const SPOrbitalNlj& b)
  {

    // short circuit check on triangularity of each one-body factor
    bool triangle_allowed = ( AllowedTriangle(a.Getj(),c.Getj(),1) && AllowedTriangle(b.Getj(),d.Getj(),1));
    if (!triangle_allowed )
      return 0.;

    // evaluate matrix element
    int phase = ParitySign(d.Getj() + a.Getj() + J);
    double matrix_element = phase * Wigner6J(c.Getj(),d.Getj(),J,b.Getj(),a.Getj(),1)
      * ShellAngularMomentumVectorOBRME(angular_momentum_operator,c,a) * ShellAngularMomentumVectorOBRME(angular_momentum_operator,d,b);
  	
    return matrix_element;
  }

  // ShellAngularMomentumVectorDotTBME uses <b||T||a> [via
  // (cd|T.T|ab)] to compute the TBMEs <cd|T1.T2|ab> of the "two-body
  // part" of a squared angular momentum operator.
  //
  // Computes <cd|V|ab>_AS for pp/nn states, or (cd|V|ab)_pn for pn
  //    states.
  //
  // After csbasis (55) and (56).  See mac "spin operator" notes page 2.

  double ShellAngularMomentumVectorDotTBME (AngularMomentumType angular_momentum_operator,
					    TwoSpeciesStateType operator_species, 
					    TwoSpeciesStateType state_type, int J, 
					    const TwoBodyStateNlj& s2, const TwoBodyStateNlj& s1)
  {
    
    SPOrbitalNlj a = s1.a1;
    SPOrbitalNlj b = s1.a2;
    SPOrbitalNlj c = s2.a1;
    SPOrbitalNlj d = s2.a2;
    
    //std::cout << "DEBUG: " << c << d << a << b << J << " " << state_type << std::endl;
    
    double matrix_element;

    // short circuit if this sector does not contribute to the operator of interest
    // 
    // For V_(T.T) two-body term, proton operator only has pp sector
    // and neutron operator only has nn sector.
    if (!((operator_species == state_type) || (operator_species == kPN)))
      return 0.;
			
    if ((state_type == kPP) || (state_type == kNN))
      {
	// like nucleon case
	matrix_element = ShellAngularMomentumVectorDotTBMEProduct(angular_momentum_operator, J, c, d, a, b);
	int phase = - ParitySign(J - a.Getj() - b.Getj());
	matrix_element += phase * ShellAngularMomentumVectorDotTBMEProduct(angular_momentum_operator, J, c, d, b, a);
      }
    else if (state_type == kPN)
      {
	// proton-neutron case
	matrix_element = ShellAngularMomentumVectorDotTBMEProduct(angular_momentum_operator, J, c, d, a, b);
      }
    
    return matrix_element;
  }

  // TwoBodyMatrixSectorAddAngularMomentum adds a multiple of a
  // squared angular momentum operator to a TBME matrix sector.
  //
  // Arguments
  //   scale: overall scale factor
  //   angular_momentum_operator: identifies momentum operator type (kOrbital, kSpin, kTotal)
  //   operator_species: whether operator is total or restricted (kPP for pure proton, nNN for 
  //     pure neutron, or kPN for *total* (pp+pn+nn) operator
  //   A: atomic mass
  //   destination_matrix: matrix to add to
  //   state_type, J, g: sector to add to

  void TwoBodyMatrixSectorAddAngularMomentum (double scale, 
					      AngularMomentumType angular_momentum_operator, 
					      TwoSpeciesStateType operator_species, 
					      int A,
					      TwoBodyMatrixNljTzJP& destination_matrix,
					      const SectorNljTzJP& sector)
  {
    // recover sector properties
    const TwoSpeciesStateType state_type = sector.GetStateType();
    const int J = sector.GetJ();
    const int g = sector.GetGrade();
    const TwoBodyBasisNljTzJP basis = destination_matrix.GetTwoBodyBasis();
    const int dimension = basis.GetDimension(state_type,J,g);

    // for canonical pairs of states in destination two-body space (in lexicographical order)
    for (int k1 = 0; k1 < dimension; ++k1)
      for (int k2 = k1; k2 < dimension; ++k2)
	{
	  // identify states
	  TwoBodyStateNlj s1 = basis.GetState(state_type,J,g,k1);
	  TwoBodyStateNlj s2 = basis.GetState(state_type,J,g,k2);

	  // calculate matrix element (pn or AS)
	  double matrix_element_t2 
	    = ShellAngularMomentumScalarTBME(angular_momentum_operator,operator_species,state_type,J,s2,s1);
	  double matrix_element_t1t2 
	    = ShellAngularMomentumVectorDotTBME(angular_momentum_operator,operator_species,state_type,J,s2,s1);
	  double matrix_element = 1./(A-1)*matrix_element_t2 + 2*matrix_element_t1t2;
	  matrix_element *= scale;

	  // store matrix element
	  double destination_matrix_element = destination_matrix.GetMatrixElementUNAS(state_type, s1, s2);
	  destination_matrix_element += matrix_element;
	  destination_matrix.SetMatrixElementUNAS(state_type, s1, s2, destination_matrix_element); 
	}
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
