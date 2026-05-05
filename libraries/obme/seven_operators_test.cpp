/****************************************************************
  seven_operators.cpp

  Vicor Duménil
  University of Notre Dame & LPC Caen



  [1] : Multipole operators in semileptonic weak and electromagnetic interactions 
        with nuclei: Harmonic oscillator single-particle matrix elements, 
        TW Donnelly, WC Haxton - Atomic Data and Nuclear Data Tables, 1979
  [2] : SevenOperators, a Mathematica script for harmonic oscillator nuclear 
        matrix elements arising in semileptonic electroweak interactions
        W Haxton, C Lunardini, Computer Physics Communications, 2008
  [3] : Semileptonic weak and electromagnetic interactions with nuclei: Nuclear 
        current operators through order (v/c) nucleon2
        BD Serot - Nuclear Physics A, 1978
        
****************************************************************/

#include <cmath>
#include <fstream>
#include <iostream>

#include "spline/spline.h"
#include "spline/spline_me.h"

#include "am/wigner_gsl_twice.h"
#include "am/halfint.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>


////////////////////////////////////////////////////////////////
// test code
////////////////////////////////////////////////////////////////

// Spherical Bessel function j_L(x)
struct SphericalBesselParams
// Parameters for a spherical Bessel function j_L(q*r).
//
// Fields:
//   order (int): angular momentum order L of j_L
//   q (double): wave number q (units inverse of the length parameter b)
{
  int order;
  double q;
};


double SphericalBesselFunctionEval(double r, void * p)
  {
    SphericalBesselParams* params = static_cast< SphericalBesselParams*>(p);
    int order = (params->order);
    double q = (params->q);

    return gsl_sf_bessel_jl(order, q*r);
 }


double y_var(double q, double b)
{
	return std::pow(q * b / 2, 2);
}

/*double Hat(int j)
{
	return std::sqrt(2 * j + 1);
}*/

// Test BesselMatrixElement
double Test_BesselMatrixElement(int ni, int li, int bi, int nf, int lf, int bf, int L, double q) 
// Calculate the matrix element :
//   <n' l' | j_L(qr) | n l>
{
    double me_tot ;
    
	std::cout << "Test BesselMatrixElement" << std::endl;
	double me1 ; 

	gsl_function f;
	struct SphericalBesselParams params = {L, q};
	f.function = &SphericalBesselFunctionEval;
	f.params = &params;

	me1 = spline::RadialMatrixElementOfFunction(
	    ni, li, bi, spline::BasisType::kOscillator,
	    nf, lf, bf, spline::BasisType::kOscillator,
	    spline::OperatorType::kR, &f
	  );
	std::cout << "<" << nf << "lf:" << lf << "| " << "j_" << L << "(qr)" << " |" << ni << "li:" << li << ">" << std::endl; 
	std::cout << "Using RadialMatrixElementOfFunction: " << me1 << std::endl;

	me_tot = me1 ;
	std::cout << "Total matrix Element : " << me_tot << "\n" << std::endl;

	return me_tot;
}

// Test BesselMatrixElement_Minus
double Test_BesselMatrixElement_Minus(int ni, int li, int bi, int nf, int lf, int bf, int L, double q)
// Calculate the matrix element :
//   <n' l' | j_L(qr)(d_r - l/r) | n l>
{
    double me_tot ;
    
	std::cout << "Test BesselMatrixElement_Minus" << std::endl;
	if (ni==0) {
		double me1 ;

		double y = y_var(q, bi);
		double prefactor = -std::pow(8*y, -0.5);
		double prefactor_me1 = am::Hat2(2*(li+1)); 

		gsl_function f;
		struct SphericalBesselParams params = {L, q};
		f.function = &SphericalBesselFunctionEval;
		f.params = &params;

		me1 = spline::RadialMatrixElementOfFunction(
			    ni, li+1, bi, spline::BasisType::kOscillator,
			    nf, lf, bf, spline::BasisType::kOscillator,
			    spline::OperatorType::kR, &f
			  );
		std::cout << "<" << nf << "lf:" << lf << "| " << "j_" << L << "(qr)" << " |" << ni << "li:" << li << ">" << std::endl; 
		std::cout << "Using RadialMatrixElementOfFunction: " << me1 << std::endl;

		me_tot = prefactor * prefactor_me1 * me1;
		std::cout << "Total matrix Element : " << me_tot << "\n" << std::endl;
		
	} 
	else if (ni==1) {
		double me1, me2 ;

		double y = y_var(q, bi);
		double prefactor = -std::pow(8*y, -0.5) / std::sqrt(2);
		double prefactor_me1 = am::Hat2(2*(li+1)) * am::Hat2(2*(li+1)) ;
		double prefactor_me2 = am::Hat2(2*(li+2)) * am::Hat2(2*(li+3));

		gsl_function f;
		struct SphericalBesselParams params = {L, q};
		f.function = &SphericalBesselFunctionEval;
		f.params = &params;

		me1 = spline::RadialMatrixElementOfFunction(
					    ni, li+1, bi, spline::BasisType::kOscillator,
					    nf, lf, bf, spline::BasisType::kOscillator,
					    spline::OperatorType::kR, &f
					  );

		me2 = spline::RadialMatrixElementOfFunction(
							    ni, li+3, bi, spline::BasisType::kOscillator,
							    nf, lf, bf, spline::BasisType::kOscillator,
							    spline::OperatorType::kR, &f
							  );
							  
		std::cout << "<" << nf << "lf:" << lf << "| " << "j_" << L << "(qr) (d_r - l/r)" << " |" << ni << "li:" << li << ">" << std::endl; 
		std::cout << "Using RadialMatrixElementOfFunction: " << me1 << std::endl;
		std::cout << "Using RadialMatrixElementOfFunction: " << me2 << std::endl;

		me_tot = prefactor * (
							  prefactor_me1 * me1 
							  - prefactor_me2 * me2
							  );
		std::cout << "Total matrix Element : " << me_tot << "\n" << std::endl;
	}
	else if (ni==2) {
		double me1, me2, me3 ;

		double y = y_var(q, bi);
		double prefactor = -std::pow(8*y, -0.5) / std::sqrt(8);
		double prefactor_me1 = am::Hat2(2*(li+1)) * am::Hat2(2*(li+2)) * am::Hat2(2*(li+1));
		double prefactor_me2 = 2 * am::Hat2(2*(li+2)) * am::Hat2(2*(li+2)) * am::Hat2(2*(li+3));
		double prefactor_me3 = am::Hat2(2*(li+3)) * am::Hat2(2*(li+4)) * am::Hat2(2*(li+5));

		gsl_function f;
		struct SphericalBesselParams params = {L, q};
		f.function = &SphericalBesselFunctionEval;
		f.params = &params;

		me1 = spline::RadialMatrixElementOfFunction(
					    ni, li+1, bi, spline::BasisType::kOscillator,
					    nf, lf, bf, spline::BasisType::kOscillator,
					    spline::OperatorType::kR, &f
					  );

		me2 = spline::RadialMatrixElementOfFunction(
					    ni, li+3, bi, spline::BasisType::kOscillator,
					    nf, lf, bf, spline::BasisType::kOscillator,
					    spline::OperatorType::kR, &f
					  );

		me3 = spline::RadialMatrixElementOfFunction(
					    ni, li+5, bi, spline::BasisType::kOscillator,
					    nf, lf, bf, spline::BasisType::kOscillator,
					    spline::OperatorType::kR, &f
					  );

		std::cout << "<" << nf << "lf:" << lf << "| " << "j_" << L << "(qr) (d_r - l/r)" << " |" << ni << "li:" << li << ">" << std::endl; 
		std::cout << "Using RadialMatrixElementOfFunction: " << me1 << std::endl;
		std::cout << "Using RadialMatrixElementOfFunction: " << me2 << std::endl;
		std::cout << "Using RadialMatrixElementOfFunction: " << me3 << std::endl;	

		me_tot = prefactor * (
							  prefactor_me1 * me1 
							  - prefactor_me2 * me2
							  + prefactor_me3 * me3
							 );
		std::cout << "Total matrix Element : " << me_tot << "\n" << std::endl;
	}

	return me_tot;
}


// Test BesselMatrixElement_Plus
double Test_BesselMatrixElement_Plus(int ni, int li, int bi, int nf, int lf, int bf, int L, double q)
// Calculate the matrix element :
//   <n' l' | j_L(qr)(d_r + (l+1)/r) | n l>
{
	double me_tot;
	
	std::cout << "Test BesselMatrixElement_Minus" << std::endl;
	if (ni==0) {
		double me1, me2 ;

		double y = y_var(q, bi);
		double prefactor = std::pow(8*y, -0.5);
		double prefactor_me1 = 2 * am::Hat2(2*li); 
		double prefactor_me2 = am::Hat2(2*(li+1));

		gsl_function f;
		struct SphericalBesselParams params = {L, q};
		f.function = &SphericalBesselFunctionEval;
		f.params = &params;

		me1 = spline::RadialMatrixElementOfFunction(
			    ni, li-1, bi, spline::BasisType::kOscillator,
			    nf, lf, bf, spline::BasisType::kOscillator,
			    spline::OperatorType::kR, &f
			  );

		me2 = spline::RadialMatrixElementOfFunction(
			    ni, li+1, bi, spline::BasisType::kOscillator,
			    nf, lf, bf, spline::BasisType::kOscillator,
			    spline::OperatorType::kR, &f
			  );
					  
		std::cout << "<" << nf << "lf:" << lf << "| " << "j_" << L << "(qr)" << " |" << ni << "li:" << li << ">" << std::endl; 
		std::cout << "Using RadialMatrixElementOfFunction: " << me1 << std::endl;
		std::cout << "Using RadialMatrixElementOfFunction: " << me2 << std::endl;

		me_tot = prefactor * (
							  prefactor_me1 * me1
							  - prefactor_me2 * me2
							 );
		std::cout << "Total matrix Element : " << me_tot << "\n" << std::endl;
		
	} 
	else if (ni==1) {
		double me1, me2, me3, me4 ;

		double y = y_var(q, bi);
		double prefactor = std::pow(8*y, -0.5) / std::sqrt(2);
		double prefactor_me1 = am::Hat2(2*(li+1)) * 2*am::Hat2(2*li) ;
		double prefactor_me2 = am::Hat2(2*(li+1)) * am::Hat2(2*(li+1)) ;
		double prefactor_me3 = am::Hat2(2*(li+2)) * 2*am::Hat2(2*(li+2));
		double prefactor_me4 = am::Hat2(2*(li+2)) * am::Hat2(2*(li+3));

		gsl_function f;
		struct SphericalBesselParams params = {L, q};
		f.function = &SphericalBesselFunctionEval;
		f.params = &params;

		me1 = spline::RadialMatrixElementOfFunction(
					    ni, li-1, bi, spline::BasisType::kOscillator,
					    nf, lf, bf, spline::BasisType::kOscillator,
					    spline::OperatorType::kR, &f
					  );

		me2 = spline::RadialMatrixElementOfFunction(
					    ni, li+1, bi, spline::BasisType::kOscillator,
					    nf, lf, bf, spline::BasisType::kOscillator,
					    spline::OperatorType::kR, &f
					  );

		me3 = spline::RadialMatrixElementOfFunction(
					    ni, li+1, bi, spline::BasisType::kOscillator,
					    nf, lf, bf, spline::BasisType::kOscillator,
					    spline::OperatorType::kR, &f
					  );

		me4 = spline::RadialMatrixElementOfFunction(
					    ni, li+3, bi, spline::BasisType::kOscillator,
					    nf, lf, bf, spline::BasisType::kOscillator,
					    spline::OperatorType::kR, &f
					  );
							  
		std::cout << "<" << nf << "lf:" << lf << "| " << "j_" << L << "(qr) (d_r - l/r)" << " |" << ni << "li:" << li << ">" << std::endl; 
		std::cout << "Using RadialMatrixElementOfFunction: " << me1 << std::endl;
		std::cout << "Using RadialMatrixElementOfFunction: " << me2 << std::endl;
		std::cout << "Using RadialMatrixElementOfFunction: " << me3 << std::endl;
		std::cout << "Using RadialMatrixElementOfFunction: " << me4 << std::endl;

		me_tot = prefactor * (
			prefactor_me1 * me1
			- prefactor_me2 * me2
			- prefactor_me3 * me3
			+ prefactor_me4 * me4
		);
		std::cout << "Total matrix Element : " << me_tot << "\n" << std::endl;
	}
	else if (ni==2) {
		double me1, me2, me3, me4, me5, me6 ;

		double y = y_var(q, bi);
		double prefactor = std::pow(8*y, -0.5) / std::sqrt(8);
		double prefactor_me1 = am::Hat2(2*(li+1)) * am::Hat2(2*(li+2)) * 2*am::Hat2(2*li);
		double prefactor_me2 = am::Hat2(2*(li+1)) * am::Hat2(2*(li+2)) * am::Hat2(2*(li+1));
		double prefactor_me3 = 2*am::Hat2(2*(li+2)) * am::Hat2(2*(li+2)) * 2*am::Hat2(2*(li+2));
		double prefactor_me4 = 2*am::Hat2(2*(li+2)) * am::Hat2(2*(li+2)) * am::Hat2(2*(li+3));
		double prefactor_me5 = am::Hat2(2*(li+3)) * am::Hat2(2*(li+4)) * 2*am::Hat2(2*(li+4));
		double prefactor_me6 = am::Hat2(2*(li+3)) * am::Hat2(2*(li+4)) * am::Hat2(2*(li+5));

		gsl_function f;
		struct SphericalBesselParams params = {L, q};
		f.function = &SphericalBesselFunctionEval;
		f.params = &params;

		me1 = spline::RadialMatrixElementOfFunction(
					    ni, li-1, bi, spline::BasisType::kOscillator,
					    nf, lf, bf, spline::BasisType::kOscillator,
					    spline::OperatorType::kR, &f
					  );

		me2 = spline::RadialMatrixElementOfFunction(
					    ni, li+1, bi, spline::BasisType::kOscillator,
					    nf, lf, bf, spline::BasisType::kOscillator,
					    spline::OperatorType::kR, &f
					  );

		me3 = spline::RadialMatrixElementOfFunction(
					    ni, li+1, bi, spline::BasisType::kOscillator,
					    nf, lf, bf, spline::BasisType::kOscillator,
					    spline::OperatorType::kR, &f
					  );

		me4 = spline::RadialMatrixElementOfFunction(
					    ni, li+3, bi, spline::BasisType::kOscillator,
					    nf, lf, bf, spline::BasisType::kOscillator,
					    spline::OperatorType::kR, &f
					  );

		me5 = spline::RadialMatrixElementOfFunction(
					    ni, li+3, bi, spline::BasisType::kOscillator,
					    nf, lf, bf, spline::BasisType::kOscillator,
					    spline::OperatorType::kR, &f
					  );

		me6 = spline::RadialMatrixElementOfFunction(
					    ni, li+5, bi, spline::BasisType::kOscillator,
					    nf, lf, bf, spline::BasisType::kOscillator,
					    spline::OperatorType::kR, &f
					  );

		std::cout << "<" << nf << "lf:" << lf << "| " << "j_" << L << "(qr) (d_r - l/r)" << " |" << ni << "li:" << li << ">" << std::endl; 
		std::cout << "Using RadialMatrixElementOfFunction: " << me1 << std::endl;
		std::cout << "Using RadialMatrixElementOfFunction: " << me2 << std::endl;
		std::cout << "Using RadialMatrixElementOfFunction: " << me3 << std::endl;
		std::cout << "Using RadialMatrixElementOfFunction: " << me4 << std::endl;
		std::cout << "Using RadialMatrixElementOfFunction: " << me5 << std::endl;
		std::cout << "Using RadialMatrixElementOfFunction: " << me6 << std::endl;

		me_tot = prefactor * (
			prefactor_me1 * me1
			- prefactor_me2 * me2
			- prefactor_me3 * me3
			+ prefactor_me4 * me4
			+ prefactor_me5 * me5
			- prefactor_me6 * me6
		);
		std::cout << "Total Matrix Element : " << me_tot << "\n" << std::endl;
	}

	return me_tot ;
}


// Test MJ_MatrixElement
double Test_MJ_MatrixElement(int ni, int li, int bi, double ji, int nf, int lf, int bf, double jf, int J, double q)
// Calculate the (reduced) matrix element :
//    <n' l' j' || MJ(qr) || n l j>
{

	std::cout << "Test MJ Matrix Element" << std::endl;	
	double j6_symbol = am::Wigner6J2(2*lf, 2*jf, 2*0.5, 2*ji, 2*li, 2*J) ;
	double j3_symbol = am::Wigner3J2(2*lf, 2*J, 2*li, 0, 0, 0) ;
	double prefactor = (1/std::sqrt(4 * M_PI) * std::pow(-1, J+ji+0.5) 
					   * am::Hat2(2*lf) * am::Hat2(2*li) * am::Hat2(2*jf) * am::Hat2(2*ji) * am::Hat2(2*J)
					   * j6_symbol * j3_symbol) ;

	double MJ = prefactor * Test_BesselMatrixElement(ni, li, bi, nf, lf, bf, J, q);
	std::cout << "1/sqrt(4 pi) : " << 1/std::sqrt(4 * M_PI) << std::endl;
	std::cout << "hat(ji) : " << am::Hat2(2*ji) << std::endl;
	std::cout << "prefactor : " << 1/std::sqrt(4 * M_PI) * std::pow(-1, J+ji+0.5) * Hat(lf) * Hat(li) * Hat(jf) * Hat(ji) * Hat(J) << std::endl;
	std::cout << "(-1)^ : " << std::pow(-1, J+ji+0.5) << std::endl;
	std::cout << "Wigner6J : " << j6_symbol << std::endl;
	std::cout << "Wigner3J : " << j3_symbol << std::endl;
	std::cout << "MJ : " << MJ << "\n" << std::endl;
	
	return MJ;
}


// Test MJLSigma_MatrixElement
double Test_MJLSigma_MatrixElement(int ni, int li, int bi, double ji, int nf, int lf, int bf, double jf, int J, int L, double q)
// Calculate the (reduced) matrix element :
//    <n' l' j' || MJL(qr) \sigma || n l j>
{

	double j9_symbol = am::Wigner9J2(2*lf, 2*li, 2*L, 2*0.5, 2*0.5, 2*1, 2*jf, 2*ji, 2*J) ;
	double j3_symbol = am::Wigner3J2(2*lf, 2*L, 2*li, 0, 0, 0) ;
	double prefactor = (1/std::sqrt(4 * M_PI) * std::pow(-1, lf) * std::sqrt(6)
					   * am::Hat2(2*lf) * am::Hat2(2*li) * am::Hat2(2*jf) * am::Hat2(2*ji) * am::Hat2(2*L) * am::Hat2(2*J)
					   * j9_symbol * j3_symbol) ; 

	double MJLSigma = prefactor * Test_BesselMatrixElement(ni, li, bi, nf, lf, bf, J, q);

	return MJLSigma;
}


// Test MJLNabla_MatrixElement
double Test_MJLNabla_MatrixElement(int ni, int li, int bi, double ji, int nf, int lf, int bf, double jf, int J, int L, double q)
// Calculate the (reduced) matrix element :
//    <n' l' j' || MJL(qr) \nabla / q || n l j>
{

	double j6_symbol = am::Wigner6J2(2*lf, 2*jf, 2*0.5, 2*ji, 2*li, 2*J) ;
	double prefactor = (1/std::sqrt(4 * M_PI) * std::pow(-1, L+ji+0.5)
					   * am::Hat2(2*lf) * am::Hat2(2*jf) * am::Hat2(2*ji) * am::Hat2(2*L) * am::Hat2(2*J)
					   * j6_symbol) ;

	double j6_symbol_term1 = am::Wigner6J2(2*L, 2*1, 2*J, 2*li, 2*lf, 2*(li+1)) ;
	double j3_symbol_term1 = am::Wigner3J2(2*lf, 2*L, 2*(li+1), 0, 0, 0) ;
	double prefactor_term1 = std::sqrt(li + 1) * am::Hat2(2*(li+1)) * j6_symbol_term1 * j3_symbol_term1 ; 
	double term1 = Test_BesselMatrixElement_Minus(ni, li, bi, nf, lf, bf, J, q);

	// add test if li=0 --> issue 
	double j6_symbol_term2 = am::Wigner6J2(2*L, 2*1, 2*J, 2*li, 2*lf, 2*(li-1)) ;
	double j3_symbol_term2 = am::Wigner3J2(2*lf, 2*L, 2*(li-1), 0, 0, 0) ;
	double prefactor_term2 = std::sqrt(li) * am::Hat2(2*(li-1)) * j6_symbol_term2 * j3_symbol_term2 ; 
	double term2 = Test_BesselMatrixElement_Plus(ni, li, bi, nf, lf, bf, J, q);

	double MJLNabla = prefactor * (
								   - prefactor_term1 * term1
								   + prefactor_term2 * term2
								   );

	return MJLNabla;

}


// Test MJSigmaNabla_MatrixElement
double Test_MJSigmaNabla_MatrixElement(int ni, int li, int bi, double ji, int nf, int lf, int bf, double jf, int J, double q)
// Calculate the (reduced) matrix element :
//    <n' l' j' || MJ(qr) \sigma \nabla / q|| n l j>
{

	double j6_symbol = am::Wigner6J2(2*lf, 2*jf, 2*0.5, 2*ji, 2*(2*ji-li), 2*J) ;
	double j3_symbol = am::Wigner3J2(2*lf, 2*J, 2*(2*ji-1), 0, 0, 0) ;
	double prefactor = (1/std::sqrt(4 * M_PI) * std::pow(-1, lf)
					   * am::Hat2(2*lf) * am::Hat2(2*jf) * am::Hat2(2*ji) * am::Hat2(2*(2*ji-1)) * am::Hat2(2*J)
					   * j6_symbol * j3_symbol) ;

	double term1 = 0.0 ;
	double term2 = 0.0 ;
	
	if (ji==(li+0.5)) {
		term1 = Test_BesselMatrixElement_Minus(ni, li, bi, nf, lf, bf, J, q);
	}

	if (ji==(li-0.5)) {
		term2 = Test_BesselMatrixElement_Plus(ni, li, bi, nf, lf, bf, J, q);
	}

	double MJSigmaNabla = prefactor * (-term1 + term2);

	return MJSigmaNabla;
	
}


////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{


  Test_BesselMatrixElement(0, 0, 1.0, 0, 0, 1.0, 0, 1.0);
  Test_BesselMatrixElement(0, 0, 1.0, 0, 0, 1.0, 1, 1.0);
  Test_BesselMatrixElement(0, 0, 1.0, 0, 0, 1.0, 2, 1.0);

  Test_BesselMatrixElement_Minus(0, 0, 1.0, 0, 0, 1.0, 0, 1.0);
  Test_BesselMatrixElement_Minus(0, 0, 1.0, 0, 0, 1.0, 1, 1.0);
  Test_BesselMatrixElement_Minus(0, 0, 1.0, 0, 0, 1.0, 2, 1.0);

  Test_BesselMatrixElement_Minus(1, 0, 1.0, 0, 0, 1.0, 0, 1.0);

  Test_MJ_MatrixElement(0, 0, 1.0, 0.5, 0, 0, 1.0, 0.5, 0, 1.0);

  // termination
  return 0;
}
