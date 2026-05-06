/****************************************************************
  seven_operators_test.cpp

  Victor Duménil
  University of Notre Dame & LPC Caen



  [1] : Multipole operators in semileptonic weak and electromagnetic interactions 
        with nuclei: Harmonic oscillator single-particle matrix elements, 
        TW Donnelly, WC Haxton - Atomic Data and Nuclear Data Tables, 1979
  [2] : SevenOperators, a Mathematica script for harmonic oscillator nuclear 
        matrix elements arising in semileptonic electroweak interactions
        W Haxton, C Lunardini - Computer Physics Communications, 2008
  [3] : Semileptonic weak and electromagnetic interactions with nuclei: Nuclear 
        current operators through order (v/c) nucleon2
        BD Serot - Nuclear Physics A, 1978


  TO DO : 
    - Add the other 'seven operators' in ref [3]
        
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
//   q (double): impulsion q (units inverse of the length parameter b)
{
  int order;
  double q;
};


double SphericalBesselFunctionEval(double r, void * p)
// Arguments :
//   order : 
//   q : 
//   r : radius
// Return :
//   spherical Bessel function : j_order(q*r)
{
    SphericalBesselParams* params = static_cast< SphericalBesselParams*>(p);
    int order = (params->order);
    double q = (params->q);

    return gsl_sf_bessel_jl(order, q*r);
 }


// y variable
double y_var(double q, double b)
// Arguments :
//   q : impulsion
//   b : harmonic oscillator parameter
// Return :
//   (q*b/2)^2
{
	return std::pow(q * b / 2, 2);
}


// Nodal quantum number N 
// Inverse notation with ref [1] & [2] !!
// Be careful with these numbers (n & N), 
// the principal quantum number (N) in ref [1] & [2], doesn't seem to be the 'principal' quantum number
// here n is the true 'principal' quantum number (n=0, 1, 2 by convention, different convention use in ref [1] & [2]) 
// where n=0 correponds to --> 0s, 0p, ... 
int NodalQuantumNumber(int n, int l)
// Arguments :
//   n : principal quantum number (>=0)
//   l : angular momentum
// Return :
//   N = (n - l)/2 /1 > 0
{
	return (n - l)/2 + 1;
}


// Parity & Physical Condition
// Normal parity : (-1)^li x (-1)^lf x (-1)^J == 1
// Abnormal parity : (-1)^li x (-1)^lf x (-1)^(J+1) == 1

// NormalPhysicalCondition
bool NormalPhysicalCondition(int ni, int li, double ji, int nf, int lf, double jf, int J)
// Arguments :
//   ni, nf : initial/final principal quantum number 
//   li, lf : initial/final angular momentum 
//   ji, jf : initial/final total angular momentum
//   J : rank
// Return :
//   true / false if the normal physical condition is respected
{
	bool parity_ok = std::pow(-1, li + lf + J) == 1;

	int Ni = NodalQuantumNumber(ni, li);
	int Nf = NodalQuantumNumber(nf, lf);

	return (parity_ok && (std::abs(ji-jf) <= J) && (J <= (ji+jf))); //&& (Ni>0) && (Nf>0) 
}


// AbnormalPhysicalCondition
bool AbnormalPhysicalCondition(int ni, int li, double ji, int nf, int lf, double jf, int J)
// Arguments :
//   ni, nf : initial/final principal quantum number 
//   li, lf : initial/final angular momentum 
//   ji, jf : initial/final total angular momentum
//   J : rank
// Return :
//   true / false if the abnormal physical condition is respected
{
	bool parity_ok = std::pow(-1, li + lf + J + 1) == 1; 

	int Ni = NodalQuantumNumber(ni, li);
	int Nf = NodalQuantumNumber(nf, lf);

	return (parity_ok && (std::abs(ji-jf) <= J) && (J <= (ji+jf))); //&& (Ni>0) && (Nf>0) 
}


// Test BesselMatrixElement
double Test_BesselMatrixElement(int ni, int li, int bi, int nf, int lf, int bf, int L, double q) 
// Arguments :
//   ni, nf : initial/final principal quantum number 
//   li, lf : initial/final angular momentum 
//   bi, bf : initial/final harmonic oscillator parameter
//   L : rank
//   q : impulsion
// Return the matrix element :
//   <n' l' | j_L(qr) | n l>
{
	std::cout << "Test BesselMatrixElement" << std::endl;
	double me_tot ;
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


// We give the relation to calculate the 3 'Bessel' matrix elements which appear in eq.3 in ref [1]
// <n' l' j' | j_L(rho) | n l j> ; <n' l' j' | j_L(rho)(d_rho - l/rho) | n l j> ; <n' l' j' | j_L(rho)(d_rho + (l+1)/rho) | n l j> 
// Test BesselMatrixElement_Minus
double Test_BesselMatrixElement_Minus(int ni, int li, int bi, int nf, int lf, int bf, int L, double q)
// Return the matrix element :
//   <n' l' | j_L(qr)(d_r - l/r) | n l>
{
	std::cout << "Test BesselMatrixElement_Minus" << std::endl;
	double me_tot ;
	int Ni = NodalQuantumNumber(ni, li); // the reccurence relations are deduced in function of the Nodal quantum number
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
// Return the matrix element :
//   <n' l' | j_L(qr)(d_r + (l+1)/r) | n l>
{
	std::cout << "Test BesselMatrixElement_Minus" << std::endl;
	double me_tot ;
	int Ni = NodalQuantumNumber(ni, li); // the reccurence relations are deduced in function of the Nodal quantum number
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


// We give the 4 reduced matrix elements (eq.3 in ref [1])
// <n' l' j' || MJ(qr) || n l j> ; <n' l' j' || MJL(qr) sigma || n l j>
// <n' l' j' || MJL(qr) nabla/q || n l j> ; <n' l' j' || MJ(qr) sigma nabla/q || n l j>
// Test MJ_MatrixElement
double Test_MJ_MatrixElement(int ni, int li, int bi, double ji, int nf, int lf, int bf, double jf, int J, double q)
// Calculate the (reduced) matrix element :
//    <n' l' j' || MJ(qr) || n l j>
{
	std::cout << "Test MJ Matrix Element" << std::endl;	
	double MJ = 0.0 ;
	if (J>=0) { 
		double j6_symbol = am::Wigner6J2(2*lf, 2*jf, 2*0.5, 2*ji, 2*li, 2*J) ;
		double j3_symbol = am::Wigner3J2(2*lf, 2*J, 2*li, 0, 0, 0) ;
		double prefactor = (1/std::sqrt(4 * M_PI) * std::pow(-1, J+ji+0.5) 
						   * am::Hat2(2*lf) * am::Hat2(2*li) * am::Hat2(2*jf) * am::Hat2(2*ji) * am::Hat2(2*J)
						   * j6_symbol * j3_symbol) ;

		// MJ matrix element calculation
		MJ = prefactor * Test_BesselMatrixElement(ni, li, bi, nf, lf, bf, J, q);
		/*std::cout << "1/sqrt(4 pi) : " << 1/std::sqrt(4 * M_PI) << std::endl;
		std::cout << "hat(ji) : " << am::Hat2(2*ji) << std::endl;
		std::cout << "prefactor : " << 1/std::sqrt(4 * M_PI) * std::pow(-1, J+ji+0.5) * Hat(lf) * Hat(li) * Hat(jf) * Hat(ji) * Hat(J) << std::endl;
		std::cout << "(-1)^ : " << std::pow(-1, J+ji+0.5) << std::endl;
		std::cout << "Wigner6J : " << j6_symbol << std::endl;
		std::cout << "Wigner3J : " << j3_symbol << std::endl;*/
		std::cout << "MJ : " << MJ << "\n" << std::endl;
	}
	
	return MJ;
}


// Test MJLSigma_MatrixElement
double Test_MJLSigma_MatrixElement(int ni, int li, int bi, double ji, int nf, int lf, int bf, double jf, int J, int L, double q)
// Calculate the (reduced) matrix element :
//    <n' l' j' || MJL(qr) \sigma || n l j>
{
	std::cout << "MJLSigma Matrix Element" << std::endl;
	double MJLSigma = 0.0 ;
	if (J>=0 && L>=0) {
		double j9_symbol = am::Wigner9J2(2*lf, 2*li, 2*L, 2*0.5, 2*0.5, 2*1, 2*jf, 2*ji, 2*J) ;
		double j3_symbol = am::Wigner3J2(2*lf, 2*L, 2*li, 0, 0, 0) ;
		double prefactor = (1/std::sqrt(4 * M_PI) * std::pow(-1, lf) * std::sqrt(6)
						   * am::Hat2(2*lf) * am::Hat2(2*li) * am::Hat2(2*jf) * am::Hat2(2*ji) * am::Hat2(2*L) * am::Hat2(2*J)
						   * j9_symbol * j3_symbol) ; 
		
		// MJLSigma matrix element Calculation
		MJLSigma = prefactor * Test_BesselMatrixElement(ni, li, bi, nf, lf, bf, L, q);
	}
	return MJLSigma;
}


// Test MJLNabla_MatrixElement
double Test_MJLNabla_MatrixElement(int ni, int li, int bi, double ji, int nf, int lf, int bf, double jf, int J, int L, double q)
// Calculate the (reduced) matrix element :
//    <n' l' j' || MJL(qr) \nabla / q || n l j>
{
	std::cout << "MJLNabla Matrix Element" << std::endl;
	double MJLNabla = 0.0 ;
	if (J>=0 && L>=0) {
		double j6_symbol = am::Wigner6J2(2*lf, 2*jf, 2*0.5, 2*ji, 2*li, 2*J) ;
		double prefactor = (1/std::sqrt(4 * M_PI) * std::pow(-1, L+ji+0.5)
						   * am::Hat2(2*lf) * am::Hat2(2*jf) * am::Hat2(2*ji) * am::Hat2(2*L) * am::Hat2(2*J)
						   * j6_symbol) ;

		// Calculation of 1st term
		double j6_symbol_term1 = am::Wigner6J2(2*L, 2*1, 2*J, 2*li, 2*lf, 2*(li+1)) ;
		double j3_symbol_term1 = am::Wigner3J2(2*lf, 2*L, 2*(li+1), 0, 0, 0) ;
		double prefactor_term1 = std::sqrt(li + 1) * am::Hat2(2*(li+1)) * j6_symbol_term1 * j3_symbol_term1 ; 
		double term1 = Test_BesselMatrixElement_Minus(ni, li, bi, nf, lf, bf, L, q);

		// Calculation of 2nd term
		// li > 0 
		double prefactor_term2 = 0.0 ;
		double term2 = 0.0 ;
		if (li>0) {		 
			double j6_symbol_term2 = am::Wigner6J2(2*L, 2*1, 2*J, 2*li, 2*lf, 2*(li-1)) ;
			double j3_symbol_term2 = am::Wigner3J2(2*lf, 2*L, 2*(li-1), 0, 0, 0) ;
			prefactor_term2 = std::sqrt(li) * am::Hat2(2*(li-1)) * j6_symbol_term2 * j3_symbol_term2 ; 
			term2 = Test_BesselMatrixElement_Plus(ni, li, bi, nf, lf, bf, L, q);
		}

		// MJLNabla matrix element calculation
		MJLNabla = prefactor * (
							   - prefactor_term1 * term1
							   + prefactor_term2 * term2
							   );
	}
	return MJLNabla;

}


// Test MJSigmaNabla_MatrixElement
double Test_MJSigmaNabla_MatrixElement(int ni, int li, int bi, double ji, int nf, int lf, int bf, double jf, int J, double q)
// Calculate the (reduced) matrix element :
//    <n' l' j' || MJ(qr) \sigma \nabla / q|| n l j>
{
	std::cout << "MJSigmaNabla Matrix Element" << std::endl;
	double MJSigmaNabla = 0.0 ;
	if (J>=0) {
		double j6_symbol = am::Wigner6J2(2*lf, 2*jf, 2*0.5, 2*ji, 2*(2*ji-li), 2*J) ;
		double j3_symbol = am::Wigner3J2(2*lf, 2*J, 2*(2*ji-li), 0, 0, 0) ;
		double prefactor = (1/std::sqrt(4 * M_PI) * std::pow(-1, lf)
						   * am::Hat2(2*lf) * am::Hat2(2*jf) * am::Hat2(2*ji) * am::Hat2(2*(2*ji-li)) * am::Hat2(2*J)
						   * j6_symbol * j3_symbol) ;

		double term1 = 0.0 ;
		double term2 = 0.0 ;

		// Calculation of 1st term
		if (ji==(li+0.5)) {
			term1 = Test_BesselMatrixElement_Minus(ni, li, bi, nf, lf, bf, J, q);
		}
		
		// Calculation of 2nd term
		if (ji==(li-0.5)) {
			term2 = Test_BesselMatrixElement_Plus(ni, li, bi, nf, lf, bf, J, q);
		}

		// MJSigmaNabla matrix element calculation
		MJSigmaNabla = prefactor * (-term1 + term2);
	}
	return MJSigmaNabla;
	
}


// Seven basis single-particle operators
// Here we calculate <n' l' j' || \hat{O}_J(qr) || n l j>
// where the operator \hat{O}_J(qr) corresponds to eq.(1) in ref [1] :
// MUST satisfy the Normal parity : M_J(qr) ; Δ'_J(qr) ; Σ_J(qr)
// MUST satisfy the Abnormal parity : Δ_J(qr) ; Σ'_J(qr) ; Σ''_J(qr) ; Ω_J(qr) ; Ω'_J(qr)

// Test MJ_SevenOprator
double Test_MJ_SevenOperator(int ni, int li, int bi, double ji, int nf, int lf, int bf, double jf, int J, double q)
// Calculate the (reduced) matrix element : 
//    <n' l' j' || M_J(qr) || n l j>
{
    std::cout << "MJ Seven Operator" << std::endl ;
    double MJ = 0.0;
	if (NormalPhysicalCondition(ni, li, ji, nf, lf, jf, J)) {
		MJ = Test_MJ_MatrixElement(ni, li, bi, ji, nf, lf, bf, jf, J, q) ;
	}
	std::cout << "MJ : " << MJ << "\n" << std::endl;

	return MJ ;
}


// Test DeltaJ_SevenOperator
double Test_DeltaJ_SevenOperator(int ni, int li, int bi, double ji, int nf, int lf, int bf, double jf, int J, double q)
// Calculate the (reduced) matrix element : 
//    <n' l' j' || Δ_J(qr) || n l j>
{
	std::cout << "DeltaJ Seven Operator" << std::endl ;
	double DeltaJ = 0.0 ;
	if (AbnormalPhysicalCondition(ni, li, ji, nf, lf, jf, J)) {
		DeltaJ = Test_MJLNabla_MatrixElement(ni, li, bi, ji, nf, lf, bf, jf, J, J, q) ;
	}
	std::cout << "DeltaJ : " << DeltaJ << "\n" << std::endl;
	
	return DeltaJ ;
}


// Test DeltaJP_SevenOperator
double Test_DeltaJP_SevenOperator(int ni, int li, int bi, double ji, int nf, int lf, int bf, double jf, int J, double q)
// Calculate the (reduced) matrix element : 
//    <n' l' j' || Δ'_J(qr) || n l j>
{
	std::cout << "DeltaJP Seven Operator" << std::endl ;
	double DeltaJP = 0.0 ;
	if (NormalPhysicalCondition(ni, li, ji, nf, lf, jf, J)) {
		DeltaJP = 1/am::Hat2(2*J) * (- std::sqrt(J) * Test_MJLNabla_MatrixElement(ni, li, bi, ji, nf, lf, bf, jf, J, J+1, q) 
									 + std::sqrt(J + 1) * Test_MJLNabla_MatrixElement(ni, li, bi, ji, nf, lf, bf, jf, J, J-1, q));
	}
	std::cout << "DeltaJP : " << DeltaJP << "\n" << std::endl;
	
	return DeltaJP ;
}


// Test SigmaJ_SevenOperator
double Test_SigmaJ_SevenOperator(int ni, int li, int bi, double ji, int nf, int lf, int bf, double jf, int J, double q)
// Calculate the (reduced) matrix element : 
//    <n' l' j' || Σ_J(qr) || n l j>
{
	std::cout << "SigmaJ Seven Operator" << std::endl ;
	double SigmaJ = 0.0 ;
	if (NormalPhysicalCondition(ni, li, ji, nf, lf, jf, J)) {
		SigmaJ = Test_MJLSigma_MatrixElement(ni, li, bi, ji, nf, lf, bf, jf, J, J, q) ;
	}
	std::cout << "SigmaJ : " << SigmaJ << "\n" << std::endl;
	
	return SigmaJ ;
}


// Test SigmaJP_SevenOperator
double Test_SigmaJP_SevenOperator(int ni, int li, int bi, double ji, int nf, int lf, int bf, double jf, int J, double q)
// Calculate the (reduced) matrix element : 
//    <n' l' j' || Σ'_J(qr) || n l j>
{
	std::cout << "SigmaJP Seven Operator" << std::endl ;
	double SigmaJP = 0.0 ;
	if (AbnormalPhysicalCondition(ni, li, ji, nf, lf, jf, J)) {
		SigmaJP = 1/am::Hat2(2*J) * (- std::sqrt(J) * Test_MJLSigma_MatrixElement(ni, li, bi, ji, nf, lf, bf, jf, J, J+1, q) 
										+ std::sqrt(J + 1) * Test_MJLSigma_MatrixElement(ni, li, bi, ji, nf, lf, bf, jf, J, J-1, q));
	}
	std::cout << "SigmaJP : " << SigmaJP << "\n" << std::endl;
	
	return SigmaJP ;
}


// Test SigmaJPP_SevenOperator
double Test_SigmaJPP_SevenOperator(int ni, int li, int bi, double ji, int nf, int lf, int bf, double jf, int J, double q)
// Calculate the (reduced) matrix element : 
//    <n' l' j' || Σ''_J(qr) || n l j>
{
	std::cout << "SigmaJPP Seven Operator" << std::endl ;
	double SigmaJPP = 0.0 ;
	if (AbnormalPhysicalCondition(ni, li, ji, nf, lf, jf, J)) {
		SigmaJPP = 1/am::Hat2(2*J) * (std::sqrt(J + 1) * Test_MJLSigma_MatrixElement(ni, li, bi, ji, nf, lf, bf, jf, J, J+1, q) 
										+ std::sqrt(J) * Test_MJLSigma_MatrixElement(ni, li, bi, ji, nf, lf, bf, jf, J, J-1, q));
	}
	std::cout << "SigmaJPP : " << SigmaJPP << "\n" << std::endl;
	
	return SigmaJPP ;
}


// Test OmegaJ_SevenOperator
double Test_OmegaJ_SevenOperator(int ni, int li, int bi, double ji, int nf, int lf, int bf, double jf, int J, double q)
// Calculate the (reduced) matrix element : 
//    <n' l' j' || Ω_J(qr) || n l j>
{
	std::cout << "OmegaJ Seven Operator" << std::endl ;
	double OmegaJ = 0.0 ;
	if (AbnormalPhysicalCondition(ni, li, ji, nf, lf, jf, J)) {
		OmegaJ = Test_MJSigmaNabla_MatrixElement(ni, li, bi, ji, nf, lf, bf, jf, J, q) ;
	}
	std::cout << "OmegaJ : " << OmegaJ << "\n" << std::endl;
	
	return OmegaJ ;
}


// Test OmegaJP_SevenOperator
double Test_OmegaJP_SevenOperator(int ni, int li, int bi, double ji, int nf, int lf, int bf, double jf, int J, double q)
// Calculate the (reduced) matrix element : 
//    <n' l' j' || Ω'_J(qr) || n l j>
{
	std::cout << "OmegaJP Seven Operator" << std::endl ;
	double OmegaJP = 0.0 ;
	if (AbnormalPhysicalCondition(ni, li, ji, nf, lf, jf, J)) {
		OmegaJP = (Test_OmegaJ_SevenOperator(ni, li, bi, ji, nf, lf, bf, jf, J, q) 
					  + 0.5 * Test_SigmaJPP_SevenOperator(ni, li, bi, ji, nf, lf, bf, jf, J, q)) ;
	}
	std::cout << "OmegaJP : " << OmegaJP << "\n" << std::endl;
	
	return OmegaJP ;
}


////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  // Bessel Matrix Element
  //Test_BesselMatrixElement(0, 0, 1.0, 0, 0, 1.0, 0, 1.0);
  //Test_BesselMatrixElement(0, 0, 1.0, 0, 0, 1.0, 1, 1.0);
  //Test_BesselMatrixElement(0, 0, 1.0, 0, 0, 1.0, 2, 1.0);

  // Bessel Matrix Element Minus
  //Test_BesselMatrixElement_Minus(0, 0, 1.0, 0, 0, 1.0, 0, 1.0);
  //Test_BesselMatrixElement_Minus(0, 0, 1.0, 0, 0, 1.0, 1, 1.0);
  //Test_BesselMatrixElement_Minus(0, 0, 1.0, 0, 0, 1.0, 2, 1.0);
  
  //Test_BesselMatrixElement_Minus(1, 0, 1.0, 0, 0, 1.0, 0, 1.0);

  // MJ Matrix Element
  //Test_MJ_MatrixElement(0, 0, 1.0, 0.5, 0, 0, 1.0, 0.5, 0, 1.0);

  // MJ Seven Operator
  Test_MJ_SevenOperator(0, 0, 1.0, 0.5, 0, 0, 1.0, 0.5, 0, 1.0);
  Test_MJ_SevenOperator(0, 0, 1.0, 0.5, 0, 1, 1.0, 0.5, 1, 1.0);

  // DeltaJ Seven Operator
  Test_DeltaJ_SevenOperator(0, 1, 1.0, 0.5, 0, 1, 1.0, 0.5, 1, 1.0);
  Test_DeltaJ_SevenOperator(0, 1, 1.0, 0.5, 0, 2, 1.0, 1.5, 2, 1.0);

  // DeltaJP Seven Operator
  Test_DeltaJP_SevenOperator(0, 0, 1.0, 0.5, 0, 2, 1.0, 1.5, 2, 1.0);
  Test_DeltaJP_SevenOperator(0, 0, 1.0, 0.5, 0, 1, 1.0, 0.5, 1, 1.0);

  // SigmaJ Seven Operator
  Test_SigmaJ_SevenOperator(0, 0, 1.0, 0.5, 0, 2, 1.0, 1.5, 2, 1.0);
  Test_SigmaJ_SevenOperator(0, 0, 1.0, 0.5, 0, 1, 1.0, 0.5, 1, 1.0);

  // SigmaJP Seven Operator
  Test_SigmaJP_SevenOperator(0, 1, 1.0, 0.5, 0, 1, 1.0, 0.5, 1, 1.0);
  Test_SigmaJP_SevenOperator(0, 1, 1.0, 0.5, 0, 2, 1.0, 1.5, 2, 1.0);

  // SigmaJPP Seven Operator
  Test_SigmaJPP_SevenOperator(0, 1, 1.0, 0.5, 0, 1, 1.0, 0.5, 1, 1.0);
  Test_SigmaJPP_SevenOperator(0, 1, 1.0, 0.5, 0, 2, 1.0, 1.5, 2, 1.0);

  // OmegaJ Seven Operator
  Test_OmegaJ_SevenOperator(0, 1, 1.0, 0.5, 0, 1, 1.0, 0.5, 1, 1.0);
  Test_OmegaJ_SevenOperator(0, 1, 1.0, 0.5, 0, 2, 1.0, 1.5, 2, 1.0);

  // OmegaJP Seven Operator
  Test_OmegaJP_SevenOperator(0, 1, 1.0, 0.5, 0, 1, 1.0, 0.5, 1, 1.0);
  Test_OmegaJP_SevenOperator(0, 1, 1.0, 0.5, 0, 2, 1.0, 1.5, 2, 1.0);
  
  // termination
  return 0;
}
