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
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>


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
	double jl ;
	if (q < 0) {
		if (order % 2 == 0) {
			jl = gsl_sf_bessel_jl(order, std::abs(q)*r);
		}
		else {
			jl = -gsl_sf_bessel_jl(order, std::abs(q)*r);
		}
	}  
	else {
		jl = gsl_sf_bessel_jl(order, q*r);
	}

    return jl;
};


// y variable
double y_var(double q, double b)
// Arguments :
//   q : impulsion
//   b : harmonic oscillator parameter
// Return :
//   (q*b/2)^2
{
	return std::pow(q * b / 2.0, 2);
};


// Nodal quantum number N 
// NOT use here
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
};


// Parity & Physical Condition
// Normal parity : (-1)^li x (-1)^lf x (-1)^J == 1
// Abnormal parity : (-1)^li x (-1)^lf x (-1)^(J+1) == 1

bool NormalParity(int li, int lf, int J)
//
{
	return std::pow(-1, li + lf + J) == 1;
};


bool AbnormalParity(int li, int lf, int J)
//
{
	return std::pow(-1, li + lf + J + 1) == 1;
};


bool TriangularCondition(double ji, double jf, int J)
//
{
	return (std::abs(ji-jf) <= J) && (J <= (ji+jf)) ;
};


// NormalPhysicalCondition
bool NormalPhysicalCondition(int li, double ji, int lf, double jf, int J)
// Arguments :
//   li, lf : initial/final angular momentum 
//   ji, jf : initial/final total angular momentum
//   J : rank
// Return :
//   true / false if the normal physical condition is respected
{
	bool parity_ok = NormalParity(li, lf, J);

	return (parity_ok && TriangularCondition(ji, jf, J));
};


// AbnormalPhysicalCondition
bool AbnormalPhysicalCondition(int li, double ji, int lf, double jf, int J)
// Arguments :
//   li, lf : initial/final angular momentum 
//   ji, jf : initial/final total angular momentum
//   J : rank
// Return :
//   true / false if the abnormal physical condition is respected
{
	bool parity_ok = AbnormalParity(li, lf, J);

	return (parity_ok && TriangularCondition(ji, jf, J));
};


// Translationally invariant term for one-body matrix elements
// Translationally invariant matrix elements of general one-body operators, Petr Navrátil, 2021

double TranslationallyInvariantTerm(int A)
//
{
	/*if (A <= 1) {
		return 1.0;
	}
	else {
		return -std::sqrt((A-1.0)/A);
	}*/
	return 1.0;
}

// We give the relation to calculate the 3 'Bessel' matrix elements which appear in eq.3 in ref [1]
// <n' l' j' | j_L(rho) | n l j> ; <n' l' j' | j_L(rho)(d_rho - l/rho) | n l j> ; <n' l' j' | j_L(rho)(d_rho + (l+1)/rho) | n l j> 
// BesselMatrixElement
/*double BesselMatrixElement(int ni, int li, int bi, int nf, int lf, int bf, int L, double q) 
// Arguments :
//   ni, nf : initial/final principal quantum number 
//   li, lf : initial/final angular momentum 
//   bi, bf : initial/final harmonic oscillator parameter
//   L : rank
//   q : impulsion
// Return the matrix element :
//   <n' l' | j_L(qr) | n l>
{
	double me_tot = 0.0;
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

	me_tot = me1 ;

	return me_tot;
};


// BesselMatrixElement_Minus
double BesselMatrixElement_Minus(int ni, int li, int bi, int nf, int lf, int bf, int L, double q)
// Return the matrix element :
//   <n' l' | j_L(qr)(d_r - l/r) | n l>
{
	double me_tot = 0.0;
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

		me_tot = prefactor * prefactor_me1 * me1;
		
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

		me_tot = prefactor * (
							  prefactor_me1 * me1 
							  - prefactor_me2 * me2
							  );
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

		me_tot = prefactor * (
							  prefactor_me1 * me1 
							  - prefactor_me2 * me2
							  + prefactor_me3 * me3
							 );
	}

	return me_tot;
};


// BesselMatrixElement_Plus
double BesselMatrixElement_Plus(int ni, int li, int bi, int nf, int lf, int bf, int L, double q)
// Return the matrix element :
//   <n' l' | j_L(qr)(d_r + (l+1)/r) | n l>
{
	double me_tot = 0.0;
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

		me_tot = prefactor * (
							  prefactor_me1 * me1
							  - prefactor_me2 * me2
							 );
		
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

		me_tot = prefactor * (
			prefactor_me1 * me1
			- prefactor_me2 * me2
			- prefactor_me3 * me3
			+ prefactor_me4 * me4
		);
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

		me_tot = prefactor * (
			prefactor_me1 * me1
			- prefactor_me2 * me2
			- prefactor_me3 * me3
			+ prefactor_me4 * me4
			+ prefactor_me5 * me5
			- prefactor_me6 * me6
		);
	}

	return me_tot ;
};*/


double factorial2(int n) {
    if (n <= 0) return 1.0;
    double result = 1.0;
    for (int k = 1; k <= n; ++k) {
        result *= (2 * k + 1);
    }
    return result;
}

// We give the relation to calculate the 3 'Bessel' matrix elements which appear in eq.3 in ref [1]
// <n' l' j' | j_L(rho) | n l j> ; <n' l' j' | j_L(rho)(d_rho - l/rho) | n l j> ; <n' l' j' | j_L(rho)(d_rho + (l+1)/rho) | n l j> 
// --- Implementation of Basic Functions ---
double BF1(double y, int ni, int li, int nf, int lf, int L) {
    return (std::pow(2.0, L) / factorial2(L) *
            std::pow(y, L / 2.0) * std::exp(-y) *
            std::sqrt(gsl_sf_fact(ni - 1) * gsl_sf_fact(nf - 1)));
}

double BF2(double y, int ni, int li, int nf, int lf, int L) {
    return std::sqrt(gsl_sf_gamma(nf + lf + 0.5) * gsl_sf_gamma(ni + li + 0.5));
}

double S1(int ni, int li, int mi, int nf, int lf, int mf) {
    return (std::pow(-1.0, mi + mf) /
            (gsl_sf_fact(mi) * gsl_sf_fact(mf) *
             gsl_sf_fact(ni - 1 - mi) * gsl_sf_fact(nf - 1 - mf)));
}

double S2(int ni, int li, int mi, int nf, int lf, int mf, int L) {
    double numerator = gsl_sf_gamma((li + lf + L + 2 * mi + 2 * mf + 3) / 2.0);
    double denominator = gsl_sf_gamma(li + mi + 1.5) * gsl_sf_gamma(lf + mf + 1.5);
    return numerator / denominator;
}

double S3(double y, int ni, int li, int mi, int nf, int lf, int mf, int L) {
    double a = (L - li - lf - 2 * mi - 2 * mf) / 2.0;
    double b = L + 1.5;
    return gsl_sf_hyperg_1F1(a, b, y);
}

double BF3(double y, int ni, int li, int nf, int lf, int L) {
    double total = 0.0;
    for (int mi = 0; mi < ni; ++mi) {
        for (int mf = 0; mf < nf; ++mf) {
            total += S1(ni, li, mi, nf, lf, mf) *
                     S2(ni, li, mi, nf, lf, mf, L) *
                     S3(y, ni, li, mi, nf, lf, mf, L);
        }
    }
    return total;
}

// --- BesselElement ---
double BesselElement(double y, int ni, int li, int nf, int lf, int L) {
    return BF1(y, ni, li, nf, lf, L) * BF2(y, ni, li, nf, lf, L) * BF3(y, ni, li, nf, lf, L);
}

// --- Gradiant Bessel Elements (Minus) ---
double BF1A(double y, int ni, int li, int nf, int lf, int L) {
    return (std::pow(2.0, L - 1) / factorial2(L) *
            std::pow(y, (L - 1) / 2.0) * std::exp(-y) *
            std::sqrt(gsl_sf_fact(ni - 1) * gsl_sf_fact(nf - 1)));
}

double S2A(int ni, int li, int mi, int nf, int lf, int mf, int L) {
    double numerator = gsl_sf_gamma((L + li + lf + 2 * mi + 2 * mf + 2) / 2.0);
    double denominator = gsl_sf_gamma(li + mi + 1.5) * gsl_sf_gamma(lf + mf + 1.5);
    return numerator / denominator;
}

double S3A(double y, int ni, int li, int mi, int nf, int lf, int mf, int L) {
    double a1 = (L - li - lf - 2 * mi - 2 * mf - 1) / 2.0;
    double a2 = (L - li - lf - 2 * mi - 2 * mf + 1) / 2.0;
    double b = L + 1.5;
    return (-(li + lf + L + 2 * mi + 2 * mf + 2) / 2.0 * gsl_sf_hyperg_1F1(a1, b, y) +
            2 * mi * gsl_sf_hyperg_1F1(a2, b, y));
}

double BF3A(double y, int ni, int li, int nf, int lf, int L) {
    double total = 0.0;
    for (int mi = 0; mi < ni; ++mi) {
        for (int mf = 0; mf < nf; ++mf) {
            total += S1(ni, li, mi, nf, lf, mf) *
                     S2A(ni, li, mi, nf, lf, mf, L) *
                     S3A(y, ni, li, mi, nf, lf, mf, L);
        }
    }
    return total;
}

double BesselElementMinus(double y, int ni, int li, int nf, int lf, int L) {
    return BF1A(y, ni, li, nf, lf, L) * BF2(y, ni, li, nf, lf, L) * BF3A(y, ni, li, nf, lf, L);
}

// --- Gradiant Bessel Elements (Plus) ---
double S4A(double y, int ni, int li, int mi, int nf, int lf, int mf, int L) {
    double a1 = (L - li - lf - 2 * mi - 2 * mf - 1) / 2.0;
    double a2 = (L - li - lf - 2 * mi - 2 * mf + 1) / 2.0;
    double b = L + 1.5;
    return (-(li + lf + L + 2 * mi + 2 * mf + 2) / 2.0 * gsl_sf_hyperg_1F1(a1, b, y) +
            (2 * li + 2 * mi + 1) * gsl_sf_hyperg_1F1(a2, b, y));
}

double BF4A(double y, int ni, int li, int nf, int lf, int L) {
    double total = 0.0;
    for (int mi = 0; mi < ni; ++mi) {
        for (int mf = 0; mf < nf; ++mf) {
            total += S1(ni, li, mi, nf, lf, mf) *
                     S2A(ni, li, mi, nf, lf, mf, L) *
                     S4A(y, ni, li, mi, nf, lf, mf, L);
        }
    }
    return total;
}

double BesselElementPlus(double y, int ni, int li, int nf, int lf, int L) {
    return BF1A(y, ni, li, nf, lf, L) * BF2(y, ni, li, nf, lf, L) * BF4A(y, ni, li, nf, lf, L);
}


// --- Wrappers ---
// BesselMatrixElement
double BesselMatrixElement(int ni, int li, int bi, int nf, int lf, int bf, int L, double q) 
// Return the matrix element :
//   <n' l' | j_L(qr) | n l>
{
    double b = bi; 
    double y = (q * b / 2.0) * (q * b / 2.0);
	ni = ni + 1;
	nf = nf + 1; 
    return BesselElement(y, ni, li, nf, lf, L);
}


// BesselMatrixElement_Minus
double BesselMatrixElement_Minus(int ni, int li, int bi, int nf, int lf, int bf, int L, double q) 
// Return the matrix element :
//   <n' l' | j_L(qr)(d_r - l/r) | n l>
{
    double b = bi; 
    double y = (q * b / 2.0) * (q * b / 2.0);
	ni = ni + 1;
	nf = nf + 1;
    return BesselElementMinus(y, ni, li, nf, lf, L);
}


// BesselMatrixElement_Plus
double BesselMatrixElement_Plus(int ni, int li, int bi, int nf, int lf, int bf, int L, double q) 
// Return the matrix element :
//   <n' l' | j_L(qr)(d_r + (l+1)/r) | n l>
{
    double b = bi; 
    double y = (q * b / 2.0) * (q * b / 2.0);
	ni = ni + 1;
	nf = nf + 1;
    return BesselElementPlus(y, ni, li, nf, lf, L);
}


// We give the 4 'Bessel' matrix elements which appear in [3]
// BesselMatrixElement_Minus_Minus = <n' l' j' | j_L(rho)(d_rho - (l+1)/rho)(d_rho - l/rho) | n l j> ; 
// BesselMatrixElement_Minus_Plus = <n' l' j' | j_L(rho)(d_rho + (l+2)/rho)(d_rho - l/rho) | n l j> ;
// BesselMatrixElement_Plus_Plus = <n' l' j' | j_L(rho)(d_rho + (l)/rho)(d_rho + (l+1)/rho) | n l j> ;
// BesselMatrixElement_Plus_Minus = <n' l' j' | j_L(rho)(d_rho - (l-1)/rho)(d_rho + (l+1)/rho) | n l j> ;

// BesselMatrixElement_Minus_Minus
double BesselMatrixElement_Minus_Minus(int ni, int li, int bi, int nf, int lf, int bf, int L, double q)
// Return the matrix element :
//   <n' l' | j_L(qr)(d_r - (l+1)/r)(d_r - l/r) | n l>
{
	double me_tot ;

	return me_tot ;
};


// BesselMatrixElement_Minus_Plus
double BesselMatrixElement_Minus_Plus(int ni, int li, int bi, int nf, int lf, int bf, int L, double q)
// Return the matrix element :
//   <n' l' | j_L(qr)(d_r + (l+2)/r)(d_r - l/r) | n l>
{
	double me_tot ;

	return me_tot ;
};


// BesselMatrixElement_Plus_Plus
double BesselMatrixElement_Plus_Plus(int ni, int li, int bi, int nf, int lf, int bf, int L, double q)
// Return the matrix element :
//   <n' l' | j_L(qr)(d_r + (l)/r)(d_r + (l+1)/r) | n l>
{
	double me_tot ;

	return me_tot ;
};


// BesselMatrixElement_Plus_Minus
double BesselMatrixElement_Plus_Minus(int ni, int li, int bi, int nf, int lf, int bf, int L, double q)
// Return the matrix element :
//   <n' l' | j_L(qr)(d_r - (l-1)/r)(d_r + (l+1)/r) | n l>
{
	double me_tot ;

	return me_tot ;
};


// We give the 4 reduced matrix elements (eq.3 in ref [1])
// <n' l' j' || MJ(qr) || n l j> ; <n' l' j' || MJL(qr) σ || n l j>
// <n' l' j' || MJL(qr) ∇/q || n l j> ; <n' l' j' || MJ(qr) σ ∇/q || n l j>
// MJ_MatrixElement
double MJ_MatrixElement(int ni, int li, double ji, double bi, int nf, int lf, double jf, double bf, int J, double q, int A)
// Calculate the (reduced) matrix element :
//    <n' l' j' || MJ(qr) || n l j>
{
	double MJ = 0.0 ;
	if (J>=0) { 
		double j6_symbol = am::Wigner6J2(2*lf, 2*jf, 2*0.5, 2*ji, 2*li, 2*J) ;
		double j3_symbol = am::Wigner3J2(2*lf, 2*J, 2*li, 0, 0, 0) ;
		double prefactor = (1/std::sqrt(4 * M_PI) * std::pow(-1, J+ji+0.5) 
						   * am::Hat2(2*lf) * am::Hat2(2*li) * am::Hat2(2*jf) * am::Hat2(2*ji) * am::Hat2(2*J)
						   * j6_symbol * j3_symbol) ;

		// MJ matrix element calculation
		MJ = prefactor * BesselMatrixElement(ni, li, bi, nf, lf, bf, J, q);
	}
	return MJ;
};


// MJLSigma_MatrixElement
double MJLSigma_MatrixElement(int ni, int li, double ji, double bi, int nf, int lf, double jf, double bf, int J, int L, double q, int A)
// Calculate the (reduced) matrix element :
//    <n' l' j' || MJL(qr) \sigma || n l j>
{
	double MJLSigma = 0.0 ;
	if (J>=0 && L>=0) {
		double j9_symbol = am::Wigner9J2(2*lf, 2*li, 2*L, 2*0.5, 2*0.5, 2*1, 2*jf, 2*ji, 2*J) ;
		double j3_symbol = am::Wigner3J2(2*lf, 2*L, 2*li, 0, 0, 0) ;
		double prefactor = (1/std::sqrt(4 * M_PI) * std::pow(-1, lf) * std::sqrt(6)
						   * am::Hat2(2*lf) * am::Hat2(2*li) * am::Hat2(2*jf) * am::Hat2(2*ji) * am::Hat2(2*L) * am::Hat2(2*J)
						   * j9_symbol * j3_symbol) ; 
		
		// MJLSigma matrix element Calculation
		MJLSigma = prefactor * BesselMatrixElement(ni, li, bi, nf, lf, bf, L, q);
	}
	return MJLSigma;
};


// MJLNabla_MatrixElement
double MJLNabla_MatrixElement(int ni, int li, double ji, double bi, int nf, int lf, double jf, double bf, int J, int L, double q, int A)
// Calculate the (reduced) matrix element :
//    <n' l' j' || MJL(qr) \nabla / q || n l j>
{
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
		double term1 = TranslationallyInvariantTerm(A) * BesselMatrixElement_Minus(ni, li, bi, nf, lf, bf, L, q);

		// Calculation of 2nd term
		// li > 0 
		double prefactor_term2 = 0.0 ;
		double term2 = 0.0 ;
		if (li>0) {		 
			double j6_symbol_term2 = am::Wigner6J2(2*L, 2*1, 2*J, 2*li, 2*lf, 2*(li-1)) ;
			double j3_symbol_term2 = am::Wigner3J2(2*lf, 2*L, 2*(li-1), 0, 0, 0) ;
			prefactor_term2 = std::sqrt(li) * am::Hat2(2*(li-1)) * j6_symbol_term2 * j3_symbol_term2 ; 
			term2 = TranslationallyInvariantTerm(A) * BesselMatrixElement_Plus(ni, li, bi, nf, lf, bf, L, q);
		}

		// MJLNabla matrix element calculation
		MJLNabla = prefactor * (
							   - prefactor_term1 * term1
							   + prefactor_term2 * term2
							   );
	}
	return MJLNabla;
};


// MJSigmaNabla_MatrixElement
double MJSigmaNabla_MatrixElement(int ni, int li, double ji, double bi, int nf, int lf, double jf, double bf, int J, double q, int A)
// Calculate the (reduced) matrix element :
//    <n' l' j' || MJ(qr) \sigma \nabla / q|| n l j>
{
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
			term1 = TranslationallyInvariantTerm(A) * BesselMatrixElement_Minus(ni, li, bi, nf, lf, bf, J, q);
		}
		
		// Calculation of 2nd term
		if (ji==(li-0.5)) {
			term2 = TranslationallyInvariantTerm(A) * BesselMatrixElement_Plus(ni, li, bi, nf, lf, bf, J, q);
		}

		// MJSigmaNabla matrix element calculation
		MJSigmaNabla = prefactor * (-term1 + term2);
	}
	return MJSigmaNabla;
};


// We give the 6 others (reduced) matrix elements which appear in ref [3]
// for calculate the other 'seven' operators in ref [3]
// <n' l' j' || MJ(qr) (∇/q)^2 || n l j>  ; <n' l' j' || (MJL(qr) σ) (∇/q)^2 || n l j>
// <n' l' j' || (MJL(qr) ∇/q) (σ ∇/q) || n l j>  ; <n' l' j' || i MJL(qr) (σ x ∇/q) || n l j>
// <n' l' j' || [MK(qr) σ)_L ∇/q]_J || n l j>  ; <n' l' j' || [MK(qr) ∇/q)_L σ]_J || n l j>
/*
// MJNablaSquare_MatrixElement
double MJNablaSquare_MatrixElement(int ni, int li, double ji, double bi, int nf, int lf, double jf, double bf, int J, double q)
// Calculate the (reduced) matrix element :
//    <n' l' j' || MJ(qr) (∇/q)^2 || n l j>
{
	double MJNablaSquare = 0.0 ;

	return MJNablaSquare;
};


// MJLSigmaNablaSquare_MatrixElement
double MJLSigmaNablaSquare_MatrixElement(int ni, int li, double ji, double bi, int nf, int lf, double jf, double bf, int J, int L, double q)
// Calculate the (reduced) matrix element :
//    <n' l' j' || (MJL(qr) σ) (∇/q)^2 || n l j>
{
	double MJLSigmaNablaSquare = 0.0 ;

	return MJLSigmaNablaSquare;
};	


// MJLSigmaNablaSigmaNabla_MatrixElement
double MJLNablaSigmaNabla_MatrixElement(int ni, int li, double ji, double bi, int nf, int lf, double jf, double bf, int J, int L, double q)
// Calculate the (reduced) matrix element :
//    <n' l' j' || (MJL(qr) ∇/q) (σ ∇/q) || n l j>
{
	double MJLSigmaNablaSigmaNabla = 0.0 ;
	double prefactor = (1/std::sqrt(4 * M_PI) * std::pow(-1, li+lf+1) * std::sqrt(6)
					   * am::Hat2(2*J) * am::Hat2(2*li) * am::Hat2(2*L) * am::Hat2(2*ji) * am::Hat2(2*jf) ) ;

	double j6_symbol1 = am::Wigner6J2(2*lf, 2*(li+1), 2*J, 2*ji, 2*jf, 2*0.5) ;
	double j6_symbol2 = am::Wigner6J2(2*li, 2*1, 2*(li+1), 2*0.5, 2*ji, 2*0.5) ;
	double prefactor1 = j6_symbol1 * j6_symbol2 * std::sqrt((li + 1) * (2*li + 3))
	
	double j6_symbol_term1 = am::Wigner6J2(2*L, 2*J, 2*1.0, 2*(li+1), 2*(li+2), 2*lf) ;
	double j3_symbol1_term1 = am::Wigner3J2(2*lf,2*L, 2*(li+2), 0, 0, 0) ;
	double j3_symbol2_term1 = am::Wigner3J2(2*(li+2), 2*1.0, 2*(li+1), 0, 0, 0) ;
	double prefactor_term1 = j6_symbol_term1 * j3_symbol1_term1 / j3_symbol2_term1 * (li + 2) / (2*li + 3) ;
	double term1 = prefactor_term1 * BesselMatrixElement_Minus_Minus(ni, li, bi, nf, lf, bf, L, q);

	double j6_symbol_term2 = am::Wigner6J2(2*L, 2*J, 2*1.0, 2*(li+1), 2*li, 2*lf) ;
	double j3_symbol1_term2 = am::Wigner3J2(2*lf,2*L, 2*li, 0, 0, 0) ;
	double j3_symbol2_term2 = am::Wigner3J2(2*li, 2*1.0, 2*(li+1), 0, 0, 0) ;
	double prefactor_term2 = j6_symbol_term2 * j3_symbol1_term2 / j3_symbol2_term2 * (li + 1) / (2*li + 3) ;
	double term2 = prefactor_term2 * BesselMatrixElement_Minus_Plus(ni, li, bi, nf, lf, bf, L, q);

	double j6_symbol3 = am::Wigner6J2(2*lf, 2*(li-1), 2*J, 2*ji, 2*jf, 2*0.5) ;
	double j6_symbol4 = am::Wigner6J2(2*li, 2*1, 2*(li-1), 2*0.5, 2*ji, 2*0.5) ;
	double prefactor2 = j6_symbol3 * j6_symbol4 * std::sqrt(li * (2*li + 1)) ;

	double j6_symbol_term3 = am::Wigner6J2(2*L, 2*J, 2*1.0, 2*(li-1), 2*li, 2*lf) ;
	double j3_symbol1_term3 = am::Wigner3J2(2*lf,2*L, 2*li, 0, 0, 0) ;
	double j3_symbol2_term3 = am::Wigner3J2(2*li, 2*1.0, 2*(li-1), 0, 0, 0) ;
	double prefactor_term3 = j6_symbol_term3 * j3_symbol1_term3 / j3_symbol2_term3 * li / (2*li - 1) ;
	double term3 = prefactor_term3 * BesselMatrixElement_Plus_Minus(ni, li, bi, nf, lf, bf, L, q);

	double j6_symbol_term4 = am::Wigner6J2(2*L, 2*J, 2*1.0, 2*(li-1), 2*(li-2), 2*lf) ;
	double j3_symbol1_term4 = am::Wigner3J2(2*lif,2*L, 2*(li-2), 0, 0, 0) ;
	double j3_symbol2_term4 = am::Wigner3J2(2*(li-2), 2*1.0, 2*(li-1), 0, 0, 0) ;
	double prefactor_term4 = j6_symbol_term4 * j3_symbol1_term4 / j3_symbol2_term4 * (li - 1) / (2*li - 1) ;
	double term4 = prefactor_term4 * BesselMatrixElement_Plus_Plus(ni, li, bi, nf, lf, bf, L, q);

	MJLSigmaNablaSigmaNabla = prefactor * (
									       prefactor1 * (term1 + term2)
									       - prefactor2 * (term3 +  term4)
									   );

	return MJLSigmaNablaSigmaNabla;
};


// MJLSigmaCrossNabla_MatrixElement
double MJLSigmaCrossNabla_MatrixElement(int ni, int li, double ji, double bi, int nf, int lf, double jf, double bf, int J, int L, double q)
// Calculate the (reduced) matrix element :
//    <n' l' j' || i MJL(qr) (σ x ∇/q) || n l j>
{
	double MJLSigmaCrossNabla = 0.0 ;

	return MJLSigmaCrossNabla;
};


// MKSigmaLNablaJ_MatrixElement
double MKSigmaLNablaJ_MatrixElement(int ni, int li, double ji, double bi, int nf, int lf, double jf, double bf, int J, int K, double q)
// Calculate the (reduced) matrix element :
//    <n' l' j' || [MK(qr) σ)_L ∇/ q]_J || n l j>
{	
	double MKSigmaLNablaJ = 0.0 ;

	return MKSigmaLNablaJ;
};


// MKNablaLSigmaJ_MatrixElement
double MKNablaLSigmaJ_MatrixElement(int ni, int li, double ji, double bi, int nf, int lf, double jf, double bf, int J, int K, int L, double q)
// Calculate the (reduced) matrix element :
//    <n' l' j' || [MK(qr) ∇/q)_L σ]_J || n l j>
{
	double MKSigmaNablaLJ = 0.0 ;

	return MKSigmaNablaLJ;
};
*/

// Seven basis single-particle operators
// Here we calculate <n' l' j' || \hat{O}_J(qr) || n l j>
// where the operator \hat{O}_J(qr) corresponds to eq.(1) in ref [1] :
// MUST satisfy the Normal parity : M_J(qr) ; Δ'_J(qr) ; Σ_J(qr)
// MUST satisfy the Abnormal parity : Δ_J(qr) ; Σ'_J(qr) ; Σ''_J(qr) ; Ω_J(qr) ; Ω'_J(qr)

// MJ_SevenOprator
double MJ_SevenOperator(int ni, int li, double ji, double bi, int nf, int lf, double jf, double bf, int J, double q, int A)
// Calculate the (reduced) matrix element : 
//    <n' l' j' || M_J(qr) || n l j>
{
    double MJ = MJ_MatrixElement(ni, li, ji, bi, nf, lf, jf, bf, J, q, A) ;

	return MJ ;
};


// DeltaJ_SevenOperator
double DeltaJ_SevenOperator(int ni, int li, double ji, double bi, int nf, int lf, double jf, double bf, int J, double q, int A)
// Calculate the (reduced) matrix element : 
//    <n' l' j' || Δ_J(qr) || n l j>
{
	double DeltaJ = TranslationallyInvariantTerm(A) * MJLNabla_MatrixElement(ni, li, ji, bi, nf, lf, jf, bf, J, J, q, A) ;
	
	return DeltaJ ;
};


// DeltaJP_SevenOperator
double DeltaJP_SevenOperator(int ni, int li, double ji, double bi, int nf, int lf, double jf, double bf, int J, double q, int A)
// Calculate the (reduced) matrix element : 
//    <n' l' j' || Δ'_J(qr) || n l j>
{
	double DeltaJP = TranslationallyInvariantTerm(A) * 1/am::Hat2(2*J) * (
										- std::sqrt(J) * MJLNabla_MatrixElement(ni, li, ji, bi, nf, lf, jf, bf, J, J+1, q, A) 
									    + std::sqrt(J + 1) * MJLNabla_MatrixElement(ni, li, ji, bi, nf, lf, jf, bf, J, J-1, q, A));
	
	return DeltaJP ;
};


// SigmaJ_SevenOperator
double SigmaJ_SevenOperator(int ni, int li, double ji, double bi, int nf, int lf, double jf, double bf, int J, double q, int A)
// Calculate the (reduced) matrix element : 
//    <n' l' j' || Σ_J(qr) || n l j>
{
	double SigmaJ = MJLSigma_MatrixElement(ni, li, ji, bi, nf, lf, jf, bf, J, J, q, A) ;
	
	return SigmaJ ;
};


// SigmaJP_SevenOperator
double SigmaJP_SevenOperator(int ni, int li, double ji, double bi, int nf, int lf, double jf, double bf, int J, double q, int A)
// Calculate the (reduced) matrix element : 
//    <n' l' j' || Σ'_J(qr) || n l j>
{
	double SigmaJP = TranslationallyInvariantTerm(A) * 1/am::Hat2(2*J) * (
										- std::sqrt(J) * MJLSigma_MatrixElement(ni, li, ji, bi, nf, lf, jf, bf, J, J+1, q, A) 
										+ std::sqrt(J + 1) * MJLSigma_MatrixElement(ni, li, ji, bi, nf, lf, jf, bf, J, J-1, q, A));
	
	return SigmaJP ;
};


// SigmaJPP_SevenOperator
double SigmaJPP_SevenOperator(int ni, int li, double ji, double bi, int nf, int lf, double jf, double bf, int J, double q, int A)
// Calculate the (reduced) matrix element : 
//    <n' l' j' || Σ''_J(qr) || n l j>
{
	double SigmaJPP = TranslationallyInvariantTerm(A) * 1/am::Hat2(2*J) * (
										 std::sqrt(J + 1) * MJLSigma_MatrixElement(ni, li, ji, bi, nf, lf, jf, bf, J, J+1, q, A) 
										 + std::sqrt(J) * MJLSigma_MatrixElement(ni, li, ji, bi, nf, lf, jf, bf, J, J-1, q, A));
	
	return SigmaJPP ;
};


// OmegaJ_SevenOperator
double OmegaJ_SevenOperator(int ni, int li, double ji, double bi, int nf, int lf, double jf, double bf, int J, double q, int A)
// Calculate the (reduced) matrix element : 
//    <n' l' j' || Ω_J(qr) || n l j>
{
	double OmegaJ = MJSigmaNabla_MatrixElement(ni, li, ji, bi, nf, lf, jf, bf, J, q, A) ;
		
	return OmegaJ ;
};


// OmegaJP_SevenOperator
double OmegaJP_SevenOperator(int ni, int li, double ji, double bi, int nf, int lf, double jf, double bf, int J, double q, int A)
// Calculate the (reduced) matrix element : 
//    <n' l' j' || Ω'_J(qr) || n l j>
{
	double OmegaJP = (OmegaJ_SevenOperator(ni, li, ji, bi, nf, lf, jf, bf, J, q, A) 
					  + 0.5 * SigmaJPP_SevenOperator(ni, li, ji, bi, nf, lf, jf, bf, J, q, A)) ;
	
	return OmegaJP ;
};

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  double y = 0.25 ;
  double b = 1.0 ;
  double q = 1.0 ;
  int A = 6 ;
  // Bessel Matrix Element
  //Test_BesselMatrixElement(0, 0, 1.0, 0, 0, 1.0, 0, q);
  //Test_BesselMatrixElement(0, 0, 1.0, 0, 0, 1.0, 1, q);
  //Test_BesselMatrixElement(0, 0, 1.0, 0, 0, 1.0, 2, q);

  // Bessel Matrix Element Minus
  //Test_BesselMatrixElement_Minus(0, 0, 1.0, 0, 0, 1.0, 0, q);
  //Test_BesselMatrixElement_Minus(0, 0, 1.0, 0, 0, 1.0, 1, q);
  //Test_BesselMatrixElement_Minus(0, 0, 1.0, 0, 0, 1.0, 2, q);
  
  //Test_BesselMatrixElement_Minus(1, 0, 1.0, 0, 0, 1.0, 0, q);

  // MJ Matrix Element
  //Test_MJ_MatrixElement(0, 0, 0.5, b, 0, 0, 0.5, b, 0, q);

  /*
  // MJ Seven Operator
  Test_MJ_SevenOperator(0, 0, 0.5, b, 0, 0, 0.5, b, 0, q);
  Test_MJ_SevenOperator(0, 0, 0.5, b, 0, 1, 0.5, b, 1, q);

  // DeltaJ Seven Operator
  Test_DeltaJ_SevenOperator(0, 1, 0.5, b, 0, 1, 0.5, b, 1, q);
  Test_DeltaJ_SevenOperator(0, 1, 0.5, b, 0, 2, 1.5, b, 2, q);

  // DeltaJP Seven Operator
  Test_DeltaJP_SevenOperator(0, 0, 0.5, b, 0, 2, 1.5, b, 2, q);
  Test_DeltaJP_SevenOperator(0, 0, 0.5, b, 0, 1, 0.5, b, 1, q);

  // SigmaJ Seven Operator
  Test_SigmaJ_SevenOperator(0, 0, 0.5, b, 0, 2, 1.5, b, 2, q);
  Test_SigmaJ_SevenOperator(0, 0, 0.5, b, 0, 1, 0.5, b, 1, q);

  // SigmaJP Seven Operator
  Test_SigmaJP_SevenOperator(0, 1, 0.5, b, 0, 1, 0.5, b, 1, q);
  Test_SigmaJP_SevenOperator(0, 1, 0.5, b, 0, 2, 1.5, b, 2, q);

  // SigmaJPP Seven Operator
  Test_SigmaJPP_SevenOperator(0, 1, 0.5, b, 0, 1, 0.5, b, 1, q);
  Test_SigmaJPP_SevenOperator(0, 1, 0.5, b, 0, 2, 1.5, b, 2, q);

  // OmegaJ Seven Operator
  Test_OmegaJ_SevenOperator(0, 1, 0.5, b, 0, 1, 0.5, b, 1, q);
  Test_OmegaJ_SevenOperator(0, 1, 0.5, b, 0, 2, 1.5, b, 2, q);

  // OmegaJP Seven Operator
  Test_OmegaJP_SevenOperator(0, 1, 0.5, b, 0, 1, 0.5, b, 1, q);
  Test_OmegaJP_SevenOperator(0, 1, 0.5, b, 0, 2, 1.5, b, 2, q);
  */

  std::cout << "Example p111 Donnelly & Haxton, 1979" << std::endl;
  std::cout << "DeltaJ_SevenOperator : " << DeltaJ_SevenOperator(0, 2, 2.5, b, 0, 2, 2.5, b, 1, q, A) << std::endl;
  std::cout << "Expected value : " << 1.0/std::sqrt(4.0*M_PI) * std::pow(y, (1.0-1.0)/2.0) * std::exp(-y) * 2.0/5.0 * std::sqrt(35) * (-1.0 + 2.0/5.0*y) << "\n" << std::endl;

  std::cout << "SigmaJP_SevenOperator : " << SigmaJP_SevenOperator(0, 2, 2.5, b, 0, 2, 2.5, b, 1, q, A) << std::endl;
  std::cout << "Expected value : " << 1/std::sqrt(4.0*M_PI) * std::pow(y, (1.0-1.0)/2.0) * std::exp(-y) * 2.0/5.0 * std::sqrt(35) * (1.0 - 8.0/5.0*y + 12.0/35.0*std::pow(y, 2)) << "\n" << std::endl;


  std::cout << "Example Table III" << std::endl;
  std::cout << "MJ_SevenOperator : " << MJ_SevenOperator(0, 0, 0.5, b, 0, 1, 0.5, b, 1, q, A) << std::endl;
  std::cout << "Expected value : " << 1/std::sqrt(4.0*M_PI) * std::pow(y, (1.0-2.0)/2.0) * std::exp(-y) * std::sqrt(3.0) * 2.0/3.0 * (0.0 - 1.0*y) << "\n" << std::endl;
  
  std::cout << "DeltaJP_SevenOperator : " << DeltaJP_SevenOperator(0, 0, 0.5, b, 0, 1, 0.5, b, 1, q, A) << std::endl;
  std::cout << "Expected value : " << 1/std::sqrt(4.0*M_PI) * std::pow(y, (1.0-2.0)/2.0) * std::exp(-y) * 1.0/6.0 * std::sqrt(6.0) * (1.0 + 0.0*y) << "\n" << std::endl;
  
  std::cout << "SigmaJ_SevenOperator : " << SigmaJ_SevenOperator(0, 0, 0.5, b, 0, 1, 0.5, b, 1, q, A) << std::endl;
  std::cout << "Expected value : " << 1/std::sqrt(4.0*M_PI) * std::pow(y, (1.0-2.0)/2.0) * std::exp(-y) * 2.0/3.0 * std::sqrt(6.0) * (0.0 + 1.0*y) << "\n" << std::endl;

  
  std::cout << "Example Table VI" << std::endl;
  std::cout << "DeltaJ_SevenOperator : " << DeltaJ_SevenOperator(0, 1, 0.5, b, 0, 1, 1.5, b, 1, q, A) << std::endl;
  std::cout << "Expected value : " << 1/std::sqrt(4.0*M_PI) * std::pow(y, (1.0-1.0)/2.0) * std::exp(-y) * 1.0/3.0 * std::sqrt(2.0) * (-1.0 + 0.0*y) << "\n" << std::endl;
  
  std::cout << "SigmaJP_SevenOperator : " << SigmaJP_SevenOperator(0, 1, 0.5, b, 0, 1, 1.5, b, 1, q, A) << std::endl;
  std::cout << "Expected value : " << 1/std::sqrt(4.0*M_PI) * std::pow(y, (1.0-1.0)/2.0) * std::exp(-y) * 4.0/3.0 * std::sqrt(2.0) * (-1.0 + 0.5*y) << "\n" << std::endl;
  
  std::cout << "SigmaJPP_SevenOperator : " << SigmaJPP_SevenOperator(0, 1, 0.5, b, 0, 1, 1.5, b, 1, q, A) << std::endl;
  std::cout << "Expected value : " << 1/std::sqrt(4.0*M_PI) * std::pow(y, (1.0-1.0)/2.0) * std::exp(-y) * 4.0/3.0 * std::sqrt(1.0) * (-1.0 + 1.0*y) << "\n" << std::endl;

  std::cout << "OmegaJP_SevenOperator : " << OmegaJP_SevenOperator(0, 1, 0.5, b, 0, 1, 1.5, b, 1, q, A) << std::endl;
  std::cout << "Expected value : " << 1/std::sqrt(4.0*M_PI) * std::pow(y, (1.0-1.0)/2.0) * std::exp(-y) * 1.0 * std::sqrt(1.0) * (-1.0 + 0.0*y) << "\n" << std::endl;

  
  std::cout << "Example Table IX" << std::endl;
  std::cout << "MJ_SevenOperator : " << MJ_SevenOperator(0, 1, 0.5, b, 0, 2, 2.5, b, 3, q, A) << std::endl;
  std::cout << "Expected value : " << 1/std::sqrt(4.0*M_PI) * std::pow(y, (3.0-2.0)/2.0) * std::exp(-y) * std::sqrt(15.0) * 4.0/15.0 * (0.0 - 1.0*y) << "\n" << std::endl;
  
  std::cout << "DeltaJP_SevenOperator : " << DeltaJP_SevenOperator(0, 1, 0.5, b, 0, 2, 2.5, b, 3, q, A) << std::endl;
  std::cout << "Expected value : " << 1/std::sqrt(4.0*M_PI) * std::pow(y, (3.0-2.0)/2.0) * std::exp(-y) * 2.0/15.0 * std::sqrt(5.0) * (1.0 + 0.0*y) << "\n" << std::endl;

  std::cout << "SigmaJ_SevenOperator : " << SigmaJ_SevenOperator(0, 1, 0.5, b, 0, 2, 2.5, b, 3, q, A) << std::endl;
  std::cout << "Expected value : " << 1/std::sqrt(4.0*M_PI) * std::pow(y, (3.0-2.0)/2.0) * std::exp(-y) * 8.0/15.0 * std::sqrt(5.0) * (0 - 1.0*y) << "\n" << std::endl;

  
  std::cout << "Example Table X" << std::endl;
  std::cout << "DeltaJ_SevenOperator : " << DeltaJ_SevenOperator(0, 1, 0.5, b, 0, 2, 2.5, b, 2, q, A) << std::endl;
  std::cout << "Expected value : " << 1/std::sqrt(4.0*M_PI) * std::pow(y, (2.0-1.0)/2.0) * std::exp(-y) * 2.0/15.0 * std::sqrt(10.0) * (-1.0 + 0.0*y) << "\n" << std::endl;

  std::cout << "SigmaJP_SevenOperator : " << SigmaJP_SevenOperator(0, 1, 0.5, b, 0, 2, 2.5, b, 2, q, A) << std::endl;
  std::cout << "Expected value : " << 1/std::sqrt(4.0*M_PI) * std::pow(y, (2.0-1.0)/2.0) * std::exp(-y) * 4.0/5.0 * std::sqrt(10.0) * (-1.0 + 1.0/3.0*y) << "\n" << std::endl;

  std::cout << "SigmaJPP_SevenOperator : " << SigmaJPP_SevenOperator(0, 1, 0.5, b, 0, 2, 2.5, b, 2, q, A) << std::endl;
  std::cout << "Expected value : " << 1/std::sqrt(4.0*M_PI) * std::pow(y, (2.0-1.0)/2.0) * std::exp(-y) * 8.0/15.0 * std::sqrt(15.0) * (-1.0 + 0.5*y) << "\n" << std::endl;

  std::cout << "OmegaJP_SevenOperator : " << OmegaJP_SevenOperator(0, 1, 0.5, b, 0, 2, 2.5, b, 2, q, A) << std::endl;
  std::cout << "Expected value : " << 1/std::sqrt(4.0*M_PI) * std::pow(y, (2.0-1.0)/2.0) * std::exp(-y) * 1.0/5.0 * std::sqrt(15.0) * (-1.0 + 0.0*y) << "\n" << std::endl;
    
  
  // termination
  return 0;
}
