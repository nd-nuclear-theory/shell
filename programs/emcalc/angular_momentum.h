/**************************************
  angular_momentum.h
 
  Angular momentum algebra utility functions.
                                  
  Created by Mark Caprio on 12/02/10.
  Code extracted from halfint originally by Ke Cai on 6/18/10:
    - Renamed functions addj() to ProductAngularMomenta()
      and validTriple() to AllowedTriangle()
    - Simplified implementation of ProductAngularMomenta() 
      and AllowedTriangle()
  2/20/11 (mac): Revisions to integration with halfint.

**************************************/

#ifndef angular_momentum_h
#define angular_momentum_h

#include <cmath>
#include <iostream>
#include <vector>

#include "halfint.h"
#include "halfint_bound.h"


////////////////////////////////
// angular momentum arithmetic
////////////////////////////////

// Overloading: Versions for both halfint and int arguments are
// provided, though the int version is strictly unnecessary due to
// automatic type conversion from halfint to int.

// angular momentum hat symbol

inline
double Hat(const HalfInt& j)
{
	return sqrt(static_cast<double>(TwiceValue(j)+1));
}

inline
double Hat(int j)
{
	return sqrt(static_cast<double>(2*j+1));
}

// phase sign (-)^sum

// DEBUGGING: Note that C++ % operator is not guaranteed positive
// definite unless both arguments are nonnegative.  Although taking
// abs() does not in general preserve modular equivalency class, it
// does suffice for the present purpose (checking evenness or
// oddness).  Without abs(), e.g., ParitySign(-1) can result in
// failure, from a remainder result other than 0 or 1.

inline 
int ParitySign(const HalfInt& sum)
{

	int remainder = abs(IValue(sum)) % 2;
	return (remainder == 0) ? +1 : -1;
}

inline 
int ParitySign(int sum)
{

	int remainder = abs(sum) % 2;
	return (remainder == 0) ? +1 : -1;
}


////////////////////////////////
// triangle inequality
////////////////////////////////

// Checks if the three HalfInts are coupled legally, i.e. they form a closed triangle and the parity is valid.
// Returns true if both triangle condition and parity condition are met, false if either is not met
bool AllowedTriangle(const HalfInt&, const HalfInt&, const HalfInt&);

// Returns a vector of angular momenta that j1 and j2 can be coupled to under the triangle inequality
// CONSIDER: generating from TriangleBound
// TODO: reverse order of iteration so j increasing
std::vector<HalfInt> ProductAngularMomenta(const HalfInt&, const HalfInt&);

// Returns interval of angular momenta that j1 and j2 can be coupled to under the triangle inequality
HalfIntBound TriangleBound(const HalfInt&, const HalfInt&);

#endif
