/**************************************
 angular_momentum.cpp
 
 Angular momentum algebra utility functions.
                                 
 Created by Mark Caprio on 12/02/10
 Last changes by mac 2/16/11
                                   
 **************************************/

#include "angular_momentum.h"

using namespace std;

bool AllowedTriangle(const HalfInt& h1, const HalfInt& h2, const HalfInt& h3)
{
	bool triangular, proper_integrity;
	triangular = ((abs(h1-h2) <= h3) && (h3 <= (h1+h2)));
        proper_integrity = (h1.TwiceValue() + h2.TwiceValue() + h3.TwiceValue())%2 
		== 0;
	return triangular && proper_integrity;
}

vector<HalfInt> ProductAngularMomenta (const HalfInt& j1, const HalfInt& j2)
{
	int i = 0;
	vector<HalfInt> result;
	HalfInt jmax = j1+j2;
	HalfInt jmin = abs(j1-j2);
	for (HalfInt jsum = jmax; jsum >=jmin; --jsum) 
		result.push_back(jsum);
	return result;
}

HalfIntBound TriangleBound(const HalfInt& j1, const HalfInt& j2)
{
	return HalfIntBound(abs(j1-j2),j1+j2);
}
