/******************************************************************************
  
  Created by M. A. Caprio
  University of Notre Dame
  November 27, 2010

  Last modified 12/24/10 (mac).
  Patch header path 7/4/16.

******************************************************************************/

#include "vector_tuple.h"

#include "am/halfint.h"

using namespace std;


int main(int argc, char **argv)
{

	// VectorTuple tests
	VectorTuple<int,3> z1;
	z1[0] = 1; z1[1] = 2; z1[2] = 3;
	cout << z1[0] << z1[1] << z1[2] << endl;
	VectorTuple<int,3> z2(5);
	cout << z2[0] << endl;

	VectorTuple<HalfInt,3> v1(HalfInt(1,2));
	VectorTuple<HalfInt,3> v3(HalfInt(3,2));
	VectorTuple<HalfInt,3> vs;
	cout << v1[0] << endl;
	vs = +v1;
	cout << vs[0] << endl;
	vs = -v1;
	cout << vs[0] << endl;
	vs = v3;
	vs -= v1;
	cout << vs[0] << endl;
	vs = v1 + v3;
	cout << vs[0] << endl;
	cout << (v1<v3) << (v1==v3) << (v1>v3) << endl;

	// try 1-based numbering
	VectorTuple<int,3,1> z3;
	z3[1] = 1; z3[2] = 2; z3[3] = 3;
	cout << z3[1] << z3[2] << z3[3] << endl;

	// test stream output
	cout << z3 << endl;
	cout << v1 << endl;
	VectorTuple<HalfInt,3>::SetDelimiters("< "," : "," >");
	cout << v1 << endl;

	cout << "****" << endl;

	cout << v1.size()  << endl;

	cout << "****" << endl;


	// demo of transform
	// vector<HalfInt> va(4,9), vb(4,11), vc(4);
	// transform(va.begin(),va.end(),vb.begin(),vc.begin(),BinaryPlus<HalfInt>);
	// cout << vc[0] << endl;


	// termination
	return 0;
}
