/******************************************************************************
  
  Created by M. A. Caprio, University of Notre Dame.

  4/23/11 (mac): Originated.
  10/9/16 (pjf): Rename mcpp -> mcutils.

******************************************************************************/


#include <iostream>
#include <iomanip>

#include "am/halfint.h"
#include "mcutils/profiling.h"

#include "legacy/pair_indexing.h"

using namespace std;
using namespace legacy;

int main(int argc, char **argv)
{

	cout << "****************" << endl;
	cout << "lookup table" << endl;

	PairLookupArray<float> pla(kUpperTriangular,4,3.14159);
	cout << pla.MinorDimension() << " " << pla.MajorDimension() << endl;
	cout << pla(0,0) << endl;
	pla(0,0)=2.7;
	cout << pla(0,0) << endl;
	// cout << pla(4,99999) << endl;  // note -- no range checking

	cout << "****************" << endl;
	cout << "assertions" << endl;
	// cout << pla(4,4) << endl;
	// cout << pla(3,2) << endl;

	// termination
	return 0;
}
