/******************************************************************************
  
  Created by M. A. Caprio, University of Notre Dame.

  4/23/11 (mac): Originated.
  Last modified 4/23/11.

******************************************************************************/


#include <iostream>
#include <iomanip>

#include <halfint/halfint.h>
#include <mcpp/profiling.h>

#include <shell/pair_indexing.h>

using namespace std;
using namespace shell;

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
