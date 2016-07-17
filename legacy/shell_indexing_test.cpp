/******************************************************************************

  test_indexing.cpp
  
  M. A. Caprio, University of Notre Dame.

  2/15/11 (mac): Originated.
  4/25/25 (mac): Tests of lj space indexing.

******************************************************************************/


#include <shell/shell_indexing_nlj.h>
#include <shell/shell_indexing_lj.h>
#include <shell/shell_indexing_nl.h>

#include <iostream>
#include <iomanip>

#include <halfint/halfint.h>
#include <mcpp/profiling.h>

using namespace std;
using namespace shell;

int main(int argc, char **argv)
{

  ////////////////////////////////
  // dimension checks
  ////////////////////////////////

  // basic tests

  cout << "Nlj orbital counts" << endl;
  for (int N = 0; N<6; ++N)
    cout << LevelCountNlj(N) << " ";
  cout << endl;

  ////////////////////////////////
  // Nlj indexing
  ////////////////////////////////

  cout << "****************" << endl;
  cout << "Nlj basic indexing tests" << endl;

  SPOrbitalNlj x(0);
  SPOrbitalNlj y(0,HalfInt(1,2));
  cout << x.GetIndex() << " " << x.GetIndex1() << " " << x.GetN() 
       << " " << x.Getl() << " " << x.Getj() 
       << " " << x 
       << endl;
  cout << y.GetIndex() << y << endl;
  ++y;
  cout << y.GetIndex() << y << endl;
  SPOrbitalNlj y2(2,HalfInt(1,2));
  cout << y2 << endl;

  for (int k = 0; k<=11; ++k)
    cout << SPOrbitalNlj(k) << " " << k << endl; 
  cout << endl;

  // example use of temporary instance for calculation
  cout << SPOrbitalNlj(10) << endl;

  cout << "****************" << endl;
  cout << "Nl basic indexing tests" << endl;
			
  for (int N = 0; N<6; ++N)
    cout << " " << LevelCountNl(N);
  cout << endl;

  for (int k = 1; k<=12; ++k)
    cout << SPOrbitalNl(k) << " " << k << endl; 
  cout << endl;

  for (SPOrbitalNl orbital(1); orbital.GetIndex() <=12; ++orbital)
    cout << orbital << " " << orbital.GetIndex() << endl; // << " " << orbital.GetIndex1()
  cout << endl;

  SPOrbitalNl z(5,3);
  cout << z << " " << z.GetIndex() << endl;


  ////////////////////////////////
  // lj indexing
  ////////////////////////////////

  cout << "****************" << endl;
  cout << "lj basic indexing tests" << endl;

  SPSpacelj sx;
  SPSpacelj sy(2);
  SPSpacelj sz(1,HalfInt(3,2));
  cout << sx << " " << sy  << " " << sz
       << endl;

  for (SPSpacelj space; space <= SPSpaceljMax(4); ++space)
    cout << "  " << space << endl; 
  cout << endl;

  // termination
  return 0;
}
