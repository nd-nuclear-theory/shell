#ifndef radial_integrals_h
#define radial_integrals_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm> 
#include <cmath>

using namespace std;

class RadialIntegrals
{
 public:
  void Initialize(istream& );
  double GetElement(int& ,int& ,int& ,int&);
  void Free();

 private:
  int lmax_, dl_, nmax_, mmax_;
  vector<vector<vector<vector<double> > > > matrix_elements_;
};

#endif
