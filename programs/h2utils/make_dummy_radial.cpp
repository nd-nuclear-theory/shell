/****************************************************************
  make_dummy_radial.cpp

  Make dummy identity radial matrices for h2xform development.

  Mark A. Caprio
  University of Notre Dame

  10/30/16 (mac): Created, based on radial_io_test.

****************************************************************/

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "radial/radial_io.h"

#include "eigen3/Eigen/Core"
#include "basis/nlj_orbital.h"

////////////////////////////////////////////////////////////////
// test code
////////////////////////////////////////////////////////////////

void MakeRadialOut(
    const std::string& filename,
    shell::RadialOperatorType operator_type, int power,
    int Nmax
  ) 
{
  // set up bra/ket space
  basis::OrbitalSpaceLJPN bra_subspace(Nmax);
  basis::OrbitalSpaceLJPN ket_subspace(Nmax);
  int l0max = power;
  int Tz0 = 0;
  basis::OrbitalSectorsLJPN sectors(bra_subspace,ket_subspace,l0max,Tz0);
  shell::OutRadialStream os(filename,bra_subspace,ket_subspace,sectors,operator_type);
  
  // generate matrices
  basis::MatrixVector matrices;
  basis::SetOperatorToIdentity(sectors,matrices);

  // write to file
  os.Write(matrices);
  os.Close();
}


int main(int argc, char **argv) {
  int Nmax=10;
  MakeRadialOut("test/radial-me-dummy-r1-Nmax10.dat",shell::RadialOperatorType::kR,1,Nmax);
  MakeRadialOut("test/radial-me-dummy-r2-Nmax10.dat",shell::RadialOperatorType::kR,2,Nmax);
  MakeRadialOut("test/radial-me-dummy-k1-Nmax10.dat",shell::RadialOperatorType::kK,1,Nmax);
  MakeRadialOut("test/radial-me-dummy-k2-Nmax10.dat",shell::RadialOperatorType::kK,2,Nmax);
  MakeRadialOut("test/radial-xform-dummy-Nmax10.dat",shell::RadialOperatorType::kR,0,Nmax);
}
