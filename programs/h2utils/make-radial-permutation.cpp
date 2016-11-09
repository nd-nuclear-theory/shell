/****************************************************************
  make_radial_permutation.cpp

  Make radial overlaps for a basis reordering.

  Mark A. Caprio
  University of Notre Dame

  11/6/16 (mac): Created, based on make_dummy_radial.

****************************************************************/

#include <iomanip>
#include <iostream>
#include <string>

#include "eigen3/Eigen/Dense"

#include "basis/nlj_orbital.h"
#include "radial/radial_io.h"

////////////////////////////////////////////////////////////////
// test code
////////////////////////////////////////////////////////////////

void MakeRadialOut(
    const std::string& filename,
    int Nmax
  ) 
{
  // set up bra/ket space
  basis::OrbitalSpaceLJPN bra_subspace(Nmax);
  basis::OrbitalSpaceLJPN ket_subspace(Nmax);
  int l0max = 0;
  int Tz0 = 0;
  shell::RadialOperatorType operator_type = shell::RadialOperatorType::kO;
  basis::OrbitalSectorsLJPN sectors(bra_subspace,ket_subspace,l0max,Tz0);
  shell::OutRadialStream os(filename,bra_subspace,ket_subspace,sectors,operator_type);
  
  // generate matrices
  basis::MatrixVector matrices;
  basis::SetOperatorToZero(sectors,matrices);
  for (int sector_index=0; sector_index<matrices.size(); ++sector_index)
    {
      Eigen::MatrixXd& matrix = matrices[sector_index];
      assert(matrix.rows()==matrix.cols());
      int dimension = matrix.rows();

      for (int row_index=0; row_index<dimension; ++row_index)
        {
          matrix(row_index,dimension-row_index-1) = 1;  // make "slash" pattern matrix
        }

    }

  // write to file
  os.Write(matrices);
  os.Close();
}


int main(int argc, char **argv) {
  /// MakeRadialOut("test/radial-olap-permute-Nmax02.dat",2);
  MakeRadialOut("test/radial-olap-permute-Nmax06.dat",6);
}
