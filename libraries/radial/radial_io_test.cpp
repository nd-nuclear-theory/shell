/****************************************************************
  radial_io_test.cpp

  Patrick J. Fasano
  University of Notre Dame

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

int element_index = 0;
Eigen::MatrixXd GenerateDummyMatrix(int rows, int columns) {
  Eigen::MatrixXd dummy(rows, columns);

  for (int j=0; j < rows; ++j) {
    for (int k=0; k < columns; ++k) {
      dummy(j, k) = ++element_index;
    }
  }
  return dummy;
}

basis::MatrixVector TestRadialOut(const std::string& filename, bool verbose = false) {
  std::cout << "Radial matrix elements -- write test" << std::endl;

  // set up bra space
  std::cout << "Bra space" << std::endl;
  int bra_Nmax = 20;
  basis::OrbitalSpaceLJPN bra_space(bra_Nmax);
  std::cout << bra_space.DebugStr();

  // set up ket space
  std::cout << "Ket space" << std::endl;
  int ket_Nmax = 20;
  basis::OrbitalSpaceLJPN ket_space(ket_Nmax);
  std::cout << ket_space.DebugStr();

  // set up operator
  std::cout << "Operator" << std::endl;
  shell::RadialOperator rSquared('r', 2, 0);
  std::cout << rSquared.DebugStr();

  // set up output stream
  std::cout << "Output stream" << std::endl;
  shell::OutRadialStream os(filename, bra_space, ket_space, rSquared);

  if (verbose) {
    // print sectors
    std::cout << "Sectors" << std::endl;
    std::cout << os.sectors().DebugStr();
  }

  // generate matrices
  basis::MatrixVector matrices;
  for (int sector_index=0; sector_index < os.sectors().size(); ++sector_index) {
    // get sector
    const basis::BaseSector<basis::OrbitalSubspaceLJPN> sector =
        os.sectors().GetSector(sector_index);
    if (verbose) std::cout << "generating sector " << sector_index << std::endl;
    // std::cout << std::get<0>(sector.Key()) << " " << std::get<1>(sector.Key()) << std::endl;

    // generate stupid numbers
    Eigen::MatrixXd dummy_matrix = GenerateDummyMatrix(
      sector.bra_subspace().size(), sector.ket_subspace().size());

    matrices.push_back(dummy_matrix);
  }

  // write to file
  std::cout << "Writing to file" << std::endl;
  os.Write(matrices);

  return matrices;
}

basis::MatrixVector TestRadialIn(const std::string& filename) {
  std::cout << "Radial matrix elements -- read test" << std::endl;

  // set up input stream
  std::cout << "Input stream" << std::endl;
  shell::InRadialStream is(filename);

  // show bra space
  std::cout << "Bra space" << std::endl;
  std::cout << is.bra_orbital_space().DebugStr();

  // show ket space
  std::cout << "Ket space" << std::endl;
  std::cout << is.ket_orbital_space().DebugStr();

  // show sectors
  std::cout << "Sectors" << std::endl;
  std::cout << is.sectors().DebugStr();

  return is.Read();
}

int main(int argc, char **argv) {
  std::string filename("test/radial_out_test.dat");
  basis::MatrixVector matrices_out = TestRadialOut(filename, false);
  basis::MatrixVector matrices_in = TestRadialIn(filename);
  bool status = (matrices_out == matrices_in);

  // termination
  return !status;
}
