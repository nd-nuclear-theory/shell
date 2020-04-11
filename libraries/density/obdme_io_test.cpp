/****************************************************************
  obdme_io_test.cpp

  Patrick J. Fasano
  University of Notre Dame

****************************************************************/

#include <fstream>
#include <iomanip>
#include <iostream>

#include "density/obdme_io.h"
#include "mcutils/eigen.h"

////////////////////////////////////////////////////////////////
// test code
////////////////////////////////////////////////////////////////

void TestMultiFileRead(
    const std::string& data_filename, const std::string& info_filename,
    const basis::OrbitalSpaceLJPN& orbital_space,
    const HalfInt& J_bra, const int& g_bra, const int& n_bra,
    const HalfInt& J_ket, const int& g_ket, const int& n_ket
  )
{
  std::cout << "Reading obdme data" << std::endl;

  std::cout << "  constructing reader" << std::endl;
  shell::InOBDMEStreamMulti reader(
      info_filename, data_filename, orbital_space,
      J_bra, g_bra, n_bra,
      J_ket, g_ket, n_ket
    );

  // for (auto& info_line : reader.obdme_info()) {
  //   std::cout << "bra labels: " << std::endl;
  //   std::cout << "  sp: " << int(std::get<0>(info_line.bra_labels))
  //             << "  n:  " << std::get<1>(info_line.bra_labels)
  //             << "  l:  " << std::get<2>(info_line.bra_labels)
  //             << "  j:  " << std::get<3>(info_line.bra_labels).TwiceValue()
  //             << std::endl;

  //   std::cout << "ket labels: " << std::endl;
  //   std::cout << "  sp: " << int(std::get<0>(info_line.ket_labels))
  //             << "  n:  " << std::get<1>(info_line.ket_labels)
  //             << "  l:  " << std::get<2>(info_line.ket_labels)
  //             << "  j:  " << std::get<3>(info_line.ket_labels).TwiceValue()
  //             << std::endl;

  //   std::cout << "Multipole: " << info_line.multipole << std::endl << std::endl;
  // }
  std::cout << "  reading multipole 0" << std::endl;
  basis::OperatorBlocks<double> matrices;
  basis::OrbitalSectorsLJPN sectors;
  reader.GetMultipole(0, sectors, matrices);

  for (auto& matrix : matrices) {
    std::cout << mcutils::FormatMatrix(matrix, "16.8e") << std::endl;
  }
}

void TestSingleFileRead(
    const std::string& filename,
    const basis::OrbitalSpaceLJPN& orbital_space,
    const HalfInt& J_bra, const int& g_bra, const int& n_bra,
    const HalfInt& J_ket, const int& g_ket, const int& n_ket
  )
{
  std::cout << "Reading single-file OBDMEs" << std::endl;

  shell::InOBDMEStreamSingle stream(
      filename, orbital_space,
      J_bra, g_bra, n_bra,
      J_ket, g_ket, n_ket
    );

  std::cout << "  reading multipole 2" << std::endl;
  basis::OperatorBlocks<double> matrices;
  basis::OrbitalSectorsLJPN sectors;
  stream.GetMultipole(2, sectors, matrices);

  for (auto& matrix : matrices) {
    std::cout << mcutils::FormatMatrix(matrix, "16.8e") << std::endl << std::endl;
  }

}

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  // Read orbitals
  std::string orbital_filename("test/orbitals.dat");
  std::ifstream is(orbital_filename);
  std::vector<basis::OrbitalPNInfo> input_orbitals =
    basis::ParseOrbitalPNStream(is, true);
  basis::OrbitalSpaceLJPN space(input_orbitals);

  std::string info_filename("test/mfdn.rppobdme.info");
  std::string data_filename("test/mfdn.statrobdme.seq001.2J00.n01.2T00");
  TestMultiFileRead(data_filename, info_filename, space, 0, 0, 1, 0, 0, 1);

  orbital_filename = std::string("test/orbitals-single.dat");
  is = std::ifstream(orbital_filename);
  input_orbitals =
    basis::ParseOrbitalPNStream(is, true);
  std::string filename("test/mfdn.robdme");
  space = basis::OrbitalSpaceLJPN(input_orbitals);
  TestSingleFileRead(filename, space, HalfInt(3,2), 1, 1, HalfInt(3,2), 1, 1);

  // termination
  return 0;
}
