/****************************************************************
  obme_io_test.cpp

  Patrick J. Fasano
  University of Notre Dame

****************************************************************/

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "basis/nlj_orbital.h"
#include "obme/obme_io.h"
#include "obme/radial.h"
#include "obme/obme.h"


////////////////////////////////////////////////////////////////
// test code
////////////////////////////////////////////////////////////////

basis::OperatorBlocks<double> TestRadial(const std::string& filename, int J0, int g0, int order, bool verbose = false) {
  std::cout << "Radial matrix elements test" << std::endl;

  // set up space
  std::cout << "Space" << std::endl;
  int Nmax = 12;
  basis::OrbitalSpaceLJPN space(Nmax);
  std::cout << space.DebugStr();

  // print sectors
  std::cout << "Sectors" << std::endl;
  basis::OrbitalSectorsLJPN sectors(space, space, J0, g0, 0);
  std::cout << "J0: "   << sectors.J0()
            << "g0: "   << sectors.g0()
            << " Tz0: " << sectors.Tz0()
            << std::endl;
  std::cout << sectors.DebugStr();

  // set up output stream
  std::cout << "Output stream" << std::endl;
  shell::OutOBMEStream os(
    filename, space, space, sectors,
    basis::OneBodyOperatorType::kRadial
  );

  // generate matrices
  basis::OperatorBlocks<double> matrices;
  shell::GenerateRadialOperator(
    shell::RadialBasisType::kOscillator,
    shell::RadialOperatorType::kR,
    order,
    space,
    sectors,
    matrices
  );

  // write to file
  std::cout << "Writing to file" << std::endl;
  os.Write(matrices);

  return matrices;
}

basis::OperatorBlocks<double> TestOperator(const std::string& filename, int J0, int g0, int order, bool verbose = false) {
  std::cout << "One-body operator matrix elements test" << std::endl;

  // set up space
  std::cout << "Space" << std::endl;
  int Nmax = 12;
  basis::OrbitalSpaceLJPN space(Nmax);
  std::cout << space.DebugStr();

  // print sectors
  std::cout << "Sectors" << std::endl;
  basis::OrbitalSectorsLJPN sectors(space, space, J0, g0, 0);
  std::cout << "J0: "   << sectors.J0()
            << "g0: "   << sectors.g0()
            << " Tz0: " << sectors.Tz0()
            << std::endl;
  std::cout << sectors.DebugStr();

  // set up output stream
  std::cout << "Output stream" << std::endl;
  shell::OutOBMEStream os(
    filename, space, space, sectors,
    basis::OneBodyOperatorType::kRadial
  );

  // generate matrices
  basis::OperatorBlocks<double> matrices;
  shell::SolidHarmonicOneBodyOperator(
    shell::RadialBasisType::kOscillator,
    shell::RadialOperatorType::kR,
    order,
    space,
    sectors,
    matrices
  );

  // write to file
  std::cout << "Writing to file" << std::endl;
  os.Write(matrices);

  return matrices;
}

basis::OperatorBlocks<double> TestProduct(const std::string& filename, bool verbose = false) {
  std::cout << "One-body operator matrix elements test" << std::endl;

  // set up space
  std::cout << "Space" << std::endl;
  int Nmax = 12;
  basis::OrbitalSpaceLJPN space(Nmax);
  std::cout << space.DebugStr();

  // print sectors
  std::cout << "Sectors" << std::endl;
  basis::OrbitalSectorsLJPN sectors(space, space, 1, 1, 0);
  std::cout << "J0: "   << sectors.J0()
            << "g0: "   << sectors.g0()
            << " Tz0: " << sectors.Tz0()
            << std::endl;
  std::cout << sectors.DebugStr();

  // generate matrices
  basis::OperatorBlocks<double> matrices;
  shell::SolidHarmonicOneBodyOperator(
    shell::RadialBasisType::kOscillator,
    shell::RadialOperatorType::kR,
    1,
    space,
    sectors,
    matrices
  );

  // print output sectors
  std::cout << "Sectors" << std::endl;
  basis::OrbitalSectorsLJPN output_sectors(space, space, 0, 0, 0);
  std::cout << "J0: "   << output_sectors.J0()
            << "g0: "   << output_sectors.g0()
            << " Tz0: " << output_sectors.Tz0()
            << std::endl;
  std::cout << output_sectors.DebugStr();

  // generate matrices
  basis::OperatorBlocks<double> output_matrices;
  shell::OneBodyOperatorTensorProduct(space, sectors, matrices, sectors, matrices, output_sectors, output_matrices);

  // set up output stream
  std::cout << "Output stream" << std::endl;
  shell::OutOBMEStream os(
    filename, space, space, output_sectors,
    basis::OneBodyOperatorType::kRadial
  );

  // write to file
  std::cout << "Writing to file" << std::endl;
  os.Write(output_matrices);

  return output_matrices;
}


int main(int argc, char **argv) {
  std::string radial_filename("test/radial_r1_test.dat");
  TestRadial(radial_filename, 1, 1, 1, false);
  radial_filename = ("test/radial_r2_test.dat");
  TestRadial(radial_filename, 0, 0, 2, false);
  std::string operator_filename("test/obme_r1_test.dat");
  TestOperator(operator_filename, 1, 1, 1, false);
  operator_filename = "test/obme_r2_test.dat";
  TestOperator(operator_filename, 0, 0, 2, false);
  std::string product_filename("test/obme_r.r_test.dat");
  TestProduct(product_filename, false);

  // termination
  return 0;
}
