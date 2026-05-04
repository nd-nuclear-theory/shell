/****************************************************************
  obme_io_test.cpp

  Patrick J. Fasano
  University of Notre Dame

****************************************************************/

#include <cassert>
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
            << " g0: "   << sectors.g0()
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

basis::OperatorBlocks<double> TestProduct(const std::string& filename, bool verbose = false)
// Generate r.r as product operator.
{
  std::cout << "One-body operator matrix elements test (r.r)" << std::endl;

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


void TestProductNonscalar(bool verbose = false)
// Generate r^2 Y2 as product operator.
//
// Spot check:
//
// <0s_1/2 || Q_2 || 0d_5/2>_Edmonds = 1.338 by Suhonen Table 6.4
//
// <0s_1/2 || Q_2 || 0d_5/2>_Rose = 1/sqrt(2)*<0s_1/2 || Q_2 || 0d_5/2>_Edmonds
//   = 0.9461
//
// obme_r2xY2_test.dat
//
// # bra subspace labels: l = 0 2j = 1 2Tz = 1
// # ket subspace labels: l = 2 2j = 5 2Tz = 1
//   9.46174696e-01   0.00000000e+00
//  -1.54509681e+00   1.44530572e+00
//   6.90988299e-01  -2.58544147e+00
//
// obme_r2Y2_test.dat
//
// # bra subspace labels: l = 0 2j = 1 2Tz = 1
// # ket subspace labels: l = 2 2j = 5 2Tz = 1
//   9.46174696e-01   0.00000000e+00
//  -1.54509681e+00   1.44530572e+00
//   6.90988299e-01  -2.58544147e+00
  
{
  std::cout << "One-body operator matrix elements test (r^2 Y2)" << std::endl;

  // set up space
  std::cout << "Space" << std::endl;
  int Nmax = 4;
  basis::OrbitalSpaceLJPN space(Nmax);
  std::cout << space.DebugStr();

  // generate r^2 operator
  std::cout << "Operator A (r^2)" << std::endl;
  basis::OrbitalSectorsLJPN sectors_a(space, space, 0, 0, 0);
  std::cout << " J0: " << sectors_a.J0()
            << " g0: " << sectors_a.g0()
            << " Tz0: " << sectors_a.Tz0()
            << std::endl;
  std::cout << sectors_a.DebugStr();
  basis::OperatorBlocks<double> matrices_a;
  shell::SolidHarmonicOneBodyOperator(
    shell::RadialBasisType::kOscillator,
    shell::RadialOperatorType::kR,
    2,
    space,
    sectors_a,
    matrices_a
  );
  // write to file
  std::cout << "Writing to file" << std::endl;
  {
    shell::OutOBMEStream os(
        "test/obme_r2xY2_r2-factor_test.dat", space, space, sectors_a,
        basis::OneBodyOperatorType::kRadial
      );
    os.Write(matrices_a);
  }

  
  // generate Y2 operator
  std::cout << "Operator A (Y2)" << std::endl;
  basis::OrbitalSectorsLJPN sectors_b(space, space, 2, 0, 0);
  std::cout << " J0: " << sectors_b.J0()
            << " g0: " << sectors_b.g0()
            << " Tz0: " << sectors_b.Tz0()
            << std::endl;
  std::cout << sectors_b.DebugStr();
  basis::OperatorBlocks<double> matrices_b;
  shell::SolidHarmonicOneBodyOperator(
    shell::RadialBasisType::kOscillator,
    shell::RadialOperatorType::kR,
    0,
    space,
    sectors_b,
    matrices_b
  );
  basis::ScalarMultiplyOperator(sectors_b, matrices_b, std::sqrt(5/(4*M_PI)));
  {
    shell::OutOBMEStream os(
        "test/obme_r2xY2_Y2-factor_test.dat", space, space, sectors_b,
        basis::OneBodyOperatorType::kRadial
      );
    os.Write(matrices_b);
  }
    
  // print output sectors
  std::cout << "Sectors" << std::endl;
  basis::OrbitalSectorsLJPN output_sectors(space, space, 2, 0, 0);
  std::cout << "J0: "   << output_sectors.J0()
            << "g0: "   << output_sectors.g0()
            << " Tz0: " << output_sectors.Tz0()
            << std::endl;
  std::cout << output_sectors.DebugStr();

  // generate matrices
  basis::OperatorBlocks<double> output_matrices;
  shell::OneBodyOperatorTensorProduct(space, sectors_a, matrices_a, sectors_b, matrices_b, output_sectors, output_matrices);
  
  // write to file
  std::cout << "Writing to file" << std::endl;
  {
    shell::OutOBMEStream os(
        "test/obme_r2xY2_test.dat", space, space, output_sectors,
        basis::OneBodyOperatorType::kRadial
      );
    os.Write(output_matrices);
  }

  // generate r^2*Y2 directly for comparison
  std::cout << "Comparison operator (r^2 Y2)" << std::endl;
  basis::OrbitalSectorsLJPN sectors_comparison(space, space, 2, 0, 0);
  std::cout << " J0: " << sectors_comparison.J0()
            << " g0: " << sectors_comparison.g0()
            << " Tz0: " << sectors_comparison.Tz0()
            << std::endl;
  // std::cout << sectors_comparison.DebugStr();
  basis::OperatorBlocks<double> matrices_comparison;
  shell::SolidHarmonicOneBodyOperator(
    shell::RadialBasisType::kOscillator,
    shell::RadialOperatorType::kR,
    2,
    space,
    sectors_comparison,
    matrices_comparison
  );
  basis::ScalarMultiplyOperator(sectors_comparison, matrices_comparison, std::sqrt(5/(4*M_PI)));

  // write to file
  std::cout << "Writing to file" << std::endl;
  {
    shell::OutOBMEStream os(
        "test/obme_r2Y2_test.dat", space, space, sectors_comparison,
        basis::OneBodyOperatorType::kRadial
      );
    os.Write(matrices_comparison);
  }
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
  TestProductNonscalar(false);

  // termination
  return 0;
}
