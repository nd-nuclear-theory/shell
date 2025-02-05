/****************************************************************
  me2j_io.cpp

  Zhou Zhou
  University of Notre Dame

****************************************************************/

////////////////////////////////////////////////////////////////
// two-body matrix elements ordering: translate between me2j and jjjttz orderings
////////////////////////////////////////////////////////////////
// ABCD (jjjt/jjjttz)
// NA+NB<NC+ND
// and if NA+NB==NC+ND
// A<=B
// C<=D
// A<=C
// if A==C, B<=D
// ABCD are single particle indexes for labels of each single particle state
//
// ABCD (me2j)
// B<=A
// C<=A
// if C==A, D<=B
// else D<=C
//
// To match me2j into jjjttz ordering, the shortcut is to do:
// if ND+NC<NB+NA or (ND+NC==NB+NA and D<=B), return DCBA
// else, return BADC
////////////////////////////////////////////////////////////////

// Notes on known flavors of binary me2j file:
//
//   - TUD (me2j-f2): header consists of 255 bytes; matrix elements are float (based on menj 2.1.0)
//
//   - PN: header consists of 40 bytes; matrix elements are double (based on example from Livermore)
//     + 02/04/25 (mac): should confirm whether or not there might be FORTRAN record delimiters
//
//   - miyagi: no header; matrix elements can be either float or double

#include "tbme/me2j_io.h"

#include <cstddef>
#include <cstring>
#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <memory>

#include "am/halfint_fmt.h"  // for diagnostics
#include "basis/nlj_orbital.h"
#include "fmt/format.h"
#include "mcutils/io.h"
#include "mcutils/parsing.h"

namespace shell {

  ////////////////////////////////////////////////////////////////
  // file text/binary I/O mode identification
  ////////////////////////////////////////////////////////////////

  const std::array<const char*,2> kMe2jModeDescription({"text","binary"});
  
  Me2jMode DeducedIOModeMe2j(const std::string& filename)
  {
    if (filename.length() < 3 )
      {
        // prevent compare on underlength string
        std::cerr << "ERROR: Me2j file I/O: No extension found (too short) in filename " << filename << std::endl;
        exit(EXIT_FAILURE);
      }
    else if ( ! filename.compare(filename.length()-3,3,"bin") )
      return Me2jMode::kBinary;
    else
      return Me2jMode::kText;
  }

  void ReadMe2jFile(
      const basis::TwoBodySpaceJJJTTz& space,
      const basis::TwoBodySectorsJJJTTz& sectors,
      basis::OperatorBlocks<double>& matrices,
      const std::string filename
    )
  {
    // binary file parameters
    const std::size_t header_length = 255;
    const std::size_t float_size = 4;
    assert((float_size == 4) || (float_size == 8));

    // validate operator labels
    // TODO (mac): instead, provide initialization of sectors?
    if (sectors.J0()!=0 || sectors.g0()!=0 || sectors.Tz0()!=0) {
      std::cerr << "ERROR: Provided operator has unsupported (J0,g0,Tz0)!=(0,0,0)." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    
    // extract space truncation
    int N1max = space.N1max();
    int N2max = space.N2max();

    // choose file format (text or binary)
    Me2jMode me2j_mode = DeducedIOModeMe2j(filename);
    std::ios_base::openmode mode_argument;
    if (me2j_mode == Me2jMode::kText) {
      mode_argument = std::ios_base::in;
    } else {
      mode_argument = (std::ios_base::in | std::ios_base::binary);
    }

    // write diagnostics
    std::cout << fmt::format("  File: {}", filename) << std::endl;
    std::cout << fmt::format("  Format: {}", kMe2jModeDescription[int(me2j_mode)]) << std::endl;
    std::cout << fmt::format("  Truncation: N1max {} N2max {}", N1max, N2max) << std::endl;
    
    // open input file
    std::ifstream is(filename.c_str(), mode_argument);

    // skip file header
    if (me2j_mode == Me2jMode::kText) { // only text files have a header line
      is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    } else {
      // presently hard coded for TUD header
      assert(header_length==255);
      char header[header_length];
      mcutils::ReadBinary<char>(is, header, header_length);
      if(std::strstr(header,"me2j-f2-bin")==NULL) {
           std::cerr << "ERROR: unrecognized me2j file header" << std::endl;
           std::exit(EXIT_FAILURE);
         }
    }

    // set up orbital indexing
    const auto orbitals = basis::OrbitalSubspacePN(basis::OrbitalSpeciesPN::kP, N1max);  // set species label to "proton" arbitrarily as dummy
    int num_orbitals = orbitals.size();
    
    // iterate over matrix elements to read
    int tbme_count=0;
    for (int a = 0; a < num_orbitals; a++) {
      for (int b = 0; b <= a; b++) {
        // iterate over bra orbitals

        // extract bra orbital labels
        const auto orbital_a = orbitals.GetState(a);
        int Na = orbital_a.N();
        int ga = orbital_a.g();
        HalfInt ja = orbital_a.j();
        const auto orbital_b = orbitals.GetState(b);
        int Nb = orbital_b.N();
        int gb = orbital_b.g();
        HalfInt jb = orbital_b.j();

        // apply truncation condition on bra
        if (Na+Nb > N2max)
          continue;
        
        for (int c = 0; c <= a; c++) {
          int dmax = (a == c) ? b : c;
          for (int d = 0; d <= dmax; d++) {
            // iterate over ket orbitals

            // extract ket orbital labels
            const auto orbital_c = orbitals.GetState(c);
            int Nc = orbital_c.N();
            int gc = orbital_c.g();
            HalfInt jc = orbital_c.j();
            const auto orbital_d = orbitals.GetState(d);
            int Nd = orbital_d.N();
            int gd = orbital_d.g();
            HalfInt jd = orbital_d.j();

            // apply truncation condition on ket
            if (Nc+Nd > N2max)
              continue;

            // apply parity selection rule (for g0=0 operator)
            int gab = (ga+gb)%2;
            int gcd = (gc+gd)%2;
            if (gab!=gcd)
              continue;

            // apply angular momentum coupling constraint on orbitals
            int Jmin = std::max(std::abs(int(ja-jb)),std::abs(int(jc-jd)));
            int Jmax = std::min(int(ja+jb),int(jc+jd));
            if (Jmin > Jmax)
              continue;
            
            for (int J = Jmin; J <= Jmax; J++) {
              for (int T = 0; T <= 1; T++) {
                for (int Tz = -T; Tz <= T; Tz++) {
                  // iterate over two-body state labels

                  // read matrix element
                  tbme_count++;
                  double matrix_element;
                  if (me2j_mode == Me2jMode::kText) {
                    is >> matrix_element;
                  } else {
                    if (float_size == 4) {
                      float temp_matrix_element;
                      mcutils::ReadBinary<float>(is, temp_matrix_element);
                      matrix_element = double(temp_matrix_element);
                    } else {
                      mcutils::ReadBinary<double>(is, matrix_element);
                    }
                  }

                  // store matrix element (with canonicalization)
                  if (int((ja+jb+jc+jd))%2==1) {
                    matrix_element *= -1;
                  }
                  if ((Nc+Nd)<(Na+Nb) || ((Nc+Nd)==(Na+Nb) && d <= b)) {
                    // to make the indexes ordered the same way as jjjttz_operator
                    basis::SetTwoBodyOperatorMatrixElementJJJTTz( // save as dcba in matrices
                      space,
                      basis::TwoBodySubspaceJJJTTz::SubspaceLabelsType(J,T,gcd,Tz),
                      basis::TwoBodySubspaceJJJTTz::SubspaceLabelsType(J,T,gab,Tz), // note that gcd == gab
                      basis::TwoBodyStateJJJTTz::StateLabelsType(Nd,jd,Nc,jc),
                      basis::TwoBodyStateJJJTTz::StateLabelsType(Nb,jb,Na,ja),
                      sectors,
                      matrices,
                      matrix_element
                      // 0
                    );
                  } else {
                    basis::SetTwoBodyOperatorMatrixElementJJJTTz( // save as dcba in matrices
                      space,
                      basis::TwoBodySubspaceJJJTTz::SubspaceLabelsType(J,T,gab,Tz),
                      basis::TwoBodySubspaceJJJTTz::SubspaceLabelsType(J,T,gcd,Tz), // note that gcd == gab
                      basis::TwoBodyStateJJJTTz::StateLabelsType(Nb,jb,Na,ja),
                      basis::TwoBodyStateJJJTTz::StateLabelsType(Nd,jd,Nc,jc),
                      sectors,
                      matrices,
                      matrix_element
                      // test
                    );
                  }
                }
              }
            }
          }
        }
      }
    }

    // write diagnostics
    std::cout << fmt::format("  Matrix elements: {}", tbme_count) << std::endl;

    // check for unexpected file length (at least in text mode)
    if (me2j_mode == Me2jMode::kText) {
      double dummy;
      is >> dummy;
      if (is) {
        std::cerr << "ERROR: more matrix elements available in file than expected" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    
  }

  void WriteMe2jFile(
      const basis::TwoBodySpaceJJJTTz& space,
      const basis::TwoBodySectorsJJJTTz& sectors,
      const basis::OperatorBlocks<double>& matrices,
      const std::string filename
    )
  {
    // binary file parameters
    const std::size_t header_length = 255;
    const std::size_t float_size = 4;
    assert((float_size == 4) || (float_size == 8));

    // validate operator labels
    if (sectors.J0()!=0 || sectors.g0()!=0 || sectors.Tz0()!=0) {
      std::cerr << "ERROR: Provided operator has unsupported (J0,g0,Tz0)!=(0,0,0)." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    
    // extract space truncation
    int N1max = space.N1max();
    int N2max = space.N2max();
    
    // choose file format (text or binary)
    Me2jMode me2j_mode = DeducedIOModeMe2j(filename);
    std::ios_base::openmode mode_argument;
    if (me2j_mode == Me2jMode::kText) {
      mode_argument = std::ios_base::out;
    } else {
      mode_argument = (std::ios_base::out | std::ios_base::binary);
    }

    // write diagnostics
    std::cout << fmt::format("  File: {}", filename) << std::endl;
    std::cout << fmt::format("  Format: {}", kMe2jModeDescription[int(me2j_mode)]) << std::endl;
    std::cout << fmt::format("  Truncation: N1max {} N2max {}", N1max, N2max) << std::endl;
    
    // open output file
    std::ofstream os(filename.c_str(), mode_argument);
    if (me2j_mode == Me2jMode::kText) {
      os.precision(7);
      os << std::fixed;
      os << std::setw(12);
    }
    
    // write file header
    if (me2j_mode == Me2jMode::kText) { // only text files have a header line
      os << "(*** written by shell (https://github.com/nd-nuclear-theory/shell) ***)" << std::endl;
    } else {
      // presently hard coded for TUD header
      // 02/04/24 (mac): can implement more flexibly using std::string and c_str()
      assert(header_length==255);
      char header[header_length]="me2j-f2-bin";
      memset(&header[sizeof("me2j-f2-bin")], '\0', header_length-sizeof("me2j-f2-bin"));
      header[header_length-1] = '\0';
      mcutils::WriteBinary<char>(os, header, header_length);
    }

    // set up orbital indexing
    const auto orbitals = basis::OrbitalSubspacePN(basis::OrbitalSpeciesPN::kP, N1max);  // set species label to "proton" arbitrarily as dummy
    int num_orbitals = orbitals.size();
    
    // iterate over matrix elements to write
    int tbme_count=0;
    for (int a = 0; a < num_orbitals; a++) {
      for (int b = 0; b <= a; b++) {
        // iterate over bra orbitals

        // extract bra orbital labels
        const auto orbital_a = orbitals.GetState(a);
        int Na = orbital_a.N();
        int ga = orbital_a.g();
        HalfInt ja = orbital_a.j();
        const auto orbital_b = orbitals.GetState(b);
        int Nb = orbital_b.N();
        int gb = orbital_b.g();
        HalfInt jb = orbital_b.j();

        // apply truncation condition on bra
        if (Na+Nb > N2max)
          continue;
        
        for (int c = 0; c <= a; c++) {
          int dmax = (a == c) ? b : c;
          for (int d = 0; d <= dmax; d++) {
            // iterate over ket orbitals

            // extract ket orbital labels
            const auto orbital_c = orbitals.GetState(c);
            int Nc = orbital_c.N();
            int gc = orbital_c.g();
            HalfInt jc = orbital_c.j();
            const auto orbital_d = orbitals.GetState(d);
            int Nd = orbital_d.N();
            int gd = orbital_d.g();
            HalfInt jd = orbital_d.j();
            
            // apply truncation condition on ket
            if (Nc+Nd > N2max)
              continue;

            // apply parity selection rule (for g0=0 operator)
            int gab = (ga+gb)%2;
            int gcd = (gc+gd)%2;
            if (gab!=gcd)
              continue;

            // apply angular momentum coupling constraint on orbitals
            int Jmin = std::max(std::abs(int(ja-jb)),std::abs(int(jc-jd)));
            int Jmax = std::min(int(ja+jb),int(jc+jd));
            if (Jmin > Jmax)
              continue;

            for (int J = Jmin; J <= Jmax; J++) {
              for (int T = 0; T <= 1; T++) {
                for (int Tz = -T; Tz <= T; Tz++) {
                  // iterate over two-body state labels

                  // retrieve matrix element (with canonicalization)
                  double matrix_element;
                  if ((Nc+Nd)<(Na+Nb) || ((Nc+Nd)==(Na+Nb) && d <= b)) {
                    // std::cout
                    //   << fmt::format(
                    //     "A: {} {} {} {}   {} {} {} {} {}",
                    //     J,T,gcd,Tz,
                    //     J,T,gab,Tz,
                    //     basis::TwoBodySubspaceJJJTTz::SubspaceLabelsType(J,T,gcd,Tz)<=basis::TwoBodySubspaceJJJTTz::SubspaceLabelsType(J,T,gab,Tz)
                    //     )
                    //   <<std::endl;
                    matrix_element = basis::GetTwoBodyOperatorMatrixElementJJJTTz(
                      space,
                      basis::TwoBodySubspaceJJJTTz::SubspaceLabelsType(J,T,gcd,Tz),
                      basis::TwoBodySubspaceJJJTTz::SubspaceLabelsType(J,T,gab,Tz), // note that gcd == gab
                      basis::TwoBodyStateJJJTTz::StateLabelsType(Nd,jd,Nc,jc),
                      basis::TwoBodyStateJJJTTz::StateLabelsType(Nb,jb,Na,ja),
                      sectors,
                      matrices
                    );
                  } else {
                    // std::cout
                    //   << fmt::format(
                    //     "B: {} {} {} {}   {} {} {} {} {}",
                    //     J,T,gab,Tz,
                    //     J,T,gcd,Tz,
                    //     basis::TwoBodySubspaceJJJTTz::SubspaceLabelsType(J,T,gab,Tz)<=basis::TwoBodySubspaceJJJTTz::SubspaceLabelsType(J,T,gcd,Tz)
                    //     )
                    //   <<std::endl;
                    matrix_element = basis::GetTwoBodyOperatorMatrixElementJJJTTz(
                      space,
                      basis::TwoBodySubspaceJJJTTz::SubspaceLabelsType(J,T,gab,Tz),
                      basis::TwoBodySubspaceJJJTTz::SubspaceLabelsType(J,T,gcd,Tz), // note that gcd == gab
                      basis::TwoBodyStateJJJTTz::StateLabelsType(Nb,jb,Na,ja),
                      basis::TwoBodyStateJJJTTz::StateLabelsType(Nd,jd,Nc,jc),
                      sectors,
                      matrices
                    );
                  }
                  if ((int(ja+jb+jc+jd))%2==1) {
                    matrix_element *= -1;
                  }

                  // write matrix element
                  tbme_count++;
                  if (me2j_mode == Me2jMode::kText) {
                    os << " " << std::setw(12) << matrix_element;
                    if (tbme_count%10==0) {
                      os << std::endl;
                    }
                  } else {
                    if (float_size == 4) {
                      mcutils::WriteBinary<float>(os, float(matrix_element));
                    } else {
                      mcutils::WriteBinary<double>(os, matrix_element);
                    }
                  }
                  
                }
              }
            }
          }
        }
      }
    }

    // write diagnostics
    std::cout << fmt::format("  Matrix elements: {}", tbme_count) << std::endl;

  }

  ////////////////////////////////////////////////////////////////
} // namespace
