/****************************************************************
  me2j_io.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "tbme/me2j_io.h"

#include <cstddef>
#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <memory>

#include "fmt/format.h"
#include "mcutils/io.h"
#include "mcutils/parsing.h"

namespace shell {

  ////////////////////////////////////////////////////////////////
  // file text/binary I/O mode identification
  ////////////////////////////////////////////////////////////////


  // const std::array<const char*,2> kMe2jModeDescription({"text","binary"});
  // const std::array<const char*,2> kMe2jModeExtension({".dat",".bin"});

  Me2jMode DeducedIOModeMe2j(const std::string& filename)
  {
    if (filename.length() < 3 )
      {
        // prevent compare on underlength string
        std::cerr << "Me2j file I/O: No extension found (too short) in filename " << filename << std::endl;
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
    // check space truncation
    if (space.N1max()!=space.N2max()) {
      std::cout << "Space is not constructed with a 2-body Nmax truncation." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    int Nmax = space.N2max();

    // choose file format (text or binary)
    Me2jMode me2j_mode = DeducedIOModeMe2j(filename);
    std::ios_base::openmode mode_argument;
    if (me2j_mode == Me2jMode::kText) {
      mode_argument = std::ios_base::in;
    } else {
      mode_argument = (std::ios_base::in | std::ios_base::binary);
    }
    std::ifstream is(filename.c_str(), mode_argument);
    // skip header
    if (me2j_mode == Me2jMode::kText) { // only text files have a header line
      is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    } else {
      double temp;
      for (size_t i = 0; i < 5; i++) {
        mcutils::ReadBinary<double>(is, temp);
      }
    }

    // find array size for given Nmax
    int array_size = 0;
    for (int N = 0; N <= Nmax; N++) {
      for (int l = 0; l <= N; l++) {
        if ((N-l)%2 == 1) {
          continue;
        }
        int n = (N-l)/2;
        for (HalfInt s = HalfInt(-1,2); s <= HalfInt(1,2); s++) {
          HalfInt j=HalfInt(l)+s;
          if (j < HalfInt(1,2)) {
            continue;
          }
          array_size+=1;
        }
      }
    }

    // save arrays for (N,l,j) indexing
    int N_array[array_size];
    int l_array[array_size];
    HalfInt j_array[array_size];
    int array_index = 0;
    for (int N = 0; N <= Nmax; N++) {
      for (int l = 0; l <= N; l++) {
        if ((N-l)%2 == 1) {
          continue;
        }
        int n = (N-l)/2;
        for (HalfInt s = HalfInt(-1,2); s <= HalfInt(1,2); s++) {
          HalfInt j = HalfInt(l)+s;
          if (j < HalfInt(1,2)) {
            continue;
          }
          N_array[array_index] = N;
          l_array[array_index] = l;
          j_array[array_index] = j;
          // std::cout << "N,l,j" << N << " " << l << " " << j << std::endl;
          array_index += 1;
        }
      }
    }

    // enumerate and go through input matrix elements, then save result to matrices
    double matrix_element=0;
    int count=0;
    for (int a = 0; a < array_size; a++) {
      int Na = N_array[a];
      int la = l_array[a];
      HalfInt ja = j_array[a];
      for (int b = 0; b <= a; b++) {
        int Nb = N_array[b];
        if (Na+Nb > Nmax) {
          continue;
        }
        int lb = l_array[b];
        HalfInt jb = j_array[b];
        for (int c = 0; c <= a; c++) {
          int Nc = N_array[c];
          int lc = l_array[c];
          HalfInt jc = j_array[c];
          int dmax = c;
          if (a == c) {
            dmax = b;
          }
          for (int d = 0; d <= dmax; d++) {
            int Nd = N_array[d];
            if (Nc+Nd > Nmax) {
              continue;
            }
            int ld = l_array[d];
            HalfInt jd = j_array[d];
            int gab = (la+lb)%2;
            int gcd = (lc+ld)%2;
            if (gab!=gcd) { // parity selection rule for g0=0 operators
              continue;
            }
            int Jmin = std::max(std::abs(int(ja-jb)),std::abs(int(jc-jd)));
            int Jmax = std::min(int(ja+jb),int(jc+jd));
            if (Jmin > Jmax) {
              continue;
            }
            for (int J = Jmin; J <= Jmax; J++) {
              for (int T = 0; T <= 1; T++) {
                for (int Tz = -T; Tz <= T; Tz++) {
                  count++;
                  if (me2j_mode == Me2jMode::kText) {
                    is >> matrix_element;
                  } else {
                    float temp_matrix_element;
                    // mcutils::ReadBinary<float>(is, temp_matrix_element);
                    // matrix_element = double(temp_matrix_element);
                    mcutils::ReadBinary<double>(is, matrix_element);
                    // std::cout << matrix_element << std::endl;
                  }
                  // std::cout << matrix_element << std::endl;
                  // if (!is) {
                  //   std::cout << "reading more than there are in the file" << std::endl;
                  //   exit(EXIT_FAILURE);
                  // }
                  if (int((ja+jb+jc+jd))%2==1) {
                    matrix_element *= -1;
                  }
                  if (d <= b) {
                    basis::SetTwoBodyOperatorMatrixElementJJJTTz( // save as dcba in matrices
                      space,
                      basis::TwoBodySubspaceJJJTTz::SubspaceLabelsType(J,T,gcd,Tz),
                      basis::TwoBodySubspaceJJJTTz::SubspaceLabelsType(J,T,gab,Tz), // note that gcd == gab
                      basis::TwoBodyStateJJJTTz::StateLabelsType(Nd,jd,Nc,jc),
                      basis::TwoBodyStateJJJTTz::StateLabelsType(Nb,jb,Na,ja),
                      sectors,
                      matrices,
                      matrix_element
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
                    );
                  }
                }
              }
            }
          }
        }
      }
    }
    std::cout << count << std::endl;
    is >> matrix_element;
    if (is) {
      std::cout << "reading less than there are in the file" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  void WriteMe2jFile(
      const basis::TwoBodySpaceJJJTTz& space,
      const basis::TwoBodySectorsJJJTTz& sectors,
      const basis::OperatorBlocks<double>& matrices,
      const std::string filename
    )
  {
    // check space truncation
    if (space.N1max()!=space.N2max()) {
      std::cout << "Space is not constructed with a 2-body Nmax truncation." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    int Nmax = space.N2max();

    // choose file format (text or binary)
    Me2jMode me2j_mode = DeducedIOModeMe2j(filename);
    std::ios_base::openmode mode_argument;
    if (me2j_mode == Me2jMode::kText) {
      mode_argument = std::ios_base::out;
    } else {
      mode_argument = (std::ios_base::out | std::ios_base::binary);
    }
    std::ofstream os(filename.c_str(), mode_argument);
    os.precision(7);
    os << std::fixed;
    os << std::setw(12);
    if (me2j_mode == Me2jMode::kText) { // only text files have a header line
      os << "written by me2j_io(basis)" << std::endl;
    }

    // find array size for given Nmax
    int array_size = 0;
    for (int N = 0; N <= Nmax; N++) {
      for (int l = 0; l <= N; l++) {
        if ((N-l)%2 == 1) {
          continue;
        }
        int n = (N-l)/2;
        for (HalfInt s = HalfInt(-1,2); s <= HalfInt(1,2); s++) {
          HalfInt j=HalfInt(l)+s;
          if (j < HalfInt(1,2)) {
            continue;
          }
          array_size+=1;
        }
      }
    }

    // save arrays for (N,l,j) indexing
    int N_array[array_size];
    int l_array[array_size];
    HalfInt j_array[array_size];
    int array_index = 0;
    for (int N = 0; N <= Nmax; N++) {
      for (int l = 0; l <= N; l++) {
        if ((N-l)%2 == 1) {
          continue;
        }
        int n = (N-l)/2;
        for (HalfInt s = HalfInt(-1,2); s <= HalfInt(1,2); s++) {
          HalfInt j = HalfInt(l)+s;
          if (j < HalfInt(1,2)) {
            continue;
          }
          N_array[array_index] = N;
          l_array[array_index] = l;
          j_array[array_index] = j;
          array_index += 1;
        }
      }
    }


    // enumerate and go through all matrix elements related to me2j, then write to the file
    int count=0;
    double matrix_element;
    for (int a = 0; a < array_size; a++) {
      int Na = N_array[a];
      int la = l_array[a];
      HalfInt ja = j_array[a];
      for (int b = 0; b <= a; b++) {
        int Nb = N_array[b];
        if (Na+Nb > Nmax) {
          continue;
        }
        int lb = l_array[b];
        HalfInt jb = j_array[b];
        for (int c = 0; c <= a; c++) {
          int Nc = N_array[c];
          int lc = l_array[c];
          HalfInt jc = j_array[c];
          int dmax = c;
          if (a == c) {
            dmax = b;
          }
          for (int d = 0; d <= dmax; d++) {
            int Nd = N_array[d];
            if (Nc+Nd > Nmax) {
              continue;
            }
            int ld = l_array[d];
            HalfInt jd = j_array[d];
            int gab = (la+lb)%2;
            int gcd = (lc+ld)%2;
            if (gab!=gcd) { // parity selection rule for g0=0 operators
              continue;
            }
            int Jmin = std::max(std::abs(int(ja-jb)),std::abs(int(jc-jd)));
            int Jmax = std::min(int(ja+jb),int(jc+jd));
            if (Jmin > Jmax) {
              continue;
            }
            for (int J = Jmin; J <= Jmax; J++) {
              for (int T = 0; T <= 1; T++) {
                for (int Tz = -T; Tz <= T; Tz++) {
                  count++;
                  if (d <= b) {
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
                  os << matrix_element << " ";
                  if (count%10==0) {
                    os << std::endl;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
} // namespace
