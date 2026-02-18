/****************************************************************
  me2j_io.h

  Defines I/O class and functions for MFDn me2j interaction file formats.

  Zhou Zhou
  University of Notre Dame

  + Created by zz ~12/17/23.
  + 02/04/25 (mac):
    - Generalize me2j I/O routines to handle general (N1max,N2max) truncation.
    - Add diagnostic output.
    - Add precondition tests on operator labels.
    - Use basis::OrbitalSubspacePN for orbital indexing.
    - Streamline loop organization for readability.
    - Define parameters for binary file header length and float size (implementation pending).
 

****************************************************************/

#ifndef ME2J_IO_H_
#define ME2J_IO_H_

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <memory>

#include <Eigen/Core>

#include "basis/jjjttz_scheme.h"
#include "basis/operator.h"
#include "basis/jjjttz_operator.h"

namespace shell {

  enum class Me2jMode {kText,kBinary};
  // text/binary mode

  // notational definitions for me2j file modes
  //
  // Use of these arrays requires conversion of the Me2jMode to int.
  extern const std::array<const char*,2> kMe2jModeDescription; // ({"text","binary"});

  Me2jMode DeducedIOModeMe2j(const std::string& filename);
  // Deduce me2j file mode from filename extension.
  //
  //   .dat: ascii format
  //   .bin: binary format
  //
  // Arguments:
  //   filename (string) : filename from which to deduce mode
  //
  // Returns:
  //   (Me2jMode) : the mode

  void ReadMe2jFile(
      const basis::TwoBodySpaceJJJTTz& space,
      const basis::TwoBodySectorsJJJTTz& sectors,
      basis::OperatorBlocks<double>& matrices,
      const std::string filename,
      std::size_t float_size = 4
    );
  // Read me2j files and save as jjjttz format in memory (AS).

  void WriteMe2jFile(
      const basis::TwoBodySpaceJJJTTz& space,
      const basis::TwoBodySectorsJJJTTz& sectors,
      const basis::OperatorBlocks<double>& matrices,
      const std::string filename,
      std::size_t float_size = 4
    );
  // Write me2j files from jjjttz format in memory (AS).

  ////////////////////////////////////////////////////////////////
} // namespace


#endif
