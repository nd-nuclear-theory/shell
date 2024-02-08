/****************************************************************
  me2j_io.h

  Defines I/O class and functions for MFDn me2j interaction file formats.

  Zhou Zhou
  University of Notre Dame
****************************************************************/

#ifndef ME2J_IO_H_
#define ME2J_IO_H_

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <memory>

#include "eigen3/Eigen/Core"

#include "basis/jjjttz_scheme.h"
#include "basis/operator.h"
#include "basis/jjjttz_operator.h"

namespace shell {

  enum class Me2jMode {kText,kBinary};
  // text/binary mode
  //
  // Note: Not currently implementing matrix mode (for interchange
  // with applied math group).

  // notational definitions for me2j file modes
  //
  // Use of these arrays requires conversion of the Me2jMode to int.

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
      const std::string filename
    );
  // Read me2j files and save as jjjttz format in memory (AS).

  void WriteMe2jFile(
      const basis::TwoBodySpaceJJJTTz& space,
      const basis::TwoBodySectorsJJJTTz& sectors,
      const basis::OperatorBlocks<double>& matrices,
      const std::string filename
    );
  // Write me2j files from jjjttz format in memory (AS).

  ////////////////////////////////////////////////////////////////
} // namespace


#endif
