/****************************************************************
  me2j_io.h

  Defines I/O class for MFDn me2j interaction file formats.

  Normalization convention: All matrix elements are stored as NAS
  RMEs.  These RMEs are stored under the group theory Wigner-Eckart
  normalization convention (i.e., "no dimension factor out front, just
  the Clebsch"), but, for scalar operators, note that this RME is
  equivalently, and more simply, the branched ME (with M'=M).

  Symmetrization convention: The full square matrix is *not* populated
  on diagonal sectors.  For these sectors, the lower triangle is
  zero-initialized on input and ignored on output.

  Mark A. Caprio
  University of Notre Dame

  + 08/31/12 (mac): Adapted from mfdn_io.h (as mfdn_me2j).
     - Converted I/O to class, from functions acting on stream
       argument.
     - Added support for binary I/O and format code.
  + 01/28/13 (mac): Complete implementation.
  + 07/25/14 (mac): Add matrix (".mat") format for Me2j files.
  + 04/25/15 (mac): Reformat source file.
  + 10/11/16 (mac,pjf):
    - Rename to me2j_io.
    - Integrate into shell project build.
  + 10/19/16 (mac): Complete implementation for Me2j Version0.
  + 10/25/16 (mac): Add InMe2jStream::SeekToSector.
  + 11/01/16 (mac):
    - Convert from AS to NAS storage.
  + 11/13/16 (mac): Implement Me2j Version15099 binary output.
  + 11/28/16 (mac): Add Tz0 to me2j format 15099 output header.
  + 10/19/17 (mac): Add optional on-the-fly conversion from AS to
    NAS matrix elements on output.
  + 01/22/18 (mac): Begin implementing nonzero Tz0.
  + 02/12/19 (pjf): Finish implementing nonzero Tz0.
  + 02/21/19 (pjf):
    - Remove Tz0!=0 support from me2jv15099.
    - Implement Me2j Version15200.
  + 05/09/19 (pjf): Use std::size_t for basis indices and sizes.
  + 05/10/19 (pjf): Ensure Fortran records aren't larger than allowed.
  + 05/30/19 (pjf): Use long int/integer(kind=8) for number of matrix elements
    in header of binary file.
  + 06/03/19 (pjf): Implement version 15200 binary I/O.
  + 08/28/19 (pjf): Correctly extract implicit one-body truncation from
    orbitals listed in version 15200 file.
  + 09/06/19 (pjf): Ensure (for version 15200) that orbital space weight maxes
    match two-body space truncation.
  + 10/10/20 (pjf): Dramatically improve binary I/O performance by using
    buffered reads and writes (sector-at-a-time I/O).
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
  extern const std::array<const char*,2> kMe2jModeDescription; // ({"text","binary"});
  extern const std::array<const char*,2> kMe2jModeExtension; // ({"dat","bin","int"}); or "bin" vs not "bin"

  Me2jMode DeducedIOMode(const std::string& filename);
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

  void WriteMe2jFile(
      const basis::TwoBodySpaceJJJTTz& space,
      const basis::TwoBodySectorsJJJTTz& sectors,
      const basis::OperatorBlocks<double>& matrices,
      const std::string filename
    );
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace


#endif
