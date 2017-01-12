/************************************************************//**
  @file obdme_io.h

  Defines I/O classes for one-body density matrix element access.

  Language: C++11

  Patrick J. Fasano
  University of Notre Dame

  + 11/22/16 (pjf): Created, based on radial_io.
  + 1/11/17 (pjf): Implemented MFDn v15 OBDME input.
    - InOBDMEReader parses info files at construction time.
    - The same reader can parse multiple data files.
    - Quantum numbers stored in data files are currently ignored.

****************************************************************/

#ifndef OBDME_IO_H_
#define OBDME_IO_H_

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include "eigen3/Eigen/Core"

#include "basis/operator.h"
#include "basis/nlj_orbital.h"

namespace shell {

/**
 * Class for reading one-body density matrix elements
 */
class InOBDMEReader {
 public:
  /**
   * Default constructor -- provided since required for certain
   * purposes by STL container classes (e.g., std::vector::resize)
   */
  InOBDMEReader() = default;

  InOBDMEReader(
      const std::string& info_filename,
      const basis::OrbitalSpaceLJPN& orbital_space
    );
  // Construct a reader by parsing an info file.

  // destructor
  ~InOBDMEReader() {
    delete info_stream_ptr_;
  }

  // Get sectors associated with a particular multipole order
  void SetToIndexing(int order, basis::OrbitalSectorsLJPN& sectors);
  // extract a single multiple order from a data file
  void ReadMultipole(const std::string& data_filename, int order, basis::MatrixVector& matrices);

  // info file data structure
  struct InfoLine {
    basis::FullOrbitalLabels bra_labels;
    basis::FullOrbitalLabels ket_labels;
    int multipole;
    InfoLine(basis::FullOrbitalLabels bra, basis::FullOrbitalLabels ket, int mp)
    : bra_labels(bra), ket_labels(ket), multipole(mp) {}
  };

  // accessors
  const std::vector<InfoLine>& obdme_info() const { return obdme_info_; }

 private:

  // read info header
  void ReadInfoHeader();

  // read info
  void ReadInfo();

  // read data header
  void ReadDataHeader(std::ifstream& data_stream, int& data_line_count);

  // indexing accessors
  const basis::OrbitalSpaceLJPN& orbital_space() const {
    return orbital_space_;
  }

  // indexing information
  basis::OrbitalSpaceLJPN orbital_space_;

  // filename
  std::string info_filename_;

  // info file stream
  std::ifstream& info_stream() const {return *info_stream_ptr_;}  // alias for convenience
  std::ifstream* info_stream_ptr_;
  int line_count_;

  // info file header
  int max_K_;
  int num_proton_obdme_;
  int num_neutron_obdme_;

  // info container
  std::vector<InfoLine> obdme_info_;
};

};  // namespace shell
#endif  // RADIAL_IO_H_
