/************************************************************//**
  @file obdme_io.h

  Defines I/O classes for one-body density matrix element access.

  Language: C++11

  Patrick J. Fasano
  University of Notre Dame

  + 11/22/16 (pjf): Created, based on radial_io.
  + 1/11/17 (pjf): Implemented MFDn v15 OBDME input.
    - InOBDMEStreamMulti parses info files at construction time.
    - The same reader can parse multiple data files.
    - Quantum numbers stored in data files are currently ignored.
  + 10/16/17 (pjf):
    - Define reader with g0 and Tz0.
    - TODO: check these values against those stored in data files.
  + 07/25/17 (pjf):
    - Add support for version 1600 OBDMEs from postprocessor.
    - Reorganize class structure so that multi-file OBDME
      formats (v1405/v1500) and single-file formats (v1600)
      share a common storage scheme.
    - Store all OBDMEs for a single data file.


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
 * Base class for unified one-body density matrix element retrieval
 */
class InOBDMEStream {
 public:
  InOBDMEStream() = default;
  void GetMultipole(int K, basis::OrbitalSectorsLJPN& sectors, basis::OperatorBlocks<double>& matrices) const;
  int min_K() const { return min_K_; }
  int max_K() const { return max_K_; }
  int g0()    const { return g0_; }
  int Tz0()   const { return Tz0_; }

 protected:
  InOBDMEStream(
      const basis::OrbitalSpaceLJPN& orbital_space,
      int g0 = 0,
      int Tz0 = 0
    ) : orbital_space_(orbital_space), g0_(g0), Tz0_(Tz0) {};

  // allocate and zero indexing and matrices
  void InitStorage();

  // indexing information
  basis::OrbitalSpaceLJPN orbital_space_;
  int g0_, Tz0_;
  int min_K_, max_K_;

  // indexing accessors
  const basis::OrbitalSpaceLJPN& orbital_space() const {
   return orbital_space_;
  }

  // matrix element storage
  std::vector<basis::OrbitalSectorsLJPN> sectors_;
  std::vector<basis::OperatorBlocks<double>> matrices_;
};


/**
 * Class for reading old-style, multi-file one-body density matrix elements
 */
class InOBDMEStreamMulti : public InOBDMEStream {
 public:
  /**
   * Default constructor -- provided since required for certain
   * purposes by STL container classes (e.g., std::vector::resize)
   */
  InOBDMEStreamMulti() = default;

  InOBDMEStreamMulti(
      const std::string& info_filename,
      const std::string& data_filename,
      const basis::OrbitalSpaceLJPN& orbital_space,
      int g0 = 0,
      int Tz0 = 0
    );
  // Construct a reader by parsing an info file.

  // destructor
  ~InOBDMEStreamMulti() {
    if (info_stream_ptr_)
      delete info_stream_ptr_;
  }

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
  void ReadInfoHeader1405();
  void ReadInfoHeader1500();

  // read info
  void ReadInfo();
  void ReadInfo1405();
  void ReadInfo1500();

  // read data header
  void ReadDataHeader(std::ifstream& data_stream, int& data_line_count) const;
  // void ReadDataHeader1405(std::ifstream& data_stream, int& data_line_count) const;
  // void ReadDataHeader1500(std::ifstream& data_stream, int& data_line_count) const;

  // read data
  void ReadData();

  // filename
  std::string info_filename_;
  std::string data_filename_;

  // info file stream
  std::ifstream& info_stream() const {return *info_stream_ptr_;}  // alias for convenience
  std::ifstream* info_stream_ptr_;
  int line_count_;

  // info file header
  int version_number_;
  int num_proton_obdme_;
  int num_neutron_obdme_;

  // info container
  std::vector<InfoLine> obdme_info_;
};


/**
 * Class for reading new-style, single-file one-body density matrix elements
 */
class InOBDMEStreamSingle : public InOBDMEStream {
 public:
  /**
   * Default constructor -- provided since required for certain
   * purposes by STL container classes (e.g., std::vector::resize)
   */
  InOBDMEStreamSingle() = default;

  InOBDMEStreamSingle(
      const std::string& filename,
      const basis::OrbitalSpaceLJPN& orbital_space
    );
  // Construct a reader by parsing an info file.

  // destructor
  ~InOBDMEStreamSingle() {
    if (stream_ptr_)
      delete stream_ptr_;
  }

private:

  // read info header
  void ReadHeader();
  void ReadHeader1600();

  // read data
  void ReadData();
  void ReadData1600();

  // filename
  std::string filename_;

  // file stream
  std::ifstream& stream() const {return *stream_ptr_;}  // alias for convenience
  std::ifstream* stream_ptr_;
  int line_count_;

  // file header
  int version_number_;
  int num_proton_obdme_;
  int num_neutron_obdme_;

};

};  // namespace shell
#endif  // RADIAL_IO_H_
