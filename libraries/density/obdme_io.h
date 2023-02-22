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
  + 04/03/19 (pjf):
    - Use mcutils::GetLine for input.
    - Modify reading from version 1520 OBDME files (formerly known as 1600).
    - Convert to Rose convention on input, for consistency with other
      one-body operators.
  + 05/09/19 (pjf): Use std::size_t for basis indices and sizes.
  + 05/28/19 (pjf):
    + Rename K->j0, max_K->j0_max, and min_K->j0_min.
    + Deprecate max_K() and min_K() accessors.
    + Modify access specifications and provide accessors for matrices and
      sectors as a function of j0.
    + Fix indexing problems for j0_min != 0.
  + 08/17/19 (pjf): Fix conversion to Rose convention.
  + 02/27/20 (pjf): Restructure class hierarchy:
    - Make InOBDMEStream a virtual base class.
    - Move common accessors into InOBDMEStream.
    - Require quantum numbers to construct stream.

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
  virtual ~InOBDMEStream() = default;

  void GetMultipole(
      int j0,
      basis::OrbitalSectorsLJPN& sectors,
      basis::OperatorBlocks<double>& matrices
    ) const;

  // accessors
  int Z_bra() const { return Z_bra_; }
  int N_bra() const { return N_bra_; }
  int A_bra() const { return Z_bra_ + N_bra_; }
  HalfInt Tz_bra() const { return HalfInt(Z_bra_-N_bra_, 2); }
  HalfInt J_bra() const { return J_bra_; }
  HalfInt M_bra() const { return M_bra_; }
  int g_bra() const { return g_bra_; }
  int n_bra() const { return n_bra_; }
  double T_bra() const { return T_bra_; }

  int Z_ket() const { return Z_ket_; }
  int N_ket() const { return N_ket_; }
  int A_ket() const { return Z_ket_ + N_ket_; }
  HalfInt Tz_ket() const { return HalfInt(Z_ket_-N_ket_, 2); }
  HalfInt J_ket() const { return J_ket_; }
  HalfInt M_ket() const { return M_ket_; }
  int g_ket() const { return g_ket_; }
  int n_ket() const { return n_ket_; }
  double T_ket() const { return T_ket_; }
  int j0_min() const { return j0_min_; }
  int j0_max() const { return j0_max_; }
  int g0()     const { return (g_bra_+g_ket_)%2; }
  int Tz0()    const { return int(Tz_bra()-Tz_ket());}

  // indexing accessors
  const basis::OrbitalSpaceLJPN& orbital_space() const { return orbital_space_; }

  const basis::OrbitalSectorsLJPN& sectors(int j0) const {
    assert((j0 >= j0_min()) && (j0 <= j0_max()));
    return sectors_.at(j0-j0_min());
  }

  const basis::OperatorBlocks<double>& matrices(int j0) const {
    assert((j0 >= j0_min()) && (j0 <= j0_max()));
    return matrices_.at(j0-j0_min());
  }

 protected:
  InOBDMEStream(
      const basis::OrbitalSpaceLJPN& orbital_space,
      HalfInt J_bra, int g_bra, int n_bra,
      HalfInt J_ket, int g_ket, int n_ket
    ) : orbital_space_(orbital_space),
        J_bra_(J_bra), g_bra_(g_bra), n_bra_(n_bra),
        J_ket_(J_ket), g_ket_(g_ket), n_ket_(n_ket)
    {};

  // allocate and zero indexing and matrices
  void InitStorage();

  // header quantum numbers (extracted from file)
  int Z_bra_, Z_ket_;
  int N_bra_, N_ket_;
  HalfInt M_bra_, M_ket_;
  float T_bra_, T_ket_;
  // indexing information
  basis::OrbitalSpaceLJPN orbital_space_;
  int j0_min_, j0_max_;

  // indexing accessors
  basis::OrbitalSectorsLJPN& sectors(int j0) {
    assert((j0 >= j0_min()) && (j0 <= j0_max()));
    return sectors_.at(j0-j0_min());
  }

  basis::OperatorBlocks<double>& matrices(int j0) {
    assert((j0 >= j0_min()) && (j0 <= j0_max()));
    return matrices_.at(j0-j0_min());
  }

 private:
  // state quantum numbers
  HalfInt J_bra_, J_ket_;
  int g_bra_, g_ket_;
  int n_bra_, n_ket_;

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
      HalfInt J_bra, int g_bra, int n_bra,
      HalfInt J_ket, int g_ket, int n_ket
    );
  // Construct a reader by parsing an info file.

  // destructor
  ~InOBDMEStreamMulti() {
    if (info_stream_ptr_)
      delete info_stream_ptr_;
  }

 private:
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

  // read info header
  void ReadInfoHeader();
  void ReadInfoHeader1405();
  void ReadInfoHeader1500();

  // read info
  void ReadInfo();
  void ReadInfo1405();
  void ReadInfo1500();

  // read data header
  void ReadDataHeader(std::ifstream& data_stream, int& data_line_count);
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
  std::size_t num_proton_obdme_;
  std::size_t num_neutron_obdme_;

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
      const basis::OrbitalSpaceLJPN& orbital_space,
      HalfInt J_bra, int g_bra, int n_bra,
      HalfInt J_ket, int g_ket, int n_ket
    );
  // Construct a reader by parsing an info file.

  // accessors
  double E_bra() const { return E_bra_; }
  double E_ket() const { return E_ket_; }

  // destructor
  ~InOBDMEStreamSingle() {
    if (stream_ptr_)
      delete stream_ptr_;
  }

private:

  // read info header
  void ReadHeader();
  void ReadHeader1520();

  // read data
  void ReadData();
  void ReadData1520();

  // filename
  std::string filename_;

  // file stream
  std::ifstream& stream() const {return *stream_ptr_;}  // alias for convenience
  std::ifstream* stream_ptr_;
  int line_count_;

  // file header
  int version_number_;
  int seq_bra_, seq_ket_;
  float E_bra_, E_ket_;

  std::size_t num_proton_obdme_;
  std::size_t num_neutron_obdme_;
  basis::OrbitalPNList orbital_list_;

};

};  // namespace shell
#endif  // RADIAL_IO_H_
