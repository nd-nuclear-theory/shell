/// @file
/****************************************************************
  radial_io.h

  Defines I/O classes for radial matrix element storage.

  Language: C++11

  Patrick J. Fasano
  University of Notre Dame

  + 10/25/16 (pjf): Created, modeling after h2_io.h.
  + 10/28/16 (pjf): Updated interface
   - RadialOperator is a simple enum class.
   - Don't store matrices internally to streams.
   - Read() takes an empty MatrixVector&.
   - Constructors read/write headers but not data.
   - OutRadialStream takes an OrbitalSectorsLJPN at construction.
  + 10/29/16 (mac):
   - Add error checking on stream open.
   - Fix double close on stream destruction.
  + 10/31/16 (mac):
   - Rename RadialOperator to RadialOperatorType and similarly rename
     stream accessor to radial_operator_type().
   - Move OutRadialStream initializations into initializer list.
  + 11/2/16 (pjf): Added RadialOperatorType::kO.

****************************************************************/

#ifndef RADIAL_IO_H_
#define RADIAL_IO_H_

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include "eigen3/Eigen/Core"

#include "basis/operator.h"
#include "basis/nlj_orbital.h"

namespace shell {

/**
 * Radial IDs
 */
enum class RadialOperatorType : char {
  kR = 'r', kK = 'k', kO = 'o'
};

/**
 * Base stream case with common attributes for input and output radial streams.
 */
class RadialStreamBase {
 public:
  /**
   * Default constructor -- provided since required for certain
   * purposes by STL container classes (e.g., std::vector::resize)
   */
  RadialStreamBase() = default;

  explicit RadialStreamBase(const std::string& filename)
    : filename_(filename), sector_index_(0) {}

  RadialStreamBase(
      const std::string& filename,
      const basis::OrbitalSpaceLJPN& bra_space,
      const basis::OrbitalSpaceLJPN& ket_space,
      const basis::OrbitalSectorsLJPN& sectors,
      const RadialOperatorType radial_operator_type
    )
    : filename_(filename), sector_index_(0),
    bra_orbital_space_(bra_space),
    ket_orbital_space_(ket_space),
    radial_operator_type_(radial_operator_type),
    sectors_(sectors) {}

  // indexing accessors
  const basis::OrbitalSpaceLJPN& bra_orbital_space() const {
    return bra_orbital_space_;
  }
  const basis::OrbitalSpaceLJPN& ket_orbital_space() const {
    return ket_orbital_space_;
  }
  const basis::OrbitalSectorsLJPN& sectors() const { return sectors_; }
  const RadialOperatorType& radial_operator_type() const {
    return radial_operator_type_;
  }

 protected:
  // indexing information
  basis::OrbitalSpaceLJPN bra_orbital_space_;
  basis::OrbitalSpaceLJPN ket_orbital_space_;
  basis::OrbitalSectorsLJPN sectors_;
  RadialOperatorType radial_operator_type_;

  // current pointer
  int sector_index_;

  // filename
  std::string filename_;
};

/**
 * Input stream for radial matrix element file.
 */
class InRadialStream : public RadialStreamBase {
 public:
  /**
   * Default constructor -- provided since required for certain
   * purposes by STL container classes (e.g., std::vector::resize)
   */
  InRadialStream() : stream_ptr_(NULL), line_count_(0) {}

  explicit InRadialStream(const std::string& filename);

  // destructor
  ~InRadialStream() {
    delete stream_ptr_;
  }

  // I/O
  void Read(basis::MatrixVector& matrices);
  void Close();

 private:
  void ReadHeader();
  Eigen::MatrixXd ReadNextSector();

  // file stream
  std::ifstream& stream() const {return *stream_ptr_;}  // alias for convenience
  std::ifstream* stream_ptr_;
  int line_count_;
};

/**
 * Output stream for radial matrix element file.
 */
class OutRadialStream : public RadialStreamBase {
 public:
  /**
  * Default constructor -- provided since required for certain
  * purposes by STL container classes (e.g., std::vector::resize)
  */
  OutRadialStream() : stream_ptr_(NULL) {}

  explicit OutRadialStream(
      const std::string& filename,
      const basis::OrbitalSpaceLJPN& bra_space,
      const basis::OrbitalSpaceLJPN& ket_space,
      const basis::OrbitalSectorsLJPN& sectors,
      const RadialOperatorType radial_operator_type);

  // destructor
  ~OutRadialStream() {
    delete stream_ptr_;
  }

  // I/O
  void Write(const basis::MatrixVector& matrices);
  void Close();

 private:
  void WriteHeader();
  void WriteNextSector(const Eigen::MatrixXd& matrix);

  // file stream
  std::ofstream& stream() const {return *stream_ptr_;}  // alias for convenience
  std::ofstream* stream_ptr_;
  int line_count_;
};

};  // namespace shell
#endif  // RADIAL_IO_H_
