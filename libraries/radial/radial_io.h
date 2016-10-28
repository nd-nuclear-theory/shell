/// @file
/****************************************************************
  radial_io.h

  Defines I/O classes for radial matrix element storage.

  Language: C++11

  Patrick J. Fasano
  University of Notre Dame

  + 10/25/16 (pjf): Created, modeling after h2_io.h.

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
 * Container struct describing a radial operator.
 */
struct RadialOperator {
  char name;
  int l0max;
  int Tz0;

  RadialOperator() = default;

  RadialOperator(const char name, int l0max, int Tz0)
    : name(name), l0max(l0max), Tz0(Tz0) {}

  std::string DebugStr() const;
};

/**
 * Common kinematic operators
 */
// enum class KinematicOperator : RadialOperator {
//   kR = RadialOperator('r', 1, 0),
//   kRSquared = RadialOperator('r', 2, 0),
//   kP = RadialOperator('p', 1, 0),
//   kPSquared = RadialOperator('p', 2, 0)
//   kOverlap = RadialOperator('i', 0, 1)
// };

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

  // indexing accessors
  const basis::OrbitalSpaceLJPN& bra_orbital_space() const {
    return bra_orbital_space_;
  }
  const basis::OrbitalSpaceLJPN& ket_orbital_space() const {
    return ket_orbital_space_;
  }
  const basis::OrbitalSectorsLJPN& sectors() const { return sectors_; }
  const RadialOperator& radial_operator() const {
    return radial_operator_;
  }

 protected:
  // indexing information
  basis::OrbitalSpaceLJPN bra_orbital_space_;
  basis::OrbitalSpaceLJPN ket_orbital_space_;
  basis::OrbitalSectorsLJPN sectors_;
  RadialOperator radial_operator_;

  // current pointer
  int sector_index_;

  // matrix storage
  basis::MatrixVector matrices_;

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
    Close();
  }

  // I/O
  const basis::MatrixVector& Read() const {return matrices_;}
  void Close();

 private:
  void ReadHeader();
  void ReadNextSector();

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
      const RadialOperator& radial_operator);

  // destructor
  ~OutRadialStream() {
    Close();
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
