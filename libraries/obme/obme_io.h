/****************************************************************
  @file obme_io.h

  Defines I/O classes for one-body matrix element storage.

  Language: C++11

  Patrick J. Fasano
  University of Notre Dame

  + 10/25/16 (pjf): Created, modeling after h2_io.h.
  + 10/28/16 (pjf): Updated interface
   - RadialOperator is a simple enum class.
   - Don't store matrices internally to streams.
   - Read() takes an empty OperatorBlocks<double>&.
   - Constructors read/write headers but not data.
   - OutOBMEStream takes an OrbitalSectorsLJPN at construction.
  + 10/29/16 (mac):
   - Add error checking on stream open.
   - Fix double close on stream destruction.
  + 10/31/16 (mac):
   - Rename RadialOperator to RadialOperatorType and similarly rename
     stream accessor to radial_operator_type().
   - Move OutOBMEStream initializations into initializer list.
  + 11/2/16 (pjf): Added RadialOperatorType::kO.
  + 11/3/16 (mac): Provide InOBMEStream::SetToIndexing and hide direct
     indexing accessors (copies are easily invalidated).
  + 08/11/17 (pjf): Add verbose_mode option to OutOBMEStream.
  + 09/20/17 (pjf): Output Tz labels in verbose_mode.
  + 10/12/17 (pjf): Add support for generic and spherically-constrained operators:
    - Create new format (version 1) of radial file.
    - Separate operator type from operator truncation; add generic RadialOperatorType.
    - Store operator power separately for monomial radial operator.
  + 8/10/18 (pjf): Use mcutils::GetLine for more robust I/O.
  + 02/20/19 (pjf): Use ParseOrbitalStream() and OrbitalDefinitionStr() for
      orbital I/O.
  + 04/03/19 (pjf):
    - Use correct orbital formats for older files.
    - Prepend orbital index for orbitals in version 0 and 1 files.
****************************************************************/

#ifndef OBME_IO_H_
#define OBME_IO_H_

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "eigen3/Eigen/Core"

#include "basis/nlj_orbital.h"
#include "obme/obme_operator.h"

namespace shell {
/**
 * Base stream case with common attributes for input and output radial
 * streams.
 */
class OBMEStreamBase {
 public:
  /**
   * Default constructor -- provided since required for certain
   * purposes by STL container classes (e.g., std::vector::resize)
   */
  OBMEStreamBase() = default;

  explicit OBMEStreamBase(const std::string& filename)
      : filename_(filename), sector_index_(0) {}

  OBMEStreamBase(const std::string& filename,
                 const basis::OrbitalSpaceLJPN& bra_space,
                 const basis::OrbitalSpaceLJPN& ket_space,
                 const basis::OrbitalSectorsLJPN& sectors,
                 const basis::OneBodyOperatorType operator_type,
                 const RadialOperatorType radial_operator_type,
                 int radial_operator_power)
      : filename_(filename),
        sector_index_(0),
        bra_orbital_space_(bra_space),
        ket_orbital_space_(ket_space),
        operator_type_(operator_type),
        radial_operator_type_(radial_operator_type),
        radial_operator_power_(radial_operator_power),
        sectors_(sectors) {}

  // operator type accessor
  const basis::OneBodyOperatorType& operator_type() const { return operator_type_; }
  const RadialOperatorType& radial_operator_type() const { return radial_operator_type_; }

  // operator power accessor
  int radial_operator_power() const { return radial_operator_power_; }

 protected:
  // indexing accessors
  const basis::OrbitalSpaceLJPN& bra_orbital_space() const { return bra_orbital_space_; }
  const basis::OrbitalSpaceLJPN& ket_orbital_space() const { return ket_orbital_space_; }
  const basis::OrbitalSectorsLJPN& sectors() const { return sectors_; }
  // Warning: If you copy sectors(), caveat emptor.  This object
  // contains references to subspaces in bra_orbital_space_ and
  // ket_orbital_space_.  If the present RadialStream goes out of
  // scope and is destroyed, these references will be invalidated and
  // will point to "garbage" subspaces.

  // operator information
  basis::OneBodyOperatorType operator_type_;
  RadialOperatorType radial_operator_type_;
  int radial_operator_power_;

  // indexing information
  basis::OrbitalSpaceLJPN bra_orbital_space_;
  basis::OrbitalSpaceLJPN ket_orbital_space_;
  basis::OrbitalSectorsLJPN sectors_;

  // current pointer
  int sector_index_;

  // filename
  std::string filename_;
};

/**
 * Input stream for radial matrix element file.
 */
class InOBMEStream : public OBMEStreamBase {
 public:
  /**
   * Default constructor -- provided since required for certain
   * purposes by STL container classes (e.g., std::vector::resize)
   */
  InOBMEStream() : stream_ptr_(NULL), line_count_(0) {}

  explicit InOBMEStream(const std::string& filename);

  // destructor
  ~InOBMEStream() { delete stream_ptr_; }

  // I/O
  void SetToIndexing(basis::OrbitalSpaceLJPN& bra_orbital_space__,
                     basis::OrbitalSpaceLJPN& ket_orbital_space__,
                     basis::OrbitalSectorsLJPN& sectors__);
  // Set space and sectors indexing variables to point to fresh copies
  // which will not be invalidated when stream is destroyed.
  //
  // If the internal bra and ket orbital spaces are equal, it is
  // safe to use the same target variable to hold both.

  void Read(basis::OperatorBlocks<double>& matrices);
  void Close();

 private:
  void ReadHeader();
  Eigen::MatrixXd ReadNextSector();

  // file stream -- alias for convenience
  std::ifstream& stream() const { return *stream_ptr_; }
  std::ifstream* stream_ptr_;
  int line_count_;
};

/**
 * Output stream for radial matrix element file.
 */
class OutOBMEStream : public OBMEStreamBase {
 public:
  /**
   * Default constructor -- provided since required for certain
   * purposes by STL container classes (e.g., std::vector::resize)
   */
  OutOBMEStream() : stream_ptr_(NULL) {}

  explicit OutOBMEStream(const std::string& filename,
                         const basis::OrbitalSpaceLJPN& bra_space,
                         const basis::OrbitalSpaceLJPN& ket_space,
                         const basis::OrbitalSectorsLJPN& sectors,
                         const basis::OneBodyOperatorType operator_type,
                         const RadialOperatorType radial_operator_type,
                         int radial_operator_power,
                         const std::string& format_str,
                         bool verbose_mode = true);

  explicit OutOBMEStream(const std::string& filename,
                         const basis::OrbitalSpaceLJPN& bra_space,
                         const basis::OrbitalSpaceLJPN& ket_space,
                         const basis::OrbitalSectorsLJPN& sectors,
                         const basis::OneBodyOperatorType operator_type,
                         const RadialOperatorType radial_operator_type,
                         int radial_operator_power)
      : OutOBMEStream(filename, bra_space, ket_space, sectors, operator_type,
                      radial_operator_type, radial_operator_power, "16.8e", true) {}

  explicit OutOBMEStream(const std::string& filename, const basis::OrbitalSpaceLJPN& bra_space,
                         const basis::OrbitalSpaceLJPN& ket_space,
                         const basis::OrbitalSectorsLJPN& sectors,
                         const basis::OneBodyOperatorType operator_type,
                         const RadialOperatorType radial_operator_type,
                         int radial_operator_power, bool verbose_mode)
      : OutOBMEStream(filename, bra_space, ket_space, sectors, operator_type,
                      radial_operator_type, radial_operator_power, "16.8e", verbose_mode) {}

  // destructor
  ~OutOBMEStream() { delete stream_ptr_; }

  // I/O
  void Write(const basis::OperatorBlocks<double>& matrices);
  void Close();

 private:
  void WriteHeader();
  void WriteNextSector(const Eigen::MatrixXd& matrix);

  // formatting
  std::string format_str_;
  bool verbose_mode_;

  // file stream -- alias for convenience
  std::ofstream& stream() const { return *stream_ptr_; }
  std::ofstream* stream_ptr_;
  int line_count_;
};

};      // namespace shell
#endif  // OBME_IO_H_
