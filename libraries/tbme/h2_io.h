/****************************************************************
  h2_io.h

  Defines I/O class for MFDn H2 interaction file formats.

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

  + 08/31/12 (mac): Adapted from mfdn_io.h (as mfdn_h2).
     - Converted I/O to class, from functions acting on stream
       argument.
     - Added support for binary I/O and format code.
  + 01/28/13 (mac): Complete implementation.
  + 07/25/14 (mac): Add matrix (".mat") format for H2 files.
  + 04/25/15 (mac): Reformat source file.
  + 10/11/16 (mac,pjf):
    - Rename to h2_io.
    - Integrate into shell project build.
  + 10/19/16 (mac): Complete implementation for H2 Version0.
  + 10/25/16 (mac): Add InH2Stream::SeekToSector.
  + 11/01/16 (mac):
    - Convert from AS to NAS storage.
  + 11/13/16 (mac): Implement H2 Version15099 binary output.
  + 11/28/16 (mac): Add Tz0 to h2 format 15099 output header.
  + 10/19/17 (mac): Add optional on-the-fly conversion from AS to
    NAS matrix elements on output.
  + 01/22/18 (mac): Begin implementing nonzero Tz0.
  + 02/12/19 (pjf): Finish implementing nonzero Tz0.
  + 02/21/19 (pjf):
    - Remove Tz0!=0 support from h2v15099.
    - Implement H2 Version15200.
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

#ifndef H2_IO_H_
#define H2_IO_H_

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <memory>

#include <Eigen/Core>

#include "basis/jjjpn_scheme.h"

namespace shell {


  ////////////////////////////////////////////////////////////////
  // H2 stream mode control
  ////////////////////////////////////////////////////////////////

  // H2 format code
  //
  //   0: oscillator-like sp indexing, used with MFDn v14
  //   15000: Pieter's preliminary general orbitals test format for MFDn v15 beta00 (not supported)
  //   15099: general orbitals specification (June 2016)
  //   15200: add support for nonzero Tz0 (February 2019) -- WIP

  typedef int H2Format;
  constexpr H2Format kVersion0=0;
  constexpr H2Format kVersion15099=15099;
  constexpr H2Format kVersion15200=15200;
  // enum class H2Format
  // {kVersion0=0,kVersion15000=15000,kVersion15099=15099,kVersion15200=15200};

  enum class H2Mode {kText,kBinary,kMatrix};
  // text/binary mode
  //
  // Note: Not currently implementing matrix mode (for interchange
  // with applied math group).

  // notational definitions for h2 file modes
  //
  // Use of these arrays requires conversion of the H2Mode to int.
  extern const std::array<const char*,3> kH2ModeDescription; // ({"text","binary","matrix"});
  extern const std::array<const char*,3> kH2ModeExtension; // ({".dat",".bin",".mat"});

  extern const std::unordered_map<H2Format,basis::TwoBodySpaceJJJPNOrdering> kH2SpaceOrdering;
    // {{kVersion0,kPN},{kVersion15099,kPN},{kVersion15200,kTz}}
  // space ordering

  H2Mode DeducedIOMode(const std::string& filename);
  // Deduce h2 file mode from filename extension.
  //
  //   .dat: ascii format
  //   .bin: binary format
  //   .mat: ad hoc matrix format
  //     for interface with applied mathematicians
  //
  // Arguments:
  //   filename (string) : filename from which to deduce mode
  //
  // Returns:
  //   (H2Mode) : the mode

  ////////////////////////////////////////////////////////////////
  // H2StreamBase
  ////////////////////////////////////////////////////////////////

  // sector index is kept along with stream -- this helps to properly
  // keep track of whether or not sector exists and (especially) whether
  // or not it is at a file record boundary at the truncation defined
  // for *this* stream, which is not necessarily the truncation defined
  // for the sector iterator of the control loop which will be
  // reading/writing this stream

  class H2StreamBase
  // Base class with common attributes for input and output H2
  // streams.
  {
  public:

    // constructors

    H2StreamBase() = default;
    // default constructor -- provided since needed by derived type
    // default constructors

    H2StreamBase(const std::string& filename)
      // Default constructor to zero-initialize some of the POD
      // fields.
      : filename_(filename), sector_index_(0),
      num_sectors_by_type_({0,0,0}), size_by_type_({0,0,0}), Jmax_by_type_({0,0,0})
      {}

    // indexing accessors
    const basis::OrbitalSpacePN& orbital_space() const
    {
      return orbital_space_;
    }
    const basis::TwoBodySpaceJJJPN& space() const
    {
      return space_;
    }
    const basis::TwoBodySectorsJJJPN& sectors() const
    {
      return sectors_;
    }

    // stream status accessors
    H2Format h2_format() const
    {
      return h2_format_;
    }
    H2Mode h2_mode() const
    {
      return h2_mode_;
    }

    // size accessors
    const std::array<std::size_t,3>& num_sectors_by_type() const
    {
      return num_sectors_by_type_;
    }
    std::size_t num_sectors() const
    {
      // std::size_t total = num_sectors_by_type_[0]+num_sectors_by_type_[1]+num_sectors_by_type_[2];  // okay, but...
      std::size_t total = sectors().size();  // simpler...
      return total;
    }
    const std::array<std::size_t,3>& size_by_type() const
    {
      return size_by_type_;
    }
    std::size_t size() const
    {
      std::size_t total = size_by_type_[0]+size_by_type_[1]+size_by_type_[2];
      return total;
    }
    const std::array<int,3>& Jmax_by_type() const
    {
      return Jmax_by_type_;
    }
    // diagnostics
    std::string DiagnosticStr() const;

  protected:

    // indexing information
    basis::OrbitalSpacePN orbital_space_;
    basis::TwoBodySpaceJJJPN space_;
    basis::TwoBodySectorsJJJPN sectors_;

    // size information (for pp/nn/pn)
    std::array<std::size_t,3> num_sectors_by_type_;  // number of sectors
    std::array<std::size_t,3> size_by_type_;  // number of matrix elements
    std::array<int,3> Jmax_by_type_;  // twice Jmax

    // mode information
    H2Format h2_format_;
    H2Mode h2_mode_;

    // filename
    std::string filename_;

    // current pointer
    std::size_t sector_index_;

    bool SectorIsFirstOfType() const;
    // Determine if sector_index is first of its pp/nn/pn type.
    bool SectorIsLastOfType() const;
    // Determine if sector_index is last of its pp/nn/pn type.

  };

  ////////////////////////////////////////////////////////////////
  // InH2Stream
  ////////////////////////////////////////////////////////////////

  class InH2Stream
    : public H2StreamBase
  // Input stream for H2 file
  {
  public:

    // constructors

    InH2Stream() : stream_ptr_(nullptr) {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    explicit InH2Stream(const std::string& filename);

    // destructor
    ~InH2Stream ()
      {
        Close();
      };

    // I/O

    void ReadSector(std::size_t sector_index, Eigen::MatrixXd& matrix);
    // Read current sector.

    void Close();

  private:

    // format-specific implementation methods
    void ReadVersion();
    void SkipSector();
    // Skip (read but do not store) current sector.

    void SeekToSector(std::size_t seek_index);
    // Skip (read but do not store) through sectors until arriving at
    // given sector.

    // ... Version0
    void ReadHeader_Version0();
    void ReadSector_Version0(Eigen::MatrixXd& matrix);
    void SkipSector_Version0();

    // ... Version15099
    void ReadHeader_Version15099();
    void ReadSector_Version15099(Eigen::MatrixXd& matrix);
    void SkipSector_Version15099();

    // ... Version15200
    void ReadHeader_Version15200();
    void ReadSector_Version15200(Eigen::MatrixXd& matrix);
    void SkipSector_Version15200();

    // file stream
    std::ifstream& stream() const {return *stream_ptr_;}  // alias for convenience
    std::unique_ptr<std::ifstream> stream_ptr_;
    int line_count_;  // for text mode input
    std::vector<std::ifstream::pos_type> stream_sector_positions_;

  };

  ////////////////////////////////////////////////////////////////
  // OutH2Stream
  ////////////////////////////////////////////////////////////////

  class OutH2Stream
    : public H2StreamBase
  // Output stream for H2 file
  {
  public:

    // constructors

    OutH2Stream() : stream_ptr_(nullptr) {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    OutH2Stream(
        const std::string& filename,
        const basis::OrbitalSpacePN& orbital_space,
        const basis::TwoBodySpaceJJJPN& space,
        const basis::TwoBodySectorsJJJPN& sectors,
        H2Format h2_format
      );

    // destructor
    ~OutH2Stream()
      {
        Close();
      };

    // I/O
    void WriteSector(
        std::size_t sector_index,
        const Eigen::MatrixXd& matrix,
        basis::NormalizationConversion conversion_mode = basis::NormalizationConversion::kNone
      );
    void Close();

    // debugging
    // std::ofstream* stream_ptr() const {return stream_ptr_;}

  private:

    // convenience accessor (for internal use)
    std::ofstream& stream() const {return *stream_ptr_;}

    // format-specific implementation methods
    void WriteVersion();

    // ... Version0
    void WriteHeader_Version0();
    void WriteSector_Version0(const Eigen::MatrixXd& matrix, basis::NormalizationConversion conversion_mode);

    // ... Version15099
    void WriteHeader_Version15099();
    void WriteSector_Version15099(const Eigen::MatrixXd& matrix, basis::NormalizationConversion conversion_mode);

    // ... Version15100
    void WriteHeader_Version15200();
    void WriteSector_Version15200(const Eigen::MatrixXd& matrix, basis::NormalizationConversion conversion_mode);

    // file stream
    // std::ofstream& stream() const {return *stream_ptr_;}  // alias for convenience
    std::unique_ptr<std::ofstream> stream_ptr_;

  };

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
