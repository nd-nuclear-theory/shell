/****************************************************************
  h2_io.h                       

  Defines I/O class for MFDn H2 interaction file formats.
                                  
  Mark A. Caprio
  University of Notre Dame

  8/31/12 (mac): Adapted from mfdn_io.h (as mfdn_h2).  
     -- Converted I/O to class, from functions acting on stream
     argument.
     -- Added support for binary I/O and format code.
  1/28/13 (mac): Complete implementation.
  7/25/14 (mac): Add matrix (".mat") format for H2 files.
  4/25/15 (mac): Reformat source file.
  10/11/16 (mac,pjf):
    -- Rename to h2_io.
    -- Integrate into shell project build.
  10/19/16 (mac): Complete implementation for H2 Version0.
  10/25/16 (mac): Add InH2Stream::SeekToSector.
     
****************************************************************/

#ifndef H2_IO_H_
#define H2_IO_H_

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include "eigen3/Eigen/Core"

#include "basis/jjjpn_scheme.h"

namespace shell {


  ////////////////////////////////////////////////////////////////
  // H2 stream mode control
  ////////////////////////////////////////////////////////////////

  // H2 format code
  //
  //   0: oscillator-like sp indexing, used with MFDn v14
  //   15000: Pieter's preliminary general orbitals test format for MFDn v15 beta00
  //   15099: general orbitals specification June 2016

  typedef int H2Format;
  const int kVersion0=0;
  const int kVersion15000=15000;
  const int kVersion15099=15099;
  // enum class H2Format
  // {kVersion0=0,kVersion15000=15000,kVersion15099=15099};

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
    const std::array<int,3>& num_sectors_by_type() const
    {
      return num_sectors_by_type_;
    }
    int num_sectors() const
    {
      // int total = num_sectors_by_type_[0]+num_sectors_by_type_[1]+num_sectors_by_type_[2];  // okay, but...
      int total = sectors().size();  // simpler...
      return total;
    }
    const std::array<int,3>& size_by_type() const
    {
      return size_by_type_;
    }
    int size() const
    {
      int total = size_by_type_[0]+size_by_type_[1]+size_by_type_[2];
      return total;
    }
    const std::array<int,3>& Jmax_by_type() const
    {
      return Jmax_by_type_;
    }
    bool SectorIsFirstOfType() const;
    // Determine if current sector is first of its pp/nn/pn type.
    bool SectorIsLastOfType() const;
    // Determine if current sector is last of its pp/nn/pn type.

    // diagnostics
    std::string DiagnosticStr() const;

  protected:

    // indexing information
    basis::OrbitalSpacePN orbital_space_;
    basis::TwoBodySpaceJJJPN space_;
    basis::TwoBodySectorsJJJPN sectors_;

    // size information (for pp/nn/pn)
    std::array<int,3> num_sectors_by_type_;  // number of sectors
    std::array<int,3> size_by_type_;  // number of matrix elements
    std::array<int,3> Jmax_by_type_;  // twice Jmax

    // mode information
    H2Format h2_format_;
    H2Mode h2_mode_;

    // filename
    std::string filename_;

    // current pointer
    int sector_index_;
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

    InH2Stream() : stream_ptr_(NULL) {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    explicit InH2Stream(const std::string& filename);

    // destructor
    ~InH2Stream ()
      {
        Close();
        delete stream_ptr_;
      };

    // I/O

    void ReadSector(Eigen::MatrixXd& matrix);
    // Read current sector.

    void SkipSector();
    // Skip (read but do not store) current sector.

    void SeekToSector(int seek_index);
    // Skip (read but do not store) through sectors until arriving at
    // given sector.

    void Close();

  private:

    // format-specific implementation methods
    void ReadVersion();
    void ReadOrSkipSector(Eigen::MatrixXd& matrix,bool store);

    // ... Version0
    void ReadHeader_Version0();
    void ReadSector_Version0(Eigen::MatrixXd& matrix, bool store);

    // ... Version15099

    // file stream
    std::ifstream& stream() const {return *stream_ptr_;}  // alias for convenience
    std::ifstream* stream_ptr_;
    int line_count_;  // for text mode input

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

    OutH2Stream() : stream_ptr_(NULL) {};
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
        delete stream_ptr_;
      };

    // I/O
    void WriteSector(const Eigen::MatrixXd& matrix);
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
    void WriteSector_Version0(const Eigen::MatrixXd& matrix);

    // ... Version15099
    void WriteHeader_Version15099();
    void WriteSector_Version15099(const Eigen::MatrixXd& matrix);

    // file stream
    // std::ofstream& stream() const {return *stream_ptr_;}  // alias for convenience
    std::ofstream* stream_ptr_;

  };

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
