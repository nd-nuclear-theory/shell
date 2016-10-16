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
     
****************************************************************/

#ifndef H2_IO_H_
#define H2_IO_H_

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include "eigen3/Eigen/Core"

#include "basis/jjjpnorb_scheme.h"

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
  // H2 streams
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

    // constructor
    H2StreamBase(const std::string& filename)
      // Default constructor to zero-initialize some of the POD
      // fields.
      : filename_(filename), sector_index_(0), num_sectors_by_type_({0,0,0}), size_by_type_({0,0,0})
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
      return sectors_by_type_;
    }
    int num_sectors() const
    {
      int total = num_sectors_by_type_[0]+num_sectors_by_type_[1]+num_sectors_by_type_[2];
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

    // mode information
    H2Format h2_format_;
    H2Mode h2_mode_;

    // filename
    std::string filename_;

    // current pointer
    int sector_index_;
  };

  class InH2Stream
    : public H2StreamBase
  // Input stream for H2 file
  {
  public:

    // constructor
    InH2Stream(const std::string& filename);

    // destructor
    ~InH2Stream ()
      {
        Close();
        delete stream_ptr_;
      };

    // I/O
    void ReadSector (Eigen::MatrixXd& matrix);
    void SkipSector (Eigen::MatrixXd& matrix);
    void Close ();

  private:

    // format-specific implementation methods
    void ReadVersion_Text();
    void ReadVersion_Binary();
    void ReadHeader_Version0_Text();
    void ReadHeader_Version0_Binary();
    void ReadSector_Version0_Text(Eigen::MatrixXd& matrix, bool store);
    void ReadSector_Version0_Binary(Eigen::MatrixXd& matrix, bool store);
    void ConstructIndexing_Version0(int N1max, int N2max, int size_pp_nn, int size_pn);

    // file stream
    std::ifstream* stream_ptr_;
    int line_count_;  // for text mode input

  };

  class OutH2Stream
    : public H2StreamBase
  // Output stream for H2 file
  {
  public:

    // constructor
    OutH2Stream(
        const std::string& filename,
        const basis::OrbitalSpacePN& orbital_space,
        const basis::TwoBodySpaceJJJPN& space,
        const basis::TwoBodySectorsJJJPN& sectors,
        H2Format h2_format, H2Mode h2_mode
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

  private:
    
    // specialized methods
    void WriteVersion_Text();
    void WriteVersion_Binary();
    void WriteHeader_Version0_Text();
    void WriteHeader_Version0_Binary();
    void WriteSector_Version0_Text(const Eigen::MatrixXd& matrix);
    void WriteSector_Version0_Binary(const Eigen::MatrixXd& matrix);

    // file stream
    std::ofstream* stream_ptr_;

  };



  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif