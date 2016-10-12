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

#include "legacy/shell_2body.h"

namespace shell {

  ////////////////////////////////////////////////////////////////
  // H2 header data
  ////////////////////////////////////////////////////////////////

  // Note: may ultimately hide header as internal to stream object

  struct MFDnH2Header
  {
    int format;
    int num_types;
    int N1b, N2b;
    int matrix_size_PPNN, matrix_size_PN;

    void Initialize(const legacy::TwoBodyBasisNljTzJP& basis);
  };

  ////////////////////////////////////////////////////////////////
  // H2 streams
  ////////////////////////////////////////////////////////////////

  // sector index is kept along with stream -- this helps to properly
  // keep track of whether or not sector exists and (especially) whether
  // or not it is at a file record boundary at the truncation defined
  // for *this* stream, which is not necessarily the truncation defined
  // for the sector iterator of the control loop which will be
  // reading/writing this stream

  // text/binary mode switching

  enum class H2TextBinaryMode { kText, kMatrix, kBinary };

  H2TextBinaryMode H2IOMode (const std::string& filename);
  // Deduce h2 file mode from extension.


  class OutH2Stream
  // Output stream for H2 file
  {
  public:

    // constructor
    OutH2Stream(
        const std::string& filename,
        const basis::OrbitalSpacePN& orbital_space,
        const basis::TwoBodySpaceJJJPN& space,
        const basis::TwoBodySectorsJJJPN& sectors
      );
    : filename_(filename),
        orbital_space_(orbital_space), space_(space), sectors_(sectors),
        stream_(0), sector_index_(0), size_by_type_({0,0,0})
        {};  // TODO compute size_by_type_ properly

    // destructor
    ~OutH2Stream()
      {
        Close();
        delete stream_;
      };

    // diagnostics
    std::string DiagnosticStr() const;

    // I/O
    void WriteSector(const Eigen::MatrixXd& matrix);
    void Close();

  private:
    // specialized methods
    void WriteHeaderText() const;
    void WriteHeaderBin() const;
    void WriteSectorText(const Eigen::MatrixXd& matrix);
    void WriteSectorMatrix(const Eigen::MatrixXd& matrix);
    void WriteSectorBin(const Eigen::MatrixXd& matrix);

    // file stream data
    H2TextBinaryMode text_binary_mode_;
    std::ofstream* stream_;
    std::string filename_;

    // mode information
    int format;

    // indexing information
    const basis::OrbitalSpacePN& orbital_space_;
    const basis::TwoBodySpaceJJJPN& space_;
    const basis::TwoBodySectorsJJJPN& sectors_;
    int std::array<3,int> size_by_type_;

    // current pointer
    int sector_index_;
  };



  class InH2Stream
  // Input stream for H2 file
  {
  public:

    // constructor
    InH2Stream(const std::string& filename);
    : filename_(filename)
    {};  // TODO compute everything properly

    // destructor
    ~InMFDnH2Stream ()
      {
        Close();
        delete stream_;
      };

    // diagnostics
    std::string DiagnosticStr() const;

    // I/O
    void ReadSector (Eigen::MatrixXd& matrix);
    void SkipSector (Eigen::MatrixXd& matrix);
    void Close ();

    // accessors
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

  private:
    // specialized methods
    void ReadHeaderText();
    void ReadHeaderBin();
    void ReadSectorText(Eigen::MatrixXd& matrix, bool store);
    void ReadSectorMatrix (Eigen::MatrixXd& matrix, bool store);
    void ReadSectorBin(Eigen::MatrixXd& matrix, bool store);

    // file stream data
    H2TextBinaryMode text_binary_mode_;
    std::ifstream* stream_;
    std::string filename_;

    // mode information
    int format;

    // indexing information
    basis::OrbitalSpacePN orbital_space_;
    basis::TwoBodySpaceJJJPN space_;
    basis::TwoBodySectorsJJJPN sectors_;
    int std::array<3,int> size_by_type_;

    // current pointer
    int sector_index_;
  };




  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
