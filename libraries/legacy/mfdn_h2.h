/****************************************************************
  mfdn_h2.h                       

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
  10/11/16 (mac): Integrate into shell project build.
     
****************************************************************/

#ifndef MFDN_H2_H_
#define MFDN_H2_H_

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include "legacy/shell_2body.h"

namespace legacy {

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

  // sector iterator is kept along with stream -- this helps to properly
  // keep track of whether or not sector exists and (especially) whether
  // or not it is at a file record boundary at the truncation defined
  // for *this* stream, which is not necessarily the truncation defined
  // for the sector iterator of the control loop which will be
  // reading/writing this stream

  // text/binary mode switching

  enum H2TextBinaryMode { kText, kMatrix, kBinary };
  H2TextBinaryMode H2IOMode (const std::string& filename);

  // output stream for H2 file

  class OutMFDnH2Stream
  {
  public:
    OutMFDnH2Stream () : stream_(0) {
    };
    ~OutMFDnH2Stream ()  {delete stream_;};
    void Open (const std::string& filename, const MFDnH2Header& header);
    void PrintDiagnostic (std::ostream& log_stream = std::cout);
    void WriteSector (const legacy::TwoBodyMatrixNljTzJP& matrix);
    void Close ();

  private:
    // specialized methods
    void WriteHeaderText() const;
    void WriteHeaderBin() const;
    void WriteSectorText (const legacy::TwoBodyMatrixNljTzJP& matrix);
    void WriteSectorMatrix (const legacy::TwoBodyMatrixNljTzJP& matrix);
    void WriteSectorBin (const legacy::TwoBodyMatrixNljTzJP& matrix);

    // data
    std::ofstream* stream_;
    std::string filename_;
    H2TextBinaryMode text_binary_mode_;  
    MFDnH2Header header_;
    legacy::SectorNljTzJP sector_;
  };


  // input stream for H2 file

  class InMFDnH2Stream
  {
  public:
    InMFDnH2Stream () : stream_(0) {
    };
    ~InMFDnH2Stream ()  {delete stream_;};
    void Open (const std::string& filename, MFDnH2Header& header);
    void PrintDiagnostic (std::ostream& log_stream = std::cout);
    void ReadSector (legacy::TwoBodyMatrixNljTzJP& matrix, bool store = true);
    void SkipSector (legacy::TwoBodyMatrixNljTzJP& matrix); // wrapper for ReadSector with no storage
    void Close ();

  private:
    // specialized methods
    void ReadHeaderText();
    void ReadHeaderBin();
    void ReadSectorText (legacy::TwoBodyMatrixNljTzJP& matrix, bool store);
    void ReadSectorMatrix (legacy::TwoBodyMatrixNljTzJP& matrix, bool store);
    void ReadSectorBin (legacy::TwoBodyMatrixNljTzJP& matrix, bool store);

    // data
    std::ifstream* stream_;
    std::string filename_;
    H2TextBinaryMode text_binary_mode_;  
    MFDnH2Header header_;
    legacy::SectorNljTzJP sector_;
  };




  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif