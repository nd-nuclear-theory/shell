/****************************************************************
  h2_io.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "tbme/h2_io.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include "fmt/format.h"
#include "mcutils/io.h"
#include "mcutils/parsing.h"

namespace shell {

  ////////////////////////////////////////////////////////////////
  // binary I/O parameters
  ////////////////////////////////////////////////////////////////

  // code assumes int and float are both 4-byte
  const int kIntegerSize = 4;
  const int kFloatSize = 4;

  ////////////////////////////////////////////////////////////////
  // pp/nn/pn matrix size counting
  ////////////////////////////////////////////////////////////////

  void EvaluateJmaxByType(
      const basis::TwoBodySpaceJJJPN& space,
      std::array<int,3>& Jmax_by_type
    )
  // Extract Jmax value by two-body species type for JJJPN space.
  //
  // Arguments:
  //   space (..., input) : space
  //   Jmax_by_type (std::array<int,3>, output) : maximum J values
  {

    // initialize
    Jmax_by_type = std::array<int,3>({0,0,0});

    // scan subspaces for Jmax
    for (int subspace_index=0; subspace_index<space.size(); ++subspace_index)
      {
        // extract subspace
        const basis::TwoBodySubspaceJJJPN& subspace = space.GetSubspace(subspace_index);
        basis::TwoBodySpeciesPN two_body_species = subspace.two_body_species();

        // incorporate sector Jmax
        Jmax_by_type[int(two_body_species)] = std::max(
            Jmax_by_type[int(two_body_species)],
            subspace.J()
          );
      }
  }

  void EvaluateCountsByType(
      const basis::TwoBodySectorsJJJPN& sectors,
      std::array<int,3>& num_sectors_by_type,
      std::array<int,3>& size_by_type
    )
  // Count sectors and matrix elements in the upper triangular portion
  // of a set of pp/nn/pn sectors.
  //
  // Based on basis::UpperTriangularEntries.
  //
  // Lower triangular sectors are ignored.  In diagonal sectors,
  // only upper-triangular entries are counted.
  //
  // Arguments:
  //   sectors (..., input) : container for sectors
  //   num_sectors_by_type (std::array<int,3>, output) : number of sectors
  //   size_by_type (std::array<int,3>, output) : number of matrix elements
  {

    // initialize
    num_sectors_by_type = std::array<int,3>({0,0,0});
    size_by_type = std::array<int,3>({0,0,0});

    // count over sectors
    for (int sector_index=0; sector_index<sectors.size(); ++sector_index)
      {

        // extract sector
        const typename basis::TwoBodySectorsJJJPN::SectorType& sector = sectors.GetSector(sector_index);
        const typename basis::TwoBodySectorsJJJPN::SubspaceType& bra_subspace = sector.bra_subspace();
        const typename basis::TwoBodySectorsJJJPN::SubspaceType& ket_subspace = sector.ket_subspace();

        // characterize ket pp/nn/pn type for sector
        basis::TwoBodySpeciesPN two_body_species = ket_subspace.two_body_species();

        // count sector
        ++num_sectors_by_type[int(two_body_species)];

        // count sector entries
        int sector_entries = 0;
        if (sector.IsDiagonal())
          // diagonal sector
          {
            int dimension = ket_subspace.size();
            sector_entries = dimension*(dimension+1)/2;
          }
        else if (sector.IsUpperTriangle())
          // upper triangle sector (but not diagonal)
          {
            int bra_dimension = bra_subspace.size();
            int ket_dimension = ket_subspace.size();
            sector_entries = bra_dimension*ket_dimension;
          }
        size_by_type[int(two_body_species)] += sector_entries;
      }
  }

  ////////////////////////////////////////////////////////////////
  // file text/binary I/O mode identification
  ////////////////////////////////////////////////////////////////


  const std::array<const char*,3> kH2ModeDescription({"text","binary","matrix"});
  const std::array<const char*,3> kH2ModeExtension({".dat",".bin",".mat"});

  H2Mode DeducedIOMode(const std::string& filename)
  {
    if (filename.length() < 4 )
      {
        // prevent compare on underlength string
        std::cerr << "H2 file I/O: No extension found (too short) in filename " << filename << std::endl;
        exit(EXIT_FAILURE);
      }
    else if ( ! filename.compare(filename.length()-4,4,".dat") )
      return H2Mode::kText;
    else if ( ! filename.compare(filename.length()-4,4,".mat") )
      return H2Mode::kMatrix;
    else if ( ! filename.compare(filename.length()-4,4,".bin") )
      return H2Mode::kBinary;
    else
      {
        std::cerr << "H2 file I/O: Extension unrecognized in filename " << filename << std::endl;
        exit(EXIT_FAILURE);
      }
  }

  ////////////////////////////////////////////////////////////////
  // H2StreamBase
  ////////////////////////////////////////////////////////////////

  bool H2StreamBase::SectorIsFirstOfType() const
  {
    const typename basis::TwoBodySectorsJJJPN::SectorType& sector = sectors().GetSector(sector_index_);
    int current_type_index = int(sector.ket_subspace().two_body_species());

    // count sectors of prior types (exclusive)
    int prior_type_sectors = 0;
    for (int type_index=0; type_index<current_type_index;++type_index)
      prior_type_sectors+=num_sectors_by_type()[type_index];

    bool is_first_of_type = (sector_index_==prior_type_sectors);
    return is_first_of_type;
  }

  bool H2StreamBase::SectorIsLastOfType() const
  {
    const typename basis::TwoBodySectorsJJJPN::SectorType& sector = sectors().GetSector(sector_index_);
    int current_type_index = int(sector.ket_subspace().two_body_species());

    // count sectors of prior types (inclusive)
    int prior_type_sectors = 0;
    for (int type_index=0; type_index<=current_type_index;++type_index)
      prior_type_sectors+=num_sectors_by_type()[type_index];

    bool is_last_of_type = (sector_index_+1==prior_type_sectors);
    return is_last_of_type;
  }

  std::string H2StreamBase::DiagnosticStr() const
  {

    // TODO: make use of oscillator_like_ check for two-body space when available

    std::string str = fmt::format(
        "  File: {}\n"
        "  Format: {} ({})\n"
        "  Operator properties: J0 {} g0 {} Tz0 {}\n"
        "  Orbitals: p {} n {} (oscillator-like {})\n"
        // "  One-body truncation: p {:.4f} n {:.4f}\n"
        // "  Two-body truncation: pp {:.4f} nn {:.4f} pn {:.4f}\n"
        "  Truncation: p {:.4f} n {:.4f} pp {:.4f} nn {:.4f} pn {:.4f}\n"
        "  Sectors: pp {} nn {} pn {} => total {}\n"
        "  Matrix elements: pp {} nn {} pn {} => total {}\n",
        filename_,
        int(h2_format()),kH2ModeDescription[int(h2_mode())],
        sectors().J0(),sectors().g0(),sectors().Tz0(),
        orbital_space().GetSubspace(0).size(),orbital_space().GetSubspace(1).size(),orbital_space().is_oscillator_like(),
        space().weight_max().one_body[0],space().weight_max().one_body[1],
        space().weight_max().two_body[0],space().weight_max().two_body[1],space().weight_max().two_body[2],
        num_sectors_by_type()[0],num_sectors_by_type()[1],num_sectors_by_type()[2],num_sectors(),
        size_by_type()[0],size_by_type()[1],size_by_type()[2],size()
      );
    return str;
  }

  ////////////////////////////////////////////////////////////////
  // format/mode-specific I/O: version number
  ////////////////////////////////////////////////////////////////

  // read version: text mode
  void InH2Stream::ReadVersion()
  {

    if (h2_mode()==H2Mode::kText)
      {
        // set up line input
        std::string line;

        // line: format code
        {
          ++line_count_;
          std::getline(stream(),line);
          std::istringstream line_stream(line);
          line_stream >> h2_format_;
          ParsingCheck(line_stream,line_count_,line);
        }
      }
    else if (h2_mode()==H2Mode::kBinary)
      {
        int num_fields = 1;
        int bytes = num_fields * kIntegerSize;
        mcutils::VerifyBinary<int>(stream(),bytes,"Encountered unexpected value in H2 file","record delimiter");
        mcutils::ReadBinary<int>(stream(),h2_format_);
        mcutils::VerifyBinary<int>(stream(),bytes,"Encountered unexpected value in H2 file","record delimiter");
      }

    StreamCheck(bool(stream()),filename_,"Failure while reading H2 file version code");

  };

  void OutH2Stream::WriteVersion()
  {
    if (h2_mode()==H2Mode::kText)
      {
        stream()
          << fmt::format("{:10d}",int(h2_format())) << std::endl;
      }
    else if (h2_mode()==H2Mode::kBinary)
      {
        int num_fields = 1;
        int bytes = num_fields * kIntegerSize;
        mcutils::WriteBinary<int>(stream(),bytes);
        mcutils::WriteBinary<int>(stream(),h2_format());
        mcutils::WriteBinary<int>(stream(),bytes);
      }

    StreamCheck(bool(stream()),filename_,"Failure while writing H2 file version code");

  }


  ////////////////////////////////////////////////////////////////
  // InH2Stream
  ////////////////////////////////////////////////////////////////

  InH2Stream::InH2Stream(
      const std::string& filename
    )
    : H2StreamBase(filename), line_count_(0)
  {

    // determine mode
    h2_mode_ = DeducedIOMode(filename_);

    // open stream
    std::ios_base::openmode mode_argument = std::ios_base::in;
    if (h2_mode_==H2Mode::kBinary)
      mode_argument |= std::ios_base::binary;
    stream_ptr_ = new std::ifstream(filename_,mode_argument);
    StreamCheck(bool(stream()),filename_,"Failure opening H2 file for input");

    // read format version and header and set up indexing
    ReadVersion();
    if (h2_format()==kVersion0)
      ReadHeader_Version0();
    else if (h2_format()==kVersion15099)
      ReadHeader_Version15099();
    else
      {
        std::cerr << "Unsupported version encountered when opening H2 file " << filename_ << std::endl;
        std::exit(EXIT_FAILURE);
      }
  }

  void InH2Stream::Close()
  {
    stream().close();
  };


  void InH2Stream::ReadOrSkipSector(Eigen::MatrixXd& matrix, bool store)
  {
    // read sector
    if (h2_format()==kVersion0)
      ReadSector_Version0(matrix,store);
    else if (h2_format()==kVersion15099)
      ReadSector_Version15099(matrix,store);
    else
      // format version was already checked when reading header, so we should
      // never get here, unless perhaps someday we implement header-only support
      // for some file format version
      {
        std::cerr << "Unsupported version encountered when reading H2 file " << filename_ << std::endl;
        std::exit(EXIT_FAILURE);
      }

    // validate status
    StreamCheck(bool(stream()),filename_,"Failure while reading H2 file sector");

    // increment to next sector
    sector_index_++;

  }

  void InH2Stream::ReadSector(Eigen::MatrixXd& matrix)
  {
    bool store = true;
    ReadOrSkipSector(matrix,store);
  }

  void InH2Stream::SkipSector()
  {
    // Note that matrix object is still needed as dummy, to use
    // ReadSector machinery, even though matrix is not used
    Eigen::MatrixXd matrix;

    bool store = false;
    ReadOrSkipSector(matrix,store);
  }

  void InH2Stream::SeekToSector(int seek_index)
  {
    assert(seek_index!=basis::kNone);
    assert(seek_index<sectors().size());
    while (sector_index_<seek_index)
      SkipSector();
  }


  ////////////////////////////////////////////////////////////////
  // OutH2Stream
  ////////////////////////////////////////////////////////////////

  OutH2Stream::OutH2Stream(
      const std::string& filename,
      const basis::OrbitalSpacePN& orbital_space,
      const basis::TwoBodySpaceJJJPN& space,
      const basis::TwoBodySectorsJJJPN& sectors,
      H2Format h2_format
    )
    : H2StreamBase(filename)
  {

    // copy in indexing
    orbital_space_ = orbital_space;
    space_ = space;
    sectors_ = sectors;

    // store counts by type
    EvaluateJmaxByType(space_,Jmax_by_type_);
    EvaluateCountsByType(sectors_,num_sectors_by_type_,size_by_type_);

    // set format and mode
    h2_format_ = h2_format;
    h2_mode_ = DeducedIOMode(filename_);

    // open stream
    std::ios_base::openmode mode_argument = std::ios_base::out;
    if (h2_mode_==H2Mode::kBinary)
      mode_argument |= std::ios_base::binary;
    stream_ptr_ = new std::ofstream(filename_,mode_argument);
    StreamCheck(bool(stream()),filename_,"Failure opening H2 file for output");

    // write version
    WriteVersion();

    // write header
    if (h2_format_==kVersion0)
      WriteHeader_Version0();
    else if (h2_format_==kVersion15099)
      WriteHeader_Version15099();
    else
      {
        std::cerr << "Unsupported version encountered when opening H2 file " << filename_ << std::endl;
        std::exit(EXIT_FAILURE);
      }
    StreamCheck(bool(stream()),filename_,"Failure while writing H2 file header");
  }

  void OutH2Stream::WriteSector(
      const Eigen::MatrixXd& matrix,
      basis::NormalizationConversion conversion_mode
    )
  {
    // validate matrix dimensions
    const typename basis::TwoBodySectorsJJJPN::SectorType& sector = sectors().GetSector(sector_index_);
    const typename basis::TwoBodySectorsJJJPN::SubspaceType& bra_subspace = sector.bra_subspace();
    const typename basis::TwoBodySectorsJJJPN::SubspaceType& ket_subspace = sector.ket_subspace();
    assert((matrix.rows()==sector.bra_subspace().size())&&(matrix.cols()==sector.ket_subspace().size()));

    // validate output as NAS
    assert(
        (conversion_mode == basis::NormalizationConversion::kNone)
        || (conversion_mode == basis::NormalizationConversion::kASToNAS)
      );

    // write sector
    if (h2_format()==kVersion0)
      WriteSector_Version0(matrix,conversion_mode);
    else if (h2_format()==kVersion15099)
      WriteSector_Version15099(matrix,conversion_mode);
    else
      assert(false);  // format version was already checked when writing header

    // validate status
    StreamCheck(bool(stream()),filename_,"Failure while writing H2 file sector");

    // increment to next sector
    sector_index_++;

  }

  void OutH2Stream::Close()
  {
    stream().close();
  };
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

////////////////////////////////////////////////////////////////
// version-specific implementation
////////////////////////////////////////////////////////////////
//
// When adding a new include file, also update module.mk:
//   - Add to module_extras file list.
//   - Add to dependency for h2_io.o.

#include "h2_io_v0.cpp"
#include "h2_io_v15099.cpp"


