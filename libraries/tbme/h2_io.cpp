/****************************************************************
  h2_io.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "tbme/h2_io.h"

#include <cstddef>
#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <memory>

#include "fmt/format.h"
#include "mcutils/io.h"
#include "mcutils/parsing.h"

namespace shell {

  ////////////////////////////////////////////////////////////////
  // binary I/O parameters
  ////////////////////////////////////////////////////////////////

  // code assumes int and float are both 4-byte
  constexpr int kIntegerSize = 4;
  static_assert(sizeof(int)==kIntegerSize, "Integers are not 4 bytes.");
  constexpr int kLongIntegerSize = 8;
  static_assert(sizeof(long int)==kLongIntegerSize, "Long integers are not 8 bytes.");
  constexpr int kFloatSize = 4;
  static_assert(sizeof(float)==kFloatSize, "Floats are not 4 bytes.");

  // Fortran allows a maximum record length of (2**31-1) minus the bytes
  // for record overhead
  constexpr long int kMaxRecordLength = (long int)2147483647 - 2*kIntegerSize;

  ////////////////////////////////////////////////////////////////
  // space ordering
  ////////////////////////////////////////////////////////////////

  const std::unordered_map<H2Format,basis::TwoBodySpaceJJJPNOrdering> kH2SpaceOrdering({
      {kVersion0,     basis::TwoBodySpaceJJJPNOrdering::kPN},
      {kVersion15099, basis::TwoBodySpaceJJJPNOrdering::kPN},
      {kVersion15200, basis::TwoBodySpaceJJJPNOrdering::kTz},
    });

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
    for (std::size_t subspace_index=0; subspace_index<space.size(); ++subspace_index)
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

  std::size_t GetSectorCount(const basis::TwoBodySectorsJJJPN::SectorType& sector)
  {
    const auto& bra_subspace = sector.bra_subspace();
    const auto& ket_subspace = sector.ket_subspace();
    std::size_t sector_entries = 0;
    if (sector.IsDiagonal())
      // diagonal sector
      {
        std::size_t dimension = ket_subspace.size();
        sector_entries = dimension*(dimension+1)/2;
      }
    else if (sector.IsUpperTriangle())
      // upper triangle sector (but not diagonal)
      {
        std::size_t bra_dimension = bra_subspace.size();
        std::size_t ket_dimension = ket_subspace.size();
        sector_entries = bra_dimension*ket_dimension;
      }
    return sector_entries;
  }

  void EvaluateCountsByType(
      const basis::TwoBodySectorsJJJPN& sectors,
      std::array<std::size_t,3>& num_sectors_by_type,
      std::array<std::size_t,3>& size_by_type
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
  //   num_sectors_by_type (std::array<std::size_t,3>, output) : number of sectors
  //   size_by_type (std::array<std::size_t,3>, output) : number of matrix elements
  {

    // initialize
    num_sectors_by_type = {0,0,0};
    size_by_type = {0,0,0};

    // count over sectors
    for (std::size_t sector_index=0; sector_index<sectors.size(); ++sector_index)
      {

        // extract sector
        const auto& sector = sectors.GetSector(sector_index);
        const auto& bra_subspace = sector.bra_subspace();
        const auto& ket_subspace = sector.ket_subspace();

        // characterize bra pp/nn/pn type for sector
        basis::TwoBodySpeciesPN two_body_species = bra_subspace.two_body_species();

        // count sector
        ++num_sectors_by_type[int(two_body_species)];

        // count sector entries
        size_by_type[int(two_body_species)] += GetSectorCount(sector);
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
    assert((sector_index_>=0) && (sector_index_<sectors().size()));
    // short-circuit on first sector overall
    if (sector_index_ == 0)
      return true;
    const auto& current_sector = sectors().GetSector(sector_index_);
    auto current_type = current_sector.bra_subspace().two_body_species();

    const auto& prev_sector = sectors().GetSector(sector_index_-1);
    auto prev_type = prev_sector.bra_subspace().two_body_species();

    bool is_first_of_type = (current_type!=prev_type);
    return is_first_of_type;
  }

  bool H2StreamBase::SectorIsLastOfType() const
  {
    assert((sector_index_>=0) && (sector_index_<sectors().size()));
    // short-circuit on first sector overall
    if (sector_index_ == sectors().size()-1)
      return true;
    const auto& current_sector = sectors().GetSector(sector_index_);
    auto current_type = current_sector.bra_subspace().two_body_species();

    const auto& next_sector = sectors().GetSector(sector_index_+1);
    auto next_type = next_sector.bra_subspace().two_body_species();

    bool is_last_of_type = (current_type!=next_type);
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
    : H2StreamBase(filename), line_count_(0), stream_sector_positions_()
  {

    // determine mode
    h2_mode_ = DeducedIOMode(filename_);

    // open stream
    std::ios_base::openmode mode_argument = std::ios_base::in;
    if (h2_mode_==H2Mode::kBinary)
      mode_argument |= std::ios_base::binary;
    stream_ptr_ = std::unique_ptr<std::ifstream>(new std::ifstream(filename_,mode_argument));
    StreamCheck(bool(stream()),filename_,"Failure opening H2 file for input");

    // read format version and header and set up indexing
    ReadVersion();
    if (h2_format()==kVersion0)
      ReadHeader_Version0();
    else if (h2_format()==kVersion15099)
      ReadHeader_Version15099();
    else if (h2_format()==kVersion15200)
      ReadHeader_Version15200();
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


  void InH2Stream::ReadSector(std::size_t sector_index, Eigen::MatrixXd& matrix)
  {
    // jump to correct sector
    SeekToSector(sector_index);

    // read sector
    if (h2_format()==kVersion0)
      ReadSector_Version0(matrix);
    else if (h2_format()==kVersion15099)
      ReadSector_Version15099(matrix);
    else if (h2_format()==kVersion15200)
      ReadSector_Version15200(matrix);
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

    // store next sector location if not seen before
    if (sector_index_ > stream_sector_positions_.size()-1)
      stream_sector_positions_.push_back(stream().tellg());
  }

  void InH2Stream::SkipSector()
  {
    // skip sector
    if (h2_format()==kVersion0)
      SkipSector_Version0();
    else if (h2_format()==kVersion15099)
      SkipSector_Version15099();
    else if (h2_format()==kVersion15200)
      SkipSector_Version15200();
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

    // store next sector location if not seen before
    if (sector_index_ > stream_sector_positions_.size()-1)
      stream_sector_positions_.push_back(stream().tellg());
  }

  void InH2Stream::SeekToSector(std::size_t seek_index)
  {
    assert(seek_index!=basis::kNone);
    assert(seek_index<sectors().size());
    // check if we can seek directly to file pointer position
    if (seek_index <= stream_sector_positions_.size()-1)
      // sector has already been encountered in file scan
      {
        stream().seekg(stream_sector_positions_.at(seek_index));
        sector_index_ = seek_index;
      }
    else
      // need to seek forward to new part of file
      {
        while (sector_index_<seek_index)
          SkipSector();
      }
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
    stream_ptr_ = std::unique_ptr<std::ofstream>(new std::ofstream(filename_,mode_argument));
    StreamCheck(bool(stream()),filename_,"Failure opening H2 file for output");

    // write version
    WriteVersion();

    // write header
    if (h2_format_==kVersion0)
      WriteHeader_Version0();
    else if (h2_format_==kVersion15099)
      WriteHeader_Version15099();
    else if (h2_format_==kVersion15200)
      WriteHeader_Version15200();
    else
      {
        std::cerr << "Unsupported version encountered when opening H2 file " << filename_ << std::endl;
        std::exit(EXIT_FAILURE);
      }
    StreamCheck(bool(stream()),filename_,"Failure while writing H2 file header");
  }

  void OutH2Stream::WriteSector(
      std::size_t sector_index,
      const Eigen::MatrixXd& matrix,
      basis::NormalizationConversion conversion_mode
    )
  {
    // validate matrix dimensions
    const auto& sector = sectors().GetSector(sector_index_);
    const auto& bra_subspace = sector.bra_subspace();
    const auto& ket_subspace = sector.ket_subspace();
    assert(sector_index==sector_index_);
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
    else if (h2_format()==kVersion15200)
      WriteSector_Version15200(matrix,conversion_mode);
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
#include "h2_io_v15200.cpp"


