
  ////////////////////////////////////////////////////////////////
  // Me2j stream mode control
  ////////////////////////////////////////////////////////////////

  typedef int Me2jFormat;
  constexpr Me2jFormat kVersion0=0;

  class Me2jStreamBase
  // Base class with common attributes for input and output Me2j
  // streams.
  {
  public:

    // constructors

    Me2jStreamBase() = default;
    // default constructor -- provided since needed by derived type
    // default constructors

    Me2jStreamBase(const std::string& filename)
      // Default constructor to zero-initialize some of the POD
      // fields.
      : filename_(filename), sector_index_(0),
      num_sectors_by_type_({0,0,0}), size_by_type_({0,0,0}), Jmax_by_type_({0,0,0})
      {}

    // indexing accessors
    // const basis::OrbitalSpaceTTz& orbital_space() const
    // {
    //   return orbital_space_;
    // }
    const basis::TwoBodySpaceJJJTTz& space() const
    {
      return space_;
    }
    const basis::TwoBodySectorsJJJTTz& sectors() const
    {
      return sectors_;
    }

    // stream status accessors
    Me2jFormat me2j_format() const
    {
      return me2j_format_;
    }
    Me2jMode me2j_mode() const
    {
      return me2j_mode_;
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
    // basis::OrbitalSpaceTTz orbital_space_;
    basis::TwoBodySpaceJJJTTz space_;
    basis::TwoBodySectorsJJJTTz sectors_;

    // size information (for pp/nn/pn) TODO: these three may not matter anymore, and only appears in DiagnosticStr
    std::array<std::size_t,3> num_sectors_by_type_;  // number of sectors TODO: replace with lenth of sectors_
    std::array<std::size_t,3> size_by_type_;  // number of matrix elements
    std::array<int,3> Jmax_by_type_;  // twice Jmax

    // mode information
    Me2jFormat me2j_format_;
    Me2jMode me2j_mode_;

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
  // InMe2jStream
  ////////////////////////////////////////////////////////////////

  class InMe2jStream
    : public Me2jStreamBase
  // Input stream for Me2j file
  {
  public:

    // constructors

    InMe2jStream() : stream_ptr_(nullptr) {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    explicit InMe2jStream(const std::string& filename);

    // destructor
    ~InMe2jStream ()
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

    // file stream
    std::ifstream& stream() const {return *stream_ptr_;}  // alias for convenience
    std::unique_ptr<std::ifstream> stream_ptr_;
    int line_count_;  // for text mode input
    std::vector<std::ifstream::pos_type> stream_sector_positions_;

  };



  // ////////////////////////////////////////////////////////////////
  // // OutMe2jStream
  // ////////////////////////////////////////////////////////////////
  //
  // class OutMe2jStream
  //   : public Me2jStreamBase
  // // Output stream for Me2j file
  // {
  // public:
  //
  //   // constructors
  //
  //   OutMe2jStream() : stream_ptr_(nullptr) {};
  //   // default constructor -- provided since required for certain
  //   // purposes by STL container classes (e.g., std::vector::resize)
  //
  //   OutMe2jStream(
  //       const std::string& filename,
  //       const basis::OrbitalSpaceTTz& orbital_space,
  //       const basis::TwoBodySpaceJJJTTz& space,
  //       const basis::TwoBodySectorsJJJTTz& sectors,
  //       Me2jFormat me2j_format
  //     );
  //
  //   // destructor
  //   ~OutMe2jStream()
  //     {
  //       Close();
  //     };
  //
  //   // I/O
  //   void WriteSector(
  //       std::size_t sector_index,
  //       const Eigen::MatrixXd& matrix,
  //       basis::NormalizationConversion conversion_mode = basis::NormalizationConversion::kNone
  //     );
  //   void Close();
  //
  //   // debugging
  //   // std::ofstream* stream_ptr() const {return stream_ptr_;}
  //
  // private:
  //
  //   // convenience accessor (for internal use)
  //   std::ofstream& stream() const {return *stream_ptr_;}
  //
  //   // format-specific implementation methods
  //   void WriteVersion();
  //
  //   // ... Version0
  //   void WriteHeader_Version0();
  //   void WriteSector_Version0(const Eigen::MatrixXd& matrix, basis::NormalizationConversion conversion_mode);
  //
  //   // file stream
  //   // std::ofstream& stream() const {return *stream_ptr_;}  // alias for convenience
  //   std::unique_ptr<std::ofstream> stream_ptr_;
  //
  // };


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
    // pp/nn/pn matrix size counting
    ////////////////////////////////////////////////////////////////

    // void EvaluateJmaxByType(
    //     const basis::TwoBodySpaceJJJTTz& space,
    //     std::array<int,3>& Jmax_by_type
    //   )
    // // Extract Jmax value by two-body species type for JJJPN space.
    // //
    // // Arguments:
    // //   space (..., input) : space
    // //   Jmax_by_type (std::array<int,3>, output) : maximum J values
    // {
    //
    //   // initialize
    //   Jmax_by_type = std::array<int,3>({0,0,0});
    //
    //   // scan subspaces for Jmax
    //   for (std::size_t subspace_index=0; subspace_index<space.size(); ++subspace_index)
    //     {
    //       // extract subspace
    //       const basis::TwoBodySubspaceJJJTTz& subspace = space.GetSubspace(subspace_index);
    //       basis::TwoBodySpeciesTTz two_body_species = subspace.two_body_species();
    //
    //       // incorporate sector Jmax
    //       Jmax_by_type[int(two_body_species)] = std::max(
    //           Jmax_by_type[int(two_body_species)],
    //           subspace.J()
    //         );
    //     }
    // }

    std::size_t GetSectorCount(const basis::TwoBodySectorsJJJTTz::SectorType& sector)
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

    // void EvaluateCountsByType(
    //     const basis::TwoBodySectorsJJJTTz& sectors,
    //     std::array<std::size_t,3>& num_sectors_by_type,
    //     std::array<std::size_t,3>& size_by_type
    //   )
    // // Count sectors and matrix elements in the upper triangular portion
    // // of a set of pp/nn/pn sectors.
    // //
    // // Based on basis::UpperTriangularEntries.
    // //
    // // Lower triangular sectors are ignored.  In diagonal sectors,
    // // only upper-triangular entries are counted.
    // //
    // // Arguments:
    // //   sectors (..., input) : container for sectors
    // //   num_sectors_by_type (std::array<std::size_t,3>, output) : number of sectors
    // //   size_by_type (std::array<std::size_t,3>, output) : number of matrix elements
    // {
    //
    //   // initialize
    //   num_sectors_by_type = {0,0,0};
    //   size_by_type = {0,0,0};
    //
    //   // count over sectors
    //   for (std::size_t sector_index=0; sector_index<sectors.size(); ++sector_index)
    //     {
    //
    //       // extract sector
    //       const auto& sector = sectors.GetSector(sector_index);
    //       const auto& bra_subspace = sector.bra_subspace();
    //       const auto& ket_subspace = sector.ket_subspace();
    //
    //       // characterize bra pp/nn/pn type for sector
    //       basis::TwoBodySpeciesTTz two_body_species = bra_subspace.two_body_species();
    //
    //       // count sector
    //       ++num_sectors_by_type[int(two_body_species)];
    //
    //       // count sector entries
    //       size_by_type[int(two_body_species)] += GetSectorCount(sector);
    //     }
    // }
    ////////////////////////////////////////////////////////////////
    // Me2jStreamBase
    ////////////////////////////////////////////////////////////////

    bool Me2jStreamBase::SectorIsFirstOfType() const
    {
      assert((sector_index_!=basis::kNone) && (sector_index_<sectors().size()));
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

    bool Me2jStreamBase::SectorIsLastOfType() const
    {
      assert((sector_index_!=basis::kNone) && (sector_index_<sectors().size()));
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

    std::string Me2jStreamBase::DiagnosticStr() const
    {
      // TODO
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
          int(me2j_format()),kMe2jModeDescription[int(me2j_mode())],
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
    void InMe2jStream::ReadVersion()
    {

      if (me2j_mode()==Me2jMode::kText)
        {
          // set up line input
          std::string line;

          // line: format code
          {
            ++line_count_;
            std::getline(stream(),line);
            std::istringstream line_stream(line);
            line_stream >> me2j_format_;
            mcutils::ParsingCheck(line_stream,line_count_,line);
          }
        }
      else if (me2j_mode()==Me2jMode::kBinary)
        {
          int num_fields = 1;
          int bytes = num_fields * kIntegerSize;
          mcutils::VerifyBinary<int>(stream(),bytes,"Encountered unexpected value in Me2j file","record delimiter");
          mcutils::ReadBinary<int>(stream(),me2j_format_);
          mcutils::VerifyBinary<int>(stream(),bytes,"Encountered unexpected value in Me2j file","record delimiter");
        }

      mcutils::StreamCheck(bool(stream()),filename_,"Failure while reading Me2j file version code");

    };

    // void OutMe2jStream::WriteVersion()
    // {
    //   if (me2j_mode()==Me2jMode::kText)
    //     {
    //       stream()
    //         << fmt::format("{:10d}",int(me2j_format())) << std::endl;
    //     }
    //   else if (me2j_mode()==Me2jMode::kBinary)
    //     {
    //       int num_fields = 1;
    //       int bytes = num_fields * kIntegerSize;
    //       mcutils::WriteBinary<int>(stream(),bytes);
    //       mcutils::WriteBinary<int>(stream(),me2j_format());
    //       mcutils::WriteBinary<int>(stream(),bytes);
    //     }
    //
    //   mcutils::StreamCheck(bool(stream()),filename_,"Failure while writing Me2j file version code");
    //
    // }


      ////////////////////////////////////////////////////////////////
      // InMe2jStream
      ////////////////////////////////////////////////////////////////

      InMe2jStream::InMe2jStream(
          const std::string& filename
        )
        : Me2jStreamBase(filename), line_count_(0), stream_sector_positions_()
      {

        // determine mode
        me2j_mode_ = DeducedIOMode(filename_);

        // open stream
        std::ios_base::openmode mode_argument = std::ios_base::in;
        if (me2j_mode_==Me2jMode::kBinary)
          mode_argument |= std::ios_base::binary;
        stream_ptr_ = std::unique_ptr<std::ifstream>(new std::ifstream(filename_,mode_argument));
        mcutils::StreamCheck(bool(stream()),filename_,"Failure opening Me2j file for input");

        // read format version and header and set up indexing
        ReadVersion();
        ReadHeader_Version0();
      }

      void InMe2jStream::Close()
      {
        stream().close();
      };


      void InMe2jStream::ReadSector(std::size_t sector_index, Eigen::MatrixXd& matrix)
      {
        // jump to correct sector
        SeekToSector(sector_index);

        // read sector
        ReadSector_Version0(matrix);

        // validate status
        mcutils::StreamCheck(bool(stream()),filename_,"Failure while reading Me2j file sector");

        // increment to next sector
        sector_index_++;

        // store next sector location if not seen before
        if (sector_index_ > stream_sector_positions_.size()-1)
          stream_sector_positions_.push_back(stream().tellg());
      }

      void InMe2jStream::SkipSector()
      {
        // skip sector
        SkipSector_Version0();

        // validate status
        mcutils::StreamCheck(bool(stream()),filename_,"Failure while reading Me2j file sector");

        // increment to next sector
        sector_index_++;

        // store next sector location if not seen before
        if (sector_index_ > stream_sector_positions_.size()-1)
          stream_sector_positions_.push_back(stream().tellg());
      }

      void InMe2jStream::SeekToSector(std::size_t seek_index)
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

      // ////////////////////////////////////////////////////////////////
      // // OutMe2jStream
      // ////////////////////////////////////////////////////////////////
      //
      // OutMe2jStream::OutMe2jStream(
      //     const std::string& filename,
      //     const basis::OrbitalSpaceTTz& orbital_space,
      //     const basis::TwoBodySpaceJJJTTz& space,
      //     const basis::TwoBodySectorsJJJTTz& sectors,
      //     Me2jFormat me2j_format
      //   )
      //   : Me2jStreamBase(filename)
      // {
      //
      //   // copy in indexing
      //   orbital_space_ = orbital_space;
      //   space_ = space;
      //   sectors_ = sectors;
      //
      //   // store counts by type
      //   EvaluateJmaxByType(space_,Jmax_by_type_);
      //   EvaluateCountsByType(sectors_,num_sectors_by_type_,size_by_type_);
      //
      //   // set format and mode
      //   me2j_format_ = me2j_format;
      //   me2j_mode_ = DeducedIOMode(filename_);
      //
      //   // open stream
      //   std::ios_base::openmode mode_argument = std::ios_base::out;
      //   if (me2j_mode_==Me2jMode::kBinary)
      //     mode_argument |= std::ios_base::binary;
      //   stream_ptr_ = std::unique_ptr<std::ofstream>(new std::ofstream(filename_,mode_argument));
      //   mcutils::StreamCheck(bool(stream()),filename_,"Failure opening Me2j file for output");
      //
      //   // write version
      //   WriteVersion();
      //
      //   // write header
      //   if (me2j_format_==kVersion0)
      //     WriteHeader_Version0();
      //   else
      //     {
      //       std::cerr << "Unsupported version encountered when opening Me2j file " << filename_ << std::endl;
      //       std::exit(EXIT_FAILURE);
      //     }
      //   mcutils::StreamCheck(bool(stream()),filename_,"Failure while writing Me2j file header");
      // }
      //
      // void OutMe2jStream::WriteSector(
      //     std::size_t sector_index,
      //     const Eigen::MatrixXd& matrix,
      //     basis::NormalizationConversion conversion_mode
      //   )
      // {
      //   // validate matrix dimensions
      //   const auto& sector = sectors().GetSector(sector_index_);
      //   const auto& bra_subspace = sector.bra_subspace();
      //   const auto& ket_subspace = sector.ket_subspace();
      //   assert(sector_index==sector_index_);
      //   assert((matrix.rows()==sector.bra_subspace().size())&&(matrix.cols()==sector.ket_subspace().size()));
      //
      //   // validate output as NAS
      //   assert(
      //       (conversion_mode == basis::NormalizationConversion::kNone)
      //       || (conversion_mode == basis::NormalizationConversion::kASToNAS)
      //     );
      //
      //   // write sector
      //   if (me2j_format()==kVersion0)
      //     WriteSector_Version0(matrix,conversion_mode);
      //   else
      //     assert(false);  // format version was already checked when writing header
      //
      //   // validate status
      //   mcutils::StreamCheck(bool(stream()),filename_,"Failure while writing Me2j file sector");
      //
      //   // increment to next sector
      //   sector_index_++;
      //
      // }
      //
      // void OutMe2jStream::Close()
      // {
      //   stream().close();
      // };
      ////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////
    } // namespace

    ////////////////////////////////////////////////////////////////
    // version-specific implementation
    ////////////////////////////////////////////////////////////////
    //
    // When adding a new include file, also update module.mk:
    //   - Add to module_extras file list.
    //   - Add to dependency for me2j_io.o.

    #include "me2j_io_v0.cpp"

    ////////////////////////////////////////////////////////////////
    // format/mode-specific I/O: v0
    ////////////////////////////////////////////////////////////////

    #include "tbme/me2j_io.h"  // include for IDE tools

    #include "mcutils/parsing.h"
    #include "mcutils/io.h"

    namespace shell {

      void InMe2jStream::Read(int emax, Eigen::MatrixXd& matrix)
      {
        // Note: Although Me2j format Version0 only suports *diagonal* sectors
        // (operator conserves Tz, J, g), we treat the bra and ket
        // subspaces as distinct below to maintain direct analogy of code
        // with that for later Me2j formats.  We only modify this by
        // asserting that the sector is diagonal (and using the ket
        // subspace labels as "the" labels).

        int lmax=emax;


        // extract sector
        const auto& sector = sectors().GetSector(sector_index_);
        const auto& bra_subspace = sector.bra_subspace();
        const auto& ket_subspace = sector.ket_subspace();

        // std::cout << fmt::format(
        //     "sector {} type {} first {} last {}",
        //     sector_index_,basis::kTwoBodySpeciesTTzCodeDecimal[int(ket_subspace.two_body_species())],
        //     SectorIsFirstOfType(),SectorIsLastOfType()
        //   )
        //           << std::endl;

        // verify that sector is canonical
        //
        // This is a check that the caller's sector construction
        // followed the specification that only "upper triangle"
        // sectors are stored.
        assert(sector.IsUpperTriangle());

        // allocate matrix
        matrix = Eigen::MatrixXd::Zero(bra_subspace.size(),ket_subspace.size());

        // Version0 restrictions -- sector is actually diagonal (in two_body_species, J, g)
        assert(sector.IsDiagonal());
        // const int two_body_species_code = basis::kTwoBodySpeciesTTzCodeDecimal[int(ket_subspace.two_body_species())];
        // const int J = sector.ket_subspace().J();
        // const int g = sector.ket_subspace().g();

        // read FORTRAN record beginning delimiter
        if ((me2j_mode()==Me2jMode::kBinary) && SectorIsFirstOfType())
          {
            int entries = size_by_type()[int(ket_subspace.two_body_species())];
            mcutils::VerifyBinary<int>(
                stream(),entries*kIntegerSize,
                "Encountered unexpected value in Me2j file","record delimiter"
              );
          }

        // temporary buffer for binary read
        std::vector<float> buffer;
        std::size_t sector_entries = 0;
        if (me2j_mode()==Me2jMode::kBinary)
          {
            // calculate number of matrix elements in sector
            const std::size_t& dimension = ket_subspace.size();
            sector_entries = dimension*(dimension+1)/2;

            // read entire sector to temporary buffer
            buffer.resize(sector_entries, 0.);
            mcutils::ReadBinary<float>(stream(),buffer.data(),sector_entries);
          }

        // iterate over matrix elements
        std::size_t i = 0;
        for (std::size_t bra_index=0; bra_index<bra_subspace.size(); ++bra_index)
          for (std::size_t ket_index=0; ket_index<ket_subspace.size(); ++ket_index)
            {

              // diagonal sector: restrict to upper triangle
              if (sector.IsDiagonal())
                if (!(bra_index<=ket_index))
                  continue;

              // read input matrix element
              float input_matrix_element;
              if (me2j_mode()==Me2jMode::kText)
                {
                  // retrieve states
                  const basis::TwoBodyStateJJJTTz bra(bra_subspace,bra_index);
                  const basis::TwoBodyStateJJJTTz ket(ket_subspace,ket_index);

                  // input fields for text mode only
                  int input_i1, input_i2, input_i3, input_i4, input_twice_J;
                  int input_two_body_species_code;

                  // read line
                  std::string line;
                  mcutils::GetLine(stream(), line, line_count_);
                  std::istringstream line_stream(line);
                  line_stream
                    >> input_i1 >> input_i2 >> input_i3 >> input_i4
                    >> input_twice_J
                    >> input_two_body_species_code
                    >> input_matrix_element;
                  mcutils::ParsingCheck(line_stream,line_count_,line);

                  // validate input fields against expected values
                  bool inputs_as_expected = (
                      (input_i1==bra.index1()+1) && (input_i2==bra.index2()+1)
                      && (input_i3==ket.index1()+1) && (input_i4==ket.index2()+1)
                      && (input_twice_J == 2*ket.J())
                      && (input_two_body_species_code==basis::kTwoBodySpeciesTTzCodeDecimal[int(ket.two_body_species())])
                    );
                  if (!inputs_as_expected)
                    mcutils::ParsingError(line_count_,line,"Unexpected matrix element labels in input data");
                }
              else if (me2j_mode()==Me2jMode::kBinary)
                {
                  input_matrix_element = buffer[i++];
                }

              matrix(bra_index,ket_index) = input_matrix_element;
            }

        // read FORTRAN record ending delimiter
        if ((me2j_mode()==Me2jMode::kBinary) && SectorIsLastOfType())
          {
            int entries = size_by_type()[int(ket_subspace.two_body_species())];
            mcutils::VerifyBinary<int>(
                stream(),entries*kIntegerSize,
                "Encountered unexpected value in Me2j file","record delimiter"
              );
          }
      }

      // void OutMe2jStream::WriteSector_Version0(
      //     const Eigen::MatrixXd& matrix,
      //     basis::NormalizationConversion conversion_mode
      //   )
      // {
      //   // Note: Although Me2j format Version0 only suports *diagonal* sectors
      //   // (operator conserves Tz, J, g), we treat the bra and ket
      //   // subspaces as distinct below to maintain direct analogy of code
      //   // with that for later Me2j formats.  We only modify this by
      //   // asserting that the sector is diagonal (and using the ket
      //   // subspace labels as "the" labels).
      //
      //   // extract sector
      //   const auto& sector = sectors().GetSector(sector_index_);
      //   const auto& bra_subspace = sector.bra_subspace();
      //   const auto& ket_subspace = sector.ket_subspace();
      //
      //   // verify that sector is canonical
      //   //
      //   // This is a check that the caller's sector construction
      //   // followed the specification that only "upper triangle"
      //   // sectors are stored.
      //   assert(sector.IsUpperTriangle());
      //
      //   // Version0 restrictions -- sector is actually diagonal (in two_body_species, J, g)
      //   assert(sector.IsDiagonal());
      //
      //   // write FORTRAN record beginning delimiter
      //   if ((me2j_mode()==Me2jMode::kBinary) && SectorIsFirstOfType())
      //     {
      //       std::size_t entries = size_by_type()[int(ket_subspace.two_body_species())];
      //       assert(static_cast<long int>(entries*kIntegerSize) < kMaxRecordLength);
      //       mcutils::WriteBinary<int>(stream(),entries*kIntegerSize);
      //     }
      //
      //   // temporary buffer for binary write
      //   std::vector<float> buffer;
      //   std::size_t sector_entries = 0;
      //   if (me2j_mode()==Me2jMode::kBinary)
      //     {
      //       // calculate number of matrix elements in sector
      //       std::size_t dimension = ket_subspace.size();
      //       sector_entries = dimension*(dimension+1)/2;
      //
      //       // allocate buffer
      //       buffer.resize(sector_entries, 0.);
      //     }
      //
      //   // iterate over matrix elements
      //   std::size_t i = 0;
      //   for (std::size_t bra_index=0; bra_index<bra_subspace.size(); ++bra_index)
      //     for (std::size_t ket_index=0; ket_index<ket_subspace.size(); ++ket_index)
      //       {
      //
      //         // diagonal sector: restrict to upper triangle
      //         // if (sector.IsDiagonal()) // Version0 restriction -- sector is actually diagonal
      //           if (!(bra_index<=ket_index))
      //             continue;
      //
      //         // retrieve states
      //         const basis::TwoBodyStateJJJTTz bra(bra_subspace,bra_index);
      //         const basis::TwoBodyStateJJJTTz ket(ket_subspace,ket_index);
      //
      //         // determine matrix element normalization factor
      //         double conversion_factor = 1.;
      //         if (conversion_mode == basis::NormalizationConversion::kASToNAS)
      //           {
      //             if (bra.two_body_species()!=basis::TwoBodySpeciesTTz::kTTz)
      //               if (bra.index1()==bra.index2())
      //                 conversion_factor *= (1/sqrt(2.));
      //             if (ket.two_body_species()!=basis::TwoBodySpeciesTTz::kTTz)
      //               if (ket.index1()==ket.index2())
      //                 conversion_factor *= (1/sqrt(2.));
      //           }
      //
      //         // retrieve matrix element for output
      //         float output_matrix_element = conversion_factor * matrix(bra_index,ket_index);
      //
      //         if (me2j_mode()==Me2jMode::kText)
      //           {
      //             // write line: output matrix element
      //             int output_i1 = bra.index1()+1;
      //             int output_i2 = bra.index2()+1;
      //             int output_i3 = ket.index1()+1;
      //             int output_i4 = ket.index2()+1;
      //             int output_twice_J = 2*ket.J();
      //             int output_two_body_species_code=basis::kTwoBodySpeciesTTzCodeDecimal[int(ket.two_body_species())];
      //             stream()
      //               << fmt::format(
      //                   "{:3d} {:3d} {:3d} {:3d} {:3d} {:2d} {:+16.7e}",
      //                   output_i1,output_i2,output_i3,output_i4,
      //                   output_twice_J,
      //                   output_two_body_species_code,
      //                   output_matrix_element
      //                 )
      //               << std::endl;
      //           }
      //         else if (me2j_mode()==Me2jMode::kBinary)
      //           {
      //             // store to buffer
      //             buffer[i++] = output_matrix_element;
      //           }
      //       }
      //
      //   // write temporary buffer to file
      //   if (me2j_mode()==Me2jMode::kBinary)
      //     mcutils::WriteBinary<float>(stream(),buffer.data(),sector_entries);
      //
      //   // write FORTRAN record ending delimiter
      //   if ((me2j_mode()==Me2jMode::kBinary) && SectorIsLastOfType())
      //     {
      //       std::size_t entries = size_by_type()[int(ket_subspace.two_body_species())];
      //       mcutils::WriteBinary<int>(stream(),entries*kIntegerSize);
      //     }
      // }
      //
      // void TransformOperatorTwoBodyJJJTTzToTwoBodyJJJTTz(
      //     const basis::OperatorLabelsJT& operator_labels,
      //     const basis::TwoBodySpaceJJJT& two_body_jjjt_space,
      //     const std::array<basis::TwoBodySectorsJJJT,3>& two_body_jjjt_component_sectors,
      //     const std::array<basis::OperatorBlocks<double>,3>& two_body_jjjt_component_matrices,
      //     const basis::TwoBodySpaceJJJTTz& two_body_jjjpn_space,
      //     basis::TwoBodySectorsJJJTTz& two_body_jjjpn_sectors,
      //     basis::OperatorBlocks<double>& two_body_jjjpn_matrices,
      //     int Tz0
      //   )
      // {
      //   // enumerate target sectors
      //   two_body_jjjpn_sectors
      //     = basis::TwoBodySectorsJJJTTz(two_body_jjjpn_space,operator_labels.J0,operator_labels.g0,Tz0);
      //
      //   // populate matrices
      //   two_body_jjjpn_matrices.resize(two_body_jjjpn_sectors.size());
      //   for (std::size_t sector_index=0; sector_index<two_body_jjjpn_sectors.size(); ++sector_index)
      //     {
      //       // make reference to target sector
      //       const basis::TwoBodySectorsJJJTTz::SectorType& two_body_jjjpn_sector
      //         = two_body_jjjpn_sectors.GetSector(sector_index);
      //
      //       // transform
      //       Eigen::MatrixXd& matrix = two_body_jjjpn_matrices[sector_index];
      //       matrix = TwoBodyMatrixJJJTTzfromJJJTTz(
      //           operator_labels,
      //           two_body_jjjttz_space,
      //           two_body_jjjttz_component_sectors,
      //           two_body_jjjttz_component_matrices,
      //           two_body_jjjpn_sector
      //         );
      //     }
      //
      // }


    }
