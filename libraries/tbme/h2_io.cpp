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

#include "cppformat/format.h"

#include "mcutils/parsing.h"

namespace shell {

  ////////////////////////////////////////////////////////////////
  // indexing
  ////////////////////////////////////////////////////////////////

  std::array<int,3> UpperTriangularEntriesByType(
      const basis::TwoBodySectorsJJJPN& sectors
    )
  // Count entries in the upper triangular portion of a set of
  // pp/nn/pn sectors.
  //
  // Based on basis::UpperTriangularEntries.
  //
  // Lower triangular sectors are ignored.  In diagonal sectors,
  // only upper-triangular entries are counted.
  //
  // Arguments:
  //   sectors (...) : container for sectors
  //
  // Returns:
  //   (std::array<int,3>) : number of upper triangular matrix entries of each type
  {
    std::array<int,3> total_entries({0,0,0});

    for (int sector_index=0; sector_index<sectors.size(); ++sector_index)
      {
        // make reference to sector for convenience
        const typename basis::TwoBodySectorsJJJPN::SectorType& sector
          = sectors.GetSector(sector_index);

        // characterize pp/nn/pn type for sector
        assert(sector.bra_subspace().two_body_species()==sector.ket_subspace().two_body_species());  // operator Tz=0
        basis::TwoBodySpeciesPN two_body_species = sector.bra_subspace().two_body_species();

        // count sector entries
        int sector_entries = 0;
        if (sector.IsDiagonal())
          // diagonal sector
          {
            int dimension = sector.ket_subspace().size();
            sector_entries = dimension*(dimension+1)/2;
          }
        else if (sector.IsUpperTriangle())
          // upper triangle sector (but not diagonal)
          {
            int bra_dimension = sector.bra_subspace().size();
            int ket_dimension = sector.ket_subspace().size();
            sector_entries = bra_dimension*ket_dimension;
          }
          
        total_entries[int(two_body_species)] += sector_entries;
      }

    return total_entries;
  }


  ////////////////////////////////////////////////////////////////
  // file binary I/O
  ////////////////////////////////////////////////////////////////
  // code assumes int and float are both 4-byte

  const int kIntegerSize = 4;
  void WriteI4 (std::ofstream& os, int i)
  
  {
    os.write(reinterpret_cast<const char*> (&i),sizeof(int));
  }
  void ReadI4 (std::ifstream& is, int &i)
  {
    is.read(reinterpret_cast<char*> (&i),sizeof(int));
  }
  void VerifyI4 (std::ifstream& is, const int i0)
  {
    int i;
    is.read(reinterpret_cast<char*> (&i),sizeof(int));
    if ( i != i0 )
      {
	// flag mismatch
	std::cerr << "H2 file I/O: Encountered input value " << i << " when expecting " << i0 << std::endl;
	exit(EXIT_FAILURE);
      }
  }

  const int kFloatSize = 4;
  void WriteFloat (std::ofstream& os, float x)
  {
    os.write(reinterpret_cast<const char*> (&x),sizeof(float));
  }

  void ReadFloat (std::ifstream& is, float &x)
  {
    is.read(reinterpret_cast<char*> (&x),sizeof(float));
  }

  // ////////////////////////////////////////////////////////////////
  // // species ordering and labeling conventions
  // ////////////////////////////////////////////////////////////////
  // // order of state types in file
  //       
  // legacy::TwoSpeciesStateType StateTypeOrderingMFDnH2 (int i)
  // {
  //   static const legacy::TwoSpeciesStateType state_type[] = {legacy::kPP, legacy::kNN, legacy::kPN};
  //   return state_type[i];
  // }
  // 
  // // state type labels for text file
  // //   NEATER: reimplement translation using static map
  // 
  // // Note: named after Lars Onsager, whose great accomplishment was
  // // (supposedly) reported in the press as having proven "H-twelve
  // // equals H-twenty-one"
  // 
  // typedef int TwoSpeciesStateTypeOnsager;
  // const TwoSpeciesStateTypeOnsager kPPOnsager = 11;
  // const TwoSpeciesStateTypeOnsager kPNOnsager = 12;
  // const TwoSpeciesStateTypeOnsager kNNOnsager = 22;
  // 
  // legacy::TwoSpeciesStateType TwoSpeciesStateTypeFromOnsager (TwoSpeciesStateTypeOnsager state_type_onsager)
  // {
  //   legacy::TwoSpeciesStateType state_type;
  //   switch (state_type_onsager)
  //     {
  //     case kPPOnsager : 
  //       state_type = legacy::kPP;
  //       break;
  //     case kPNOnsager : 
  //       state_type = legacy::kPN;
  //       break;
  //     case kNNOnsager : 
  //       state_type = legacy::kNN;
  //       break;
  //     }
  //   return state_type;
  // }
  //
  //TwoSpeciesStateTypeOnsager TwoSpeciesStateTypeToOnsager (legacy::TwoSpeciesStateType state_type)
  //{
  //  TwoSpeciesStateTypeOnsager state_type_onsager;
  //  switch (state_type)
  //    {
  //    case legacy::kPP : 
  //      state_type_onsager = kPPOnsager;
  //      break;
  //    case legacy::kPN : 
  //      state_type_onsager = kPNOnsager;
  //      break;
  //    case legacy::kNN : 
  //      state_type_onsager = kNNOnsager;
  //      break;
  //    }
  //  return state_type_onsager;
  //}

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

  std::string H2StreamBase::DiagnosticStr() const
  {

    // diagnostic output
    std::string str = fmt::format(
        "  File: {}\n"
        "  Format: {} ({})\n"
        "  Orbitals: p {} n {} oscillator-like {}\n"
        "  One-body truncation: p {:4f} n {:4f}\n"
        "  Two-body truncation: pp {:4f} nn {:4f} pn {:4f}\n"
        "  Matrix elements: pp {} nn {} pn {} => total {}\n",
        filename_,
        int(h2_format()),kH2ModeDescription[int(h2_mode())],
        orbital_space().GetSubspace(0).size(),
        orbital_space().GetSubspace(1).size(),
        999,//orbital_space().IsOscillatorLike(),
        space().weight_max().one_body[0],space().weight_max().one_body[1],
        space().weight_max().two_body[0],space().weight_max().two_body[1],space().weight_max().two_body[2],
        size_by_type()[0],size_by_type()[1],size_by_type()[2],size()
      );
    return str;
  }

  ////////////////////////////////////////////////////////////////
  // format/mode-specific I/O: version number
  ////////////////////////////////////////////////////////////////

  // read version: text mode
  void InH2Stream::ReadVersion_Text()
  {
    // set up stream alias for convenience
    std::ifstream& stream = *stream_ptr_;

    // set up line input
    std::string line;
    line_count_ = 0;  // TODO move initialization to constructor

    // line: format code
    {
      ++line_count_;
      std::getline(stream,line);
      std::istringstream line_stream(line);
      line_stream >> h2_format_; 
      ParsingCheck(line_stream,line_count_,line);
    }

  };

  // read version: binary mode
  void InH2Stream::ReadVersion_Binary()
  {

    // set up stream alias for convenience
    std::ifstream& stream = *stream_ptr_;

    // read
    int num_fields = 1;
    int bytes = num_fields * kIntegerSize;
    VerifyI4(stream,bytes);
    ReadI4(stream,h2_format_);
    VerifyI4(stream,bytes);

  };


  void OutH2Stream::WriteVersion_Text()
  {
    // set up stream alias for convenience
    std::ofstream& stream = *stream_ptr_;

    stream
      << fmt::format("{:10d}",int(h2_format())) << std::endl;
  };

  void OutH2Stream::WriteVersion_Binary()
  {
    // set up stream alias for convenience
    std::ofstream& stream = *stream_ptr_;

    // write
    int num_fields = 1;
    int bytes = num_fields * kIntegerSize;
    WriteI4(stream,bytes);
    WriteI4(stream,h2_format());
    WriteI4(stream,bytes);
  };



  ////////////////////////////////////////////////////////////////
  // format/mode-specific I/O: header
  ////////////////////////////////////////////////////////////////

  void InH2Stream::ReadHeader_Version0_Text () 
  {
    // set up stream alias for convenience
    std::ifstream& stream = *stream_ptr_;

    // set up line input
    std::string line;

    // header parameters
    int num_types, N1max, N2max, size_pp_nn, size_pn;

    // header line 1: number of particle tyles
    {
      ++line_count_;
      std::getline(stream,line);
      std::istringstream line_stream(line);
      line_stream >> num_types; 
      ParsingCheck(line_stream,line_count_,line);
    }

    // header line 2: 1-body basis limit
    {
      ++line_count_;
      std::getline(stream,line);
      std::istringstream line_stream(line);
      line_stream >> N1max; 
      ParsingCheck(line_stream,line_count_,line);
    }

    // header line 3: 2-body basis limit
    {
      ++line_count_;
      std::getline(stream,line);
      std::istringstream line_stream(line);
      line_stream >> N2max; 
      ParsingCheck(line_stream,line_count_,line);
    }
		
    // header line 4: matrix size
    {
      ++line_count_;
      std::getline(stream,line);
      std::istringstream line_stream(line);
      line_stream >> size_pp_nn >> size_pn; 
      ParsingCheck(line_stream,line_count_,line);
    }

    // process parameters
    assert(num_types == 2);
    ConstructIndexing_Version0(N1max,N2max,size_pp_nn,size_pn);
  }

  void InH2Stream::ReadHeader_Version0_Binary () 
  {
    // set up stream alias for convenience
    std::ifstream& stream = *stream_ptr_;

    // header parameters
    int num_types, N1max, N2max, size_pp_nn, size_pn;

    // read
    int num_fields = 5;
    int bytes = num_fields * kIntegerSize;
    VerifyI4(stream,bytes);
    ReadI4(stream,.num_types);
    ReadI4(stream,N1b);
    ReadI4(stream,N2b);
    ReadI4(stream,size_pp_nn);
    ReadI4(stream,size_pn);
    VerifyI4(stream,bytes);

    // process parameters
    assert(num_types == 2);
    ConstructIndexing_Version0(int N1max, int N2max, int size_pp_nn, int size_pn);
  }

  void InH2Stream::ConstructIndexing_Version0(int N1max, int N2max, int size_pp_nn, int size_pn)
  {
    // save sizes
    size_by_type_ = std::array<int,3>({size_pp_nn,size_pp_nn,size_pn});

    // construct sp orbitals 
    orbital_space_ = OrbitalSpacePN(N1max);

    // construct two-body space
    space_ = TwoBodySpaceJJJPN(orbital_space_,WeightMax(N1max,N2max));

    // construct sectors
    int J0 = 0;
    int g0 = 0;
    sectors_ = TwoBodySectorsJJJPN(space_,J0,g0);
  }


  void OutH2Stream::WriteHeader_Version0_Text()
  {

    // set up stream alias for convenience
    std::ofstream& stream = *stream_ptr_;

    // extract parameters
    //
    // assumes oscillator-like weight cutoffs
    // TODO: add tests for oscillator-like orbitals and weight cutoffs
    const int num_types = 2;
    int N1max = int(space().weight_max().one_body[0]);
    int N2max = int(space().weight_max().two_body[0]);
    int size_pp_nn = size_by_type()[0];
    int size_pn = size_by_type()[2];

    stream
      // header line 1: number of particle species
      << fmt::format("{:10d}",num_types) << std::endl
      // header line 2: 1-body basis limit
      << fmt::format("{:10d}",N1max) << std::endl
      // header line 3: 2-body basis limit
      << fmt::format("{:10d}",N2max) << std::endl
      // header line 4: matrix size
      << fmt::format("{:10d} {:10d}",size_pp_nn,size_pn) << std::endl;
    
  };

  void OutH2Stream::WriteHeader_Version0_Binary()
  {
    // set up stream alias for convenience
    std::ofstream& stream = *stream_ptr_;

    // extract parameters
    //
    // assumes oscillator-like weight cutoffs
    // TODO: add tests for oscillator-like orbitals and weight cutoffs
    const int num_types = 2;
    int N1max = int(space().weight_max().one_body[0]);
    int N2max = int(space().weight_max().two_body[0]);
    int size_pp_nn = size_by_type()[0];
    int size_pn = size_by_type()[2];

    // record size
    int num_fields = 5;
    int bytes = num_fields * kIntegerSize;

    // write
    WriteI4(stream,bytes);
    WriteI4(stream,num_types);
    WriteI4(stream,N1max);
    WriteI4(stream,N2max);
    WriteI4(stream,size_pp_nn);
    WriteI4(stream,size_pn);
    WriteI4(stream,bytes);
    
  };

  ////////////////////////////////////////////////////////////////
  // InH2Stream
  ////////////////////////////////////////////////////////////////

  void InH2Stream::Open (const std::string& filename, H2Header& header)
  {
		
    // save arguments
    filename_ = filename;

    // set text/binary mode 
    h2_mode_ = H2IOMode (filename_);

    // open stream
    stream_ = new std::ifstream;
    std::ios_base::openmode mode_argument;
    if (h2_mode_ == kText)
      mode_argument = std::ios_base::in;
    else if (h2_mode_ == kMatrix)
      mode_argument = std::ios_base::in;
    else if (h2_mode_ == kBinary)
      mode_argument = std::ios_base::in | std::ios_base::binary;
    stream_->open(filename_.c_str(),mode_argument);
    if ( !*stream_) 
      {
	std::cerr << "Open failed on H2 file " << filename_ << std::endl;
	exit(EXIT_FAILURE);
      }

    // read header
    if (h2_mode_ == kText)
      ReadHeaderText();
    else if (h2_mode_ == kMatrix)
      ReadHeaderText();
    else if (h2_mode_ == kBinary)
      ReadHeaderBin();
    if ( !*stream_ ) 
      {
	std::cerr << "Header read failed opening H2 file " << filename_ << std::endl;
	exit(EXIT_FAILURE);
      }

    // initialize sector iterator based on header
    sector_.SetTruncation(header_.N1b,header_.N2b);

    // export header
    header = header_;

  }

  void InH2Stream::Close () {
    delete stream_;
    stream_ = 0;  // so deletion by destructor does not cause error 
  };

  void InH2Stream::ReadHeaderBin () 
  {
    // read format code
    int num_format_fields = 1;
    int format_bytes = num_format_fields * kIntegerSize;
    VerifyI4(*stream_,format_bytes);
    ReadI4(*stream_,header_.format);
    VerifyI4(*stream_,format_bytes);

    // check format
    if ( header_.format != kH2FileFormat )
      {
	// flag mismatch
	std::cerr << "H2 file I/O: Encountered input format " << header_.format << " when code presently supports only " << kH2FileFormat << std::endl;
	exit(EXIT_FAILURE);
      }


    // read header body
    int num_header_fields = 5;
    int header_bytes = num_header_fields * kIntegerSize;
    VerifyI4(*stream_,header_bytes);
    ReadI4(*stream_,header_.num_types);
    ReadI4(*stream_,header_.N1b);
    ReadI4(*stream_,header_.N2b);
    ReadI4(*stream_,header_.matrix_size_PPNN);
    ReadI4(*stream_,header_.matrix_size_PN);
    VerifyI4(*stream_,header_bytes);

    // check format
    if ( header_.num_types != 2 )
      {
	// flag mismatch
	std::cerr << "H2 file I/O: Code supports only case of num_types = 2" << std::endl;
	exit(EXIT_FAILURE);
      }


  }

  void InH2Stream::ReadSector (legacy::TwoBodyMatrixNljTzJP& matrix, bool store)
  // default for store is given in class prototype
  {
    // DEBUG: std::cout << "(" << sector_.GetStateType() << sector_.GetJ() << sector_.GetGrade() << ")" << std::flush;


    // invoke basic read code
    if (h2_mode_ == kText)
      ReadSectorText (matrix,store);
    else if (h2_mode_ == kMatrix)
      ReadSectorMatrix (matrix,store);
    else if (h2_mode_ == kBinary)
      ReadSectorBin (matrix,store);

    // validate status
    if ( !*stream_ ) 
      {
	std::cerr << "Read failure on H2 file " << filename_ << std::endl;
	exit(EXIT_FAILURE);
      }

    // increment to next sector
    sector_++;

  }

  void InH2Stream::SkipSector (legacy::TwoBodyMatrixNljTzJP& matrix)
  {
    // Note that matrix object is still needed as dummy, to use ReadSector machinery, even though matrix is not used
    ReadSector (matrix, false);
  }

  void InH2Stream::ReadSectorText (legacy::TwoBodyMatrixNljTzJP& matrix, bool store)
  {
    // recover sector properties
    const legacy::TwoSpeciesStateType state_type = sector_.GetStateType();
    const int J = sector_.GetJ();
    const int g = sector_.GetGrade();
    const int dimension = matrix.GetTwoBodyBasis().GetDimension(state_type,J,g);
		
    // for canonical pairs of states (in lexicographical order)
    for (int k1 = 0; k1 < dimension; ++k1)
      for (int k2 = k1; k2 < dimension; ++k2)
	{

	  // determine expected states

	  legacy::TwoBodyStateNlj s1 = matrix.GetTwoBodyBasis().GetState(state_type,J,g,k1);
	  legacy::TwoBodyStateNlj s2 = matrix.GetTwoBodyBasis().GetState(state_type,J,g,k2);
				
	  // read input data

	  int input_i1, input_i2, input_i3, input_i4, input_twice_J;
	  int input_state_type_onsager;
	  float input_me;
				
	  std::string line;
	  std::getline(*stream_, line);
	  std::stringstream(line) >> input_i1 >> input_i2 >> input_i3 >> input_i4 
				  >> input_twice_J 
				  >> input_state_type_onsager 
				  >> input_me;
				
	  // validate input fields against expected values
				
	  if ( ! ( 
                   (input_i1 == s1.a1.GetIndex1()) && (input_i2 == s1.a2.GetIndex1())
                   && (input_i3 == s2.a1.GetIndex1()) && (input_i4 == s2.a2.GetIndex1())
                   && (input_twice_J == 2*J)
                   && (TwoSpeciesStateTypeFromOnsager(input_state_type_onsager) == state_type)
                 ) )
	    {
	      std::cerr << "ReadTwoBodyMatrixSectorMFDnH2: unexpected matrix element indices in input data" 
			<< std::endl
			<< input_i1 << " " << input_i2 << " " << input_i3 << " " << input_i4 << " " << input_twice_J << " " << input_state_type_onsager
			<< std::endl;
	      std::exit(EXIT_FAILURE);
	    }
				
	  // save matrix element

	  if (store)
	    matrix.SetMatrixElementNAS(state_type, s1, s2, input_me);
	}
  }

  void InH2Stream::ReadSectorMatrix (legacy::TwoBodyMatrixNljTzJP& matrix, bool store)
  {
    // recover sector properties
    const legacy::TwoSpeciesStateType state_type = sector_.GetStateType();
    const int J = sector_.GetJ();
    const int g = sector_.GetGrade();
    const int dimension = matrix.GetTwoBodyBasis().GetDimension(state_type,J,g);

    // verify sector header
    std::string line;
    int twice_J = 2*J;
    int input_twice_J, input_g, input_state_type_onsager, input_dimension;

    std::getline(*stream_, line);
    std::stringstream(line) >> input_twice_J 
			    >> input_g
			    >> input_state_type_onsager
			    >> input_dimension;
    if ( ! ( 
             (input_twice_J == 2*J)
             && (input_g == g)
             && (TwoSpeciesStateTypeFromOnsager(input_state_type_onsager) == state_type)
             && (input_dimension == dimension)
           ) )
      {
	std::cerr << "ReadSectorMatrix: unexpected sector header" 
		  << std::endl;
	std::exit(EXIT_FAILURE);
      }

    // for canonical pairs of states (in lexicographical order)
    for (int k1 = 0; k1 < dimension; ++k1)
      for (int k2 = k1; k2 < dimension; ++k2)
	{

	  // determine expected states

	  legacy::TwoBodyStateNlj s1 = matrix.GetTwoBodyBasis().GetState(state_type,J,g,k1);
	  legacy::TwoBodyStateNlj s2 = matrix.GetTwoBodyBasis().GetState(state_type,J,g,k2);
				
	  float input_me;
				
	  // std::string line;
	  std::getline(*stream_, line);
	  std::stringstream(line) >> input_me;
				
	  // save matrix element

	  if (store)
	    matrix.SetMatrixElementNAS(state_type, s1, s2, input_me);
	}
  }

  void InH2Stream::ReadSectorBin (legacy::TwoBodyMatrixNljTzJP& matrix, bool store)
  {
    // recover sector properties
    const legacy::TwoSpeciesStateType state_type = sector_.GetStateType();
    const int J = sector_.GetJ();
    const int g = sector_.GetGrade();
    const int dimension = matrix.GetTwoBodyBasis().GetDimension(state_type,J,g);

    // read record beginning delimiter
    if (sector_.IsFirstOfType())
      {
	int record_bytes = matrix.GetSize(sector_.GetStateType()) * kIntegerSize;
	VerifyI4(*stream_,record_bytes);
      }
		
    // for canonical pairs of states (in lexicographical order)
    for (int k1 = 0; k1 < dimension; ++k1)
      for (int k2 = k1; k2 < dimension; ++k2)
	{
	  // determine expected states

	  legacy::TwoBodyStateNlj s1 = matrix.GetTwoBodyBasis().GetState(state_type,J,g,k1);
	  legacy::TwoBodyStateNlj s2 = matrix.GetTwoBodyBasis().GetState(state_type,J,g,k2);
				
	  // read matrix element

	  float input_me;
	  ReadFloat(*stream_,input_me);

	  // save matrix element

	  if (store)
	    matrix.SetMatrixElementNAS(state_type, s1, s2, input_me);


	}

    // write record ending delimiter
    if (sector_.IsLastOfType())
      {
	int record_bytes = matrix.GetSize(sector_.GetStateType()) * kIntegerSize;
	VerifyI4(*stream_,record_bytes);
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
      H2Format h2_format, H2Mode h2_mode
    )
    : H2StreamBase(filename)
  {

    // copy in indexing
    orbital_space_ = orbital_space;
    space_ = space;
    sectors_ = sectors;
      
    // store numbers of matrix elements
    size_by_type_ = UpperTriangularEntriesByType(sectors_);
    size_ = size_by_type_[0] + size_by_type_[1] + size_by_type_[2];

    // set format and mode
    h2_format_ = h2_format;
    h2_mode_ = DeducedIOMode(filename_);

    // open stream
    stream_ = new std::ofstream;
    std::ios_base::openmode mode_argument;
    if (h2_mode_ == H2Mode::kText)
      mode_argument = std::ios_base::out;
    else if (h2_mode_ == H2Mode::kMatrix)
      mode_argument = std::ios_base::out;
    else if (h2_mode_ == H2Mode::kBinary)
      mode_argument = std::ios_base::out | std::ios_base::binary;
    stream_->open(filename_.c_str(),mode_argument);
    if ( !*stream_) 
      {
        std::cerr << "Open failed on H2 file " << filename_ << std::endl;
        exit(EXIT_FAILURE);
      }

    // write header
    if ((h2_format_==H2Format::kVersion0)&&(h2_mode_==H2Mode::kText))
      WriteHeader<H2Format::kVersion0,H2Mode::kText>(*stream_,this);
    if ( !*stream_ ) 
      {
        std::cerr << "Header write failed while opening H2 file " << filename_ << std::endl;
        exit(EXIT_FAILURE);
      }

  }
  
  void OutH2Stream::Close()
  {
    delete stream_;
    stream_ = 0;  // so deletion by destructor does not cause error 
  };
  

//  void OutH2Stream::WriteHeaderText() const
//  {
//
//    // line "0": format code
//    *stream_
//      << fmt::format("{:10d}",int(h2_format()))
//      << std::endl;
//
//    if(h2_format()==H2Format::kVersion0)
//      // kVersion0
//      {
//
//        // TODO: add tests for oscillator-like orbitals and weight cutoffs
//        // when using "old" formats
//
//        // header line 1: number of particle species
//        const int num_types = 2;
//        *stream_
//          << fmt::format("{:10d}",num_types)
//          << std::endl;
//
//        // header line 2: 1-body basis limit
//        int N1max = int(space().weight_max().one_body[0]);  // assumes oscillator-like weight cutoffs
//        *stream_
//          << fmt::format("{:10d}",N1max)
//          << std::endl;
//
//        // header line 3: 2-body basis limit
//        int N2max = int(space().weight_max().two_body[0]);  // assumes oscillator-like weight cutoffs
//        *stream_
//          << fmt::format("{:10d}",N2max)
//          << std::endl;
//
//        // header line 4: matrix size
//        int size_pp_nn = size_by_type()[0];  // assumes oscillator-like weight cutoffs
//        int size_pn = size_by_type()[2];
//        *stream_
//          << fmt::format("{:10d} {:10d}",size_pp_nn,size_pn)
//          << std::endl;
//      }
//    else
//      // fallthrough to unsupported format
//      assert(false);
//
//  }

  void OutH2Stream::WriteHeaderBin() const
  {

    // write header entries as 4-byte fields
    //    NEATEN: There must be a better way to do this???  Was so easy in C...
		
    // write format code
    int num_format_fields = 1;
    int format_bytes = num_format_fields * kIntegerSize;
    WriteI4(*stream_,format_bytes);
    WriteI4(*stream_,int(h2_format()));
    WriteI4(*stream_,format_bytes);

    if(h2_format()==H2Format::kVersion0)
      // kVersion0
      {

        // TODO: add tests for oscillator-like orbitals and weight cutoffs
        // when using "old" formats

        // write header body

        //record begin
        int num_header_fields = 5;
        int header_bytes = num_header_fields * kIntegerSize;
        WriteI4(*stream_,header_bytes);

        // header line 1: number of particle species
        const int num_types = 2;
        WriteI4(*stream_,num_types);

        // header line 2: 1-body basis limit
        int N1max = int(space().weight_max().one_body[0]);  // assumes oscillator-like weight cutoffs
        WriteI4(*stream_,N1max);

        // header line 3: 2-body basis limit
        int N2max = int(space().weight_max().two_body[0]);  // assumes oscillator-like weight cutoffs
        WriteI4(*stream_,N2max);

        // header line 4: matrix size
        int size_pp_nn = size_by_type()[0];  // assumes oscillator-like weight cutoffs
        int size_pn = size_by_type()[2];
        WriteI4(*stream_,size_pp_nn);
        WriteI4(*stream_,size_pn);
        
        // record end
        WriteI4(*stream_,header_bytes);
      }
  }


  void OutH2Stream::WriteSector(const legacy::TwoBodyMatrixNljTzJP& matrix)
  {
    // invoke basic write code
    if (h2_mode_ == kText)
      WriteSectorText(matrix);
    else if (h2_mode_ == kMatrix)
      WriteSectorMatrix(matrix);
    else if (h2_mode_ == kBinary)
      WriteSectorBin(matrix);

    // validate status
    if ( !*stream_ ) 
      {
	std::cerr << "Write failure on H2 file " << filename_ << std::endl;
	exit(EXIT_FAILURE);
      }

    // increment to next sector
    sector_++;

  }


  void OutH2Stream::WriteSectorText (const legacy::TwoBodyMatrixNljTzJP& matrix)
  {
    // recover sector properties
    const legacy::TwoSpeciesStateType state_type = sector_.GetStateType();
    const int J = sector_.GetJ();
    const int g = sector_.GetGrade();
    const int dimension = matrix.GetTwoBodyBasis().GetDimension(state_type,J,g);
		
    // output floating point format setup
    const int field_width = 3;
    const int float_precision = 7; // precision 7 for 8 digits total
    *stream_ << std::scientific << std::setprecision(float_precision);
		
    // for canonical pairs of states (in lexicographical order)
    for (int k1 = 0; k1 < dimension; ++k1)
      for (int k2 = k1; k2 < dimension; ++k2)
	{
	  legacy::TwoBodyStateNlj s1 = matrix.GetTwoBodyBasis().GetState(state_type,J,g,k1);
	  legacy::TwoBodyStateNlj s2 = matrix.GetTwoBodyBasis().GetState(state_type,J,g,k2);
				
	  int twice_J = 2*J;
				
	  *stream_ << " " << std::setw(field_width) << s1.a1.GetIndex1() 
		   << " " << std::setw(field_width) << s1.a2.GetIndex1()
		   << " " << std::setw(field_width) << s2.a1.GetIndex1() 
		   << " " << std::setw(field_width) << s2.a2.GetIndex1()
		   << " " << std::setw(field_width) << twice_J 
		   << " " << std::setw(field_width) << TwoSpeciesStateTypeToOnsager(state_type) 
		   << " " << std::setw(16) << matrix.GetMatrixElementNAS(state_type, s1, s2) 
		   << std::endl;
	}
  }

  void OutH2Stream::WriteSectorMatrix (const legacy::TwoBodyMatrixNljTzJP& matrix)
  {
    // recover sector properties
    const legacy::TwoSpeciesStateType state_type = sector_.GetStateType();
    const int J = sector_.GetJ();
    const int g = sector_.GetGrade();
    const int dimension = matrix.GetTwoBodyBasis().GetDimension(state_type,J,g);
		
    // output floating point format setup
    const int field_width = 3;
    const int float_precision = 7; // precision 7 for 8 digits total
    *stream_ << std::scientific << std::setprecision(float_precision);
		
    // write sector header
    int twice_J = 2*J;
    *stream_ << twice_J 
	     << " " << g 
	     << " " << TwoSpeciesStateTypeToOnsager(state_type) 
	     << " " << dimension
	     << std::endl;
		
    // for canonical pairs of states (in lexicographical order)
    for (int k1 = 0; k1 < dimension; ++k1)
      for (int k2 = k1; k2 < dimension; ++k2)
	{
	  legacy::TwoBodyStateNlj s1 = matrix.GetTwoBodyBasis().GetState(state_type,J,g,k1);
	  legacy::TwoBodyStateNlj s2 = matrix.GetTwoBodyBasis().GetState(state_type,J,g,k2);
				
	  *stream_ << std::setw(16) << matrix.GetMatrixElementNAS(state_type, s1, s2) 
		   << std::endl;
	}
  }

  void OutH2Stream::WriteSectorBin (const legacy::TwoBodyMatrixNljTzJP& matrix)
  {
    // recover sector properties
    const legacy::TwoSpeciesStateType state_type = sector_.GetStateType();
    const int J = sector_.GetJ();
    const int g = sector_.GetGrade();
    const int dimension = matrix.GetTwoBodyBasis().GetDimension(state_type,J,g);

    // write record beginning delimiter
    if (sector_.IsFirstOfType())
      {
	int record_bytes = matrix.GetSize(sector_.GetStateType()) * kIntegerSize;
	WriteI4(*stream_,record_bytes);
      }
		
    // for canonical pairs of states (in lexicographical order)
    for (int k1 = 0; k1 < dimension; ++k1)
      for (int k2 = k1; k2 < dimension; ++k2)
	{
	  legacy::TwoBodyStateNlj s1 = matrix.GetTwoBodyBasis().GetState(state_type,J,g,k1);
	  legacy::TwoBodyStateNlj s2 = matrix.GetTwoBodyBasis().GetState(state_type,J,g,k2);
				
	  float value = matrix.GetMatrixElementNAS(state_type, s1, s2);
	  WriteFloat(*stream_,value);

	}

    // write record ending delimiter
    if (sector_.IsLastOfType())
      {
	int record_bytes = matrix.GetSize(sector_.GetStateType()) * kIntegerSize;
	WriteI4(*stream_,record_bytes);
      }
		
  }



  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
