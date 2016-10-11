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

#include "mcpp/parsing.h"

namespace shell {

  ////////////////////////////////////////////////////////////////
  // file format code
  ////////////////////////////////////////////////////////////////
  const int kH2FileFormat = 0;

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


  ////////////////////////////////////////////////////////////////
  // species ordering and labeling conventions
  ////////////////////////////////////////////////////////////////
  // order of state types in file
	
  legacy::TwoSpeciesStateType StateTypeOrderingMFDnH2 (int i)
  {
    static const legacy::TwoSpeciesStateType state_type[] = {legacy::kPP, legacy::kNN, legacy::kPN};
    return state_type[i];
  }

  // state type labels for text file
  //   NEATER: reimplement translation using static map

  // Note: named after Lars Onsager, whose great accomplishment was
  // (supposedly) reported in the press as having proven "H-twelve
  // equals H-twenty-one"

  typedef int TwoSpeciesStateTypeOnsager;
  const TwoSpeciesStateTypeOnsager kPPOnsager = 11;
  const TwoSpeciesStateTypeOnsager kPNOnsager = 12;
  const TwoSpeciesStateTypeOnsager kNNOnsager = 22;

  legacy::TwoSpeciesStateType TwoSpeciesStateTypeFromOnsager (TwoSpeciesStateTypeOnsager state_type_onsager)
  {
    legacy::TwoSpeciesStateType state_type;
    switch (state_type_onsager)
      {
      case kPPOnsager : 
	state_type = legacy::kPP;
	break;
      case kPNOnsager : 
	state_type = legacy::kPN;
	break;
      case kNNOnsager : 
	state_type = legacy::kNN;
	break;
      }
    return state_type;
  }

  TwoSpeciesStateTypeOnsager TwoSpeciesStateTypeToOnsager (legacy::TwoSpeciesStateType state_type)
  {
    TwoSpeciesStateTypeOnsager state_type_onsager;
    switch (state_type)
      {
      case legacy::kPP : 
	state_type_onsager = kPPOnsager;
	break;
      case legacy::kPN : 
	state_type_onsager = kPNOnsager;
	break;
      case legacy::kNN : 
	state_type_onsager = kNNOnsager;
	break;
      }
    return state_type_onsager;
  }

  // H2IOMode: Determines text/binary mode by filename
  //    .dat -- text format
  //    .bin -- binary format
  H2TextBinaryMode H2IOMode (const std::string& filename)
  {
    if (filename.length() < 4 )
      {
	// prevent compare on underlength string
	std::cerr << "H2 file I/O: No extension found (too short) in filename " << filename << std::endl;
	exit(EXIT_FAILURE);
      }
    else if ( ! filename.compare(filename.length()-4,4,".dat") )
      return kText;
    else if ( ! filename.compare(filename.length()-4,4,".mat") )
      return kMatrix;
    else if ( ! filename.compare(filename.length()-4,4,".bin") )
      return kBinary;
    else
      {
	std::cerr << "H2 file I/O: Extension unrecognized in filename " << filename << std::endl;
	exit(EXIT_FAILURE);
      }
  }

  ////////////////////////////////////////////////////////////////
  // OutMFDnH2Stream
  ////////////////////////////////////////////////////////////////

  // currently hard-coded for two-species (pn) work, but this could easily be adjusted

  void MFDnH2Header::Initialize(const legacy::TwoBodyBasisNljTzJP& basis)
  {
    format = kH2FileFormat;
    num_types = 2;
    N1b = basis.GetN1b(); 
    N2b = basis.GetN2b();
    matrix_size_PPNN = legacy::TwoBodyMatrixNljTzJPDimension(basis,legacy::kPP); 
    matrix_size_PN = legacy::TwoBodyMatrixNljTzJPDimension(basis,legacy::kPN);
  }

  ////////////////////////////////////////////////////////////////
  // OutMFDnH2Stream
  ////////////////////////////////////////////////////////////////

  void OutMFDnH2Stream::Open (const std::string& filename, const MFDnH2Header& header)
  {
		
    // save arguments
    filename_ = filename;
    header_ = header;

    // set text/binary mode 
    text_binary_mode_ = H2IOMode (filename_);

    // open stream
    stream_ = new std::ofstream;
    std::ios_base::openmode mode_argument;
    if (text_binary_mode_ == kText)
      mode_argument = std::ios_base::out;
    else if (text_binary_mode_ == kMatrix)
      mode_argument = std::ios_base::out;
    else if (text_binary_mode_ == kBinary)
      mode_argument = std::ios_base::out | std::ios_base::binary;
    stream_->open(filename_.c_str(),mode_argument);
    if ( !*stream_) 
      {
	std::cerr << "Open failed on H2 file " << filename_ << std::endl;
	exit(EXIT_FAILURE);
      }

    // write format and header
    if (text_binary_mode_ == kText)
      WriteHeaderText();
    else if (text_binary_mode_ == kMatrix)
      WriteHeaderText();
    else if (text_binary_mode_ == kBinary)
      WriteHeaderBin();
    if ( !*stream_ ) 
      {
	std::cerr << "Header write failed opening H2 file " << filename_ << std::endl;
	exit(EXIT_FAILURE);
      }

    // initialize sector iterator based on header
    sector_.SetTruncation(header_.N1b,header_.N2b);

  }

  void OutMFDnH2Stream::Close () {
    delete stream_;
    stream_ = 0;  // so deletion by destructor does not cause error 
  };

  void OutMFDnH2Stream::PrintDiagnostic (std::ostream& log_stream)
  {
    std::string text_binary_mode_string;
    if (text_binary_mode_ == kText)
      text_binary_mode_string = "text";
    else if (text_binary_mode_ == kMatrix)
      text_binary_mode_string = "matrix";
    else if (text_binary_mode_ == kBinary)
      text_binary_mode_string = "binary";
	
    // diagnostic output
    log_stream << "  Output file: " << filename_ << " (" << text_binary_mode_string << ")" << std::endl;
    log_stream << "  Mode: format " << header_.format << ", species " << header_.num_types << std::endl;
    log_stream << "  Truncation:" << " N1b " << header_.N1b << ", N2b " << header_.N2b << std::endl;
    log_stream << "  Matrix size:" << " PP/NN "  << header_.matrix_size_PPNN << ", PN " << header_.matrix_size_PN << std::endl;
	
  }

  void OutMFDnH2Stream::WriteHeaderText () const
  {
    const int field_width = 10;

    // format line
    *stream_ << " " << std::setw(field_width) << header_.format
	     << std::endl;

    // header line 1: particle species
    *stream_ << " " << std::setw(field_width) << header_.num_types 
	     << std::endl;

    // header line 2: 1-body basis limit
    *stream_ << " " << std::setw(field_width) << header_.N1b 
	     << std::endl;
	
    // header line 3: 2-body basis limit
    *stream_ << " " << std::setw(field_width) << header_.N2b
	     << std::endl;
	
    // header line 4: matrix size
    *stream_ << " " << std::setw(field_width) << header_.matrix_size_PPNN
	     << " " << std::setw(field_width) << header_.matrix_size_PN
	     << std::endl;
	
  }

  void OutMFDnH2Stream::WriteHeaderBin () const
  {
    // write header entries as 4-byte fields
    //    NEATEN: There must be a better way to do this???  Was so easy in C...
		
    // write format code
    int num_format_fields = 1;
    int format_bytes = num_format_fields * kIntegerSize;
    WriteI4(*stream_,format_bytes);
    WriteI4(*stream_,header_.format);
    WriteI4(*stream_,format_bytes);

    // write header body
    int num_header_fields = 5;
    int header_bytes = num_header_fields * kIntegerSize;
    WriteI4(*stream_,header_bytes);
    WriteI4(*stream_,header_.num_types);
    WriteI4(*stream_,header_.N1b);
    WriteI4(*stream_,header_.N2b);
    WriteI4(*stream_,header_.matrix_size_PPNN);
    WriteI4(*stream_,header_.matrix_size_PN);
    WriteI4(*stream_,header_bytes);

  }

  void OutMFDnH2Stream::WriteSector (const legacy::TwoBodyMatrixNljTzJP& matrix)
  {
    // invoke basic write code
    if (text_binary_mode_ == kText)
      WriteSectorText (matrix);
    else if (text_binary_mode_ == kMatrix)
      WriteSectorMatrix (matrix);
    else if (text_binary_mode_ == kBinary)
      WriteSectorBin (matrix);

    // validate status
    if ( !*stream_ ) 
      {
	std::cerr << "Write failure on H2 file " << filename_ << std::endl;
	exit(EXIT_FAILURE);
      }

    // increment to next sector
    sector_++;

  }

  void OutMFDnH2Stream::WriteSectorText (const legacy::TwoBodyMatrixNljTzJP& matrix)
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

  void OutMFDnH2Stream::WriteSectorMatrix (const legacy::TwoBodyMatrixNljTzJP& matrix)
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

  void OutMFDnH2Stream::WriteSectorBin (const legacy::TwoBodyMatrixNljTzJP& matrix)
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
  // InMFDnH2Stream
  ////////////////////////////////////////////////////////////////

  void InMFDnH2Stream::Open (const std::string& filename, MFDnH2Header& header)
  {
		
    // save arguments
    filename_ = filename;

    // set text/binary mode 
    text_binary_mode_ = H2IOMode (filename_);

    // open stream
    stream_ = new std::ifstream;
    std::ios_base::openmode mode_argument;
    if (text_binary_mode_ == kText)
      mode_argument = std::ios_base::in;
    else if (text_binary_mode_ == kMatrix)
      mode_argument = std::ios_base::in;
    else if (text_binary_mode_ == kBinary)
      mode_argument = std::ios_base::in | std::ios_base::binary;
    stream_->open(filename_.c_str(),mode_argument);
    if ( !*stream_) 
      {
	std::cerr << "Open failed on H2 file " << filename_ << std::endl;
	exit(EXIT_FAILURE);
      }

    // read header
    if (text_binary_mode_ == kText)
      ReadHeaderText();
    else if (text_binary_mode_ == kMatrix)
      ReadHeaderText();
    else if (text_binary_mode_ == kBinary)
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

  void InMFDnH2Stream::Close () {
    delete stream_;
    stream_ = 0;  // so deletion by destructor does not cause error 
  };

  void InMFDnH2Stream::PrintDiagnostic (std::ostream& log_stream)
  {
    std::string text_binary_mode_string;
    if (text_binary_mode_ == kText)
      text_binary_mode_string = "text";
    else if (text_binary_mode_ == kMatrix)
      text_binary_mode_string = "matrix";
    else if (text_binary_mode_ == kBinary)
      text_binary_mode_string = "binary";
	
    // diagnostic output
    log_stream << "  Input file: " << filename_ << " (" << text_binary_mode_string << ")" << std::endl;
    log_stream << "  Mode: format " << header_.format << ", species " << header_.num_types << std::endl;
    log_stream << "  Truncation:" << " N1b " << header_.N1b << ", N2b " << header_.N2b << std::endl;
    log_stream << "  Matrix size:" << " PP/NN "  << header_.matrix_size_PPNN << ", PN " << header_.matrix_size_PN << std::endl;
	
  }

  void InMFDnH2Stream::ReadHeaderText () 
  {
    std::string line;
		
    // format line
    std::getline(*stream_, line);
    std::stringstream(line) >> header_.format;

    // header line 1: particle number
    std::getline(*stream_, line);
    std::stringstream(line) >> header_.num_types;
		
    // header line 2: 1-body basis limit
    std::getline(*stream_, line);
    std::stringstream(line) >> header_.N1b;

    // header line 3: 2-body basis limit
    std::getline(*stream_, line);
    std::stringstream(line) >> header_.N2b;
		
    // header line 4: matrix size
    std::getline(*stream_, line);
    std::stringstream(line) >> header_.matrix_size_PPNN >> header_.matrix_size_PN;
	
  }

  void InMFDnH2Stream::ReadHeaderBin () 
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

  void InMFDnH2Stream::ReadSector (legacy::TwoBodyMatrixNljTzJP& matrix, bool store)
  // default for store is given in class prototype
  {
    // DEBUG: std::cout << "(" << sector_.GetStateType() << sector_.GetJ() << sector_.GetGrade() << ")" << std::flush;


    // invoke basic read code
    if (text_binary_mode_ == kText)
      ReadSectorText (matrix,store);
    else if (text_binary_mode_ == kMatrix)
      ReadSectorMatrix (matrix,store);
    else if (text_binary_mode_ == kBinary)
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

  void InMFDnH2Stream::SkipSector (legacy::TwoBodyMatrixNljTzJP& matrix)
  {
    // Note that matrix object is still needed as dummy, to use ReadSector machinery, even though matrix is not used
    ReadSector (matrix, false);
  }

  void InMFDnH2Stream::ReadSectorText (legacy::TwoBodyMatrixNljTzJP& matrix, bool store)
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

  void InMFDnH2Stream::ReadSectorMatrix (legacy::TwoBodyMatrixNljTzJP& matrix, bool store)
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

  void InMFDnH2Stream::ReadSectorBin (legacy::TwoBodyMatrixNljTzJP& matrix, bool store)
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
  ////////////////////////////////////////////////////////////////
} // namespace
