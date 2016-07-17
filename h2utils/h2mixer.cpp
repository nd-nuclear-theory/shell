/******************************************************************************

  h2mixer.cpp -- linear combination of input TBMEs to yield output TBMEs in MFDn H2 format

  Syntax: h2mixer

  Standard input:
    <N1b> <N2b> (for input)
    <N1b> <N2b> (for output)

    [Then lines of the form:]
    scale <Z> <N> <hw> (reference nucleon numbers & length scale of output basis)

    in r2k2 <basename>  (sets first four input streams as r^2, ..., using basename)
    in <inbasename1>  (implied extension .bin)
    ...

    out <outbasename1>  (implied extension .bin)
    ...

    add in <in#> <multiplier> 
    ...
    add scaled <in#> <hw_ref> <degree> <multiplier> 
    add Trel   (requires in-r2k2)
    add N <hw_tilde>  <multiplier> (requires in-r2k2)
    add NCM <hw_tilde>  <multiplier> (requires in-r2k2)
    add rrel2  (requires in-r2k2)
    add R20  (requires in-r2k2)
    add K20  (requires in-r2k2)

    add am-sqr <T>
       <T>: "L", "S", or "J"
       Appends output streams for Tp^2, Tn^2, and T^2
 
  First four input streams established with "in r2k2":
    1 = <basename>-r2.bin
    2 = <basename>-r1r2.bin
    3 = <basename>-k2.bin
    4 = <basename>-k1k2.bin

  Special negative "input stream" ids (kInID_Identity, ...) are
  defined for on-the-fly operators which may be used in linear
  combinations.

  Created by M. A. Caprio, University of Notre Dame.
  3/14/12 (mac): Originated.
  4/8/12 (mac): Add hw rescaling of homogeneous operators.
  4/11/12 (mac): Fix error in scale input for N!=Z.
  6/22/12 (mac): Add R20 and K20 operators for <NCM> calculations.
  1/30/13 (mac): Update to mfdn_h2 class-based I/O.
  5/24/15 (mac): Add squared angular momentum operator support.
  Last modified 5/25/15 (mac).

******************************************************************************/


#include <shell/mfdn_h2.h>
#include <shell/mfdn_scaling.h>
#include <shell/shell_separable.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include <mcpp/profiling.h>

using namespace std;
using namespace shell;

////////////////////////////////////////////////////////////////
// run parameter storage
////////////////////////////////////////////////////////////////

// assumed input IDs
const int kInID_r2 = 1;
const int kInID_r1r2 = 2;
const int kInID_k2 = 3;
const int kInID_k1k2 = 4;

// special input IDs
const int kInID_Identity = -1;  // identity operator
const int kInID_AngularMomentum = -2;  // squared angular momentum operators from shell_separable


// MixingParameters -- stores simple parameters for run

struct MixingParameters
{
  MixingParameters();

  int N1b_in, N2b_in;
  int N1b_out, N2b_out;
  int Z, N, A;
  double hw;
  bool scales_defined, r2k2_defined;
};

// default constructor ensures sensible initialization, especially of
// boolean flags, even though instance will have local scope in main
MixingParameters::MixingParameters() : Z(0), N(0), A(0), hw(0.), scales_defined(false), r2k2_defined(false)
{
};

// MixingAmplitudes -- vector (over out#) of vectors (for given out#) of MixingAmplitude (in#,multiplier) pairs
//    note storage is 0-based, but out# is 1-based

struct MixingAmplitude
{
  MixingAmplitude (int in_id_, double multiplier_) {
    in_id = in_id_; multiplier = multiplier_;
  };
  MixingAmplitude (int in_id_, double multiplier_, const std::string identifier_) {
    in_id = in_id_; multiplier = multiplier_; identifier = identifier_;
  };
  int in_id;
  double multiplier;  // scale factor
  std::string identifier;  // to distinguish various on-the-fly operators
};

typedef vector<MixingAmplitude> MixingAmplitudes;

// input and output matrix descriptors
//    stream is pointer and allocated dynamically to prevent horrible
//    compiler errors from default/copy/assignment issues with stream member
struct InMatrixStream
{
  InMatrixStream (const string& filename_) : matrix(), stream() {
    filename=filename_;
  };
  string filename;
  InMFDnH2Stream stream;
  TwoBodyMatrixNljTzJP matrix;
};

typedef vector< InMatrixStream > InMatrixStreams;

struct OutMatrixStream
{
  OutMatrixStream (const string& filename_) : matrix(), stream() {
    filename=filename_;
  };
  string filename;
  OutMFDnH2Stream stream;
  TwoBodyMatrixNljTzJP matrix;
  MixingAmplitudes mixing_amplitudes;
};

typedef vector< OutMatrixStream > OutMatrixStreams;
	

////////////////////////////////////////////////////////////////
// parameter input from stdin
////////////////////////////////////////////////////////////////

void ReadInstructions(MixingParameters& mixing_parameters, InMatrixStreams& in_matrix_streams, OutMatrixStreams& out_matrix_streams)
{

  std::cout << std::endl;
  std::cout << "Parameter input" << std::endl;
  std::cout << std::endl;

  string line;
  int line_count = 0;

  while ( getline(std::cin, line) )
    {
      // count line
      ++line_count;

      // set up for parsing
      istringstream line_stream(line);

      if (line_count == 1)
	{
	  // fixed header line 1
	  line_stream >> mixing_parameters.N1b_in >> mixing_parameters.N2b_in;
	  ParsingCheck (line_stream, line_count, line);

	  std::cout << "Input truncation:"  << " N1b " << mixing_parameters.N1b_in << " N2b " << mixing_parameters.N2b_in << std::endl;
	}
      else if (line_count == 2)
	{
	  // fixed header line 1
	  line_stream >> mixing_parameters.N1b_out >> mixing_parameters.N2b_out;
	  ParsingCheck (line_stream, line_count, line);

	  std::cout << "Output truncation:" << " N1b " << mixing_parameters.N1b_out << " N2b " << mixing_parameters.N2b_out << std::endl;

	  // validate as subset of input truncation
	  if ( (mixing_parameters.N1b_out > mixing_parameters.N1b_in) || (mixing_parameters.N2b_out > mixing_parameters.N2b_in) )
	    ParsingError("N1b or N2b values for output exceed corresponding values for input", line_count, line);
	}
      else
	{
	  // process keyword line
	  string keyword;
	  line_stream >> keyword;
	  // no parsing check since blank line is OK
			
	  if (keyword == "c")
	    {
	      // ignore comment
	    }
	  else if (keyword == "")
	    {
	      // ignore blank line
	    }
	  else if (keyword == "scale")
	    {
	      line_stream >> mixing_parameters.Z >> mixing_parameters.N >> mixing_parameters.hw;
	      ParsingCheck (line_stream, line_count, line);

	      mixing_parameters.scales_defined = true;
	      mixing_parameters.A = mixing_parameters.Z + mixing_parameters.N;
	      std::cout << "Scales:" << "( Z " << mixing_parameters.Z << " N " << mixing_parameters.N << ")" 
		   << " A " << mixing_parameters.A << " hw " << mixing_parameters.hw << std::endl;
	    }
	  else if (keyword == "in")
	    {
	      string filename;
	      line_stream >> filename;
	      ParsingCheck (line_stream, line_count, line);
				
	      if (filename == "r2k2")
		{
		  // special case: r2k2 input basename

		  // validate as first input streams
		  if (in_matrix_streams.size() > 0)
		    ParsingError("r2k2 basename must be defined before any other input files", line_count, line);

		  string basename;
		  line_stream >> basename;
		  ParsingCheck (line_stream, line_count, line);

		  in_matrix_streams.push_back(basename + "-r2.bin");
		  in_matrix_streams.push_back(basename + "-r1r2.bin");
		  in_matrix_streams.push_back(basename + "-k2.bin");
		  in_matrix_streams.push_back(basename + "-k1k2.bin");
		  mixing_parameters.r2k2_defined = true;
		  std::cout << "Input basename (in#1-4): " << basename << std::endl;
		}
	      else
		{
		  // generic input filename
		  in_matrix_streams.push_back(filename + ".bin");
		  std::cout << "Input filename (in#" << in_matrix_streams.size() << "): " << filename << std::endl;
		}
	    }			
	  else if (keyword == "out")
	    {
	      string filename;
	      line_stream >> filename;
	      ParsingCheck (line_stream, line_count, line);
				
	      // generic output filename
	      out_matrix_streams.push_back(filename + ".bin");
	      std::cout << "Output filename (out#" << out_matrix_streams.size() << "): " << filename << std::endl;
	    }			
	  else if (keyword == "add")
	    {
	      int out_id = out_matrix_streams.size();
				
	      string subkeyword;
	      line_stream >> subkeyword;
	      ParsingCheck (line_stream, line_count, line);

	      if (subkeyword == "in")
		{
		  // from input matrix
		  int in_id;
		  double multiplier;
					
		  line_stream >> in_id >> multiplier;
		  ParsingCheck (line_stream, line_count, line);

		  if (in_id > in_matrix_streams.size())
		    ParsingError("Addition involves undefined input matrix id", line_count, line);					
		  out_matrix_streams.back().mixing_amplitudes.push_back(MixingAmplitude(in_id,multiplier));
		  std::cout << "Adding to (out#" << out_id << "): (in#" << in_id << ") * " << multiplier << std::endl;
		}
	      else if (subkeyword == "scaled")
		{
		  // from input matrix
		  int in_id;
		  double hw_ref;
		  int degree;
		  double multiplier;
					
		  line_stream >> in_id >> hw_ref >> degree >> multiplier;
		  ParsingCheck (line_stream, line_count, line);

		  // rescaling
		  double rescaling = pow( mixing_parameters.hw/hw_ref, -0.5*degree);
		  multiplier *=rescaling;

		  if (in_id > in_matrix_streams.size())
		    ParsingError("Addition involves undefined input matrix id", line_count, line);					
		  out_matrix_streams.back().mixing_amplitudes.push_back(MixingAmplitude(in_id,multiplier));
		  std::cout << "Adding to (out#" << out_id << "): (in#" << in_id << ") * " << multiplier 
		       << " after hw rescaling " << rescaling << std::endl;
		}
	      else if (subkeyword == "Trel")
		{
		  if (!mixing_parameters.r2k2_defined)
		    ParsingError("Calculation of Trel requires prior input r2k2 definition", line_count, line);

		  int in_id;
		  double multiplier;

		  // one-body term
		  in_id = kInID_k2;
		  multiplier = mixing_parameters.hw / (2 * mixing_parameters.A);
		  out_matrix_streams.back().mixing_amplitudes.push_back(MixingAmplitude(in_id,multiplier));
		  std::cout << "Adding to (out#" << out_id << "): (in#" << in_id << ") * " << multiplier 
		       << " <== Trel" << std::endl;

		  // two-body term
		  in_id = kInID_k1k2;
		  multiplier = -2 * mixing_parameters.hw / (2 * mixing_parameters.A);
		  out_matrix_streams.back().mixing_amplitudes.push_back(MixingAmplitude(in_id,multiplier));
		  std::cout << "Adding to (out#" << out_id << "): (in#" << in_id << ") * " << multiplier 
		       << " <== Trel" << std::endl;

		}
	      else if (subkeyword == "rrel2")
		{
		  if (!mixing_parameters.r2k2_defined)
		    ParsingError("Calculation of rrel2 requires prior input r2k2 definition", line_count, line);

		  int in_id;
		  double multiplier;

		  // one-body term
		  in_id = kInID_r2;
		  multiplier = pow(OscillatorLength(mixing_parameters.hw),2) / (mixing_parameters.A * mixing_parameters.A);
		  // Note: pow(mixing_parameters.A,2) is fine under gcc but requires cast under PGI
		  out_matrix_streams.back().mixing_amplitudes.push_back(MixingAmplitude(in_id,multiplier));
		  std::cout << "Adding to (out#" << out_id << "): (in#" << in_id << ") * " << multiplier 
		       << " <== rrel2" << std::endl;

		  // two-body term
		  in_id = kInID_r1r2;
		  multiplier = -2 * pow(OscillatorLength(mixing_parameters.hw),2) / (mixing_parameters.A * mixing_parameters.A);
		  out_matrix_streams.back().mixing_amplitudes.push_back(MixingAmplitude(in_id,multiplier));
		  std::cout << "Adding to (out#" << out_id << "): (in#" << in_id << ") * " << multiplier 
		       << " <== rrel2" << std::endl;
		}
	      else if (subkeyword == "N")
		{
		  if (!mixing_parameters.r2k2_defined)
		    ParsingError("Calculation of N requires prior input r2k2 definition", line_count, line);

		  double hw_tilde, postmultiplier;
		  line_stream >> hw_tilde >> postmultiplier;
		  ParsingCheck (line_stream, line_count, line);

		  int in_id;
		  double multiplier;

		  // one-body term -- for T
		  in_id = kInID_k2;
		  multiplier = postmultiplier * mixing_parameters.hw / (2 * (mixing_parameters.A - 1) * hw_tilde);
		  out_matrix_streams.back().mixing_amplitudes.push_back(MixingAmplitude(in_id,multiplier));
		  std::cout << "Adding to (out#" << out_id << "): (in#" << in_id << ") * " << multiplier 
		       << " <== N(" << hw_tilde << ") * " << postmultiplier << std::endl;

		  // one-body term -- for r^2
		  in_id = kInID_r2;
		  multiplier = postmultiplier * hw_tilde  / (2 * (mixing_parameters.A - 1) * mixing_parameters.hw);
		  out_matrix_streams.back().mixing_amplitudes.push_back(MixingAmplitude(in_id,multiplier));
		  std::cout << "Adding to (out#" << out_id << "): (in#" << in_id << ") * " << multiplier 
		       << " <== N(" << hw_tilde << ") * " << postmultiplier << std::endl;


		  // diagonal term
		  in_id = kInID_Identity;
		  multiplier = postmultiplier * -3 / (mixing_parameters.A - 1);
		  out_matrix_streams.back().mixing_amplitudes.push_back(MixingAmplitude(in_id,multiplier));
		  std::cout << "Adding to (out#" << out_id << "): (Identity) * " << multiplier
		       << " <== N(" << hw_tilde << ") * " << postmultiplier << std::endl;

		}
	      else if (subkeyword == "NCM")
		{
		  if (!mixing_parameters.r2k2_defined)
		    ParsingError("Calculation of NCM requires prior input r2k2 definition", line_count, line);

		  double hw_tilde, postmultiplier;
		  line_stream >> hw_tilde >> postmultiplier;
		  ParsingCheck (line_stream, line_count, line);

		  int in_id;
		  double multiplier;

		  // one-body term -- for TCM
		  in_id = kInID_k2;
		  multiplier = postmultiplier * mixing_parameters.hw / (2 * mixing_parameters.A * (mixing_parameters.A - 1) * hw_tilde);
		  out_matrix_streams.back().mixing_amplitudes.push_back(MixingAmplitude(in_id,multiplier));
		  std::cout << "Adding to (out#" << out_id << "): (in#" << in_id << ") * " << multiplier 
		       << " <== NCM(" << hw_tilde << ") * " << postmultiplier << std::endl;

		  // two-body term -- for TCM
		  in_id = kInID_k1k2;
		  multiplier = postmultiplier * 2 * mixing_parameters.hw / (2 * mixing_parameters.A * hw_tilde);
		  out_matrix_streams.back().mixing_amplitudes.push_back(MixingAmplitude(in_id,multiplier));
		  std::cout << "Adding to (out#" << out_id << "): (in#" << in_id << ") * " << multiplier 
		       << " <== NCM(" << hw_tilde << ") * " << postmultiplier << std::endl;

		  // one-body term -- for RCM^2
		  in_id = kInID_r2;
		  multiplier = postmultiplier * hw_tilde  / (2 * mixing_parameters.A * (mixing_parameters.A - 1) * mixing_parameters.hw);
		  out_matrix_streams.back().mixing_amplitudes.push_back(MixingAmplitude(in_id,multiplier));
		  std::cout << "Adding to (out#" << out_id << "): (in#" << in_id << ") * " << multiplier 
		       << " <== NCM(" << hw_tilde << ") * " << postmultiplier << std::endl;


		  // two-body term -- RCM^2
		  in_id = kInID_r1r2;
		  multiplier = postmultiplier * 2 * hw_tilde / (2 * mixing_parameters.A * mixing_parameters.hw);
		  out_matrix_streams.back().mixing_amplitudes.push_back(MixingAmplitude(in_id,multiplier));
		  std::cout << "Adding to (out#" << out_id << "): (in#" << in_id << ") * " << multiplier 
		       << " <== NCM(" << hw_tilde << ") * " << postmultiplier << std::endl;


		  // diagonal term
		  in_id = kInID_Identity;
		  multiplier = postmultiplier * -3 / (mixing_parameters.A * (mixing_parameters.A - 1));
		  out_matrix_streams.back().mixing_amplitudes.push_back(MixingAmplitude(in_id,multiplier));
		  std::cout << "Adding to (out#" << out_id << "): (Identity) * " << multiplier
		       << " <== NCM(" << hw_tilde << ") * " << postmultiplier << std::endl;

		}
	      else if (subkeyword == "R20")
		{
		  if (!mixing_parameters.r2k2_defined)
		    ParsingError("Calculation of rrel2 requires prior input r2k2 definition", line_count, line);

		  int in_id;
		  double multiplier;

		  // one-body term
		  in_id = kInID_r2;
		  multiplier = 1. / (mixing_parameters.A - 1);
		  out_matrix_streams.back().mixing_amplitudes.push_back(MixingAmplitude(in_id,multiplier));
		  std::cout << "Adding to (out#" << out_id << "): (in#" << in_id << ") * " << multiplier 
		       << " <== R20" << std::endl;

		  // two-body term
		  in_id = kInID_r1r2;
		  multiplier = +2;
		  out_matrix_streams.back().mixing_amplitudes.push_back(MixingAmplitude(in_id,multiplier));
		  std::cout << "Adding to (out#" << out_id << "): (in#" << in_id << ") * " << multiplier 
		       << " <== R20" << std::endl;
		}
	      else if (subkeyword == "K20")
		{
		  if (!mixing_parameters.r2k2_defined)
		    ParsingError("Calculation of rrel2 requires prior input r2k2 definition", line_count, line);

		  int in_id;
		  double multiplier;

		  // one-body term
		  in_id = kInID_k2;
		  multiplier = 1. / (mixing_parameters.A - 1);
		  out_matrix_streams.back().mixing_amplitudes.push_back(MixingAmplitude(in_id,multiplier));
		  std::cout << "Adding to (out#" << out_id << "): (in#" << in_id << ") * " << multiplier 
		       << " <== K20" << std::endl;

		  // two-body term
		  in_id = kInID_k1k2;
		  multiplier = +2;
		  out_matrix_streams.back().mixing_amplitudes.push_back(MixingAmplitude(in_id,multiplier));
		  std::cout << "Adding to (out#" << out_id << "): (in#" << in_id << ") * " << multiplier 
		       << " <== K20" << std::endl;
		}
	      else if (subkeyword == "am-sqr")
		{
		  // set up as angular momentup operator
		  int in_id = kInID_AngularMomentum;
		  double multiplier = 1.;

		  // identify which operator
		  std::string identifier;
		  line_stream >> identifier;
		  ParsingCheck (line_stream, line_count, line);

		  // store operator
		  out_matrix_streams.back().mixing_amplitudes.push_back(MixingAmplitude(in_id,multiplier,identifier));
		  std::cout << "Adding to (out#" << out_id << "): (angular momentum " << identifier << ") * " << multiplier 
		       << std::endl;
		}
	      else
		{
		  // unrecognized operator
		  ParsingError("Unrecognized operator in add instruction", line_count, line);			
		}
				
	    }			
	  else
	    {
	      // unrecognized keyword
	      ParsingError("Unrecognized keyword", line_count, line);
	    }

	}

    }

}



////////////////////////////////////////////////////////////////
// Input storage and stream initialization
////////////////////////////////////////////////////////////////

void InitializeInput (const MixingParameters& mixing_parameters, 
		      TwoBodyBasisNljTzJP& in_basis, 
		      InMatrixStreams& in_matrix_streams)
{
  std::cout << std::endl;
  std::cout << "Input initialization" << std::endl;

  // set up input basis
  in_basis.InitializeN1bN2b(mixing_parameters.N1b_in,mixing_parameters.N2b_in);

  // initialize streams
  for (InMatrixStreams::iterator it = in_matrix_streams.begin(); it != in_matrix_streams.end(); ++it)
    {
      // set up matrix 
      it->matrix.SetBasis(in_basis);

      // attempt to open file
      MFDnH2Header header;
      it->stream.Open(it->filename,header);
      it->stream.PrintDiagnostic();
		
      // validation of truncation
	
      if ( !((header.N1b == mixing_parameters.N1b_in) && (header.N2b == mixing_parameters.N2b_in)) )
	{
	  cerr << "Input file does not match expected input truncation" << std::endl;
	  exit(EXIT_FAILURE);
	}

    }

}

////////////////////////////////////////////////////////////////
// Output storage and stream initialization
////////////////////////////////////////////////////////////////

void InitializeOutput (const MixingParameters& mixing_parameters, 
		       TwoBodyBasisNljTzJP& out_basis, 
		       OutMatrixStreams& out_matrix_streams)
{
  std::cout << std::endl;
  std::cout << "Output initialization" << std::endl;

  // set up input basis
  out_basis.InitializeN1bN2b(mixing_parameters.N1b_out,mixing_parameters.N2b_out);

  // initialize streams
  for (OutMatrixStreams::iterator it = out_matrix_streams.begin(); it != out_matrix_streams.end(); ++it)
    {
      // set up matrix 
      it->matrix.SetBasis(out_basis);

      // populate header
      MFDnH2Header header;
      header.Initialize(it->matrix.GetTwoBodyBasis());

      // open file
      it->stream.Open(it->filename,header);
      it->stream.PrintDiagnostic();
		
    }

}


////////////////////////////////////////////////////////////////
// sector processing
////////////////////////////////////////////////////////////////



void ProcessMatrices(const MixingParameters& mixing_parameters, 
		     const TwoBodyBasisNljTzJP& in_basis, 
		     InMatrixStreams& in_matrix_streams,
		     const TwoBodyBasisNljTzJP& out_basis, 
		     OutMatrixStreams& out_matrix_streams)
{
  for (SectorNljTzJP sector(in_basis); sector.InRange(); ++sector)
    {
      // progress indicator
      std::cout << "." << flush;


      // recover sector properties
      const TwoSpeciesStateType state_type = sector.GetStateType();
      const int J = sector.GetJ();
      const int g = sector.GetGrade();

      // determine if corresponding output sector exists
      bool useful_sector = (J <= out_basis.JMax(state_type));

      // process sector
      if (useful_sector)
	{
	  // sector with correponding output

	  // allocate storage and read input for sector
	  for (InMatrixStreams::iterator it = in_matrix_streams.begin(); it != in_matrix_streams.end(); ++it)
	    {
	      it->matrix.Initialize(sector);
	      it->stream.ReadSector(it->matrix);
	    }

	  // process output for sector
	  for (OutMatrixStreams::iterator it = out_matrix_streams.begin(); it != out_matrix_streams.end(); ++it)
	    {
	      // allocate storage for output 
	      it->matrix.Initialize(sector);

	      // std::cout << it->filename << " " << it->mixing_amplitudes.size() << std::endl;
	      // populate matrix
	      for (MixingAmplitudes::iterator ampl_it = it->mixing_amplitudes.begin(); ampl_it != it->mixing_amplitudes.end(); ++ampl_it)
		{

		  // std::cout << it->filename << " in#" << ampl_it->in_id << " ampl " << ampl_it->multiplier << std::endl;
		  if (ampl_it->in_id == kInID_Identity)
		    {
		      // process identity 
 
		      TwoBodyMatrixSectorAddIdentity (
						      // multiplier
						      ampl_it->multiplier, 
						      // destination matrix
						      it->matrix,
						      //sector
						      sector);
		    }
		  else if (ampl_it->in_id == kInID_AngularMomentum)
		    {

		      // process angular momentum

		      AngularMomentumType angular_momentum_operator;
		      TwoSpeciesStateType operator_species;

		      // process operator identifier
		      const std::string identifier = ampl_it->identifier;
		      // identify operator type (orbital, spin, total)
		      if (identifier[0] == 'L')
			angular_momentum_operator = kOrbital;
		      else if (identifier[0] == 'S')
			angular_momentum_operator = kSpin;
		      else if (identifier[0] == 'J')
			angular_momentum_operator = kTotal;
		      // identify nucleon species (proton, neutron, total)
		      if (identifier.size() == 1)
			operator_species = kPN;
		      else if (identifier[1] == 'p')
			operator_species = kPP;
		      else if (identifier[1] == 'n')
			operator_species = kNN;
		      // std::cout << "DEBUG: identified am " << ampl_it->in_id  << " " << identifier << " " << angular_momentum_operator << " " << operator_species << std::endl;

		      // append operator to linear combination
		      TwoBodyMatrixSectorAddAngularMomentum (// multiplier
							     ampl_it->multiplier, 
							     // operator parameters
							     angular_momentum_operator,operator_species,
							     // mass number
							     mixing_parameters.A,
							     // destination matrix
							     it->matrix,
							     //sector
							     sector);
		    }
					
		  else
		    {
		      // process normal source 
		      TwoBodyMatrixSectorAdd (
					      // source matrix
					      in_matrix_streams[ampl_it->in_id-1].matrix,  // beware 1-based id
					      // multiplier
					      ampl_it->multiplier, 
					      // destination matrix
					      it->matrix,
					      //sector
					      sector);
		    }
		}
				
	      // write output
	      it->stream.WriteSector(it->matrix);

	      // free storage for output
	      it->matrix.Free(sector);

	    }


	  // free input for sector
	  for (InMatrixStreams::iterator it = in_matrix_streams.begin(); it != in_matrix_streams.end(); ++it)
	    {
	      it->matrix.Free(sector);
	    }

	}
      else
	{
	  // skip sector with no corresponding output
	  for (InMatrixStreams::iterator it = in_matrix_streams.begin(); it != in_matrix_streams.end(); ++it)
	    {
	      it->stream.SkipSector(it->matrix);
	    }
	}


    }

  // end line for progress indicator
  std::cout << std::endl;


}

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  ////////////////////////////////////////////////////////////////
  // initialization
  ////////////////////////////////////////////////////////////////

  // header
  std::cout << std::endl;
  std::cout << "h2mixer" << std::endl;
  std::cout << std::endl;

  // start timing
  Timer total_time;
  total_time.Start();

  // run data storage
  MixingParameters mixing_parameters;
  InMatrixStreams in_matrix_streams;
  OutMatrixStreams out_matrix_streams;
  TwoBodyBasisNljTzJP in_basis, out_basis;

  // parse input instructions
  ReadInstructions(mixing_parameters, in_matrix_streams, out_matrix_streams); 

  // initialize bases and streams
  InitializeInput (mixing_parameters, in_basis, in_matrix_streams);
  InitializeOutput (mixing_parameters, out_basis, out_matrix_streams);

  // sector-by-sector evaluation
  ProcessMatrices(mixing_parameters, in_basis, in_matrix_streams, out_basis, out_matrix_streams);
	

  ////////////////////////////////////////////////////////////////
  // termination
  ////////////////////////////////////////////////////////////////
	
  // note: output stream closure handled by automatic
  // destruction of output_streams vector, which calls
  // ~OutputStream for each entry, which calls delete on pointer
  // to ofstream

  // end timing
  total_time.Stop();
  std::cout << std::endl;
  std::cout << "(Total time: " << total_time.ElapsedTime() << ")" << std::endl;
  std::cout << std::endl;


  // termination
  return 0;
}
