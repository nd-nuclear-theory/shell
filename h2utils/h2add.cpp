/******************************************************************************

  h2add.cpp -- generate linear combination of MFDn TBME files

  Note: Entire input interactions are read into memory at once.  Then
  the matrix elements needed for output are cherry picked from these.
  
  Created by M. A. Caprio, University of Notre Dame.
  9/1/11 (mac): Originated.
  10/6/11 (mac): Minor update to tty output.
  3/12/12 (mac): Updated to new protocol for iteration over two-body state types.

******************************************************************************/


#include <shell/shell_2body.h>
#include <shell/mfdn_io.h>

#include <iostream>
#include <fstream>
#include <iomanip>

#include <mcpp/profiling.h>

using namespace std;
using namespace shell;

int main(int argc, char **argv)
{

	////////////////////////////////////////////////////////////////
	// initialization
	////////////////////////////////////////////////////////////////

	// header
	cout << endl;
	cout << "h2add" << endl;
	cout << endl;

 	// parameter processing
	//   Syntax: h2add infile1 infile2 scale2 outfile N1b_max N2b_max

	string is1_name(argv[1]);
	string is2_name(argv[2]);
	double scale2 = atof(argv[3]);
	string os_name(argv[4]);
	int N1b_max = atoi(argv[5]);
	int N2b_max = atoi(argv[6]);

	////////////////////////////////////////////////////////////////
	// input #1
	////////////////////////////////////////////////////////////////

	// input stream initialization
	ifstream is1(is1_name.c_str());
	if (!is1.is_open()) 
	{
		cerr << "open failed on " << is1_name << endl;
		exit(EXIT_FAILURE);
	}

	// read input header

	MFDnH2Header is1_header;
	ReadTwoBodyMatrixHeaderMFDnH2(is1, is1_header);

	// diagnostic output

	cout << "Input file #1: " << is1_name << endl;
	cout << "Nominal nuclear parameters:" << " species " << is1_header.num_types << " Z " << is1_header.num_P << " N " << is1_header.num_N << endl;
	cout << "Truncation:" << "  N1b_max " << is1_header.N1b_max << " N2b_max " << is1_header.N2b_max << endl;
	cout << "Matrix size:" << " PP "  << is1_header.matrix_size_PP << " NN " << is1_header.matrix_size_NN << " PN " << is1_header.matrix_size_PN << endl;

	// validation of truncation
	
	if ( (is1_header.N1b_max < N1b_max) || (is1_header.N2b_max < N2b_max) )
	{
		cerr << "input file does not support requested truncation" 
		     << " N1b_max " << N1b_max << " N2b_max " << N2b_max << endl;
		exit(EXIT_FAILURE);
	}

	// input basis/matrix initialization

	TwoBodyBasisNljTzJP is1_basis;
	is1_basis.InitializeN1bN2b(is1_header.N1b_max, is1_header.N2b_max);
	TwoBodyMatrixNljTzJP is1_matrix(is1_basis);
	is1_matrix.Initialize();

	// matrix input

	Timer read_time;
	read_time.Start();

	ReadTwoBodyMatrixBodyMFDnH2 (is1, is1_matrix);

	read_time.Stop();
	cout << "Read time: " << read_time.ElapsedTime() << endl;
	cout << endl;

	// file close

	is1.close();

	////////////////////////////////////////////////////////////////
	// input #2
	////////////////////////////////////////////////////////////////

	// input stream initialization
	ifstream is2(is2_name.c_str());
	if (!is2.is_open()) 
	{
		cerr << "open failed on " << is2_name << endl;
		exit(EXIT_FAILURE);
	}

	// read input header

	MFDnH2Header is2_header;
	ReadTwoBodyMatrixHeaderMFDnH2(is2, is2_header);

	// diagnostic output

	cout << "Input file #2: " << is2_name << endl;
	cout << "Nominal nuclear parameters:" << " species " << is2_header.num_types << " Z " << is2_header.num_P << " N " << is2_header.num_N << endl;
	cout << "Truncation:" << "  N1b_max " << is2_header.N1b_max << " N2b_max " << is2_header.N2b_max << endl;
	cout << "Matrix size:" << " PP "  << is2_header.matrix_size_PP << " NN " << is2_header.matrix_size_NN << " PN " << is2_header.matrix_size_PN << endl;

	// validation of truncation
	
	if ( (is2_header.N1b_max < N1b_max) || (is2_header.N2b_max < N2b_max) )
	{
		cerr << "input file does not support requested truncation" 
		     << " N1b_max " << N1b_max << " N2b_max " << N2b_max << endl;
		exit(EXIT_FAILURE);
	}

	// input basis/matrix initialization

	TwoBodyBasisNljTzJP is2_basis;
	is2_basis.InitializeN1bN2b(is2_header.N1b_max, is2_header.N2b_max);
	TwoBodyMatrixNljTzJP is2_matrix(is2_basis);
	is2_matrix.Initialize();

	// matrix input

	/* Timer read_time; */
	read_time.Start();

	ReadTwoBodyMatrixBodyMFDnH2 (is2, is2_matrix);

	read_time.Stop();
	cout << "Read time: " << read_time.ElapsedTime() << endl;
	cout << endl;

	// file close

	is2.close();

	////////////////////////////////////////////////////////////////
	// transfer
	////////////////////////////////////////////////////////////////

	// output basis/matrix initialization

	TwoBodyBasisNljTzJP output_basis;
	output_basis.InitializeN1bN2b(N1b_max, N2b_max);
	TwoBodyMatrixNljTzJP output_matrix(output_basis);
	output_matrix.Initialize();

	// copy matrix elements

	cout << "Processing..." << endl;
	Timer transformation_time;
	transformation_time.Start();

	for (TwoSpeciesStateType state_type = kTwoSpeciesStateTypeBegin; state_type <= kTwoSpeciesStateTypeEnd; ++state_type)
		for (int J = 0; J <= output_matrix.GetTwoBodyBasis().JMax(state_type); ++J)
			for (int g = 0; g <= 1; ++g)
			{
				TwoBodyMatrixSectorCopy(is1_matrix, output_matrix, state_type, J, g);
				TwoBodyMatrixSectorAdd(is2_matrix, scale2, output_matrix, state_type, J, g);
			}

	transformation_time.Stop();
	cout << "Processing time: " << transformation_time.ElapsedTime() << endl;
	cout << endl;

	////////////////////////////////////////////////////////////////
	// output
	////////////////////////////////////////////////////////////////

	// output stream initialization

	ofstream os(os_name.c_str());
	if (!os.is_open()) 
	{
		cerr << "open failed on " << os_name << endl;
		exit(EXIT_FAILURE);
	}

	// write output header

	MFDnH2Header output_header;
	output_header = is1_header;  // copy spectator fields
	output_header.N1b_max = N1b_max;
	output_header.N2b_max = N2b_max;
	output_header.matrix_size_PP = output_matrix.GetSize(kPP);
	output_header.matrix_size_PN = output_matrix.GetSize(kPN);
	output_header.matrix_size_NN = output_matrix.GetSize(kNN);
	WriteTwoBodyMatrixHeaderMFDnH2 (os, output_header);

	// diagnostic output

	cout << "Output file: " << os_name << endl;
	cout << "Nominal nuclear parameters:" << " species " << output_header.num_types << " Z " << output_header.num_P << " N " << output_header.num_N << endl;
	cout << "Truncation:" << "  N1b_max " << output_header.N1b_max << " N2b_max " << output_header.N2b_max << endl;
	cout << "Matrix size:" << " PP "  << output_header.matrix_size_PP << " NN " << output_header.matrix_size_NN << " PN " << output_header.matrix_size_PN << endl;

	// matrix output
	Timer write_time;
	write_time.Start();

	WriteTwoBodyMatrixBodyMFDnH2 (os, output_matrix);

	write_time.Stop();
	cout << "Write time: " << write_time.ElapsedTime() << endl;
	cout << endl;

	// file close
	os.close();

	// termination
	return 0;
}
