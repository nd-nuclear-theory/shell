/******************************************************************************

  h2write.cpp -- generate test matrix in MFDn H2 format

  Output matrix is identity matrix.

  Syntax: 

    h2write output_filename N1b_max N2b_max 
    h2write --count N1b_max N2b_max 
  
  Mark A. Caprio
  University of Notre Dame

  4/18/11 (mac): Originated.
  4/25/11 (mac): Update required headers.
  3/12/12 (mac): Update to new protocol for iteration over two-body state types.
  8/31/12 (mac): Update to mfdn_h2 class-based I/O.
  10/20/13 (mac): Add special counting-only mode.
  4/25/15 (mac): Reformat source file.
  10/11/16 (mac): Integrate into shell project.


******************************************************************************/

#include <iostream>
#include <fstream>
#include <iomanip>

#include "mcutils/profiling.h"

#include "tbme/h2_io.h"

int main(int argc, char **argv)
{

  // parameter processing
  // usage message
  if (argc == 1)
    {
      std::cout << "Syntax: h2write output_filename N1b_max N2b_max" << std::endl;
      return 0;
    }

  //   output_filename N1b_max N2b_max
  std::string os_name(argv[1]);
  int N1b_max = std::atoi(argv[2]);
  int N2b_max = std::atoi(argv[3]);

  // basis initialization
  legacy::TwoBodyBasisNljTzJP output_basis;
  output_basis.InitializeN1bN2b(N1b_max, N2b_max);
	
  // diagnostic output

  std::cout << "File: " << os_name << std::endl;
  std::cout << "Truncation: " << N1b_max << " " << N2b_max << std::endl;
  std::cout << "Matrix size: ";
  for (legacy::TwoSpeciesStateType state_type = legacy::kTwoSpeciesStateTypeBegin; state_type <= legacy::kTwoSpeciesStateTypeEnd; ++state_type)
    std::cout << legacy::TwoBodyMatrixNljTzJPDimension(output_basis, state_type)  << " ";
  std::cout << std::endl;

  // abort on special filename "--count"
  if (os_name == "--count")
    {
      std::cout << "Exiting..." << std::endl;
      std::exit(0);
    }

  // matrix initialization

  legacy::TwoBodyMatrixNljTzJP output_matrix(output_basis);

  // file header

  shell::MFDnH2Header header;
  header.Initialize(output_basis);

  // output stream initialization
  shell::OutMFDnH2Stream os;
  os.Open (os_name, header);
  os.PrintDiagnostic (std::cout);
	

  // matrix output
  Timer write_time;
  write_time.Start();
  // set matrix to identity matrix
  for (legacy::SectorNljTzJP sector(output_basis); sector.InRange(); ++sector)
    {
      output_matrix.Initialize(sector);
      legacy::TwoBodyMatrixSectorAddIdentity (1.0, output_matrix, sector);
      os.WriteSector(output_matrix);
      output_matrix.Free(sector);
	
    }
  write_time.Stop();
  std::cout << "Elapsed time: " << write_time.ElapsedTime() << std::endl;

  // file close
  os.Close();

  // termination
  return 0;
}
