/******************************************************************************

  h2cp.cpp -- MFDn H2 TBME file "copy"

  Purposes:
    -- input validation
    -- text/binary conversion filter
    -- change of truncation

  Syntax: h2cp infile [ outfile [ N1b N2b ] ]

  Outline for seqential transfer of sectors:

    for each (J,Tz,g) sector defined in maximum of input and output truncations
      if given sector is defined in input truncation
        allocate sector
        read sector
      if given sector is defined in output truncation
        allocate sector
        populate sector
        write sector
      free sector (input and output)

  Mark A. Caprio
  University of Notre Dame

  4/18/11 (mac): Originated.
  4/25/11 (mac): Update required headers.
  9/1/11 (mac): Extract sector copy function to library.
  3/12/12 (mac): Update to new protocol for iteration over two-body state types.
  1/29/23 (mac): Update to mfdn_h2 class-based I/O.  Generalized master loop to allow 
    unrestricted output truncation.
  4/25/15 (mac): Reformat source file.
  10/11/16 (mac): Integrate into shell project.

******************************************************************************/


#include <iostream>
#include <iomanip>
#include <string>

#include "tbme/h2_io.h"
#include "legacy/shell_2body.h"
#include "mcutils/profiling.h"

int main(int argc, char **argv)
{

  ////////////////////////////////////////////////////////////////
  // initialization
  ////////////////////////////////////////////////////////////////

  // header
  std::cout << std::endl;
  std::cout << "h2cp" << std::endl;
  std::cout << std::endl;

  // parameter processing

  // usage message
  if (argc == 1)
    {
      std::cout << "Syntax: h2cp infile [ outfile [ N1b N2b ] ]" << std::endl;
      return 0;
    }

  // input filename
  std::string is_name(argv[1]);

  // output filename
  bool target_defined = false;
  std::string os_name;
  if (argc >= 1 + 2)
    {
      target_defined = true;
      os_name = argv[2];
    }

  // truncation
  bool truncation_defined = false;
  int N1b, N2b;
  if (argc == 1 + 4)
    {
      truncation_defined = true;
      N1b = atoi(argv[3]);
      N2b = atoi(argv[4]);
    }

  // summarize
  // std::cout << "Source: " << is_name << std::endl;
  // if (target_defined)
  // 	std::cout << "Target: " << os_name << std::endl;
  // if (truncation_defined)
  // 	std::cout << "Output truncation:" << " N1b " << N1b << " N2b " << N2b << std::endl;
  // std::cout << std::endl;

  // start timing
  Timer total_time;
  total_time.Start();

  ////////////////////////////////////////////////////////////////
  // input initialization
  ////////////////////////////////////////////////////////////////

  // input stream initialization
  shell::InMFDnH2Stream is;
  shell::MFDnH2Header is_header;
  is.Open (is_name, is_header);
  is.PrintDiagnostic ();
  std::cout << std::endl;

  // input basis/matrix initialization
  legacy::TwoBodyBasisNljTzJP input_basis;
  input_basis.InitializeN1bN2b(is_header.N1b, is_header.N2b);
  legacy::TwoBodyMatrixNljTzJP input_matrix(input_basis);

  ////////////////////////////////////////////////////////////////
  // output initialization
  ////////////////////////////////////////////////////////////////

  // provide default output truncation based on input truncation
  if (! truncation_defined)
    {
      N1b = is_header.N1b;
      N2b = is_header.N2b;
    }

  // output basis/matrix initialization
  //   only needed if target_defined but defined in any case to avoid complicated definitions due to scoping issues
  legacy::TwoBodyBasisNljTzJP output_basis;
  output_basis.InitializeN1bN2b(N1b, N2b);
  legacy::TwoBodyMatrixNljTzJP output_matrix(output_basis);

  // output stream initialization
  shell::OutMFDnH2Stream os;
  // take header contents from input
  shell::MFDnH2Header os_header; 
  os_header.Initialize(output_basis);
	
  // open output stream (if target defined)
  if (target_defined)
    {
      os.Open (os_name, os_header);
      os.PrintDiagnostic ();
      std::cout << std::endl;
    }
  else
    {
      std::cout << "  Input validation only (no target)." << std::endl;
      std::cout << std::endl;
    }


  ////////////////////////////////////////////////////////////////
  // transfer
  ////////////////////////////////////////////////////////////////

  // determine master sector iteration truncation
  int N1b_master = std::max( is_header.N1b, os_header.N1b );
  int N2b_master = std::max( is_header.N2b, os_header.N2b );

  // iterate over sectors
  for (legacy::SectorNljTzJP sector(N1b_master,N2b_master); sector.InRange(); ++sector)
    {

      // determine whether sector exists in output
      //    hence whether input data for sector needs to be stored (if available)
      bool sector_exists_in_output = target_defined && output_matrix.HasSector(sector);

      // process input
      if (input_matrix.HasSector(sector))
	{

	  // DEBUG: std::cout << " IN " << sector_exists_in_output << " " << sector.GetJ() << std::flush;
	  // std::cout << "(" << sector.GetStateType() << sector.GetJ() << sector.GetGrade() << ")" << std::flush;

	  if (sector_exists_in_output)
	    {
	      // allocate memory for sector
	      input_matrix.Initialize(sector);

	      // read sector
	      is.ReadSector(input_matrix);
	    }
	  else
	    {
	      // skip sector
	      is.SkipSector(input_matrix);
	    }
			
	  // DEBUG: std::cout << " X " << std::flush;

	}


      // process output
      if (sector_exists_in_output)
	{

	  // DEBUG: std::cout << " OUT " << std::flush;

	  // allocate memory for sector
	  output_matrix.Initialize(sector);

	  // populate matrix elements (if input exists)
	  //   else matrix elements are left undisturbed (initialized to zero)
	  if (input_matrix.HasSector(sector))
	    legacy::TwoBodyMatrixSectorCopy(input_matrix, output_matrix, sector);
			
	  // write sector
	  os.WriteSector(output_matrix);

	  // free sector
	  output_matrix.Free(sector);

	  // DEBUG: std::cout << " X " << std::flush;

	}

      // free memory for input sector (if allocated)
      if (input_matrix.HasSector(sector) && sector_exists_in_output)
	input_matrix.Free(sector);

      // progress indicator
      std::cout << "." << std::flush;
    }
	
  std::cout << std::endl;


  ////////////////////////////////////////////////////////////////
  // termination
  ////////////////////////////////////////////////////////////////

  // explicitly force closure
  //   for neatness
  is.Close();
  os.Close();

  // end timing
  total_time.Stop();
  std::cout << "(Total time: " << total_time.ElapsedTime() << ")" << std::endl;
  std::cout << std::endl;

  // termination
  return 0;
}
