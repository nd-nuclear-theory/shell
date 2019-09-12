/******************************************************************************

  h2filter.cpp -- filter MFDn H2 TBME file by keeping only the TBMEs for which
                  the number of oscillator quanta carried by the ket or bra
		  state is less than or equal to "cutoff" (the remaining TBMEs
		  are zeroed out) thus performing the interaction truncation

  Syntax:
    h2filter format input_filename output_filename

  Mark A. Caprio and Jakub Herko
  University of Notre Dame

  + 07/19 (jh): Created, based on h2conv.

******************************************************************************/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>

#include "mcutils/profiling.h"
#include "tbme/h2_io.h"

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

struct RunParameters
// Stores simple parameters for run
{
  // filenames
  std::string input_filename;
  std::string output_filename;
  // mode
  shell::H2Format output_h2_format;
};

void ProcessArguments(int argc, char **argv, RunParameters& run_parameters)
{
  // usage message
  if (argc-1 < 3)
    {
      std::cout << "Syntax: h2conv format input_filename output_filename" << std::endl;
      std::exit(EXIT_SUCCESS);
    }

  // format
  std::istringstream parameter_stream(argv[1]);
  parameter_stream >> run_parameters.output_h2_format;
  if (!parameter_stream)
    {
      std::cerr << "Expecting numeric value for H2 format argument" << std::endl;
      std::exit(EXIT_FAILURE);
    }

  // input filename
  run_parameters.input_filename = argv[2];

  // output filename
  run_parameters.output_filename = argv[3];
}

////////////////////////////////////////////////////////////////
// main program
/////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  int nmax,jj=-2,pi=-1,i1,p,i2,l1,l2,h=-1,cutoff=18; // added by J.H.

  ////////////////////////////////////////////////////////////////
  // initialization
  ////////////////////////////////////////////////////////////////

  // header
  std::cout << std::endl;
  std::cout << "h2filter -- MFDn H2 file truncation" << std::endl;
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  // read parameters
  RunParameters run_parameters;
  ProcessArguments(argc,argv,run_parameters);

  // start timing
  mcutils::SteadyTimer total_time;
  total_time.Start();

  // stream initialization

  std::cout << "Input stream" << std::endl;
  shell::InH2Stream input_stream(run_parameters.input_filename);
  std::cout << input_stream.DiagnosticStr();
  std::cout << std::endl;

  const basis::OrbitalSpacePN& orbital_space = input_stream.orbital_space();
  const basis::TwoBodySpaceJJJPN& input_space = input_stream.space();
  const basis::TwoBodySectorsJJJPN& input_sectors = input_stream.sectors();
  basis::TwoBodySpaceJJJPNOrdering space_ordering =
    shell::kH2SpaceOrdering.at(run_parameters.output_h2_format);
  const basis::TwoBodySpaceJJJPN space = basis::TwoBodySpaceJJJPN(
      orbital_space,
      input_space.weight_max(),
      space_ordering
    );
  const basis::TwoBodySectorsJJJPN sectors = basis::TwoBodySectorsJJJPN(
      space, input_sectors.J0(), input_sectors.g0(), input_sectors.Tz0()
    );

  std::cout << "Output stream" << std::endl;
  shell::OutH2Stream output_stream(
      run_parameters.output_filename,
      orbital_space,space,sectors,
      run_parameters.output_h2_format
    );
  std::cout << output_stream.DiagnosticStr();
  std::cout << std::endl;

  ////////////////////////////////////////////////////////////////
  // copy sectors
  ////////////////////////////////////////////////////////////////

  // iterate over sectors
  for (int sector_index = 0; sector_index < input_stream.num_sectors(); ++sector_index)
    {
      Eigen::MatrixXd matrix;
      // construct target sector
      const auto& sector = sectors.GetSector(sector_index);
      // locate corresponding input sector
      int input_bra_subspace_index
        = input_space.LookUpSubspaceIndex(sector.bra_subspace().labels());
      int input_ket_subspace_index
        = input_space.LookUpSubspaceIndex(sector.ket_subspace().labels());
      // Note: We cannot simply look up by target_sector's Key, since that uses target subspace indices.
      int input_sector_index
        = input_sectors.LookUpSectorIndex(input_bra_subspace_index,input_ket_subspace_index);
      input_stream.ReadSector(input_sector_index, matrix);
// beginning of a block added by J.H. performing the filtering operation
      if (sector_index == 0){
        if (matrix.rows() == 1){nmax=0;}
	else if (matrix.rows() == 4){nmax=2;}
	else if (matrix.rows() == 10){nmax=4;}
	else if (matrix.rows() == 20){nmax=6;}
	else if (matrix.rows() == 35){nmax=8;}
	else if (matrix.rows() == 56){nmax=10;}
	else if (matrix.rows() == 84){nmax=12;}
	else if (matrix.rows() == 120){nmax=14;}
	else if (matrix.rows() == 165){nmax=16;}
	else if (matrix.rows() == 220){nmax=18;}
	else if (matrix.rows() == 286){nmax=20;}
	else{std::cout << "ERROR: failed to determine Nmax" << std::endl;}
	std::cout << "CONTROL: nmax=" << nmax << std::endl;
      }
      h=h+1;
      if (2*(h/2) == h){jj=jj+2;}
      pi=-pi;
      if ((sector_index == input_stream.num_sectors()/3)||(sector_index == 2*(input_stream.num_sectors()/3))){jj=0; pi=1; h=0;}
      p=0;
      i1=0;
      for (int nk1=0; nk1<=nmax; ++nk1){
        for (int n1=nk1/2; n1>=0; --n1){
          l1=nk1-2*n1;
          for (int jj1=abs(2*l1-1); jj1<=2*l1+1; jj1=jj1+2){
            i1=i1+1;
            i2=0;
            for (int nk2=0; nk2<=nmax; ++nk2){
              for (int n2=nk2/2; n2>=0; --n2){
                l2=nk2-2*n2;
                for (int jj2=abs(2*l2-1); jj2<=2*l2+1; jj2=jj2+2){
                  i2=i2+1;
                  if (pi==1){
                    if ((abs(jj1-jj2)<=jj)&&(jj1+jj2>=jj)&&(2*((nk1+nk2)/2)==nk1+nk2)&&(nk1+nk2<=nmax)){
		      if (sector_index>=2*(input_stream.num_sectors()/3)){
		        p=p+1;
			if (nk1+nk2>cutoff){
                          for (int k=0; k<=p-1; ++k){
			    matrix(k,p-1)=0.0;
			  }
                          for (int k=p-1; k<=matrix.cols()-1; ++k){  
                            matrix(p-1,k)=0.0;
	                  }
			}
	              }
	              else if (i2>i1){
                        p=p+1;
			if (nk1+nk2>cutoff){
                          for (int k=0; k<=p-1; ++k){
		            matrix(k,p-1)=0.0;
			  }
                          for (int k=p-1; k<=matrix.cols()-1; ++k){
                            matrix(p-1,k)=0.0;
			  }
			}
	              }
                      else if ((i2==i1)&&(4*(jj/4)==jj)){
                        p=p+1;
			if (nk1+nk2>cutoff){
                          for (int k=0; k<=p-1; ++k){
                            matrix(k,p-1)=0.0;
			  }
                          for (int k=p-1; k<=matrix.cols()-1; ++k){
                            matrix(p-1,k)=0.0;
			  }
			}
		      }
		    }
		  }
                  else{
                    if ((abs(jj1-jj2)<=jj)&&(jj1+jj2>=jj)&&(2*((nk1+nk2)/2)!=nk1+nk2)&&(nk1+nk2<=nmax)){
		      if (sector_index>=2*(input_stream.num_sectors()/3)){
                        p=p+1;
			if (nk1+nk2>cutoff){
                          for (int k=0; k<=p-1; ++k){
                            matrix(k,p-1)=0.0;
			  }
                          for (int k=p-1; k<=matrix.cols()-1; ++k){
                            matrix(p-1,k)=0.0;
			  }
			}
	              }
		      else if (i2>i1){
                        p=p+1;
			if (nk1+nk2>cutoff){
                          for (int k=0; k<=p-1; ++k){
                            matrix(k,p-1)=0.0;
			  }
                          for (int k=p-1; k<=matrix.cols()-1; ++k){
                            matrix(p-1,k)=0.0;
			  }
			}
		      }
                      else if ((i2==i1)&&(4*(jj/4)==jj)){
                        p=p+1;
			if (nk1+nk2>cutoff){
                          for (int k=0; k<=p-1; ++k){
                            matrix(k,p-1)=0.0;
			  }
                          for (int k=p-1; k<=matrix.cols()-1; ++k){
                            matrix(p-1,k)=0.0;
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }

//      std::cout << "number of 2-body states for the block " << sector_index+1 << " i.e. " << jj << pi << " : " << p  << std::endl;

// end of the block added by J.H. performing the filtering operation
      output_stream.WriteSector(sector_index, matrix);
      std::cout << "." << std::flush;
    }
  std::cout << std::endl;


  ////////////////////////////////////////////////////////////////
  // termination
  ////////////////////////////////////////////////////////////////

  // explicitly force stream closure
  //   for neatness
  input_stream.Close();
  output_stream.Close();

  // end timing
  total_time.Stop();
  std::cout << "(Total time: " << total_time.ElapsedTime() << ")" << std::endl;
  std::cout << std::endl;

  // exit
  return EXIT_SUCCESS;
}
