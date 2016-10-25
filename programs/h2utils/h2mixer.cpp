/******************************************************************************

  h2mixer.cpp -- 

  Syntax:
    h2mixer

  Mark A. Caprio
  University of Notre Dame

  10/13/16 (mac): Created, succeeding prior incarnation originated 3/14/12.

******************************************************************************/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>

#include "cppformat/format.h"
#include "mcutils/profiling.h"
#include "tbme/h2_io.h"

////////////////////////////////////////////////////////////////
// containers
////////////////////////////////////////////////////////////////

struct TwoBodyOperatorIndexing
// Indexing for a two-body operator
{
  basis::OrbitalSpacePN orbital_space;
  basis::TwoBodySpaceJJJPN space;
  basis::TwoBodySectorsJJJPN sectors;
};

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

struct RunParameters
// Structure to store input parameters for run.
{
  // mode
  shell::H2Format output_h2_format;
  // direct input stream names 
  std::map<std::string,std::string> direct_input_source_definitions;  // map from source id to input filename
  // target stream names 
  std::map<std::string,std::string> target_filenames;  // map from target id to target filename
};

void ReadParameters(RunParameters& run_parameters)
// Parse control file.
{
  std::string line;
  int line_count = 0;
  while (std::getline(std::cin,line))
    {
      ++line_count;

      // set up for line parsing
      std::istringstream line_stream(line);
      std::string keyword;
      line_stream >> keyword;

      // skip blank line or hash comment line
      if ((keyword=="") || (keyword=="#"))
        continue;

      // select action based on keyword
      if (keyword=="set-target-orbitals")
        {
          // TODO
          //line_stream >> lambda1 >> mu1 >> lambda2 >> mu2; 
          ParsingCheck(line_stream,line_count,line);
        }
      else if (keyword=="set-target-truncation")
        {
          // process truncation rank code
          // if (truncation_rank_code == "ob")
          //   parameters.truncation_rank = basis::Rank::kOneBody;
          // else if (truncation_rank_code == "tb")
          //   parameters.truncation_rank = basis::Rank::kTwoBody;
          // else
          //   ParsingError(line_count,line,"unrecognized truncation rank code");
        }
      else if (keyword=="set-target-multipolarity")
        {
        }
      else
        {
          ParsingError(line_count,line,"Unrecognized keyword");
        }

      std::cout << std::endl;
    }

}

////////////////////////////////////////////////////////////////
// control code
/////////////////////////////////////////////////////////////////

void OpenTargetStreams(
    const RunParameters& run_parameters,
    const TwoBodyOperatorIndexing& target_indexing,
    std::map<std::string,shell::OutH2Stream> target_streams
  )
// Open target streams.
{
  for (const auto& id_filename : run_parameters.target_filenames)
    {
      // extract key and value
      const std::string& id = id_filename.first;
      const std::string& filename = id_filename.second;

      // create stream
      target_streams[id] = shell::OutH2Stream(
          filename,
          target_indexing.orbital_space,target_indexing.space,target_indexing.sectors,
          run_parameters.output_h2_format
        );

      // write diagnostic
      std::cout << fmt::format("Output stream: {}",id) << std::endl;
      std::cout << target_streams[id].DiagnosticStr();
      std::cout << std::endl;

    }

}


void MasterLoop()
{
}

////////////////////////////////////////////////////////////////
// main program
/////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  ////////////////////////////////////////////////////////////////
  // initialization
  ////////////////////////////////////////////////////////////////

  // header
  std::cout << std::endl;
  std::cout << "h2conv -- MFDn H2 file conversion" << std::endl;
  std::cout << std::endl;

  // read parameters
  RunParameters run_parameters;
  ReadParameters(run_parameters);

  // start timing
  Timer total_time;
  total_time.Start();

  // set up target indexing
  // basis::OrbitalSpacePN orbital_space;
  // basis::TwoBodySpaceJJJPN space;
  // basis::TwoBodySectorsJJJPN sectors;

  // ad hoc testing setup
  run_parameters.output_h2_format = 15099;
  TwoBodyOperatorIndexing target_indexing;
  target_indexing.orbital_space = basis::OrbitalSpacePN(2);
  target_indexing.space = basis::TwoBodySpaceJJJPN(target_indexing.orbital_space,basis::WeightMax(basis::Rank::kTwoBody,2));
  target_indexing.sectors = basis::TwoBodySectorsJJJPN(target_indexing.space,0,0,0);

  // set up outputs
  std::map<std::string,shell::OutH2Stream> target_streams;
  OpenTargetStreams(run_parameters,target_indexing,target_streams);


  ////////////////////////////////////////////////////////////////
  // copy sectors
  ////////////////////////////////////////////////////////////////

  // iterate over sectors
  for (int sector_index = 0; sector_index < target_indexing.sectors.size(); ++sector_index)
    {
      // generate sources

      // PopulateDirectInputSources
      // PopulateGeneratedOperatorSources
      // PopulateXformSources

      // generate outputs
    }
  std::cout << std::endl;


  ////////////////////////////////////////////////////////////////
  // termination
  ////////////////////////////////////////////////////////////////

  // explicitly force stream closure
  //   for neatness
  // TODO

  // end timing
  total_time.Stop();
  std::cout << "(Total time: " << total_time.ElapsedTime() << ")" << std::endl;
  std::cout << std::endl;

  // exit
  return EXIT_SUCCESS;
}
