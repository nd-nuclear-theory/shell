/******************************************************************************

  h2mixer.cpp -- 

  Syntax:
    h2mixer

  Mark A. Caprio
  University of Notre Dame

  10/23/16 (mac): Created, succeeding prior incarnation originated 3/14/12.
  10/25/16 (mac): Implement direct input mixing.

******************************************************************************/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>

#include "cppformat/format.h"
#include "mcutils/profiling.h"
#include "tbme/h2_io.h"
#include "tbme/two_body_mapping.h"

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

struct InputChannel
// Direct input stream and related parameters.
{
  InputChannel(const std::string& id_, const std::string& filename_)
    : id(id_), filename(filename_), stream_ptr(NULL)
  {}

  std::string id;
  std::string filename;
  shell::InH2Stream* stream_ptr;
  shell::TwoBodyMapping two_body_mapping;
  // shell::TwoBodyMapping two_body_mask;
};


struct TargetChannel
// Target stream and related parameters.
{
  TargetChannel(const std::string& id_, const std::string& filename_)
    : id(id_), filename(filename_), stream_ptr(NULL)
  {}

  std::string id;
  std::string filename;
  shell::OutH2Stream* stream_ptr;
  std::map<std::string,double> coefficients;
};

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

struct RunParameters
// Structure to store input parameters for run.
{
  // mode
  shell::H2Format output_h2_format;
};

void ReadParameters(
    RunParameters& run_parameters,
    std::vector<TargetChannel>& target_channels
)
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
// initialization code
/////////////////////////////////////////////////////////////////

void InitializeInputChannels(
    const RunParameters& run_parameters,
    const TwoBodyOperatorIndexing& target_indexing,
    std::vector<InputChannel>& input_channels
  )
{
  for (auto& input_channel : input_channels)
    {
      // create stream
      input_channel.stream_ptr = new shell::InH2Stream(input_channel.filename);

      // write diagnostic
      std::cout << fmt::format("Input stream -- {}",input_channel.id) << std::endl;
      std::cout << input_channel.stream_ptr->DiagnosticStr();
      std::cout << std::endl;

      // set up mapping
      input_channel.two_body_mapping = shell::TwoBodyMapping(
          input_channel.stream_ptr->orbital_space(),
          input_channel.stream_ptr->space(),
          target_indexing.orbital_space,
          target_indexing.space
        );

      // std::cout << input_channel.two_body_mapping.DebugStr() << std::endl;

    }

}

void InitializeTargetChannels(
    const RunParameters& run_parameters,
    const TwoBodyOperatorIndexing& target_indexing,
    std::vector<TargetChannel>& target_channels
  )
{
  for (auto& target_channel : target_channels)
    {
      // create stream
      target_channel.stream_ptr = new shell::OutH2Stream(
          target_channel.filename,
          target_indexing.orbital_space,target_indexing.space,target_indexing.sectors,
          run_parameters.output_h2_format
        );

      // write diagnostic
      std::cout << fmt::format("Output stream -- {}",target_channel.id) << std::endl;
      std::cout << target_channel.stream_ptr->DiagnosticStr();
      std::cout << std::endl;
    }

}

////////////////////////////////////////////////////////////////
// control code
/////////////////////////////////////////////////////////////////

void GenerateInputSources(
    std::vector<InputChannel>& input_channels,
    std::map<std::string,Eigen::MatrixXd>& source_matrices,
    const typename basis::TwoBodySectorsJJJPN::SectorType& target_sector
  )
{

  // iterate over input channels
  for (auto& input_channel : input_channels)
    {
            
      // locate corresponding input sector
      int input_bra_subspace_index
        = input_channel.stream_ptr->space().LookUpSubspaceIndex(target_sector.bra_subspace().labels());
      int input_ket_subspace_index
        = input_channel.stream_ptr->space().LookUpSubspaceIndex(target_sector.ket_subspace().labels());
      // Note: We cannot simply look up by target_sector's Key, since that uses target subspace indices.
      int input_sector_index
        = input_channel.stream_ptr->sectors().LookUpSectorIndex(input_bra_subspace_index,input_ket_subspace_index);
      std::cout << fmt::format("Input channel {} index {}...",input_channel.id,input_sector_index) << std::endl;
      if (input_sector_index == basis::kNone)
        continue;
      const typename basis::TwoBodySectorsJJJPN::SectorType& input_sector
        = input_channel.stream_ptr->sectors().GetSector(input_sector_index);
      // std::cout << "input" <<std::endl
      //           << input_sector.ket_subspace().LabelStr() <<std::endl
      //           << input_sector.ket_subspace().DebugStr() <<std::endl;
      // std::cout << "target" <<std::endl
      //           << target_sector.ket_subspace().LabelStr() <<std::endl
      //           << target_sector.ket_subspace().DebugStr() <<std::endl;

      // read matrix for sector
      Eigen::MatrixXd input_matrix;
      input_channel.stream_ptr->SeekToSector(input_sector_index);
      input_channel.stream_ptr->ReadSector(input_matrix);

      // remap input matrix to target indexing
      Eigen::MatrixXd remapped_matrix
        = Eigen::MatrixXd::Zero(target_sector.bra_subspace().size(),target_sector.ket_subspace().size());
      for (int input_bra_index=0; input_bra_index<input_sector.bra_subspace().size(); ++input_bra_index)
        for (int input_ket_index=0; input_ket_index<input_sector.ket_subspace().size(); ++input_ket_index)
          {

            // look up target matrix entry
            int remapped_bra_index
              = input_channel.two_body_mapping.state_mapping[input_sector.bra_subspace_index()][input_bra_index];
            if (remapped_bra_index == basis::kNone)
              continue;
            int remapped_ket_index
              = input_channel.two_body_mapping.state_mapping[input_sector.ket_subspace_index()][input_ket_index];
            if (remapped_ket_index == basis::kNone)
              continue;
            // std::cout
            //   << fmt::format("{} {} : {} {} / {} {} -> {} {} / {} {}",
            //                  input_sector.bra_subspace_index(),input_sector.ket_subspace_index(),
            //                  input_bra_index,input_ket_index,
            //                  input_sector.bra_subspace().size(),input_sector.ket_subspace().size(),
            //                  remapped_bra_index,remapped_ket_index,
            //                  target_sector.bra_subspace().size(),target_sector.ket_subspace().size()
            //     )
            //   << std::endl;
            remapped_matrix(remapped_bra_index,remapped_ket_index)
              = input_matrix(input_bra_index,input_ket_index);
          }

      // save result
      std::cout << fmt::format("Saving {}...",input_channel.id) << std::endl;
      source_matrices[input_channel.id] = remapped_matrix;
    }

}


void GenerateTargets(
    std::vector<TargetChannel>& target_channels,
    const std::map<std::string,Eigen::MatrixXd>& source_matrices,
    const typename basis::TwoBodySectorsJJJPN::SectorType& target_sector
  )
{

  // iterate over target channels
  for (auto& target_channel : target_channels)
    {

      // accumulate target matrix
      Eigen::MatrixXd target_matrix = Eigen::MatrixXd::Zero(target_sector.bra_subspace().size(),target_sector.ket_subspace().size());
      for (const auto& id_coefficient : target_channel.coefficients)
        {
          const std::string& source_id = id_coefficient.first;
          double coefficient = id_coefficient.second;
          
          auto pos = source_matrices.find(source_id);
          std::cout << fmt::format(
              "target {} source {} coefficient {} found {}",
              target_channel.id,source_id,coefficient,pos != source_matrices.end()
            )
                    << std::endl;

          if (pos != source_matrices.end())
            {
              const Eigen::MatrixXd& source_matrix = pos->second;
              target_matrix += coefficient * source_matrix;
            }
        }

      // write target matrix
      target_channel.stream_ptr->WriteSector(target_matrix);
    }

}

////////////////////////////////////////////////////////////////
// termination code
////////////////////////////////////////////////////////////////

void CloseInputChannels(
    std::vector<InputChannel>& input_channels
  )
{
  for (auto& input_channel : input_channels)
    {
      // close stream
      input_channel.stream_ptr->Close();

      // deallocate
      delete input_channel.stream_ptr;
    }
}

void CloseTargetChannels(
    std::vector<TargetChannel>& target_channels
  )
{
  for (auto& target_channel : target_channels)
    {
      // close stream
      target_channel.stream_ptr->Close();

      // deallocate
      delete target_channel.stream_ptr;
    }
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
  std::cout << "h2mixer -- MFDn H2 file generation" << std::endl;
  std::cout << std::endl;

  // read parameters
  RunParameters run_parameters;
  std::vector<InputChannel> input_channels;
  std::vector<TargetChannel> target_channels;
  // ReadParameters(run_parameters);

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
  target_indexing.orbital_space = basis::OrbitalSpacePN(0);
  target_indexing.space = basis::TwoBodySpaceJJJPN(
      target_indexing.orbital_space,basis::WeightMax(basis::Rank::kTwoBody,0)
    );
  target_indexing.sectors = basis::TwoBodySectorsJJJPN(target_indexing.space,0,0,0);

  // set up inputs
  input_channels.emplace_back("in1","test/test-tb-2.dat");

  // set up outputs
  target_channels.emplace_back("out1","test/out1.dat");
  target_channels[0].coefficients.insert({"in1",999});
  target_channels.emplace_back("out2","test/out2.dat");


  ////////////////////////////////////////////////////////////////
  // set up channels
  ////////////////////////////////////////////////////////////////

  InitializeInputChannels(run_parameters,target_indexing,input_channels);
  InitializeTargetChannels(run_parameters,target_indexing,target_channels);


  ////////////////////////////////////////////////////////////////
  // copy sectors
  ////////////////////////////////////////////////////////////////

  // iterate over sectors
  for (int sector_index = 0; sector_index < target_indexing.sectors.size(); ++sector_index)
    {
      // alias sector
      const auto& target_sector = target_indexing.sectors.GetSector(sector_index);

      
      // generate sources
      std::map<std::string,Eigen::MatrixXd> source_matrices;  // map id->matrix
      GenerateInputSources(input_channels,source_matrices,target_sector);

      // PopulateInputSources
      // PopulateGeneratedOperatorSources
      // PopulateXformSources


      // generate targets
      GenerateTargets(target_channels,source_matrices,target_sector);

      // progress indicator
      std::cout << "." << std::flush;

    }
  std::cout << std::endl;


  ////////////////////////////////////////////////////////////////
  // termination
  ////////////////////////////////////////////////////////////////

  // explicitly force stream closure
  //   for neatness
  CloseInputChannels(input_channels);
  CloseTargetChannels(target_channels);

  // end timing
  total_time.Stop();
  std::cout << "(Total time: " << total_time.ElapsedTime() << ")" << std::endl;
  std::cout << std::endl;

  // exit
  return EXIT_SUCCESS;
}
