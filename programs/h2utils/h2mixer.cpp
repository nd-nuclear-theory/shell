/******************************************************************************

  h2mixer.cpp -- 

  Syntax:
    h2mixer

  Mark A. Caprio
  University of Notre Dame

  10/23/16 (mac): Created, succeeding prior incarnation originated 3/14/12.
  10/25/16 (mac): Implement direct input mixing.
  10/29/16 (mac): Implement radial operator input.

******************************************************************************/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>

#include "basis/operator.h"
#include "cppformat/format.h"
#include "mcutils/profiling.h"
#include "radial/radial_io.h"
#include "tbme/h2_io.h"
#include "tbme/separable.h"
#include "tbme/two_body_mapping.h"

////////////////////////////////////////////////////////////////
// radial operator containers
////////////////////////////////////////////////////////////////

typedef std::tuple<shell::RadialOperator,int> RadialOperatorLabels;
// Label for radial operator as (type,power) of r or k.
//
//

struct RadialOperatorData
// Indexing and matrix elements for a radial operator or radial overlaps.
//
// For overlaps, the "bra" space is the source basis, and the "ket"
// space is the target basis.
{
  RadialOperatorData(const RadialOperatorLabels& labels_, const std::string& filename_)
    : labels(labels_), filename(filename_)
  {}

  RadialOperatorLabels labels;
  std::string filename;
  basis::OrbitalSpaceLJPN bra_orbital_space;
  basis::OrbitalSpaceLJPN ket_orbital_space;
  basis::OrbitalSectorsLJPN sectors;
  basis::MatrixVector matrices;
};

typedef std::map<RadialOperatorLabels,RadialOperatorData> RadialOperatorMap;
// Map to hold all loaded radial operators.

////////////////////////////////////////////////////////////////
// two-body operator indexing container
////////////////////////////////////////////////////////////////

struct TwoBodyOperatorIndexing
// Indexing for a two-body operator.
{
  basis::OrbitalSpacePN orbital_space;
  basis::TwoBodySpaceJJJPN space;
  basis::TwoBodySectorsJJJPN sectors;
};

////////////////////////////////////////////////////////////////
// channel containers
////////////////////////////////////////////////////////////////

struct InputChannel
// Direct input channel parameters and stream.
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

enum class OperatorClass
// Classes of generated operators (distinguished within a class by
// their id).
{kIdentity,kKinematic,kAMSqr};
const std::array<const char*,3> kOperatorClassName({"identity","kinematic","am-sqr"});
std::map<std::string,OperatorClass> kOperatorClassLookup(
    {
      {"identity",OperatorClass::kIdentity},
        {"kinematic",OperatorClass::kKinematic},
          {"am-sqr",OperatorClass::kAMSqr}
    }
  );

struct OperatorChannel
// Operator channel parameters and stream.
{
  OperatorChannel(OperatorClass operator_class_, const std::string& id_)
    : id(id_), operator_class(operator_class_)
  {}

  std::string id;
  OperatorClass operator_class;
  // std::vector<double> parameters;
  // shell::TwoBodyMapping two_body_mask;
};

struct XformChannel
// Transformed input channel parameters and stream.
{
  XformChannel(const std::string& id_, const std::string& filename_)  // TODO
    : id(id_), filename(filename_), stream_ptr(NULL)
  {}

  std::string id;
  std::string filename;
  shell::InH2Stream* stream_ptr;
  basis::WeightMax weight_max;
  shell::TwoBodyMapping two_body_mapping;
  // shell::TwoBodyMapping two_body_mask;
  // POSSIBLE: independent overlaps for each xform (e.g., Coulomb)
};

struct TargetChannel
// Target channel parameters and stream.
{
  TargetChannel(const std::string& id_, const std::string& filename_)
    : id(id_), filename(filename_), stream_ptr(NULL)
  {}

  std::string id;
  std::string filename;
  shell::OutH2Stream* stream_ptr;
  std::vector<std::pair<std::string,double>> coefficients;
};

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

struct RunParameters
// Structure to store input parameters for run.
{
  RunParameters()
    : truncation_cutoff(-1), J0(0), g0(0), Tz0(0), output_h2_format(0)
  {};

  // target truncation parameters -- oscillator-like
  basis::Rank truncation_rank;
  int truncation_cutoff;  // nonnegative value flags oscillator-like truncation

  // target truncation parameters -- general
  std::string orbital_filename;
  basis::WeightMax weight_max;

  // target operator character
  int J0, g0, Tz0;

  // atomic mass number
  int A;

  // mode
  shell::H2Format output_h2_format;

};

void ReadParameters(
    RunParameters& run_parameters,
    RadialOperatorMap& radial_operators,
    std::vector<InputChannel>& input_channels,
    std::vector<OperatorChannel>& operator_channels,
    std::vector<XformChannel>& xform_channels,
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
          line_stream >> run_parameters.orbital_filename;
          ParsingCheck(line_stream,line_count,line);
        }
      else if (keyword=="set-target-truncation")
        {
          double wp, wn, wpp, wnn, wpn;
          line_stream >> wp >> wn >> wpp >> wnn >> wpn;
          ParsingCheck(line_stream,line_count,line);
          run_parameters.weight_max = basis::WeightMax(wp,wn,wpp,wnn,wpn);
        }
      else if (keyword=="set-target-oscillator")
        {
          std::string truncation_rank_code;
          line_stream >> truncation_rank_code
                      >> run_parameters.truncation_cutoff;
          ParsingCheck(line_stream,line_count,line);

          // process truncation rank code
          if (truncation_rank_code == "ob")
            run_parameters.truncation_rank = basis::Rank::kOneBody;
          else if (truncation_rank_code == "tb")
            run_parameters.truncation_rank = basis::Rank::kTwoBody;
          else
            ParsingError(line_count,line,"unrecognized truncation rank code");
          
          // store corresponding max weights
          run_parameters.weight_max
            = basis::WeightMax(run_parameters.truncation_rank,run_parameters.truncation_cutoff);
        }
      else if (keyword=="set-target-multipolarity")
        {
          line_stream >> run_parameters.J0 >> run_parameters.g0 >> run_parameters.Tz0;
          ParsingCheck(line_stream,line_count,line);
        }
      else if (keyword=="set-mass")
        {
          line_stream >> run_parameters.A;
          ParsingCheck(line_stream,line_count,line);
        }
      else if (keyword=="set-output-format")
        {
          line_stream >> run_parameters.output_h2_format;
          ParsingCheck(line_stream,line_count,line);
        }
      else if (keyword=="define-radial-operator")
        {
          std::string type_code, filename;
          int power;
          line_stream >> type_code >> power >> filename;
          ParsingCheck(line_stream,line_count,line);
          
          // cast operator type code to enum
          if (!((type_code=="r")||(type_code=="k")))
            ParsingError(line_count,line,"Unrecognized radial operator type code");
          shell::RadialOperator radial_operator_type = shell::RadialOperator(type_code[0]);

          // store radial operator information to map
          RadialOperatorLabels labels(radial_operator_type,power);
          radial_operators.insert({labels,RadialOperatorData(labels,filename)});
        }
      else if (keyword=="define-input-source")
        {
          std::string id, filename;
          line_stream >> id >> filename;
          ParsingCheck(line_stream,line_count,line);
          input_channels.emplace_back(id,filename);
        }
      else if (keyword=="define-operator-source")
        {
          std::string operator_class_name, id;
          line_stream >> operator_class_name >> id;
          ParsingCheck(line_stream,line_count,line);

          if (!kOperatorClassLookup.count(operator_class_name))
            ParsingError(line_count,line,"Unrecognized operator class");
          OperatorClass operator_class = kOperatorClassLookup[operator_class_name];
          operator_channels.emplace_back(operator_class,id);
        }
      else if (keyword=="define-xform-source")
        {
          // std::string id, filename;
          // line_stream >> id >> filename;
          // ParsingCheck(line_stream,line_count,line);
          // input_channels.emplace_back(id,filename);
        }
      else if (keyword=="define-target")
        {
          std::string id, filename;
          line_stream >> id >> filename;
          ParsingCheck(line_stream,line_count,line);
          target_channels.emplace_back(id,filename);
        }
      else if (keyword=="add-source")
        {
          assert(target_channels.size()>0);

          std::string id;
          double coefficient;
          line_stream >> id >> coefficient;
          ParsingCheck(line_stream,line_count,line);

          // tack coefficient onto coefficient list of last member of
          // target_channels
          (--target_channels.end())->coefficients.emplace_back(id,coefficient);
        }
      else
        {
          ParsingError(line_count,line,"Unrecognized keyword");
        }

    }

}

////////////////////////////////////////////////////////////////
// initialization code
/////////////////////////////////////////////////////////////////

void InitializeRadialOperators(RadialOperatorMap& radial_operators)
{
  for (auto& labels_data : radial_operators)
    {
      RadialOperatorData& radial_operator_data = labels_data.second;
      shell::RadialOperator radial_operator_type = std::get<0>(radial_operator_data.labels);
      int radial_operator_power = std::get<1>(radial_operator_data.labels);
      std::cout
        << fmt::format(
            "Reading radial operator {}^{} from {}...",
            char(radial_operator_type),radial_operator_power,
            radial_operator_data.filename
          )
        << std::endl;

      // open radial operator file
      shell::InRadialStream operator_stream(radial_operator_data.filename);
      assert(radial_operator_type==operator_stream.radial_operator());
      assert(radial_operator_power==operator_stream.sectors().l0max());
      assert(operator_stream.sectors().Tz0()==0);

      // copy out indexing
      radial_operator_data.bra_orbital_space = operator_stream.bra_orbital_space();
      radial_operator_data.ket_orbital_space = operator_stream.ket_orbital_space();
      radial_operator_data.sectors = operator_stream.sectors();

      // read matrices
      operator_stream.Read(radial_operator_data.matrices);
      
      // close file
      operator_stream.Close();

      // diagnostic
      std::cout
        << fmt::format(
            "  {} matrix elements",
            basis::AllocatedEntries(radial_operator_data.matrices)
          )
        << std::endl;
    }
}

void InitializeTargetIndexing(
    const RunParameters& run_parameters,
    TwoBodyOperatorIndexing& target_indexing
  )
{

  // set up basis
  if (run_parameters.truncation_cutoff >= 0)
    // oscillator-like indexing
    {
      target_indexing.orbital_space = basis::OrbitalSpacePN(run_parameters.truncation_cutoff);
      target_indexing.space = basis::TwoBodySpaceJJJPN(target_indexing.orbital_space,run_parameters.weight_max);
    }
  else
    // generic indexing
    {
      std::ifstream orbital_file(run_parameters.orbital_filename);
      std::vector<basis::OrbitalPNInfo> orbital_info = basis::ParseOrbitalPNStream(orbital_file,true);
      target_indexing.orbital_space = basis::OrbitalSpacePN(orbital_info);
      target_indexing.space = basis::TwoBodySpaceJJJPN(target_indexing.orbital_space,run_parameters.weight_max);
    }

  // set up sectors
  target_indexing.sectors = basis::TwoBodySectorsJJJPN(
      target_indexing.space,
      run_parameters.J0,run_parameters.g0,run_parameters.Tz0
    );
}

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
      // std::cout << fmt::format("Input channel {} index {}...",input_channel.id,input_sector_index) << std::endl;
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
      // std::cout << fmt::format("Saving {}...",input_channel.id) << std::endl;
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
          // std::cout << fmt::format(
          //     "target {} source {} coefficient {} found {}",
          //     target_channel.id,source_id,coefficient,pos != source_matrices.end()
          //   )
          //           << std::endl;

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
  RadialOperatorMap radial_operators;
  std::vector<InputChannel> input_channels;
  std::vector<OperatorChannel> operator_channels;
  std::vector<XformChannel> xform_channels;
  std::vector<TargetChannel> target_channels;
  ReadParameters(run_parameters,radial_operators,input_channels,operator_channels,xform_channels,target_channels);

  // start timing
  Timer total_time;
  total_time.Start();

  ////////////////////////////////////////////////////////////////
  // operator setup
  ////////////////////////////////////////////////////////////////

  // read in radial operators
  InitializeRadialOperators(radial_operators);

  // set up target indexing
  TwoBodyOperatorIndexing target_indexing;
  InitializeTargetIndexing(run_parameters,target_indexing);

  // set up channels
  InitializeInputChannels(run_parameters,target_indexing,input_channels);
  // InitializeOperatorChannels(run_parameters,target_indexing,operator_channels);
  // InitializeXformChannels(run_parameters,target_indexing,xform_channels);
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
  // CloseXformChannels(xform_channels);
  CloseTargetChannels(target_channels);

  // end timing
  total_time.Stop();
  std::cout << "(Total time: " << total_time.ElapsedTime() << ")" << std::endl;
  std::cout << std::endl;

  // exit
  return EXIT_SUCCESS;
}
