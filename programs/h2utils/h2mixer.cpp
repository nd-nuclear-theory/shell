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
#include "tbme/tbme_separable.h"
#include "tbme/tbme_xform.h"
#include "tbme/two_body_mapping.h"

////////////////////////////////////////////////////////////////
// matrix diagnostic utility
////////////////////////////////////////////////////////////////

void ChopMatrix(Eigen::MatrixXd& matrix, double tolerance=1e-10)
// Round entries in matrix to zero.
{
  for (int i=0; i<matrix.rows(); ++i)
    for (int j=0; j<matrix.cols(); ++j)
      if (abs(matrix(i,j))<tolerance)
        matrix(i,j) = 0.;
}

void SymmetrizeMatrix(Eigen::MatrixXd& matrix)
// Populate lower triangle of matrix as mirror of upper triangle.
{
  for (int i=0; i<matrix.rows(); ++i)
    for (int j=0; j<matrix.cols(); ++j)
      if (i>j)
        matrix(i,j) = matrix(j,i);
}


////////////////////////////////////////////////////////////////
// radial operator containers
////////////////////////////////////////////////////////////////

typedef std::tuple<shell::RadialOperatorType,int> RadialOperatorLabels;
// Label for radial operator as (type,power) of r or k.

struct RadialOperatorData
// Indexing and matrix elements for a radial operator or radial overlaps.
//
// For overlaps, the "bra" space is the source basis, and the "ket"
// space is the target basis.
{

  RadialOperatorData() = default;

  RadialOperatorData(const RadialOperatorLabels& labels_, const std::string& filename_)
    : labels(labels_), filename(filename_)
  {}

  // operator identification
  RadialOperatorLabels labels;

  // input filename
  std::string filename;

  // operator indexing and storage
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
// Parameters and data for a source channel providing directly-input
// TBMEs.
{
  InputChannel(const std::string& id_, const std::string& filename_)
    : id(id_), filename(filename_), stream_ptr(NULL)
  {}

  std::string id;

  // input
  std::string filename;
  shell::InH2Stream* stream_ptr;

  // target mapping
  shell::TwoBodyMapping two_body_mapping;
};

struct OperatorChannel
// Parameters and data for a source channel providing a generated
// operator.
{
  OperatorChannel(const std::string& id_)
    : id(id_)
  {}

  // OperatorChannel(const std::string& id_, OperatorType operator_type_)
  //   : id(id_), operator_type(operator_type_)
  // {}

  std::string id;
  // OperatorType operator_type;
  // std::vector<double> parameters;
};

struct XformChannel
// Parameters and data for a source channel providing two-body
// transformed TBMEs.
{
  XformChannel(
      const std::string& id_, const std::string& filename_,
      const basis::WeightMax& pre_xform_weight_max_,
      const std::string& xform_filename_
    )
    : id(id_), filename(filename_),
      pre_xform_weight_max(pre_xform_weight_max_),
      xform_filename(xform_filename_),
      stream_ptr(NULL)
  {}

  std::string id;

  // input
  std::string filename;
  shell::InH2Stream* stream_ptr;

  // pre-xform storage
  basis::WeightMax pre_xform_weight_max;
  TwoBodyOperatorIndexing pre_xform_two_body_indexing;
  shell::TwoBodyMapping pre_xform_two_body_mapping;

  // xform
  std::string xform_filename;
  RadialOperatorData radial_operator_data;

  // target mapping
  shell::TwoBodyMapping two_body_mapping;
};

struct TargetChannel
// Target channel parameters and stream.
{
  TargetChannel(const std::string& id_, const std::string& filename_)
    : id(id_), filename(filename_), stream_ptr(NULL)
  {}

  std::string id;

  // construction
  std::vector<std::pair<std::string,double>> coefficients;

  // output
  std::string filename;
  shell::OutH2Stream* stream_ptr;
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
      if (keyword=="set-target-indexing")
        {
          double wp, wn, wpp, wnn, wpn;
          line_stream
            >> run_parameters.orbital_filename
            >> wp >> wn >> wpp >> wnn >> wpn;
          ParsingCheck(line_stream,line_count,line);
          run_parameters.weight_max = basis::WeightMax(wp,wn,wpp,wnn,wpn);
        }
      else if (keyword=="set-target-indexing-oscillator")
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
          shell::RadialOperatorType radial_operator_type = shell::RadialOperatorType(type_code[0]);

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
          std::string id;
          line_stream >> id;
          ParsingCheck(line_stream,line_count,line);

          // if (kOperatorTypeLookup.count(id)!=1)
          //   ParsingError(line_count,line,"Unrecognized operator ID");
          // OperatorType operator_type = kOperatorTypeLookup[id];
          // operator_channels.emplace_back(id,operator_type);

          operator_channels.emplace_back(id);
        }
      else if (keyword=="define-xform-source")
        {
          std::string id, filename, xform_filename;
          double wp, wn, wpp, wnn, wpn;
          line_stream
            >> id >> filename
            >> wp >> wn >> wpp >> wnn >> wpn
            >> xform_filename;
          ParsingCheck(line_stream,line_count,line);

          xform_channels.emplace_back(
              id,filename,basis::WeightMax(wp,wn,wpp,wnn,wpn),xform_filename
            );
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
      // extract labeling
      RadialOperatorData& radial_operator_data = labels_data.second;
      shell::RadialOperatorType radial_operator_type = std::get<0>(radial_operator_data.labels);
      int radial_operator_power = std::get<1>(radial_operator_data.labels);
      std::cout
        << fmt::format(
            "Reading radial operator {}^{} from {}...",
            char(radial_operator_type),radial_operator_power,
            radial_operator_data.filename
          )
        << std::endl;

      // open radial operator file
      shell::InRadialStream radial_operator_stream(radial_operator_data.filename);
      assert(radial_operator_type==radial_operator_stream.radial_operator_type());
      assert(radial_operator_power==radial_operator_stream.sectors().l0max());
      assert(radial_operator_stream.sectors().Tz0()==0);

      // copy out indexing
      radial_operator_data.bra_orbital_space = radial_operator_stream.bra_orbital_space();
      radial_operator_data.ket_orbital_space = radial_operator_stream.ket_orbital_space();
      radial_operator_data.sectors = radial_operator_stream.sectors();

      // read matrices
      radial_operator_stream.Read(radial_operator_data.matrices);
      
      // close file
      radial_operator_stream.Close();

      // diagnostic
      std::cout
        << fmt::format(
            "  {} matrix elements",
            basis::AllocatedEntries(radial_operator_data.matrices)
          )
        << std::endl;
    }

  std::cout << std::endl;
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
      std::cout << fmt::format("Input stream (direct) -- {}",input_channel.id) << std::endl;
      std::cout << input_channel.stream_ptr->DiagnosticStr();
      std::cout << std::endl;

      // set up target mapping
      input_channel.two_body_mapping = shell::TwoBodyMapping(
          input_channel.stream_ptr->orbital_space(),
          input_channel.stream_ptr->space(),
          target_indexing.orbital_space,
          target_indexing.space
        );

      // std::cout << input_channel.two_body_mapping.DebugStr() << std::endl;

    }
}

void InitializeXformChannels(
    const RunParameters& run_parameters,
    const TwoBodyOperatorIndexing& target_indexing,
    std::vector<XformChannel>& xform_channels
  )
{
  for (auto& xform_channel : xform_channels)
    {
      // create stream
      xform_channel.stream_ptr = new shell::InH2Stream(xform_channel.filename);

      // write diagnostic
      std::cout << fmt::format("Input stream (xform) -- {}",xform_channel.id) << std::endl;
      std::cout << xform_channel.stream_ptr->DiagnosticStr();
      std::cout << std::endl;

      // set up pre-xform indexing and mapping
      //
      // preserves input orbitals but builds two-body space with
      // user-specified truncation
      xform_channel.pre_xform_two_body_indexing.orbital_space = xform_channel.stream_ptr->orbital_space();
      xform_channel.pre_xform_two_body_indexing.space
        = basis::TwoBodySpaceJJJPN(
            xform_channel.pre_xform_two_body_indexing.orbital_space,
            xform_channel.pre_xform_weight_max
          );
      xform_channel.pre_xform_two_body_indexing.sectors = basis::TwoBodySectorsJJJPN(
          xform_channel.pre_xform_two_body_indexing.space,
          run_parameters.J0,run_parameters.g0,run_parameters.Tz0
        );
      xform_channel.pre_xform_two_body_mapping = shell::TwoBodyMapping(
          xform_channel.stream_ptr->orbital_space(),
          xform_channel.stream_ptr->space(),
          xform_channel.pre_xform_two_body_indexing.orbital_space,
          xform_channel.pre_xform_two_body_indexing.space
        );

      // read xform coefficients
      //
      // analogous to InitializeRadialOperators

      std::cout
        << fmt::format("Reading radial xform overlaps from {}...",xform_channel.xform_filename)
        << std::endl;

      // set up alias (to keep code parallel to InitializeRadialOperators)
      RadialOperatorData& radial_operator_data = xform_channel.radial_operator_data;

      // set up radial operator data
      shell::RadialOperatorType radial_operator_type = shell::RadialOperatorType::kR;
      int radial_operator_power = 0;
      radial_operator_data = 
        RadialOperatorData(RadialOperatorLabels(radial_operator_type,radial_operator_power),xform_channel.xform_filename);

      // open radial operator file
      shell::InRadialStream radial_operator_stream(radial_operator_data.filename);
      assert(radial_operator_type==radial_operator_stream.radial_operator_type());
      assert(radial_operator_power==radial_operator_stream.sectors().l0max());
      assert(radial_operator_stream.sectors().Tz0()==0);

      // copy out indexing
      radial_operator_data.bra_orbital_space = radial_operator_stream.bra_orbital_space();
      radial_operator_data.ket_orbital_space = radial_operator_stream.ket_orbital_space();
      radial_operator_data.sectors = radial_operator_stream.sectors();

      // read matrices
      radial_operator_stream.Read(radial_operator_data.matrices);
      
      // close file
      radial_operator_stream.Close();

      // diagnostic
      std::cout
        << fmt::format("  {} matrix elements",basis::AllocatedEntries(radial_operator_data.matrices))
        << std::endl;


      // set up target mapping
      xform_channel.two_body_mapping = shell::TwoBodyMapping(
          xform_channel.pre_xform_two_body_indexing.orbital_space,
          xform_channel.pre_xform_two_body_indexing.space,
          target_indexing.orbital_space,
          target_indexing.space
        );
      // std::cout << input_channel.two_body_mapping.DebugStr() << std::endl;

      std::cout << std::endl;
        

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

  // iterate over channels
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

void GenerateOperatorSources(
    const RunParameters& run_parameters,
    RadialOperatorMap& radial_operators,
    std::vector<OperatorChannel>& operator_channels,
    std::map<std::string,Eigen::MatrixXd>& source_matrices,
    const typename basis::TwoBodySectorsJJJPN::SectorType& target_sector
  )
{

  // iterate over channels
  for (auto& operator_channel : operator_channels)
    {
            
      // read matrix for sector
      Eigen::MatrixXd operator_matrix;
      
      // identity
      if (operator_channel.id=="identity")
        operator_matrix = shell::IdentityOperatorMatrixJJJPN(target_sector,run_parameters.A);
      // kinematic
      //
      // TODO neaten by mapping id to
      // (radial_operator_labels,kinematic_operator_type,radial_operator_type)
      // and refactoring
      else if (operator_channel.id=="Ursqr")
        {
          RadialOperatorLabels radial_operator_labels = RadialOperatorLabels(shell::RadialOperatorType::kR,2);
          shell::KinematicOperatorType kinematic_operator_type = shell::KinematicOperatorType::kUTSqr;
          shell::RadialOperatorType radial_operator_type = shell::RadialOperatorType::kR;
          const RadialOperatorData& radial_operator_data = radial_operators[radial_operator_labels];
          operator_matrix = shell::KinematicMatrixJJJPN(
              radial_operator_data.ket_orbital_space,radial_operator_data.sectors,radial_operator_data.matrices,
              kinematic_operator_type,radial_operator_type,
              target_sector,run_parameters.A
            );
        }
      else if (operator_channel.id=="Vr1r2")
        {
          RadialOperatorLabels radial_operator_labels = RadialOperatorLabels(shell::RadialOperatorType::kR,1);
          shell::KinematicOperatorType kinematic_operator_type = shell::KinematicOperatorType::kVT1T2;
          shell::RadialOperatorType radial_operator_type = shell::RadialOperatorType::kR;
          const RadialOperatorData& radial_operator_data = radial_operators[radial_operator_labels];
          operator_matrix = shell::KinematicMatrixJJJPN(
              radial_operator_data.ket_orbital_space,radial_operator_data.sectors,radial_operator_data.matrices,
              kinematic_operator_type,radial_operator_type,
              target_sector,run_parameters.A
            );
        }
      else if (operator_channel.id=="Uksqr")
        {
          RadialOperatorLabels radial_operator_labels = RadialOperatorLabels(shell::RadialOperatorType::kK,2);
          shell::KinematicOperatorType kinematic_operator_type = shell::KinematicOperatorType::kUTSqr;
          shell::RadialOperatorType radial_operator_type = shell::RadialOperatorType::kK;
          const RadialOperatorData& radial_operator_data = radial_operators[radial_operator_labels];
          operator_matrix = shell::KinematicMatrixJJJPN(
              radial_operator_data.ket_orbital_space,radial_operator_data.sectors,radial_operator_data.matrices,
              kinematic_operator_type,radial_operator_type,
              target_sector,run_parameters.A
            );
        }
      else if (operator_channel.id=="Vk1k2")
        {
          RadialOperatorLabels radial_operator_labels = RadialOperatorLabels(shell::RadialOperatorType::kK,1);
          shell::KinematicOperatorType kinematic_operator_type = shell::KinematicOperatorType::kVT1T2;
          shell::RadialOperatorType radial_operator_type = shell::RadialOperatorType::kK;
          const RadialOperatorData& radial_operator_data = radial_operators[radial_operator_labels];
          operator_matrix = shell::KinematicMatrixJJJPN(
              radial_operator_data.ket_orbital_space,radial_operator_data.sectors,radial_operator_data.matrices,
              kinematic_operator_type,radial_operator_type,
              target_sector,run_parameters.A
            );
        }
      // angular momentum square
      //
      // TODO neaten by mapping id to (family,species) and refactoring
      else if (operator_channel.id=="Lp")
        operator_matrix = shell::AngularMomentumMatrixJJJPN(
            shell::AngularMomentumOperatorFamily::kOrbital,shell::AngularMomentumOperatorSpecies::kP,
            target_sector,run_parameters.A
          );
      else if (operator_channel.id=="Ln")
        operator_matrix = shell::AngularMomentumMatrixJJJPN(
            shell::AngularMomentumOperatorFamily::kOrbital,shell::AngularMomentumOperatorSpecies::kN,
            target_sector,run_parameters.A
          );
      else if (operator_channel.id=="L")
        operator_matrix = shell::AngularMomentumMatrixJJJPN(
            shell::AngularMomentumOperatorFamily::kOrbital,shell::AngularMomentumOperatorSpecies::kTotal,
            target_sector,run_parameters.A
          );
      else if (operator_channel.id=="Sp")
        operator_matrix = shell::AngularMomentumMatrixJJJPN(
            shell::AngularMomentumOperatorFamily::kSpin,shell::AngularMomentumOperatorSpecies::kP,
            target_sector,run_parameters.A
          );
      else if (operator_channel.id=="Sn")
        operator_matrix = shell::AngularMomentumMatrixJJJPN(
            shell::AngularMomentumOperatorFamily::kSpin,shell::AngularMomentumOperatorSpecies::kN,
            target_sector,run_parameters.A
          );
      else if (operator_channel.id=="S")
        operator_matrix = shell::AngularMomentumMatrixJJJPN(
            shell::AngularMomentumOperatorFamily::kSpin,shell::AngularMomentumOperatorSpecies::kTotal,
            target_sector,run_parameters.A
          );
      else if (operator_channel.id=="Jp")
        operator_matrix = shell::AngularMomentumMatrixJJJPN(
            shell::AngularMomentumOperatorFamily::kTotal,shell::AngularMomentumOperatorSpecies::kP,
            target_sector,run_parameters.A
          );
      else if (operator_channel.id=="Jn")
        operator_matrix = shell::AngularMomentumMatrixJJJPN(
            shell::AngularMomentumOperatorFamily::kTotal,shell::AngularMomentumOperatorSpecies::kN,
            target_sector,run_parameters.A
          );
      else if (operator_channel.id=="J")
        operator_matrix = shell::AngularMomentumMatrixJJJPN(
            shell::AngularMomentumOperatorFamily::kTotal,shell::AngularMomentumOperatorSpecies::kTotal,
            target_sector,run_parameters.A
          );

      // save result
      // std::cout << fmt::format("Saving {}...",operator_channel.id) << std::endl;
      source_matrices[operator_channel.id] = operator_matrix;
    }

}




void GenerateXformSources(
    std::vector<XformChannel>& xform_channels,
    std::map<std::string,Eigen::MatrixXd>& source_matrices,
    const typename basis::TwoBodySectorsJJJPN::SectorType& target_sector
  )
{
  if (target_sector.IsDiagonal())
    {
      //  SymmetrizeMatrix
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
      ChopMatrix(target_matrix);  // "neaten" output by eliminating near-zero values
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

void CloseXformChannels(
    std::vector<XformChannel>& xform_channels
  )
{
  for (auto& xform_channel : xform_channels)
    {
      // close stream
      xform_channel.stream_ptr->Close();

      // deallocate
      delete xform_channel.stream_ptr;
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
  InitializeXformChannels(run_parameters,target_indexing,xform_channels);
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
      GenerateOperatorSources(run_parameters,radial_operators,operator_channels,source_matrices,target_sector);
      // GenerateXformSources

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
  CloseXformChannels(xform_channels);
  CloseTargetChannels(target_channels);

  // end timing
  total_time.Stop();
  std::cout << "(Total time: " << total_time.ElapsedTime() << ")" << std::endl;
  std::cout << std::endl;

  // exit
  return EXIT_SUCCESS;
}
