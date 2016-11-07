/******************************************************************************

  h2mixer.cpp -- two-body matrix element input, generation, transformation, and output

  Syntax:
    h2mixer

  Mark A. Caprio
  University of Notre Dame

  + 10/23/16 (mac): Created, succeeding prior incarnation originated 3/14/12.
  + 10/25/16 (mac): Implement direct input mixing.
  + 10/29/16 (mac): Implement radial operator input.
  + 11/4/16 (mac):
    - Add warnings of incomplete two-body range coverage.
    - Templatize and fix ChopMatrix.
    - Revise control structure for generating operator source matrices.
  + 11/6/16 (mac):
    - Add OpenMP/Eigen parallel initialization.
    - Update xform stream calculation to match implementation in tbme_xform.

******************************************************************************/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <set>
#include <string>

#include "basis/operator.h"
#include "cppformat/format.h"
#include "mcutils/profiling.h"
#include "radial/radial_io.h"
#include "tbme/h2_io.h"
#include "tbme/tbme_separable.h"
#include "tbme/tbme_xform.h"
#include "tbme/tbme_mapping.h"

////////////////////////////////////////////////////////////////
// matrix utilities
////////////////////////////////////////////////////////////////

template<typename tMatrixType>
void ChopMatrix(tMatrixType& matrix, double tolerance=1e-10)
// Round near-zero entries in matrix to zero.
//
// Default tolerance value is same as Mathematica Chop function's.
//
// Template arguments:
//   tMatrixType: Eigen matrix type
//
// Arguments:
//   matrix (tMatrixType): matrix to chop
//   tolerance (double, optional): truncation tolerance
{
  for (int i=0; i<matrix.rows(); ++i)
    for (int j=0; j<matrix.cols(); ++j)
      // Caution: Important to use std::abs() rather than integer abs().
      if (std::abs(matrix(i,j))<tolerance)
        matrix(i,j) = 0.;
}

template<typename tMatrixType>
void CompleteLowerTriangle(tMatrixType& matrix)
// Populate lower triangle of matrix as mirror of upper triangle.
//
// Template arguments:
//   tMatrixType: Eigen matrix type
//
// Arguments:
//   matrix (tMatrixType): matrix to complete
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
//
// Caveat: Just "copying" indexing in here is not good enough, since
// then sectors can be left pointing to deleted temporaries for the
// orbital subspaces.
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

// collections of operator definitions
std::set<std::string> kIdentityOperatorIdSet(
    {"identity","loop-test"}
  );

typedef std::tuple<shell::KinematicOperatorType,shell::RadialOperatorType,int> KinematicOperatorDefinition;
std::map<std::string,KinematicOperatorDefinition> kKinematicOperatorDefinitions(
      {
        {"Ursqr",KinematicOperatorDefinition(shell::KinematicOperatorType::kUTSqr,shell::RadialOperatorType::kR,2)},
        {"Vr1r2",KinematicOperatorDefinition(shell::KinematicOperatorType::kVT1T2,shell::RadialOperatorType::kR,1)},
        {"Uksqr",KinematicOperatorDefinition(shell::KinematicOperatorType::kUTSqr,shell::RadialOperatorType::kK,2)},
        {"Vk1k2",KinematicOperatorDefinition(shell::KinematicOperatorType::kVT1T2,shell::RadialOperatorType::kK,1)}
      }
    );

typedef std::tuple<shell::AngularMomentumOperatorFamily,shell::AngularMomentumOperatorSpecies> AngularMomentumOperatorDefinition;
std::map<std::string,AngularMomentumOperatorDefinition> kAngularMomentumOperatorDefinitions(
      {
        // Note: Lp and Ln are disallowed, since Lp^2 and Ln^2 are
        // subject to physical misinterpretation (meaning not clear, as
        // the proton or neutron orbital angular momentum "relative"
        // to what?).  Jp and Jn are similarly disallowed.

        // {"Lp",AngularMomentumOperatorDefinition(shell::AngularMomentumOperatorFamily::kOrbital,shell::AngularMomentumOperatorSpecies::kP)},
        // {"Ln",AngularMomentumOperatorDefinition(shell::AngularMomentumOperatorFamily::kOrbital,shell::AngularMomentumOperatorSpecies::kN)},
        {"L",AngularMomentumOperatorDefinition(shell::AngularMomentumOperatorFamily::kOrbital,shell::AngularMomentumOperatorSpecies::kTotal)},
        {"Sp",AngularMomentumOperatorDefinition(shell::AngularMomentumOperatorFamily::kSpin,shell::AngularMomentumOperatorSpecies::kP)},
        {"Sn",AngularMomentumOperatorDefinition(shell::AngularMomentumOperatorFamily::kSpin,shell::AngularMomentumOperatorSpecies::kN)},
        {"S",AngularMomentumOperatorDefinition(shell::AngularMomentumOperatorFamily::kSpin,shell::AngularMomentumOperatorSpecies::kTotal)},
        // {"Jp",AngularMomentumOperatorDefinition(shell::AngularMomentumOperatorFamily::kTotal,shell::AngularMomentumOperatorSpecies::kP)},
        // {"Jn",AngularMomentumOperatorDefinition(shell::AngularMomentumOperatorFamily::kTotal,shell::AngularMomentumOperatorSpecies::kN)},
        {"J",AngularMomentumOperatorDefinition(shell::AngularMomentumOperatorFamily::kTotal,shell::AngularMomentumOperatorSpecies::kTotal)}
      }
  );

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
  // shell::TwoBodyMapping two_body_mapping;
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

          if (!(
                  kIdentityOperatorIdSet.count(id)
                || kKinematicOperatorDefinitions.count(id)
                  || kAngularMomentumOperatorDefinitions.count(id)
                ))
          ParsingError(line_count,line,"Unrecognized operator ID");
          
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
            "Reading radial matrix elements for operator {}^{} from {}...",
            char(radial_operator_type),radial_operator_power,
            radial_operator_data.filename
          )
        << std::endl;

      // open radial operator file
      shell::InRadialStream radial_operator_stream(radial_operator_data.filename);
      radial_operator_stream.SetToIndexing(
          radial_operator_data.bra_orbital_space,
          radial_operator_data.ket_orbital_space,
          radial_operator_data.sectors
        );
      assert(radial_operator_type==radial_operator_stream.radial_operator_type());
      assert(radial_operator_power==radial_operator_data.sectors.l0max());
      assert(radial_operator_data.sectors.Tz0()==0);

      // read matrices
      radial_operator_stream.Read(radial_operator_data.matrices);
      
      // close file
      radial_operator_stream.Close();

      // diagnostic
      std::cout
        << fmt::format(
            "  {} sectors; {} matrix elements",
            radial_operator_data.sectors.size(),
            basis::AllocatedEntries(radial_operator_data.matrices)
          )
        << std::endl;
    }

  if (radial_operators.size()>0)
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
      if (!input_channel.two_body_mapping.range_states_covered)
        {
          std::cout
            << fmt::format(
                "INFO: input stream {}: source does not provide full coverage of target indexing",
                input_channel.id
              )
            << std::endl
            << std::endl;
        }

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
      if (!xform_channel.pre_xform_two_body_mapping.range_states_covered)
        {
          std::cout
            << fmt::format(
                "INFO: xform stream {}: source does not provide full coverage of specified pre-xform truncation",
                xform_channel.id
              )
            << std::endl
            << std::endl;
        }

      // read xform coefficients
      //
      // analogous to InitializeRadialOperators

      std::cout
        << fmt::format("Reading radial overlaps from {}...",xform_channel.xform_filename)
        << std::endl;

      // set up alias (to keep code parallel to InitializeRadialOperators)
      RadialOperatorData& radial_operator_data = xform_channel.radial_operator_data;

      // set up radial operator data
      shell::RadialOperatorType radial_operator_type = shell::RadialOperatorType::kO;
      int radial_operator_power = 0;
      radial_operator_data = 
        RadialOperatorData(RadialOperatorLabels(radial_operator_type,radial_operator_power),xform_channel.xform_filename);

      // open radial operator file
      shell::InRadialStream radial_operator_stream(radial_operator_data.filename);
      radial_operator_stream.SetToIndexing(
          radial_operator_data.bra_orbital_space,
          radial_operator_data.ket_orbital_space,
          radial_operator_data.sectors
        );
      assert(radial_operator_type==radial_operator_stream.radial_operator_type());
      assert(radial_operator_power==radial_operator_data.sectors.l0max());
      assert(radial_operator_data.sectors.Tz0()==0);

      // read matrices
      radial_operator_stream.Read(radial_operator_data.matrices);
      
      // close file
      radial_operator_stream.Close();

      // diagnostic
      std::cout
        << fmt::format(
            "  {} sectors; {} matrix elements",
            radial_operator_data.sectors.size(),
            basis::AllocatedEntries(radial_operator_data.matrices)
          )
        << std::endl;

      // diagnostic: set up target mapping
      //
      // No actual mapping is carried out, since this is instead taken
      // care of by the radial overlaps in the similarity
      // transformation.  But it is informative to check the coverage
      // of the mapping in case the user is assuming a true identity
      // mapping.
      shell::TwoBodyMapping two_body_mapping = shell::TwoBodyMapping(
          xform_channel.pre_xform_two_body_indexing.orbital_space,
          xform_channel.pre_xform_two_body_indexing.space,
          target_indexing.orbital_space,
          target_indexing.space
        );
      // std::cout << input_channel.two_body_mapping.DebugStr() << std::endl;
      if (!two_body_mapping.range_states_covered)
        {
          std::cout
            << fmt::format(
                "INFO: xform stream {}: specified pre-xform truncation does not provide full coverage of target indexing (may be okay)",
                xform_channel.id
              )
            << std::endl
            << std::endl;
        }

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

      // short circuit if no corresponding input sector
      if (input_sector_index == basis::kNone)
        continue;

      // set up alias for input sector
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

      // fill in lower triangle of matrix
      //
      // The full square matrix must be populated before remapping.
      //
      // Caution: We naively assume a symmetric matrix.  This is
      // appropriate for scalar operators, but this may need to change
      // to a more general phase relation for two-body nonscalar
      // operators.
      assert(input_channel.stream_ptr->sectors().J0()==0);
      assert(input_channel.stream_ptr->sectors().g0()==0);
      assert(input_channel.stream_ptr->sectors().Tz0()==0);
      if (input_sector.IsDiagonal())
        CompleteLowerTriangle(input_matrix);

      // remap input matrix to target indexing
      Eigen::MatrixXd matrix
        = RemappedMatrixJJJPN(
            input_sector,
            target_sector,
            input_channel.two_body_mapping,
            input_matrix
          );

      // save result
      source_matrices[input_channel.id] = matrix;
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
      
      if (operator_channel.id=="identity")
        // identity
        {
          operator_matrix = shell::IdentityOperatorMatrixJJJPN(target_sector,run_parameters.A);
        }
      else if (operator_channel.id=="loop-test")
        // loop timing test
        {
          operator_matrix = shell::TimingTestMatrixJJJPN(target_sector,true,true);
        }
      else if (kKinematicOperatorDefinitions.count(operator_channel.id))
        // kinematic
        {
          shell::KinematicOperatorType kinematic_operator_type;
          shell::RadialOperatorType radial_operator_type;
          int radial_operator_power;
          std::tie(kinematic_operator_type,radial_operator_type,radial_operator_power)
            = kKinematicOperatorDefinitions.at(operator_channel.id);
          RadialOperatorLabels radial_operator_labels = RadialOperatorLabels(radial_operator_type,radial_operator_power);
          const RadialOperatorData& radial_operator_data = radial_operators[radial_operator_labels];

          // std::cout
          //   << fmt::format(
          //       "Sector: {} {} ",
          //       target_sector.bra_subspace().LabelStr(),target_sector.ket_subspace().LabelStr()
          //     )
          //   << std::endl;
          operator_matrix = shell::KinematicMatrixJJJPN(
              radial_operator_data.ket_orbital_space,radial_operator_data.sectors,radial_operator_data.matrices,
              kinematic_operator_type,radial_operator_type,
              target_sector,run_parameters.A
            );
        }
      else if (kAngularMomentumOperatorDefinitions.count(operator_channel.id))
        // angular momentum square
        {
          shell::AngularMomentumOperatorFamily operator_family;
          shell::AngularMomentumOperatorSpecies operator_species;
          std::tie(operator_family,operator_species)
            = kAngularMomentumOperatorDefinitions.at(operator_channel.id);
          operator_matrix = shell::AngularMomentumMatrixJJJPN(
              operator_family,operator_species,
              target_sector,run_parameters.A
            );
        }

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

  // iterate over channels
  for (auto& xform_channel : xform_channels)
    {
            
      // locate corresponding input sector
      int input_bra_subspace_index
        = xform_channel.stream_ptr->space().LookUpSubspaceIndex(target_sector.bra_subspace().labels());
      int input_ket_subspace_index
        = xform_channel.stream_ptr->space().LookUpSubspaceIndex(target_sector.ket_subspace().labels());
      int input_sector_index
        = xform_channel.stream_ptr->sectors().LookUpSectorIndex(input_bra_subspace_index,input_ket_subspace_index);

      // locate corresponding intermediate pre-xform (truncated) sector
      int pre_xform_bra_subspace_index
        = xform_channel.pre_xform_two_body_indexing.space.LookUpSubspaceIndex(target_sector.bra_subspace().labels());
      int pre_xform_ket_subspace_index
        = xform_channel.pre_xform_two_body_indexing.space.LookUpSubspaceIndex(target_sector.ket_subspace().labels());
      int pre_xform_sector_index
        = xform_channel.pre_xform_two_body_indexing.sectors.LookUpSectorIndex(pre_xform_bra_subspace_index,pre_xform_ket_subspace_index);

      // short circuit if no corresponding input or intermediate sectors
      if (input_sector_index == basis::kNone)
        continue;
      if (pre_xform_sector_index == basis::kNone)
        continue;

      // set up aliases for input and intermediate pre-xform sectors
      const typename basis::TwoBodySectorsJJJPN::SectorType& input_sector
        = xform_channel.stream_ptr->sectors().GetSector(input_sector_index);
      const typename basis::TwoBodySectorsJJJPN::SectorType& pre_xform_sector
        = xform_channel.pre_xform_two_body_indexing.sectors.GetSector(pre_xform_sector_index);

      // read matrix for sector
      Eigen::MatrixXd input_matrix;
      xform_channel.stream_ptr->SeekToSector(input_sector_index);
      xform_channel.stream_ptr->ReadSector(input_matrix);

      // fill in lower triangle of matrix
      //
      // The full square matrix must be populated before remapping.
      //
      // Caution: We naively assume a symmetric matrix.  This is
      // appropriate for scalar operators, but this may need to change
      // to a more general phase relation for two-body nonscalar
      // operators.
      assert(xform_channel.stream_ptr->sectors().J0()==0);
      assert(xform_channel.stream_ptr->sectors().g0()==0);
      assert(xform_channel.stream_ptr->sectors().Tz0()==0);
      if (input_sector.IsDiagonal())
        CompleteLowerTriangle(input_matrix);

      // remap input matrix to pre-xform indexing
      //
      // This is the "truncation cut" on the two-body transformation.
      Eigen::MatrixXd pre_xform_matrix
        = RemappedMatrixJJJPN(
            input_sector,
            pre_xform_sector,
            xform_channel.pre_xform_two_body_mapping,
            input_matrix
          );

      // do the transformation
      Eigen::MatrixXd matrix
        = shell::TwoBodyTransformedMatrix(
            xform_channel.radial_operator_data.bra_orbital_space,
            xform_channel.radial_operator_data.ket_orbital_space,
            xform_channel.radial_operator_data.sectors,
            xform_channel.radial_operator_data.matrices,
            pre_xform_sector,target_sector,
            pre_xform_matrix
          );

      // save result
      source_matrices[xform_channel.id] = matrix;

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

  // set up parallelization
  
  // for now, disable Eigen internal parallelization (but we will want it later for the matmul
  std::cout
    << fmt::format("Parallelization: max_threads {}, num_procs {}",
                   omp_get_max_threads(), omp_get_num_procs()
      )
    << std::endl
    << std::endl;
  Eigen::initParallel();
  Eigen::setNbThreads(0);

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
      GenerateXformSources(xform_channels,source_matrices,target_sector);

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
