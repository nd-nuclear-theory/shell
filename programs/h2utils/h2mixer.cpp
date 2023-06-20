/******************************************************************************

  h2mixer.cpp -- two-body matrix element input, generation, transformation, and output

  Syntax:
    h2mixer

  Input format:

    set-target-indexing <orbital_filename> <wp> <wn> <wpp> <wnn> <wpn>
    set-target-indexing-oscillator <rank> <cutoff>
      rank = ob|tb
    set-target-multipolarity <J0> <g0> <Tz0>
    set-mass <A>
    set-output-format <version>
      version = 0|15099|15200
    define-xform <id> <xform_filename>
    define-ob-source <mode> <id> ...
      define-ob-source input <id> <obme_filename> <J0> <g0> <Tz0>
      define-ob-source builtin <id> [orbital_filename]
        id = identity|l|l2|s|s2|j|j2|tz|t+|t-|c+|c
      define-ob-source linear-combination <id>
        add-ob-source <id> <coefficient>
      define-ob-source tensor-product <id> <ob_factor_a_id> <ob_factor_b_id> <J0> [scale_factor]
      define-ob-source xform <id> <ob_source_id> <xform_id>
    define-tb-source <mode> <id> ...
      define-tb-source input <id> <tbme_filename>
      define-tb-source builtin <id>
        id = identity|loop-test
      define-tb-source operatorU <id> <ob_source_id>
      define-tb-source operatorV <id> <ob_factor_a_id> <ob_factor_b_id> [scale_factor]
      define-tb-source xform <id> <tbme_filename> <wp> <wn> <wpp> <wnn> <wpn> <xform_id>
    define-target <filename>
    add-source <id> <coefficient>

  Language: C++11

  Mark A. Caprio, Patrick J. Fasano
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
  + 11/7/16 (mac):
    - Remove id on target channels.
    - Update input format.
  + 11/13/16 (mac): Move matrix utilities into mcutils/eigen.
  + 09/20/17 (pjf): Add T^2 operator:
    - Accept pn-overlaps for T^2 operator.
    - Store pn-overlaps by including Tz0 in radial_operator_data.
  + 10/12/17 (pjf): Update for changes to radial_io.
  + 11/28/17 (pjf): Include version in header.
  + 10/20/18 (pjf): Relax restrictions on nonscalar operators.
  + 02/08/19 (pjf): Almost complete rewrite:
    - Add support for generating and manipulating arbitrary one-body operators.
    - Add support for generating generic "upgraded operator" two-body operators.
    - Implement *breaking changes* to input format.
    - Use polymorphism for supporting new variety of channel types.
  + 02/12/19 (pjf): Relax constraints on non-isoscalar operators.
  + 02/21/19 (pjf): Add H2 Version15200 support.
  + 05/09/19 (pjf): Use std::size_t for basis indices and sizes.
  + 05/28/19 (mac/pjf): Allow copy construction of OneBodyOperatorData and
    XformData.
  + 07/01/19 (pjf):
    - Pass RunParameters and TwoBodyOperatorIndexing to
      ConstructOneBodyOperatorData.
    - Make scale_factor optional for tensor products.
    - Make orbital_filename optional for builtin ob source.
    - Remove all builtin two-body angular momentum operator code.
  + 08/22/19 (pjf): Remove ability to overwrite existing one-body channels.
  + 09/06/19 (pjf): Truncate target orbitals before constructing orbital space.
  + 09/13/20 (pjf): Fix xform ob source.
  + 11/25/20 (pjf): Allow Eigen multithreading.
******************************************************************************/

#include <omp.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <set>
#include <string>

#include "am/am.h"
#include "basis/nlj_orbital.h"
#include "basis/operator.h"
#include "fmt/format.h"
#include "mcutils/eigen.h"
#include "mcutils/parsing.h"
#include "mcutils/profiling.h"
#include "obme/obme.h"
#include "obme/obme_io.h"
#include "obme/radial.h"
#include "tbme/h2_io.h"
#include "tbme/tbme_separable.h"
#include "tbme/tbme_xform.h"
#include "tbme/tbme_mapping.h"

////////////////////////////////////////////////////////////////
// run control
////////////////////////////////////////////////////////////////

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
// xform container
////////////////////////////////////////////////////////////////

struct XformData
// Indexing and matrix elements for a change-of-basis transformation.
//
// The "bra" space is the source basis, and the "ket"
// space is the target basis.
{
  XformData() = default;

  XformData(const std::string& id_, const std::string& filename_)
    : id(id_), filename(filename_)
  {}

  // xform identification
  std::string id;

  // input filename
  std::string filename;

  // xform indexing and storage
  basis::OrbitalSpaceLJPN bra_orbital_space;
  basis::OrbitalSpaceLJPN ket_orbital_space;
  basis::OrbitalSectorsLJPN sectors;
  basis::OperatorBlocks<double> matrices;
};

typedef std::unordered_map<std::string,XformData> XformMap;
// Map to hold all loaded xforms.

////////////////////////////////////////////////////////////////
// one-body container
////////////////////////////////////////////////////////////////

struct OneBodyOperatorData
// Indexing and matrix elements for a one-body operator.
{

  OneBodyOperatorData() = default;

  OneBodyOperatorData(const std::string& id_)
    : id(id_)
  {}

  // operator identification
  std::string id;

  // operator indexing and storage
  basis::OrbitalSpaceLJPN orbital_space;
  basis::OrbitalSectorsLJPN sectors;
  basis::OperatorBlocks<double> matrices;
};

typedef std::unordered_map<std::string,OneBodyOperatorData> OneBodyOperatorMap;
// Map to hold all loaded one-body operators.

////////////////////////////////////////////////////////////////
// one-body channels
////////////////////////////////////////////////////////////////

struct OneBodySourceChannel
// Abstract base class for a one-body source channel.
//
// Classes which derive from this are expected to implement
// ConstructOneBodyOperatorData.
{
  OneBodySourceChannel() = default;

  explicit OneBodySourceChannel(const std::string& id_)
    : id(id_)
  {}

  virtual ~OneBodySourceChannel() = default;

  // operator identification
  std::string id;

  virtual void ConstructOneBodyOperatorData(
      const RunParameters& run_parameters,
      const XformMap& xforms,
      OneBodyOperatorMap& operators,
      const TwoBodyOperatorIndexing& target_indexing
    ) = 0;
  // Construct OneBodyOperatorData object in-place.
  //
  // This is a pure-virtual function. Derived classes must provide an
  // implementation which populates `operators` with the data described
  // by the channel.
  //
  // Arguments:
  //   xforms (input, XformMap): input transformations which may be applied
  //   operators (in/out, OneBodyOperatorMap): map of OneBodyOperatorData
  //     which may contain operators necessary for construction of current
  //     channel's operator
};

typedef std::vector<std::unique_ptr<OneBodySourceChannel>> OneBodyChannelList;
// Vector of pointers to base for one-body channels.

struct OneBodyInputChannel : public OneBodySourceChannel
// Parameters for a source channel providing input OBMEs.
{
  OneBodyInputChannel(
      const std::string& id_, const std::string& filename_,
      int J0_, int g0_, int Tz0_
    )
    : OneBodySourceChannel(id_), J0(J0_), g0(g0_), Tz0(Tz0_), filename(filename_)
  {}

  void ConstructOneBodyOperatorData(
      const RunParameters& run_parameters,
      const XformMap& xforms,
      OneBodyOperatorMap& operators,
      const TwoBodyOperatorIndexing& target_indexing
    ) override;

  // operator identification
  int J0, g0, Tz0;

  // input filename
  std::string filename;
};

struct OneBodyBuiltinChannel : public OneBodySourceChannel
// Parameters for a source channel providing an operator generated internally.
{
  OneBodyBuiltinChannel() = default;

  explicit OneBodyBuiltinChannel(
      const std::string& id_, const std::string& orbital_filename_
    )
    : OneBodySourceChannel(id_), orbital_filename(orbital_filename_)
  {}

  void ConstructOneBodyOperatorData(
      const RunParameters& run_parameters,
      const XformMap& xforms,
      OneBodyOperatorMap& operators,
      const TwoBodyOperatorIndexing& target_indexing
    ) override;

  // construction
  std::string orbital_filename;
};

std::unordered_map<std::string,std::pair<am::AngularMomentumOperatorType,int>>
kAngularMomentumOneBodyOperatorDefinitions =
  {
    {"l",  {am::AngularMomentumOperatorType::kOrbital, 1}},
    {"l2", {am::AngularMomentumOperatorType::kOrbital, 2}},
    {"s",  {am::AngularMomentumOperatorType::kSpin,    1}},
    {"s2", {am::AngularMomentumOperatorType::kSpin,    2}},
    {"j",  {am::AngularMomentumOperatorType::kTotal,   1}},
    {"j2", {am::AngularMomentumOperatorType::kTotal,   2}}
  };

std::unordered_map<std::string,int>
kIsospinOneBodyOperatorDefinitions =
  {
    {"tz",  0},
    {"t+", +1},
    {"t-", -1}
  };

std::unordered_map<std::string,shell::LadderOperatorType>
kLadderOneBodyOperatorDefinitions =
  {
    {"c+", shell::LadderOperatorType::kRaising},
    {"c", shell::LadderOperatorType::kLowering},
  };

struct OneBodyLinearCombinationChannel : public OneBodySourceChannel
// Parameters for a source channel providing a linear combination of one-body
// operators.
{
  OneBodyLinearCombinationChannel() = default;

  explicit OneBodyLinearCombinationChannel(const std::string& id_)
    : OneBodySourceChannel(id_)
  {}

  void ConstructOneBodyOperatorData(
      const RunParameters& run_parameters,
      const XformMap& xforms,
      OneBodyOperatorMap& operators,
      const TwoBodyOperatorIndexing& target_indexing
    ) override;

  // construction
  std::vector<std::pair<std::string,double>> coefficients;
};

struct OneBodyTensorProductChannel : public OneBodySourceChannel
// Parameters for a source channel providing an operator constructed as a
// tensor product on the one-body space.
{
  OneBodyTensorProductChannel() = default;

  OneBodyTensorProductChannel(
      const std::string& id_,
      const std::string& factor_a_id_, const std::string& factor_b_id_,
      int J0_, double scale_factor_
    )
    : OneBodySourceChannel(id_),
      ob_factor_a_id(factor_a_id_), ob_factor_b_id(factor_b_id_),
      J0(J0_), scale_factor(scale_factor_)
  {}

  void ConstructOneBodyOperatorData(
      const RunParameters& run_parameters,
      const XformMap& xforms,
      OneBodyOperatorMap& operators,
      const TwoBodyOperatorIndexing& target_indexing
    ) override;

  // construction
  std::string ob_factor_a_id, ob_factor_b_id;
  int J0;
  double scale_factor;
};

struct OneBodyXformChannel : public OneBodySourceChannel
// Parameters for a source channel providing an operator constructed via a
// similarity-transform.
{
  OneBodyXformChannel() = default;

  OneBodyXformChannel(
      const std::string& id_,
      const std::string& source_id_,
      const std::string& xform_id_
    )
    : OneBodySourceChannel(id_), ob_source_id(source_id_), xform_id(xform_id_)
  {}

  void ConstructOneBodyOperatorData(
      const RunParameters& run_parameters,
      const XformMap& xforms,
      OneBodyOperatorMap& operators,
      const TwoBodyOperatorIndexing& target_indexing
    ) override;

  // construction
  std::string ob_source_id, xform_id;
};

////////////////////////////////////////////////////////////////
// channel containers
////////////////////////////////////////////////////////////////

typedef std::map<std::string,basis::OperatorBlock<double>> OperatorBlockMap;
// Map to hold intermediate matrices.

struct TwoBodySourceChannel
// Abstract base class for a two-body source channel.
//
// Classes which derive from this are expected to implement
// InitializeChannel, PopulateSourceMatrix, and CloseChannel.
{
  TwoBodySourceChannel() = default;

  explicit TwoBodySourceChannel(const std::string& id_)
    : id(id_)
  {}

  virtual ~TwoBodySourceChannel() = default;

  // source identification
  std::string id;

  // channel initialization
  virtual void InitializeChannel(
      const RunParameters& run_parameters,
      const XformMap& xforms,
      const OneBodyOperatorMap& one_body_operators,
      const TwoBodyOperatorIndexing& target_indexing
    ) = 0;
  // Do one-time initialization for two-body channel.
  //
  // This is a pure-virtual function. Derived classes must provide an
  // implementation which initializes the two-body channel.
  //
  // Arguments:
  //   run_parameters (input): run parameters
  //   xforms (input): one-body transformations
  //   one_body_operators (input): one-body operator data
  //   target_indexing (input): indexing for target two-body operator

  // matrix generation
  virtual void PopulateSourceMatrix(
      const RunParameters& run_parameters,
      const XformMap& xforms,
      const OneBodyOperatorMap& one_body_operators,
      OperatorBlockMap& source_matrices,
      const typename basis::TwoBodySectorsJJJPN::SectorType& target_sector
    ) = 0;
  // Populate two-body source matrix in-place.
  //
  // This is a pure-virtual function. Derived classes must provide an
  // implementation which populates source_matrices[id] with the generated
  // two-body matrix.
  //
  // Arguments:
  //   run_parameters (input): run parameters
  //   xforms (input): one-body transformations
  //   one_body_operators (input): one-body operator data
  //   source_matrices (input/output): map of source matrices
  //   target_sector (input): sector associated with source matrix to be populated

  // cleanup
  virtual void CloseChannel() = 0;
  // Do one-time cleanup for two-body channel.
  //
  // This is a pure-virtual function. Derived classes must provide an
  // implementation.
};

typedef std::vector<std::unique_ptr<TwoBodySourceChannel>> TwoBodySourceChannelList;
// Vector of pointers to base for two-body source channels.

struct InputChannel : public TwoBodySourceChannel
// Parameters and data for a source channel providing directly-input
// TBMEs.
{
  InputChannel(const std::string& id_, const std::string& filename_)
    : TwoBodySourceChannel(id_), filename(filename_), stream_ptr(NULL)
  {}

  // input
  std::string filename;
  shell::InH2Stream* stream_ptr;

  // target mapping
  shell::TwoBodyMapping two_body_mapping;

  // channel initialization
  void InitializeChannel(
      const RunParameters& run_parameters,
      const XformMap& xforms,
      const OneBodyOperatorMap& one_body_operators,
      const TwoBodyOperatorIndexing& target_indexing
    ) override;

  // matrix generation
  void PopulateSourceMatrix(
      const RunParameters& run_parameters,
      const XformMap& xforms,
      const OneBodyOperatorMap& one_body_operators,
      OperatorBlockMap& source_matrices,
      const typename basis::TwoBodySectorsJJJPN::SectorType& target_sector
    ) override;

  // cleanup
  void CloseChannel() override;
};

// collections of operator definitions
std::set<std::string> kIdentityOperatorIdSet(
    {"identity","loop-test"}
  );

struct BuiltinChannel : public TwoBodySourceChannel
// Parameters and data for a source channel providing a generated
// operator.
{
  BuiltinChannel() = default;

  explicit BuiltinChannel(const std::string& id_)
    : TwoBodySourceChannel(id_)
  {}

  // channel initialization
  void InitializeChannel(
      const RunParameters& run_parameters,
      const XformMap& xforms,
      const OneBodyOperatorMap& one_body_operators,
      const TwoBodyOperatorIndexing& target_indexing
    ) override;

  // matrix generation
  void PopulateSourceMatrix(
      const RunParameters& run_parameters,
      const XformMap& xforms,
      const OneBodyOperatorMap& one_body_operators,
      OperatorBlockMap& source_matrices,
      const typename basis::TwoBodySectorsJJJPN::SectorType& target_sector
    ) override;

  // cleanup
  void CloseChannel() override;
};

struct OperatorUChannel : public TwoBodySourceChannel
// Parameters and data for a source channel providing an "upgraded"
// one-body operator.
{
  OperatorUChannel() = default;

  OperatorUChannel(const std::string& id_, const std::string& ob_source_id_)
    : TwoBodySourceChannel(id_), ob_source_id(ob_source_id_)
  {}

  // construction
  std::string ob_source_id;

  // channel initialization
  void InitializeChannel(
      const RunParameters& run_parameters,
      const XformMap& xforms,
      const OneBodyOperatorMap& one_body_operators,
      const TwoBodyOperatorIndexing& target_indexing
    ) override;

  // matrix generation
  void PopulateSourceMatrix(
      const RunParameters& run_parameters,
      const XformMap& xforms,
      const OneBodyOperatorMap& one_body_operators,
      OperatorBlockMap& source_matrices,
      const typename basis::TwoBodySectorsJJJPN::SectorType& target_sector
    ) override;

  // cleanup
  void CloseChannel() override;
};

struct OperatorVChannel : public TwoBodySourceChannel
// Parameters and data for a source channel providing an "upgraded"
// one-body operator.
{
  OperatorVChannel() = default;

  OperatorVChannel(
      const std::string& id_,
      const std::string& ob_factor_a_id_, const std::string& ob_factor_b_id_,
      const double& scale_factor_
    )
    : TwoBodySourceChannel(id_),
      ob_factor_a_id(ob_factor_a_id_), ob_factor_b_id(ob_factor_b_id_),
      scale_factor(scale_factor_)
  {}

  // construction
  std::string ob_factor_a_id, ob_factor_b_id;
  double scale_factor;

  // channel initialization
  void InitializeChannel(
      const RunParameters& run_parameters,
      const XformMap& xforms,
      const OneBodyOperatorMap& one_body_operators,
      const TwoBodyOperatorIndexing& target_indexing
    ) override;

  // matrix generation
  void PopulateSourceMatrix(
      const RunParameters& run_parameters,
      const XformMap& xforms,
      const OneBodyOperatorMap& one_body_operators,
      OperatorBlockMap& source_matrices,
      const typename basis::TwoBodySectorsJJJPN::SectorType& target_sector
    ) override;

  // cleanup
  void CloseChannel() override;
};

struct XformChannel : public TwoBodySourceChannel
// Parameters and data for a source channel providing two-body
// transformed TBMEs.
{
  XformChannel(
      const std::string& id_, const std::string& filename_,
      const basis::WeightMax& pre_xform_weight_max_,
      const std::string& xform_id_
    )
    : TwoBodySourceChannel(id_), filename(filename_),
      pre_xform_weight_max(pre_xform_weight_max_),
      xform_id(xform_id_),
      stream_ptr(NULL)
  {}

  // input
  std::string filename;
  shell::InH2Stream* stream_ptr;

  // pre-xform storage
  basis::WeightMax pre_xform_weight_max;
  TwoBodyOperatorIndexing pre_xform_two_body_indexing;
  shell::TwoBodyMapping pre_xform_two_body_mapping;

  // xform
  std::string xform_id;

  // target mapping
  // shell::TwoBodyMapping two_body_mapping;

  // channel initialization
  void InitializeChannel(
      const RunParameters& run_parameters,
      const XformMap& xforms,
      const OneBodyOperatorMap& one_body_operators,
      const TwoBodyOperatorIndexing& target_indexing
    ) override;

  // matrix generation
  void PopulateSourceMatrix(
      const RunParameters& run_parameters,
      const XformMap& xforms,
      const OneBodyOperatorMap& one_body_operators,
      OperatorBlockMap& source_matrices,
      const typename basis::TwoBodySectorsJJJPN::SectorType& target_sector
    ) override;

  // cleanup
  void CloseChannel() override;
};

struct TargetChannel
// Target channel parameters and stream.
{
  TargetChannel(const std::string& filename_)
    : filename(filename_), stream_ptr(NULL)
  {}

  // std::string id;  // id is actually superfluous

  // construction
  std::vector<std::pair<std::string,double>> coefficients;

  // output
  std::string filename;
  shell::OutH2Stream* stream_ptr;

  // channel initialization
  void InitializeChannel(
      const RunParameters& run_parameters,
      const TwoBodyOperatorIndexing& target_indexing
    );

  // matrix generation
  void GenerateTargetMatrix(
      const OperatorBlockMap& source_matrices,
      std::size_t target_sector_index,
      const typename basis::TwoBodySectorsJJJPN::SectorType& target_sector
    );

  // cleanup
  void CloseChannel();
};

typedef std::vector<std::unique_ptr<TargetChannel>> TwoBodyTargetChannelList;
// Vector of pointers to base for two-body source channels.

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

void ReadParameters(
    RunParameters& run_parameters,
    XformMap& xforms,
    OneBodyChannelList& one_body_channels,
    TwoBodySourceChannelList& two_body_channels,
    TwoBodyTargetChannelList& target_channels
)
// Parse control file.
{

  std::string line;
  int line_count = 0;
  while (mcutils::GetLine(std::cin, line, line_count))
    {
      // set up for line parsing
      std::istringstream line_stream(line);
      std::string keyword;
      line_stream >> keyword;

      // select action based on keyword
      if (keyword=="set-target-indexing")
        {
          double wp, wn, wpp, wnn, wpn;
          line_stream
            >> run_parameters.orbital_filename
            >> wp >> wn >> wpp >> wnn >> wpn;
          mcutils::ParsingCheck(line_stream,line_count,line);
          run_parameters.weight_max = basis::WeightMax(wp,wn,wpp,wnn,wpn);
        }
      else if (keyword=="set-target-indexing-oscillator")
        {
          std::string truncation_rank_code;
          line_stream >> truncation_rank_code
                      >> run_parameters.truncation_cutoff;
          mcutils::ParsingCheck(line_stream,line_count,line);

          // process truncation rank code
          if (truncation_rank_code == "ob")
            run_parameters.truncation_rank = basis::Rank::kOneBody;
          else if (truncation_rank_code == "tb")
            run_parameters.truncation_rank = basis::Rank::kTwoBody;
          else
            mcutils::ParsingError(line_count,line,"unrecognized truncation rank code");

          // store corresponding max weights
          run_parameters.weight_max
            = basis::WeightMax(run_parameters.truncation_rank,run_parameters.truncation_cutoff);
        }
      else if (keyword=="set-target-multipolarity")
        {
          line_stream >> run_parameters.J0 >> run_parameters.g0 >> run_parameters.Tz0;
          mcutils::ParsingCheck(line_stream,line_count,line);
        }
      else if (keyword=="set-mass")
        {
          line_stream >> run_parameters.A;
          mcutils::ParsingCheck(line_stream,line_count,line);
        }
      else if (keyword=="set-output-format")
        {
          line_stream >> run_parameters.output_h2_format;
          mcutils::ParsingCheck(line_stream,line_count,line);
        }
      else if (keyword=="define-xform")
        {
          std::string id, filename;
          line_stream >> id >> filename;

          // store xform information to map
          // TODO(pjf): understand why piecewise construction is required
          xforms.emplace(
              std::piecewise_construct,
              std::forward_as_tuple(id),
              std::forward_as_tuple(id, filename)
            );
        }
      else if (keyword=="define-ob-source")
        {
          std::string sub_keyword;
          line_stream >> sub_keyword;
          mcutils::ParsingCheck(line_stream,line_count,line);

          if (sub_keyword=="input")
            {
              std::string id, filename;
              int J0, g0, Tz0;
              line_stream >> id >> filename >> J0 >> g0 >> Tz0;
              mcutils::ParsingCheck(line_stream,line_count,line);
              one_body_channels.emplace_back(
                  new OneBodyInputChannel(id, filename, J0, g0, Tz0)
                );
            }
          else if (sub_keyword=="builtin")
            {
              std::string id, orbital_filename="";
              line_stream >> id;
              mcutils::ParsingCheck(line_stream,line_count,line);
              if (!line_stream.eof())
                {
                  line_stream >> orbital_filename;
                  mcutils::ParsingCheck(line_stream,line_count,line);
                }

              if (!(
                    (id == "identity")
                    || kAngularMomentumOneBodyOperatorDefinitions.count(id)
                    || kIsospinOneBodyOperatorDefinitions.count(id)
                    || kLadderOneBodyOperatorDefinitions.count(id)
                  ))
                mcutils::ParsingError(line_count,line,"Unrecognized builtin operator ID");

              one_body_channels.emplace_back(
                  new OneBodyBuiltinChannel(id, orbital_filename)
                );
            }
          else if (sub_keyword=="linear-combination")
            {
              std::string id;
              line_stream >> id;
              mcutils::ParsingCheck(line_stream,line_count,line);

              one_body_channels.emplace_back(
                  new OneBodyLinearCombinationChannel(id)
                );
            }
          else if (sub_keyword=="tensor-product")
            {
              std::string id, ob_factor_a_id, ob_factor_b_id;
              int J0;
              double scale_factor = 1.0;
              line_stream >> id >> ob_factor_a_id >> ob_factor_b_id >> J0;
              mcutils::ParsingCheck(line_stream,line_count,line);
              if (!line_stream.eof())
                {
                  line_stream >> scale_factor;
                  mcutils::ParsingCheck(line_stream,line_count,line);
                }

              one_body_channels.emplace_back(
                  new OneBodyTensorProductChannel(
                      id, ob_factor_a_id, ob_factor_b_id, J0, scale_factor
                    )
                );
            }
          else if (sub_keyword=="xform")
            {
              std::string id, source_id, xform_id;
              line_stream >> id >> source_id >> xform_id;
              mcutils::ParsingCheck(line_stream,line_count,line);

              one_body_channels.emplace_back(
                  new OneBodyXformChannel(id, source_id, xform_id)
                );
            }
          else
            {
              mcutils::ParsingError(line_count,line,"Unrecognized sub-keyword");
            }
        }
      else if (keyword=="add-ob-source")
        {
          assert(one_body_channels.size()>0);

          std::string id;
          double coefficient;
          line_stream >> id >> coefficient;
          mcutils::ParsingCheck(line_stream,line_count,line);

          // tack coefficient onto coefficient list of last member of
          // target_channels
          if (auto* one_body_channel = dynamic_cast<OneBodyLinearCombinationChannel*>(one_body_channels.back().get()))
            {
              one_body_channel->coefficients.emplace_back(id,coefficient);
            }
          else
            {
              mcutils::ParsingError(line_count,line,"Cannot add source to last one-body channel.");
            }
        }
      else if (keyword=="define-tb-source")
        {
          std::string sub_keyword;
          line_stream >> sub_keyword;
          mcutils::ParsingCheck(line_stream,line_count,line);

          if (sub_keyword=="input")
            {
              std::string id, filename;
              line_stream >> id >> filename;
              mcutils::ParsingCheck(line_stream,line_count,line);
              two_body_channels.emplace_back(
                  new InputChannel(id,filename)
                );
            }
          else if (sub_keyword=="builtin")
            {
              std::string id;
              line_stream >> id;
              mcutils::ParsingCheck(line_stream,line_count,line);

              if (!(kIdentityOperatorIdSet.count(id)))
                mcutils::ParsingError(line_count,line,"Unrecognized operator ID");

              two_body_channels.emplace_back(
                  new BuiltinChannel(id)
                );
            }
          else if (sub_keyword=="operatorU")
            {
              std::string id, ob_source_id;
              line_stream >> id >> ob_source_id;
              mcutils::ParsingCheck(line_stream, line_count, line);

              two_body_channels.emplace_back(
                  new OperatorUChannel(id, ob_source_id)
                );
            }
          else if (sub_keyword=="operatorV")
            {
              std::string id, ob_factor_a_id, ob_factor_b_id;
              double scale_factor = 1.0;
              line_stream >> id >> ob_factor_a_id >> ob_factor_b_id;
              mcutils::ParsingCheck(line_stream, line_count, line);
              if (!line_stream.eof())
                {
                  line_stream >> scale_factor;
                  mcutils::ParsingCheck(line_stream,line_count,line);
                }

              two_body_channels.emplace_back(
                  new OperatorVChannel(
                        id, ob_factor_a_id, ob_factor_b_id, scale_factor
                      )
                );
            }
          else if (sub_keyword=="xform")
            {
              std::string id, filename, xform_id;
              double wp, wn, wpp, wnn, wpn;
              line_stream
                >> id >> filename
                >> wp >> wn >> wpp >> wnn >> wpn
                >> xform_id;
              mcutils::ParsingCheck(line_stream,line_count,line);

              two_body_channels.emplace_back(
                  new XformChannel(
                      id, filename, basis::WeightMax(wp,wn,wpp,wnn,wpn), xform_id
                    )
                );
            }
          else
            {
              mcutils::ParsingError(line_count,line,"Unrecognized sub-keyword");
            }
        }
      else if (keyword=="define-target")
        {
          std::string filename;
          line_stream >> filename;
          mcutils::ParsingCheck(line_stream,line_count,line);
          target_channels.emplace_back(
              new TargetChannel(filename)
            );
        }
      else if (keyword=="add-source")
        {
          assert(target_channels.size()>0);

          std::string id;
          double coefficient;
          line_stream >> id >> coefficient;
          mcutils::ParsingCheck(line_stream,line_count,line);

          // tack coefficient onto coefficient list of last member of
          // target_channels
          (target_channels.back())->coefficients.emplace_back(id,coefficient);
        }
      else
        {
          mcutils::ParsingError(line_count,line,"Unrecognized keyword");
        }

    }

}

////////////////////////////////////////////////////////////////
// xform initialization code
/////////////////////////////////////////////////////////////////

void InitializeXforms(
    XformMap& xforms
  )
// Loop over xforms in XformMap, populating them with data.
{
  for (auto& id_data : xforms)
  {
    XformData& xform_data = id_data.second;
    std::cout
      << fmt::format(
          "Reading xform matrix elements for xform {} from {}...",
          xform_data.id, xform_data.filename
        )
      << std::endl;

    // open xform file
    shell::InOBMEStream xform_stream(xform_data.filename);
    xform_stream.SetToIndexing(
        xform_data.bra_orbital_space,
        xform_data.ket_orbital_space,
        xform_data.sectors
      );
    assert(xform_data.sectors.J0() == 0);
    assert(xform_data.sectors.g0() == 0);
    assert(xform_data.sectors.Tz0() == 0);

    // read matrices
    xform_stream.Read(xform_data.matrices);

    // close file
    xform_stream.Close();

    // diagnostic
    std::cout
      << fmt::format(
          "  {} sectors; {} matrix elements",
          xform_data.sectors.size(),
          basis::AllocatedEntries(xform_data.matrices)
        )
      << std::endl;
  }

  if (xforms.size()>0)
    std::cout << std::endl;
}

////////////////////////////////////////////////////////////////
// one-body operator initialization code
/////////////////////////////////////////////////////////////////

void OneBodyInputChannel::ConstructOneBodyOperatorData(
    const RunParameters& run_parameters,
    const XformMap& xforms,
    OneBodyOperatorMap& operators,
    const TwoBodyOperatorIndexing& target_indexing
  )
// Populate map with data from input file.
{
  std::cout
    << fmt::format(
        "Reading one-body matrix elements for operator {} from {}...",
        id, filename
      )
    << std::endl;

  // construct new OneBodyOperatorData in-place
  OneBodyOperatorMap::iterator it;
  bool new_operator;
  std::tie(it, new_operator) = operators.emplace(id, id);
  if (!new_operator)
  {
    std::cerr << fmt::format("ERROR: Operator {} exists!", id) << std::endl;
    std::exit(EXIT_FAILURE);
  }
  OneBodyOperatorData& operator_data = it->second;

  // open radial operator file
  shell::InOBMEStream operator_stream(filename);
  operator_stream.SetToIndexing(
      operator_data.orbital_space,
      operator_data.orbital_space,
      operator_data.sectors
    );
  assert(J0 == operator_data.sectors.J0());
  assert(g0 == operator_data.sectors.g0());
  assert(Tz0 == operator_data.sectors.Tz0());

  // read matrices
  operator_stream.Read(operator_data.matrices);

  // close file
  operator_stream.Close();

  // diagnostic
  std::cout
    << fmt::format(
        "  {} sectors; {} matrix elements",
        operator_data.sectors.size(),
        basis::AllocatedEntries(operator_data.matrices)
      )
    << std::endl;
}

void OneBodyBuiltinChannel::ConstructOneBodyOperatorData(
    const RunParameters& run_parameters,
    const XformMap& xforms,
    OneBodyOperatorMap& operators,
    const TwoBodyOperatorIndexing& target_indexing
  )
// Populate map with data from built-in generation routines.
{
  std::cout
    << fmt::format(
        "Constructing one-body matrix elements for operator {} with indexing from {}...",
        id, orbital_filename
      )
    << std::endl;

  // construct new OneBodyOperatorData in-place
  OneBodyOperatorMap::iterator it;
  bool new_operator;
  std::tie(it, new_operator) = operators.emplace(id, id);
  if (!new_operator)
  {
    std::cerr << fmt::format("ERROR: Operator {} exists!", id) << std::endl;
    std::exit(EXIT_FAILURE);
  }
  OneBodyOperatorData& operator_data = it->second;

  if (orbital_filename == "")
    {
      operator_data.orbital_space = basis::OrbitalSpaceLJPN(target_indexing.orbital_space);
    }
  else
    {
      std::ifstream orbital_file(orbital_filename);
      basis::OrbitalPNList orbital_info = basis::ParseOrbitalPNStream(orbital_file,true);
      operator_data.orbital_space = basis::OrbitalSpaceLJPN(orbital_info);
    }

  if (id == "identity")
    // identity operator
    {
      operator_data.sectors = basis::OrbitalSectorsLJPN(
          operator_data.orbital_space,
          0, 0, 0 // J0, g0, Tz0
        );
      basis::SetOperatorToIdentity(operator_data.sectors, operator_data.matrices);
    }
  else if (kAngularMomentumOneBodyOperatorDefinitions.count(id))
    // angular momentum operators
    {
      am::AngularMomentumOperatorType am_operator_type;
      int power, J0;

      // set up indexing
      std::tie(am_operator_type, power) = kAngularMomentumOneBodyOperatorDefinitions.at(id);
      J0 = power%2;

      operator_data.sectors = basis::OrbitalSectorsLJPN(
          operator_data.orbital_space,
          J0, 0, 0 // J0, g0, Tz0
        );

      // populate operator
      if (power == 1)
        {
          shell::AngularMomentumOneBodyOperator(
              am_operator_type,
              operator_data.orbital_space,
              operator_data.sectors,
              operator_data.matrices
            );
        }
      else if (power == 2)
        {
          shell::AngularMomentumSquaredOneBodyOperator(
              am_operator_type,
              operator_data.orbital_space,
              operator_data.sectors,
              operator_data.matrices
            );
        }
    }
  else if (kIsospinOneBodyOperatorDefinitions.count(id))
    {
      int J0 = 0, g0 = 0;
      int Tz0 = kIsospinOneBodyOperatorDefinitions.at(id);
      operator_data.sectors = basis::OrbitalSectorsLJPN(
          operator_data.orbital_space,
          J0, g0, Tz0
        );
      shell::IsospinOneBodyOperator(
          operator_data.orbital_space,
          operator_data.sectors,
          operator_data.matrices
        );
    }
  else if (kLadderOneBodyOperatorDefinitions.count(id))
    {
      int J0 = 1, g0 = 1, Tz0 = 0;
      operator_data.sectors = basis::OrbitalSectorsLJPN(
          operator_data.orbital_space,
          J0, g0, Tz0
        );
      shell::LadderOneBodyOperator(
          shell::RadialBasisType::kOscillator,
          kLadderOneBodyOperatorDefinitions.at(id),
          operator_data.orbital_space,
          operator_data.sectors,
          operator_data.matrices
        );
    }
  else
    {
      std::cout << fmt::format("ERROR: builtin {} not yet implemented", id)
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
}

void OneBodyLinearCombinationChannel::ConstructOneBodyOperatorData(
    const RunParameters& run_parameters,
    const XformMap& xforms,
    OneBodyOperatorMap& operators,
    const TwoBodyOperatorIndexing& target_indexing
  )
{
  std::cout
    << fmt::format("Generating one-body linear combination {}", id)
    << std::endl;
  // construct new OneBodyOperatorData in-place
  OneBodyOperatorMap::iterator it;
  bool new_operator;
  std::tie(it, new_operator) = operators.emplace(id, id);
  if (!new_operator)
  {
    std::cerr << fmt::format("ERROR: Operator {} exists!", id) << std::endl;
    std::exit(EXIT_FAILURE);
  }
  OneBodyOperatorData& operator_data = it->second;

  // loop over id-coefficient pairs
  bool first_term = true;
  for (const auto& id_coefficient : coefficients)
    {
      const std::string& source_id = id_coefficient.first;
      double coefficient = id_coefficient.second;

      // look up source data
      if (!operators.count(source_id))
      {
        std::cerr << fmt::format("ERROR {}: Operator {} not found.", __LINE__, source_id)
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      const OneBodyOperatorData& source_data = operators.at(source_id);

      if (first_term)
        // copy indexing from first source
        {
          operator_data.orbital_space = source_data.orbital_space;
          int J0 = source_data.sectors.J0();
          int g0 = source_data.sectors.g0();
          int Tz0 = source_data.sectors.Tz0();
          operator_data.sectors = basis::OrbitalSectorsLJPN(
              operator_data.orbital_space, J0, g0, Tz0
            );
          basis::SetOperatorToZero(operator_data.sectors, operator_data.matrices);
          first_term = false;
        }
      else
        // assert same indexing for later terms
        {
          assert(source_data.orbital_space.OrbitalInfo() == operator_data.orbital_space.OrbitalInfo());
          assert(operator_data.sectors.J0() == source_data.sectors.J0());
          assert(operator_data.sectors.g0() == source_data.sectors.g0());
          assert(operator_data.sectors.Tz0() == source_data.sectors.Tz0());
        }

      basis::OperatorLinearCombination(
          operator_data.sectors, operator_data.matrices,
          1., operator_data.matrices,
          coefficient, source_data.matrices
        );
    }
}

void OneBodyTensorProductChannel::ConstructOneBodyOperatorData(
    const RunParameters& run_parameters,
    const XformMap& xforms,
    OneBodyOperatorMap& operators,
    const TwoBodyOperatorIndexing& target_indexing
  )
{
  std::cout
    << fmt::format(
        "Generating one-body tensor product {} = [({})({})]_{}",
        id, ob_factor_a_id, ob_factor_b_id, J0
      )
    << std::endl;
  // construct new OneBodyOperatorData in-place
  OneBodyOperatorMap::iterator it;
  bool new_operator;
  std::tie(it, new_operator) = operators.emplace(id, id);
  if (!new_operator)
  {
    std::cerr << fmt::format("ERROR: Operator {} exists!", id) << std::endl;
    std::exit(EXIT_FAILURE);
  }
  OneBodyOperatorData& operator_data = it->second;

  // look up factor data
  if (!operators.count(ob_factor_a_id))
  {
    std::cerr << fmt::format("ERROR {}: Operator {} not found.", __LINE__, ob_factor_a_id)
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const OneBodyOperatorData& data_a = operators.at(ob_factor_a_id);
  if (!operators.count(ob_factor_b_id))
  {
    std::cerr << fmt::format("ERROR {}: Operator {} not found.", __LINE__, ob_factor_b_id)
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const OneBodyOperatorData& data_b = operators.at(ob_factor_b_id);

  // construct new indexing
  assert(data_a.orbital_space.OrbitalInfo() == data_b.orbital_space.OrbitalInfo());
  operator_data.orbital_space = data_a.orbital_space;
  assert(am::AllowedTriangle(data_a.sectors.J0(), data_b.sectors.J0(), J0));
  int g0 = (data_a.sectors.g0() + data_b.sectors.g0())%2;
  int Tz0 = data_a.sectors.Tz0() + data_b.sectors.Tz0();
  operator_data.sectors = basis::OrbitalSectorsLJPN(
      operator_data.orbital_space, J0, g0, Tz0
    );

  shell::OneBodyOperatorTensorProduct(
      operator_data.orbital_space,
      data_a.sectors,
      data_a.matrices,
      data_b.sectors,
      data_b.matrices,
      operator_data.sectors,
      operator_data.matrices
    );

  // scale product after construction
  if (scale_factor != 1.0)
    basis::ScalarMultiplyOperator(operator_data.sectors, operator_data.matrices, scale_factor);
}

void OneBodyXformChannel::ConstructOneBodyOperatorData(
    const RunParameters& run_parameters,
    const XformMap& xforms,
    OneBodyOperatorMap& operators,
    const TwoBodyOperatorIndexing& target_indexing
  )
{
  std::cout
    << fmt::format(
        "Applying xform {} to one-body matrix elements for operator {}...",
        xform_id, id
      )
    << std::endl;

  // look up original operator data
  if (!operators.count(ob_source_id))
  {
    std::cerr << fmt::format("ERROR {}: Operator {} not found.", __LINE__, ob_source_id)
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const OneBodyOperatorData& source_data = operators.at(ob_source_id);

  // look up xform data
  if (!xforms.count(xform_id))
  {
    std::cerr << fmt::format("ERROR {}: xform {} not found.", __LINE__, xform_id)
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const XformData& xform_data = xforms.at(xform_id);

  // construct new OneBodyOperatorData in-place
  OneBodyOperatorMap::iterator it;
  bool new_operator;
  std::tie(it, new_operator) = operators.emplace(id, id);
  if (!new_operator)
  {
    std::cerr << fmt::format("ERROR: Operator {} exists!", id) << std::endl;
    std::exit(EXIT_FAILURE);
  }
  OneBodyOperatorData& operator_data = it->second;

  // construct new indexing
  assert(source_data.orbital_space.OrbitalInfo() == xform_data.bra_orbital_space.OrbitalInfo());
  operator_data.orbital_space = xform_data.ket_orbital_space;
  operator_data.sectors = basis::OrbitalSectorsLJPN(
      xform_data.ket_orbital_space,
      source_data.sectors.J0(), source_data.sectors.g0(), source_data.sectors.Tz0()
    );

  // perform transformation
  shell::SimilarityTransformOperator(
      xform_data.bra_orbital_space,
      xform_data.ket_orbital_space,
      xform_data.sectors,
      xform_data.matrices,
      source_data.sectors,
      source_data.matrices,
      operator_data.sectors,
      operator_data.matrices
    );

}

////////////////////////////////////////////////////////////////
// set up two-body indexing
/////////////////////////////////////////////////////////////////

void InitializeTargetIndexing(
    const RunParameters& run_parameters,
    TwoBodyOperatorIndexing& target_indexing
  )
{

  // space ordering
  basis::TwoBodySpaceJJJPNOrdering space_ordering =
    shell::kH2SpaceOrdering.at(run_parameters.output_h2_format);

  // set up basis
  if (run_parameters.truncation_cutoff >= 0)
    // oscillator-like indexing
    {
      target_indexing.orbital_space = basis::OrbitalSpacePN(run_parameters.truncation_cutoff);
    }
  else
    // generic indexing
    {
      std::ifstream orbital_file(run_parameters.orbital_filename);
      basis::OrbitalPNList orbital_info = basis::ParseOrbitalPNStream(orbital_file,true);
      orbital_info = basis::TruncateOrbitalList(run_parameters.weight_max, orbital_info);
      target_indexing.orbital_space = basis::OrbitalSpacePN(orbital_info);
    }

  // set up two-body space
  target_indexing.space = basis::TwoBodySpaceJJJPN(
      target_indexing.orbital_space, run_parameters.weight_max, space_ordering
    );
  // set up sectors
  target_indexing.sectors = basis::TwoBodySectorsJJJPN(
      target_indexing.space,
      run_parameters.J0,run_parameters.g0,run_parameters.Tz0
    );
}

////////////////////////////////////////////////////////////////
// two-body channel initialization
/////////////////////////////////////////////////////////////////

void InputChannel::InitializeChannel(
    const RunParameters& run_parameters,
    const XformMap& xforms,
    const OneBodyOperatorMap& one_body_operators,
    const TwoBodyOperatorIndexing& target_indexing
  )
{
  // create stream
  stream_ptr = new shell::InH2Stream(filename);

  // write diagnostic
  std::cout << fmt::format("Input stream (direct) -- {}", id) << std::endl;
  std::cout << stream_ptr->DiagnosticStr();
  std::cout << std::endl;

  // set up target mapping
  two_body_mapping = shell::TwoBodyMapping(
      stream_ptr->orbital_space(),
      stream_ptr->space(),
      target_indexing.orbital_space,
      target_indexing.space
    );
  // std::cout << input_channel.two_body_mapping.DebugStr() << std::endl;
  if (!two_body_mapping.range_states_covered)
    {
      std::cout
        << fmt::format(
            "INFO: input file does not provide full coverage of target indexing",
            id
          )
        << std::endl
        << std::endl;
    }
}

void BuiltinChannel::InitializeChannel(
    const RunParameters &run_parameters,
    const XformMap &xforms,
    const OneBodyOperatorMap &one_body_operators,
    const TwoBodyOperatorIndexing &target_indexing
  )
{}

void OperatorUChannel::InitializeChannel(
    const RunParameters &run_parameters,
    const XformMap &xforms,
    const OneBodyOperatorMap &one_body_operators,
    const TwoBodyOperatorIndexing &target_indexing
  )
{
  // check for existence of one-body operator
  if (!one_body_operators.count(ob_source_id))
  {
    std::cerr << fmt::format("ERROR {}: one-body operator {} not found.", __LINE__, ob_source_id)
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // check quantum numbers of one-body operator
  const OneBodyOperatorData& ob_data = one_body_operators.at(ob_source_id);
  assert(ob_data.sectors.J0() == run_parameters.J0);
  assert(ob_data.sectors.g0() == run_parameters.g0);
  assert(ob_data.sectors.Tz0() == run_parameters.Tz0);
}

void OperatorVChannel::InitializeChannel(
    const RunParameters &run_parameters,
    const XformMap &xforms,
    const OneBodyOperatorMap &one_body_operators,
    const TwoBodyOperatorIndexing &target_indexing
  )
{
  // check for existence of one-body operator
  if (!one_body_operators.count(ob_factor_a_id))
  {
    std::cerr << fmt::format("ERROR {}: one-body operator {} not found.", __LINE__, ob_factor_a_id)
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (!one_body_operators.count(ob_factor_b_id))
  {
    std::cerr << fmt::format("ERROR {}: one-body operator {} not found.", __LINE__, ob_factor_b_id)
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // check if valid tensor product can be constructed
  const OneBodyOperatorData& ob_data_a = one_body_operators.at(ob_factor_a_id);
  const OneBodyOperatorData& ob_data_b = one_body_operators.at(ob_factor_b_id);
  assert(ob_data_a.orbital_space.OrbitalInfo() == ob_data_b.orbital_space.OrbitalInfo());
  assert(am::AllowedTriangle(ob_data_a.sectors.J0(), ob_data_b.sectors.J0(), run_parameters.J0));
  assert((ob_data_a.sectors.g0()+ob_data_b.sectors.g0()+run_parameters.g0)%2==0);
  assert(ob_data_a.sectors.Tz0()+ob_data_b.sectors.Tz0() == run_parameters.Tz0);
}

void XformChannel::InitializeChannel(
    const RunParameters& run_parameters,
    const XformMap& xforms,
    const OneBodyOperatorMap& one_body_operators,
    const TwoBodyOperatorIndexing& target_indexing
  )
{
  // create stream
  stream_ptr = new shell::InH2Stream(filename);

  // write diagnostic
  std::cout << fmt::format("Input stream (xform) -- {}", id) << std::endl;
  std::cout << stream_ptr->DiagnosticStr();
  std::cout << std::endl;

  // set up pre-xform indexing and mapping
  //
  // preserves input orbitals but builds two-body space with
  // user-specified truncation
  pre_xform_two_body_indexing.orbital_space = stream_ptr->orbital_space();
  basis::TwoBodySpaceJJJPNOrdering pre_xform_space_ordering =
    shell::kH2SpaceOrdering.at(stream_ptr->h2_format());
  pre_xform_two_body_indexing.space
    = basis::TwoBodySpaceJJJPN(
        pre_xform_two_body_indexing.orbital_space,
        pre_xform_weight_max,
        pre_xform_space_ordering
      );
  pre_xform_two_body_indexing.sectors = basis::TwoBodySectorsJJJPN(
      pre_xform_two_body_indexing.space,
      run_parameters.J0,run_parameters.g0,run_parameters.Tz0
    );
  pre_xform_two_body_mapping = shell::TwoBodyMapping(
      stream_ptr->orbital_space(),
      stream_ptr->space(),
      pre_xform_two_body_indexing.orbital_space,
      pre_xform_two_body_indexing.space
    );
  if (!pre_xform_two_body_mapping.range_states_covered)
    {
      std::cout
        << fmt::format(
            "INFO: input file does not provide full coverage of specified pre-xform truncation",
            id
          )
        << std::endl
        << std::endl;
    }

  // check for xform coefficients
  std::cout
    << fmt::format("Using xform overlaps from {}...", xform_id)
    << std::endl;

  if (!xforms.count(xform_id))
  {
    std::cerr << fmt::format("ERROR {}: xform {} not found.", __LINE__, xform_id)
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // diagnostic: set up target mapping
  //
  // No actual mapping is carried out, since this is instead taken
  // care of by the radial overlaps in the similarity
  // transformation.  But it is informative to check the coverage
  // of the mapping in case the user is assuming a true identity
  // mapping.
  shell::TwoBodyMapping two_body_mapping = shell::TwoBodyMapping(
      pre_xform_two_body_indexing.orbital_space,
      pre_xform_two_body_indexing.space,
      target_indexing.orbital_space,
      target_indexing.space
    );
  // std::cout << input_channel.two_body_mapping.DebugStr() << std::endl;
  if (!two_body_mapping.range_states_covered)
    {
      std::cout
        << fmt::format(
            "INFO: specified pre-xform truncation does not"
            " provide full coverage of target indexing",
            id
          )
        << std::endl
        << std::endl;
    }

  std::cout << std::endl;
}

void TargetChannel::InitializeChannel(
    const RunParameters& run_parameters,
    const TwoBodyOperatorIndexing& target_indexing
  )
{
  // create stream
  stream_ptr = new shell::OutH2Stream(
      filename,
      target_indexing.orbital_space,target_indexing.space,target_indexing.sectors,
      run_parameters.output_h2_format
    );

  // write diagnostic
  std::cout << fmt::format("Output stream") << std::endl;
  std::cout << stream_ptr->DiagnosticStr();
  std::cout << std::endl;
}

////////////////////////////////////////////////////////////////
// two-body matrix generation
/////////////////////////////////////////////////////////////////

void InputChannel::PopulateSourceMatrix(
    const RunParameters& run_parameters,
    const XformMap& xforms,
    const OneBodyOperatorMap& one_body_operators,
    OperatorBlockMap& source_matrices,
    const typename basis::TwoBodySectorsJJJPN::SectorType& target_sector
  )
{
  // locate corresponding input sector
  std::size_t input_bra_subspace_index
    = stream_ptr->space().LookUpSubspaceIndex(target_sector.bra_subspace().labels());
  std::size_t input_ket_subspace_index
    = stream_ptr->space().LookUpSubspaceIndex(target_sector.ket_subspace().labels());
  // Note: We cannot simply look up by target_sector's Key, since that uses target subspace indices.
  std::size_t input_sector_index
    = stream_ptr->sectors().LookUpSectorIndex(input_bra_subspace_index,input_ket_subspace_index);
  // std::cout << fmt::format("Input channel {} index {}...",input_channel.id,input_sector_index) << std::endl;

  // short circuit if no corresponding input sector
  if (input_sector_index == basis::kNone)
    return;

  // set up alias for input sector
  const typename basis::TwoBodySectorsJJJPN::SectorType& input_sector
    = stream_ptr->sectors().GetSector(input_sector_index);
  // std::cout << "input" <<std::endl
  //           << input_sector.ket_subspace().LabelStr() <<std::endl
  //           << input_sector.ket_subspace().DebugStr() <<std::endl;
  // std::cout << "target" <<std::endl
  //           << target_sector.ket_subspace().LabelStr() <<std::endl
  //           << target_sector.ket_subspace().DebugStr() <<std::endl;

  // read matrix for sector
  basis::OperatorBlock<double> input_matrix;
  stream_ptr->ReadSector(input_sector_index, input_matrix);

  // fill in lower triangle of matrix
  //
  // The full square matrix must be populated before remapping.
  if (input_sector.IsDiagonal())
    mcutils::CompleteLowerTriangle(input_matrix);

  // remap input matrix to target indexing
  basis::OperatorBlock<double> matrix
    = shell::RemappedMatrixJJJPN(
        input_sector,
        target_sector,
        two_body_mapping,
        input_matrix
      );

  // save result
  source_matrices[id] = matrix;
}

void BuiltinChannel::PopulateSourceMatrix(
    const RunParameters &run_parameters,
    const XformMap &xforms,
    const OneBodyOperatorMap &one_body_operators,
    OperatorBlockMap &source_matrices,
    const typename basis::TwoBodySectorsJJJPN::SectorType &target_sector
  )
{
  // matrix for sector
  basis::OperatorBlock<double> operator_matrix;

  if (id=="identity")
    // identity
    {
      operator_matrix = shell::IdentityOperatorMatrixJJJPN(target_sector,run_parameters.A);
    }
  else if (id=="loop-test")
    // loop timing test
    {
      operator_matrix = shell::TimingTestMatrixJJJPN(target_sector,true,true);
    }

  // save result
  // std::cout << fmt::format("Saving {}...", id) << std::endl;
  source_matrices[id] = operator_matrix;
}

void OperatorUChannel::PopulateSourceMatrix(
    const RunParameters &run_parameters,
    const XformMap &xforms,
    const OneBodyOperatorMap &one_body_operators,
    OperatorBlockMap &source_matrices,
    const typename basis::TwoBodySectorsJJJPN::SectorType &target_sector
  )
{
  // read matrix for sector
  basis::OperatorBlock<double> operator_matrix;

  const OneBodyOperatorData& ob_data = one_body_operators.at(ob_source_id);
  operator_matrix = shell::UpgradeOneBodyOperatorJJJPN(
      ob_data.orbital_space,
      ob_data.sectors,
      ob_data.matrices,
      target_sector,
      run_parameters.A
    );

  // save result
  // std::cout << fmt::format("Saving {}...", id) << std::endl;
  source_matrices[id] = operator_matrix;
}

void OperatorVChannel::PopulateSourceMatrix(
    const RunParameters &run_parameters,
    const XformMap &xforms,
    const OneBodyOperatorMap &one_body_operators,
    OperatorBlockMap &source_matrices,
    const typename basis::TwoBodySectorsJJJPN::SectorType &target_sector
  )
{
  // read matrix for sector
  basis::OperatorBlock<double> operator_matrix;

  const OneBodyOperatorData& ob_data_a = one_body_operators.at(ob_factor_a_id);
  const OneBodyOperatorData& ob_data_b = one_body_operators.at(ob_factor_b_id);
  operator_matrix = shell::RacahReduceTensorProductJJJPN(
      ob_data_a.orbital_space,
      ob_data_a.sectors,
      ob_data_a.matrices,
      ob_data_b.sectors,
      ob_data_b.matrices,
      target_sector,
      run_parameters.J0
    );
  if (scale_factor != 1.0)
    operator_matrix *= scale_factor;

  // save result
  // std::cout << fmt::format("Saving {}...", id) << std::endl;
  source_matrices[id] = operator_matrix;
}

void XformChannel::PopulateSourceMatrix(
    const RunParameters& run_parameters,
    const XformMap& xforms,
    const OneBodyOperatorMap& one_body_operators,
    OperatorBlockMap& source_matrices,
    const typename basis::TwoBodySectorsJJJPN::SectorType& target_sector
  )
{
  // locate corresponding input sector
  std::size_t input_bra_subspace_index
    = stream_ptr->space().LookUpSubspaceIndex(target_sector.bra_subspace().labels());
  std::size_t input_ket_subspace_index
    = stream_ptr->space().LookUpSubspaceIndex(target_sector.ket_subspace().labels());
  std::size_t input_sector_index
    = stream_ptr->sectors().LookUpSectorIndex(input_bra_subspace_index,input_ket_subspace_index);

  // locate corresponding intermediate pre-xform (truncated) sector
  std::size_t pre_xform_bra_subspace_index
    = pre_xform_two_body_indexing.space.LookUpSubspaceIndex(target_sector.bra_subspace().labels());
  std::size_t pre_xform_ket_subspace_index
    = pre_xform_two_body_indexing.space.LookUpSubspaceIndex(target_sector.ket_subspace().labels());
  std::size_t pre_xform_sector_index
    = pre_xform_two_body_indexing.sectors.LookUpSectorIndex(pre_xform_bra_subspace_index,pre_xform_ket_subspace_index);

  // short circuit if no corresponding input or intermediate sectors
  if (input_sector_index == basis::kNone)
    return;
  if (pre_xform_sector_index == basis::kNone)
    return;

  // set up aliases for input and intermediate pre-xform sectors
  const typename basis::TwoBodySectorsJJJPN::SectorType& input_sector
    = stream_ptr->sectors().GetSector(input_sector_index);
  const typename basis::TwoBodySectorsJJJPN::SectorType& pre_xform_sector
    = pre_xform_two_body_indexing.sectors.GetSector(pre_xform_sector_index);

  // read matrix for sector
  basis::OperatorBlock<double> input_matrix;
  stream_ptr->ReadSector(input_sector_index, input_matrix);

  // fill in lower triangle of matrix
  //
  // The full square matrix must be populated before remapping.
  if (input_sector.IsDiagonal())
    mcutils::CompleteLowerTriangle(input_matrix);

  // remap input matrix to pre-xform indexing
  //
  // This is the "truncation cut" on the two-body transformation.
  basis::OperatorBlock<double> pre_xform_matrix
    = shell::RemappedMatrixJJJPN(
        input_sector,
        pre_xform_sector,
        pre_xform_two_body_mapping,
        input_matrix
      );

  // look up xform data
  //
  // nonexistence check currently suppressed, as
  // should have already occurred

  // if (!xforms.count(xform_id))
  // {
  //   std::cerr << fmt::format("ERROR {}: xform {} not found.", __LINE__, xform_id)
  //             << std::endl;
  //   std::exit(EXIT_FAILURE);
  // }
  const XformData& xform_data = xforms.at(xform_id);

  // do the transformation
  basis::OperatorBlock<double> matrix
    = shell::TwoBodyTransformedMatrix(
        xform_data.bra_orbital_space,
        xform_data.ket_orbital_space,
        xform_data.sectors,
        xform_data.matrices,
        pre_xform_sector,target_sector,
        pre_xform_matrix
      );

  // save result
  source_matrices[id] = matrix;
}

void TargetChannel::GenerateTargetMatrix(
    const OperatorBlockMap& source_matrices,
    std::size_t target_sector_index,
    const typename basis::TwoBodySectorsJJJPN::SectorType& target_sector
  )
{
  // accumulate target matrix
  basis::OperatorBlock<double> target_matrix = Eigen::MatrixXd::Zero(target_sector.bra_subspace().size(),target_sector.ket_subspace().size());
  for (const auto& id_coefficient : coefficients)
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
          const basis::OperatorBlock<double>& source_matrix = pos->second;
          target_matrix += coefficient * source_matrix;
        }
    }

  // write target matrix
  mcutils::ChopMatrix(target_matrix);  // "neaten" output by eliminating near-zero values
  stream_ptr->WriteSector(target_sector_index, target_matrix);
}

////////////////////////////////////////////////////////////////
// termination code
////////////////////////////////////////////////////////////////

void InputChannel::CloseChannel()
{
  // close stream
  stream_ptr->Close();

  // deallocate
  delete stream_ptr;
}

void BuiltinChannel::CloseChannel()
{}

void OperatorUChannel::CloseChannel()
{}

void OperatorVChannel::CloseChannel()
{}

void XformChannel::CloseChannel()
{
  // close stream
  stream_ptr->Close();

  // deallocate
  delete stream_ptr;
}

void TargetChannel::CloseChannel()
{
  // close stream
  stream_ptr->Close();

  // deallocate
  delete stream_ptr;
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
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  // read parameters
  RunParameters run_parameters;
  XformMap xforms;
  OneBodyChannelList one_body_channels;
  OneBodyOperatorMap one_body_operators;
  TwoBodySourceChannelList source_channels;
  TwoBodyTargetChannelList target_channels;
  ReadParameters(
      run_parameters, xforms, one_body_channels, source_channels, target_channels
    );

  // set up parallelization
  std::cout
    << fmt::format("Parallelization: max_threads {}, num_procs {}",
                   omp_get_max_threads(), omp_get_num_procs()
      )
    << std::endl
    << std::endl;
  Eigen::initParallel();

  // start timing
  mcutils::SteadyTimer total_time;
  total_time.Start();

  ////////////////////////////////////////////////////////////////
  // operator setup
  ////////////////////////////////////////////////////////////////

  // set up target indexing
  TwoBodyOperatorIndexing target_indexing;
  InitializeTargetIndexing(run_parameters, target_indexing);

  InitializeXforms(xforms);
  // initialize one-body operators
  for (auto& one_body_channel_ptr : one_body_channels)
  {
    one_body_channel_ptr->ConstructOneBodyOperatorData(
        run_parameters, xforms, one_body_operators, target_indexing
      );
  }

  if (one_body_channels.size() > 0)
    std::cout << std::endl;

  // set up two-body channels
  for (auto& two_body_channel_ptr : source_channels)
  {
    two_body_channel_ptr->InitializeChannel(
        run_parameters, xforms, one_body_operators, target_indexing
      );
  }
  for (auto& two_body_channel_ptr : target_channels)
  {
    two_body_channel_ptr->InitializeChannel(run_parameters, target_indexing);
  }

  ////////////////////////////////////////////////////////////////
  // copy sectors
  ////////////////////////////////////////////////////////////////

  // iterate over sectors
  for (std::size_t sector_index = 0; sector_index < target_indexing.sectors.size(); ++sector_index)
    {
      // alias sector
      const auto& target_sector = target_indexing.sectors.GetSector(sector_index);


      // generate sources
      OperatorBlockMap source_matrices;  // map id->matrix
      for (auto& two_body_channel_ptr : source_channels)
      {
        two_body_channel_ptr->PopulateSourceMatrix(
            run_parameters, xforms, one_body_operators,
            source_matrices, target_sector
          );
      }

      // generate targets
      for (auto& two_body_channel_ptr : target_channels)
      {
        two_body_channel_ptr->GenerateTargetMatrix(
            source_matrices, sector_index, target_sector
          );
      }

      // progress indicator
      std::cout << "." << std::flush;

    }
  std::cout << std::endl;


  ////////////////////////////////////////////////////////////////
  // termination
  ////////////////////////////////////////////////////////////////

  // explicitly force stream closure
  //   for neatness
  for (auto& two_body_channel_ptr : source_channels)
  {
    two_body_channel_ptr->CloseChannel();
  }
  for (auto& two_body_channel_ptr : target_channels)
  {
    two_body_channel_ptr->CloseChannel();
  }

  // end timing
  total_time.Stop();
  std::cout << "(Total time: " << total_time.ElapsedTime() << ")" << std::endl;
  std::cout << std::endl;

  // exit
  return EXIT_SUCCESS;
}
