/******************************************************************************
  @file obmixer.cpp

  Compute one-body operator matrix elements.

  Note: All matrix elements are defined in terms for unitless operators,
  i.e., r->(b^-1 r) and ik->(b ik). Equivalently, they are evaluated between
  radial functions with length parameter 1. Any alternate scaling can be
  accomplished by providing a length_parameter.

  Syntax:
    + obmixer

  Input format:
    set-basis <basis_type> <orbital_filename>
      basis_type = oscillator|laguerre
    set-length-parameter <length_parameter>
    define-xform <id> <xform_filename>
    define-source <mode> <id> ...
      define-source kinematic <id> [orbital_filename]
        id = identity|r|k|r.r|k.k
      define-source am <id> [orbital_filename]
        id = l|l2|s|s2|j|j2
      define-source isospin <id> [orbital_filename]
        id = tz|t+|t-
      define-source ladder <id> [orbital_filename]
        id = c+|c
      define-source solid-harmonic <id> <coordinate> <order> [orbital_filename]
        coordinate = r|ik
      define-source input <id> <filename> <j0> <g0> <tz0>
      define-source linear-combination <id>
        add-source <id> <coefficient>
      define-source tensor-product <id> <ob_factor_a_id> <ob_factor_b_id> <j0> [scale_factor]
      define-source xform <id> <ob_source_id> <xform_id>
    define-target <id> <filename>
  @note Currently only computes operator/radial matrix elements between harmonic oscillator
  or Laguerre basis functions with identical bra and ket spaces.

  Patrick J. Fasano
  University of Notre Dame

  + 09/08/20 (pjf): Created, based on h2mixer and radial-gen.
  + 09/10/20 (pjf): Add length parameter option.
  + 09/13/20 (pjf): Fix xform source.
  + 10/09/20 (pjf): Fix output header line.

******************************************************************************/

#include <omp.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>
// #include <set>
#include <variant>
#include <string>

#include "am/am.h"
#include "basis/nlj_orbital.h"
#include "basis/operator.h"
#include "fmt/format.h"
#include "mcutils/eigen.h"
#include "mcutils/profiling.h"
#include "obme/obme.h"
#include "obme/obme_io.h"
#include "obme/radial.h"

////////////////////////////////////////////////////////////////
// parameter definitions
////////////////////////////////////////////////////////////////

// convenience typedef for operator quantum numbers
typedef std::tuple<int,int,int> OperatorQN;

// isospin operator type enum (for std::variant)
enum class IsospinOperatorType { kRaising, kLowering, kProjection };

// variant type for operator type
typedef std::variant<
    shell::RadialOperatorType,
    am::AngularMomentumOperatorType,
    shell::LadderOperatorType,
    IsospinOperatorType
    >
  OperatorType;

// operator maps
std::unordered_map<std::string,std::tuple<OperatorType,int,OperatorQN>>
kKinematicOneBodyOperatorDefinitions =
  {
    // {id, {operator type, order, qn}}
    {"identity", {shell::RadialOperatorType::kO, 0, {0,0,0}}},
    {"r",        {shell::RadialOperatorType::kR, 1, {1,1,0}}},
    {"ik",       {shell::RadialOperatorType::kK, 1, {1,1,0}}},
    {"r.r",      {shell::RadialOperatorType::kR, 2, {0,0,0}}},
    {"ik.ik",    {shell::RadialOperatorType::kK, 2, {0,0,0}}},
  };

std::unordered_map<std::string,std::tuple<OperatorType,int,OperatorQN>>
kAngularMomentumOneBodyOperatorDefinitions =
  {
    // {id, {operator type, order, qn}}
    {"l",  {am::AngularMomentumOperatorType::kOrbital, 1, {1,0,0}}},
    {"l2", {am::AngularMomentumOperatorType::kOrbital, 2, {0,0,0}}},
    {"s",  {am::AngularMomentumOperatorType::kSpin,    1, {1,0,0}}},
    {"s2", {am::AngularMomentumOperatorType::kSpin,    2, {0,0,0}}},
    {"j",  {am::AngularMomentumOperatorType::kTotal,   1, {1,0,0}}},
    {"j2", {am::AngularMomentumOperatorType::kTotal,   2, {0,0,0}}}
  };

std::unordered_map<std::string,std::tuple<OperatorType,int,OperatorQN>>
kIsospinOneBodyOperatorDefinitions =
  {
    {"tz", {IsospinOperatorType::kProjection, 1, {0,0, 0}}},
    {"t+", {IsospinOperatorType::kRaising,    1, {0,0,+1}}},
    {"t-", {IsospinOperatorType::kLowering,   1, {0,0,-1}}}
  };

std::unordered_map<std::string,std::tuple<OperatorType,int,OperatorQN>>
kLadderOneBodyOperatorDefinitions =
  {
    {"c+", {shell::LadderOperatorType::kRaising,  1, {1,1,0}}},
    {"c",  {shell::LadderOperatorType::kLowering, 1, {1,1,0}}},
  };

std::unordered_map<std::string, shell::RadialBasisType>
kBasisTypeDefinitions =
  {
    {"oscillator", shell::RadialBasisType::kOscillator},
    {"laguerre", shell::RadialBasisType::kLaguerre}
  };

////////////////////////////////////////////////////////////////
// run control
////////////////////////////////////////////////////////////////

struct RunParameters
// Structure to store input parameters for run.
{
  RunParameters()
    : basis_type(shell::RadialBasisType::kOscillator), length_parameter(1.0), orbital_space()
  {};

  void InitializeIndexing(
      const shell::RadialBasisType& basis_type_,
      const std::string& orbital_filename_
    )
  {
    basis_type = basis_type_;
    std::ifstream orbital_file(orbital_filename_);
    basis::OrbitalPNList orbital_info = basis::ParseOrbitalPNStream(orbital_file,true);
    orbital_space = basis::OrbitalSpaceLJPN(orbital_info);
  }

  // basis parameters
  shell::RadialBasisType basis_type;
  double length_parameter;
  basis::OrbitalSpaceLJPN orbital_space;
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
      OneBodyOperatorMap& operators
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
      const OperatorQN& operator_qn_
    )
    : OneBodySourceChannel(id_), filename(filename_), operator_qn(operator_qn_)
  {}

  void ConstructOneBodyOperatorData(
      const RunParameters& run_parameters,
      const XformMap& xforms,
      OneBodyOperatorMap& operators
    ) override;

  // input filename
  std::string filename;

  // operator identification
  OperatorQN operator_qn;
};

struct OneBodyGeneratedChannel : public OneBodySourceChannel
// Parameters for a source channel providing an operator generated internally.
{
  OneBodyGeneratedChannel() = default;

  OneBodyGeneratedChannel(
      const std::string& id_,
      const OperatorType& operator_type_,
      const int& operator_order_,
      const OperatorQN& operator_qn_,
      const std::string& orbital_filename_=""
    )
    : OneBodySourceChannel(id_), operator_type(operator_type_),
      operator_order(operator_order_), operator_qn(operator_qn_),
      orbital_filename(orbital_filename_)
  {}

  void ConstructOneBodyOperatorData(
      const RunParameters& run_parameters,
      const XformMap& xforms,
      OneBodyOperatorMap& operators
    ) override;

  // operator properties
  OperatorType operator_type;
  int operator_order;
  OperatorQN operator_qn;

  // indexing (optional)
  std::string orbital_filename;
};

struct OneBodyLinearCombinationChannel : public OneBodySourceChannel
// Parameters for a source channel providing a linear combination of one-body
// operators.
{
  OneBodyLinearCombinationChannel() = default;

  explicit OneBodyLinearCombinationChannel(const std::string& id_)
    : OneBodySourceChannel(id_), coefficients()
  {}

  void ConstructOneBodyOperatorData(
      const RunParameters& run_parameters,
      const XformMap& xforms,
      OneBodyOperatorMap& operators
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
      int j0_, double scale_factor_
    )
    : OneBodySourceChannel(id_),
      ob_factor_a_id(factor_a_id_), ob_factor_b_id(factor_b_id_),
      j0(j0_), scale_factor(scale_factor_)
  {}

  void ConstructOneBodyOperatorData(
      const RunParameters& run_parameters,
      const XformMap& xforms,
      OneBodyOperatorMap& operators
    ) override;

  // construction
  std::string ob_factor_a_id, ob_factor_b_id;
  int j0;
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
      OneBodyOperatorMap& operators
    ) override;

  // construction
  std::string ob_source_id, xform_id;
};

struct TargetChannel
// Target channel parameters and stream.
{
  TargetChannel(const std::string& id_, const std::string& filename_)
    : id(id_), filename(filename_)
  {}

  // output
  std::string id;
  std::string filename;

  // write channel to file
  void WriteChannel(
      const RunParameters& run_parameters,
      const OneBodyOperatorMap& operators
    );
};

typedef std::vector<std::unique_ptr<TargetChannel>> OneBodyTargetChannelList;
// Vector of pointers to base for two-body source channels.

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

void ReadParameters(
    RunParameters& run_parameters,
    XformMap& xforms,
    OneBodyChannelList& one_body_channels,
    OneBodyTargetChannelList& target_channels
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
    if (keyword == "set-basis")
    {
      // set-ket-basis basis_type orbital_filename
      //   basis_type = oscillator|laguerre
      std::string basis_type_str, orbital_filename;
      line_stream >> basis_type_str >> orbital_filename;
      mcutils::ParsingCheck(line_stream, line_count, line);
      // convert basis
      if (kBasisTypeDefinitions.count(basis_type_str) == 0)
      {
        mcutils::ParsingError(
            line_count, line,
            "Valid analytic basis types: oscillator|laguerre"
          );
      }
      auto basis_type = kBasisTypeDefinitions.at(basis_type_str);
      mcutils::FileExistCheck(orbital_filename, true, false);
      run_parameters.InitializeIndexing(basis_type, orbital_filename);
    }
    else if (keyword=="set-length-parameter")
    {
      line_stream >> run_parameters.length_parameter;
      mcutils::ParsingCheck(line_stream, line_count, line);
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
    else if (keyword=="define-source")
    {
      std::string mode, id, orbital_filename="";
      line_stream >> mode >> id;
      mcutils::ParsingCheck(line_stream,line_count,line);

      if (mode == "kinematic")
      {
        if (!kKinematicOneBodyOperatorDefinitions.count(id))
          mcutils::ParsingError(line_count,line,"Unrecognized kinematic operator ID");
        auto [operator_type, operator_order, operator_qn]
          = kKinematicOneBodyOperatorDefinitions.at(id);
        if (!line_stream.eof())
        {
          line_stream >> orbital_filename;
          mcutils::ParsingCheck(line_stream,line_count,line);
        }
        one_body_channels.emplace_back(
            new OneBodyGeneratedChannel(
                id, operator_type, operator_order, operator_qn, orbital_filename
              )
          );
      }
      else if (mode == "am")
      {
        std::string operator_species_s;
        if (!kAngularMomentumOneBodyOperatorDefinitions.count(id))
          mcutils::ParsingError(line_count,line,"Unrecognized angular momentum operator ID");
        auto [operator_type, operator_order, operator_qn]
          = kAngularMomentumOneBodyOperatorDefinitions.at(id);
        if (!line_stream.eof())
        {
          line_stream >> orbital_filename;
          mcutils::ParsingCheck(line_stream,line_count,line);
        }
        one_body_channels.emplace_back(
            new OneBodyGeneratedChannel(
                id, operator_type, operator_order, operator_qn, orbital_filename
              )
          );
      }
      else if (mode == "isospin")
      {
        if (!kIsospinOneBodyOperatorDefinitions.count(id))
          mcutils::ParsingError(line_count,line,"Unrecognized isospin operator ID");
        auto [operator_type, operator_order, operator_qn]
          = kIsospinOneBodyOperatorDefinitions.at(id);
        if (!line_stream.eof())
        {
          line_stream >> orbital_filename;
          mcutils::ParsingCheck(line_stream,line_count,line);
        }
        one_body_channels.emplace_back(
            new OneBodyGeneratedChannel(
                id, operator_type, operator_order, operator_qn, orbital_filename
              )
          );
      }
      else if (mode == "ladder")
      {
        if (!kLadderOneBodyOperatorDefinitions.count(id))
          mcutils::ParsingError(line_count,line,"Unrecognized isospin operator ID");
        auto [operator_type, operator_order, operator_qn]
          = kLadderOneBodyOperatorDefinitions.at(id);
        if (!line_stream.eof())
        {
          line_stream >> orbital_filename;
          mcutils::ParsingCheck(line_stream,line_count,line);
        }
        one_body_channels.emplace_back(
            new OneBodyGeneratedChannel(
                id, operator_type, operator_order, operator_qn, orbital_filename
              )
          );
      }
      else if (mode == "solid-harmonic")
      {
        std::string coordinate;
        int operator_order, j0, g0, tz0=0;
        line_stream >> coordinate >> operator_order >> j0;
        mcutils::ParsingCheck(line_stream,line_count,line);

        OperatorType operator_type = shell::kCharCodeRadialOperatorType.at(coordinate);
        g0 = j0%2;
        tz0 = 0;
        if (!line_stream.eof())
        {
          line_stream >> orbital_filename;
          mcutils::ParsingCheck(line_stream,line_count,line);
        }
        one_body_channels.emplace_back(
            new OneBodyGeneratedChannel(
                id, operator_type, operator_order, {j0,g0,tz0}, orbital_filename
              )
          );
      }
      else if (mode=="input")
      {
        std::string filename;
        int j0, g0, tz0;
        line_stream >> filename >> j0 >> g0 >> tz0;
        mcutils::ParsingCheck(line_stream,line_count,line);
        one_body_channels.emplace_back(
            new OneBodyInputChannel(id, filename, {j0, g0, tz0})
          );
      }
      else if (mode=="linear-combination")
      {
        mcutils::ParsingCheck(line_stream,line_count,line);

        one_body_channels.emplace_back(
            new OneBodyLinearCombinationChannel(id)
          );
      }
      else if (mode=="tensor-product")
      {
        std::string ob_factor_a_id, ob_factor_b_id;
        int j0;
        double scale_factor = 1.0;
        line_stream >> ob_factor_a_id >> ob_factor_b_id >> j0;
        mcutils::ParsingCheck(line_stream,line_count,line);
        if (!line_stream.eof())
          {
            line_stream >> scale_factor;
            mcutils::ParsingCheck(line_stream,line_count,line);
          }

        one_body_channels.emplace_back(
            new OneBodyTensorProductChannel(
                id, ob_factor_a_id, ob_factor_b_id, j0, scale_factor
              )
          );
      }
      else if (mode=="xform")
      {
        std::string source_id, xform_id;
        line_stream >> source_id >> xform_id;
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
    else if (keyword=="add-source")
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
    else if (keyword=="define-target")
    {
      std::string id, filename;
      line_stream >> id >> filename;
      mcutils::ParsingCheck(line_stream,line_count,line);
      target_channels.emplace_back(
          new TargetChannel(id, filename)
        );
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
    assert(xform_data.sectors.j0() == 0);
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
    OneBodyOperatorMap& operators
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
  basis::OrbitalSpaceLJPN ket_orbital_space;
  shell::InOBMEStream operator_stream(filename);
  operator_stream.SetToIndexing(
      operator_data.orbital_space, ket_orbital_space, operator_data.sectors
    );
  assert(operator_data.orbital_space.OrbitalInfo() == ket_orbital_space.OrbitalInfo());
  auto& [j0, g0, tz0] = operator_qn;
  assert(j0 == operator_data.sectors.j0());
  assert(g0 == operator_data.sectors.g0());
  assert(tz0 == operator_data.sectors.Tz0());

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

void OneBodyGeneratedChannel::ConstructOneBodyOperatorData(
    const RunParameters& run_parameters,
    const XformMap& xforms,
    OneBodyOperatorMap& operators
  )
// Populate map with data from built-in generation routines.
{
  std::cout
    << fmt::format(
        "Constructing one-body matrix elements for operator {}...", id
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

  // initialize indexing
  if (orbital_filename == "")
  {
    operator_data.orbital_space = run_parameters.orbital_space;
  }
  else
  {
    std::ifstream orbital_file(orbital_filename);
    basis::OrbitalPNList orbital_info = basis::ParseOrbitalPNStream(orbital_file,true);
    operator_data.orbital_space = basis::OrbitalSpaceLJPN(orbital_info);
  }
  auto [j0, g0, tz0] = operator_qn;
  operator_data.sectors = basis::OrbitalSectorsLJPN(
      operator_data.orbital_space,
      j0, g0, tz0
    );

  if (id == "identity")
  // identity operator
  {
    basis::SetOperatorToIdentity(operator_data.sectors, operator_data.matrices);
  }
  else if (std::holds_alternative<shell::RadialOperatorType>(operator_type))
  {
    auto& radial_operator_type = std::get<shell::RadialOperatorType>(operator_type);
    shell::SolidHarmonicOneBodyOperator(
        run_parameters.basis_type,
        radial_operator_type,
        operator_order,
        operator_data.orbital_space,
        operator_data.sectors,
        operator_data.matrices
      );

    // scaling for radial power
    double scale_factor = 1.0;
    if (radial_operator_type == shell::RadialOperatorType::kR)
      scale_factor *= std::pow(run_parameters.length_parameter, operator_order);
    else if (radial_operator_type == shell::RadialOperatorType::kK)
      scale_factor *= std::pow(run_parameters.length_parameter, -operator_order);

    // convert between C and Y if not using a special kinematic operator
    if (!kKinematicOneBodyOperatorDefinitions.count(id))
      scale_factor *= Hat(operator_order) * am::kInvSqrt4Pi;

    // apply scale factor
    basis::ScalarMultiplyOperator(
        operator_data.sectors, operator_data.matrices, scale_factor
      );
  }
  else if (std::holds_alternative<am::AngularMomentumOperatorType>(operator_type))
  // angular momentum operators
  {
    auto& am_operator_type = std::get<am::AngularMomentumOperatorType>(operator_type);

    // populate operator
    if (operator_order == 1)
    {
      shell::AngularMomentumOneBodyOperator(
          am_operator_type,
          operator_data.orbital_space,
          operator_data.sectors,
          operator_data.matrices
        );
    }
    else if (operator_order == 2)
    {
      shell::AngularMomentumSquaredOneBodyOperator(
          am_operator_type,
          operator_data.orbital_space,
          operator_data.sectors,
          operator_data.matrices
        );
    }
  }
  else if (std::holds_alternative<IsospinOperatorType>(operator_type))
  {
    shell::IsospinOneBodyOperator(
        operator_data.orbital_space,
        operator_data.sectors,
        operator_data.matrices
      );
  }
  else if (std::holds_alternative<shell::LadderOperatorType>(operator_type))
  {
    auto& ladder_operator_type = std::get<shell::LadderOperatorType>(operator_type);
    LadderOneBodyOperator(
        run_parameters.basis_type,
        ladder_operator_type,
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

  // diagnostic
  std::cout
    << fmt::format(
        "  {} sectors; {} matrix elements",
        operator_data.sectors.size(),
        basis::AllocatedEntries(operator_data.matrices)
      )
    << std::endl;
}

void OneBodyLinearCombinationChannel::ConstructOneBodyOperatorData(
    const RunParameters& run_parameters,
    const XformMap& xforms,
    OneBodyOperatorMap& operators
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
          int j0 = source_data.sectors.j0();
          int g0 = source_data.sectors.g0();
          int Tz0 = source_data.sectors.Tz0();
          operator_data.sectors = basis::OrbitalSectorsLJPN(
              operator_data.orbital_space, j0, g0, Tz0
            );
          basis::SetOperatorToZero(operator_data.sectors, operator_data.matrices);
          first_term = false;
        }
      else
        // assert same quantum numbers for later terms
        {
          assert(operator_data.sectors.j0() == source_data.sectors.j0());
          assert(operator_data.sectors.g0() == source_data.sectors.g0());
          assert(operator_data.sectors.Tz0() == source_data.sectors.Tz0());
        }

      basis::OperatorLinearCombination(
          operator_data.sectors, operator_data.matrices,
          1., operator_data.matrices,
          coefficient, source_data.matrices
        );
    }

  // diagnostic
  std::cout
    << fmt::format(
        "  {} sectors; {} matrix elements",
        operator_data.sectors.size(),
        basis::AllocatedEntries(operator_data.matrices)
      )
    << std::endl;
}

void OneBodyTensorProductChannel::ConstructOneBodyOperatorData(
    const RunParameters& run_parameters,
    const XformMap& xforms,
    OneBodyOperatorMap& operators
  )
{
  std::cout
    << fmt::format(
        "Generating one-body tensor product {} = [({})({})]_{}",
        id, ob_factor_a_id, ob_factor_b_id, j0
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
  assert(am::AllowedTriangle(data_a.sectors.j0(), data_b.sectors.j0(), j0));
  int g0 = (data_a.sectors.g0() + data_b.sectors.g0())%2;
  int tz0 = data_a.sectors.Tz0() + data_b.sectors.Tz0();
  operator_data.sectors = basis::OrbitalSectorsLJPN(
      operator_data.orbital_space, j0, g0, tz0
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

  // diagnostic
  std::cout
    << fmt::format(
        "  {} sectors; {} matrix elements",
        operator_data.sectors.size(),
        basis::AllocatedEntries(operator_data.matrices)
      )
    << std::endl;
}

void OneBodyXformChannel::ConstructOneBodyOperatorData(
    const RunParameters& run_parameters,
    const XformMap& xforms,
    OneBodyOperatorMap& operators
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
      source_data.sectors.j0(), source_data.sectors.g0(), source_data.sectors.Tz0()
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

  // diagnostic
  std::cout
    << fmt::format(
        "  {} sectors; {} matrix elements",
        operator_data.sectors.size(),
        basis::AllocatedEntries(operator_data.matrices)
      )
    << std::endl;
}

void TargetChannel::WriteChannel(
    const RunParameters& run_parameters,
    const OneBodyOperatorMap& operators
  )
// Write given operator id to file, truncating if necessary.
{
  std::cout
    << fmt::format(
        "Writing one-body matrix element for operator {} to file {}...",
        id, filename
      )
    << std::endl;

  // look up operator data
  if (!operators.count(id))
  {
    std::cerr << fmt::format("ERROR {}: Operator {} not found.", __LINE__, id)
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const OneBodyOperatorData& operator_data = operators.at(id);

  if (operator_data.orbital_space.OrbitalInfo() == run_parameters.orbital_space.OrbitalInfo())
  // just open and write if operator is already on correct space
  {
    shell::OutOBMEStream os(
        filename,
        operator_data.orbital_space, operator_data.orbital_space, operator_data.sectors,
        basis::OneBodyOperatorType::kSpherical
      );
    os.Write(operator_data.matrices);
    os.Close();
  }
  else
  // transform operator (trivially) to match target indexing
  {
    // indexing for trivial transformation
    auto xform_sectors = basis::OrbitalSectorsLJPN(
        operator_data.orbital_space, run_parameters.orbital_space, 0, 0, 0
      );
    basis::OperatorBlocks<double> xform_matrices;
    shell::GenerateRadialOverlaps(
        run_parameters.basis_type,
        run_parameters.basis_type,
        1.0,
        operator_data.orbital_space,
        run_parameters.orbital_space,
        xform_sectors,
        xform_matrices
      );

    // new sectors to be written to file
    auto sectors = basis::OrbitalSectorsLJPN(
        run_parameters.orbital_space,
        operator_data.sectors.j0(), operator_data.sectors.g0(), operator_data.sectors.Tz0()
      );
    basis::OperatorBlocks<double> matrices;

    // perform trivial transformation
    shell::SimilarityTransformOperator(
        operator_data.orbital_space,
        run_parameters.orbital_space,
        xform_sectors,
        xform_matrices,
        operator_data.sectors,
        operator_data.matrices,
        sectors,
        matrices
      );

    // write to file
    shell::OutOBMEStream os(
        filename,
        run_parameters.orbital_space, run_parameters.orbital_space, sectors,
        basis::OneBodyOperatorType::kSpherical
      );
    os.Write(matrices);
    os.Close();
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
  std::cout << "obmixer -- one-body matrix element generation" << std::endl;
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  // read parameters
  RunParameters run_parameters;
  XformMap xforms;
  OneBodyChannelList one_body_channels;
  OneBodyOperatorMap one_body_operators;
  OneBodyTargetChannelList target_channels;
  ReadParameters(
      run_parameters, xforms, one_body_channels, target_channels
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

  InitializeXforms(xforms);

  // initialize one-body operators
  for (auto& one_body_channel_ptr : one_body_channels)
  {
    one_body_channel_ptr->ConstructOneBodyOperatorData(
        run_parameters, xforms, one_body_operators
      );
  }

  if (one_body_channels.size() > 0)
    std::cout << std::endl;

  // write target channels to file
  for (auto& two_body_channel_ptr : target_channels)
  {
    two_body_channel_ptr->WriteChannel(run_parameters, one_body_operators);
  }

  ////////////////////////////////////////////////////////////////
  // termination
  ////////////////////////////////////////////////////////////////

  // end timing
  total_time.Stop();
  std::cout << "(Total time: " << total_time.ElapsedTime() << ")" << std::endl;
  std::cout << std::endl;

  // exit
  return EXIT_SUCCESS;
}
