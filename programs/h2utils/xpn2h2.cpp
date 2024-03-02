/******************************************************************************

  xpn2h2.cpp -- convert BIGSTICK xpn to h2

  Syntax:
    xpn2h2 h2_format N1max N2max sps_filename xpn_filename orbitals_filename h2_filename

    programs/h2utils/xpn2h2 15099 6 6 ${HOME}/data/xpn/runmac0493-JISP16-tb-6/ncci-tb-6.sps ${HOME}/data/xpn/runmac0493-JISP16-tb-6/JISP16-tb-6-20.int orbitals.dat JISP16-tb-6-20.dat

  Mark A. Caprio
  University of Notre Dame

  + 03/01/23 (mac): Created, based on me2j2h2.

******************************************************************************/

#include <fstream>

#include "basis/jjjt_operator.h"
#include "basis/jjjpn_scheme.h"
#include "basis/jjjpn_operator.h"
#include "mcutils/eigen.h"
#include "mcutils/parsing.h"
#include "mcutils/profiling.h"
#include "tbme/h2_io.h"
#include "tbme/me2j_io.h"
#include "tbme/tbme_scheme_xform.h"

#include "moshinsky/moshinsky_xform.h"

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

struct RunParameters
// Stores simple parameters for run
{
  // filenames
  std::string sps_filename;
  std::string xpn_filename;
  std::string orbitals_filename;
  std::string h2_filename;
  // mode
  shell::H2Format output_h2_format;
  // truncation
  basis::WeightMax weight_max;
};

void ProcessArguments(int argc, char **argv, RunParameters& run_parameters)
{
  // usage message
  if (argc-1 < 7)
    {
      std::cout << "Syntax: xpn2h2 h2_format N1max N2max sps_filename xpn_filename orbitals_filename h2_filename" << std::endl;
      std::exit(EXIT_SUCCESS);
    }

  // format
  {
    std::istringstream parameter_stream(argv[1]);
    parameter_stream >> run_parameters.output_h2_format;
    if (!parameter_stream)
      {
        std::cerr << "Expecting numeric value for H2 format argument" << std::endl;
        std::exit(EXIT_FAILURE);
      }
  }

  // truncation
  int N1max, N2max;
  {
    std::istringstream parameter_stream(argv[2]);
    parameter_stream >> N1max;
    if (!parameter_stream)
      {
        std::cerr << "Expecting numeric value for N1max argument" << std::endl;
        std::exit(EXIT_FAILURE);
      }
  }
  {
    std::istringstream parameter_stream(argv[3]);
    parameter_stream >> N2max;
    if (!parameter_stream)
      {
        std::cerr << "Expecting numeric value for N2max argument" << std::endl;
        std::exit(EXIT_FAILURE);
      }
  }
  run_parameters.weight_max = basis::WeightMax(N1max, N2max);


  // input sps filename
  run_parameters.sps_filename = argv[4];

  // input xpn filename
  run_parameters.xpn_filename = argv[5];

  // output orbitals filename
  run_parameters.orbitals_filename = argv[6];

  // output h2 filename
  run_parameters.h2_filename = argv[7];
  
}

basis::OrbitalSpacePN ReadOrbitals(const std::string& filename)
// Read BIGSTICK sps file for single-particle space.
//
// See Sec 4.2 "Defining the model space" of Johnson arXiv:1801.08432.
{

  // open sps file
  std::ifstream is(filename);
  if (!is)
    {
      std::cout << "ERROR: Failure opening sps file" << std::endl;
      std::exit(EXIT_SUCCESS);
    }

  // set up line counter for use in error messages
  std::string line;
  int line_count = 0;

  // parse header
  // line 1: mode
  // Only "iso" mode supported.
  {
    mcutils::GetLine(is,line,line_count);
    std::istringstream line_stream(line);
    std::string mode;
    line_stream >> mode;
    mcutils::ParsingCheck(line_stream,line_count,line);
    assert(mode=="iso");
  }

  // line 2: number of orbitals (same for p and n)
  int num_orbitals;
  {
    mcutils::GetLine(is,line,line_count);
    std::istringstream line_stream(line);
    std::string mode;
    line_stream >> num_orbitals;
    mcutils::ParsingCheck(line_stream,line_count,line);
    assert(num_orbitals>0);
  }

  // read in orbital info
  basis::OrbitalPNList states;
  for (int orbital_index=0; orbital_index<num_orbitals; ++orbital_index)
    {
      mcutils::GetLine(is,line,line_count);
      std::istringstream line_stream(line);
      std::string mode;
      float n_raw, l_raw, j_raw, w_raw;
      line_stream >> n_raw >> l_raw >> j_raw >> w_raw;
      mcutils::ParsingCheck(line_stream,line_count,line);
      //std::cout << fmt::format("{} {} {} {}", n_raw, l_raw, j_raw, w_raw) << std::endl;

      // convert orbital parameters
      int n = int(n_raw);
      int l = int(l_raw);
      HalfInt j = HalfInt(2*j_raw, 2);
      float weight = w_raw;

      // store as proton orbital
      basis::OrbitalPNInfo proton_orbital(basis::OrbitalSpeciesPN::kP, n, l, j, weight);
      states.push_back(proton_orbital);
    }
  
  // append identical list of neutron orbitals
  for (int orbital_index=0; orbital_index<num_orbitals; ++orbital_index)
    {
      const basis::OrbitalPNInfo& proton_orbital = states[orbital_index];
      basis::OrbitalPNInfo neutron_orbital(basis::OrbitalSpeciesPN::kN, proton_orbital.n, proton_orbital.l, proton_orbital.j, proton_orbital.weight);
      states.push_back(neutron_orbital);
    }
  basis::OrbitalSpacePN orbital_space(states);

  return orbital_space;
}


int main(int argc, char **argv)
{
  // header
  std::cout << std::endl;
  std::cout << "xpn2h2  -- xpn to h2 conversion" << std::endl;
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  // read parameters
  RunParameters run_parameters;
  ProcessArguments(argc,argv,run_parameters);

  // read orbitals from sps file
  basis::OrbitalSpacePN orbital_space = ReadOrbitals(run_parameters.sps_filename);
  
  // write orbitals
  std::ofstream os(run_parameters.orbitals_filename);
  os << basis::OrbitalDefinitionStr(orbital_space.OrbitalInfo(), true, basis::MFDnOrbitalFormat::kVersion15099);

  // set up operator storage
  // Note: Both xpn and h2 use NAS storage, so use NAS internally as well.
  int J0 = 0;
  int g0 = 0;
  int Tz0 = 0;
  const basis::TwoBodySpaceJJJPN two_body_space(
      orbital_space,
      run_parameters.weight_max,
      shell::kH2SpaceOrdering.at(run_parameters.output_h2_format)
    );
  const basis::TwoBodySectorsJJJPN two_body_sectors(
      two_body_space,
      J0, g0, Tz0
    );
  basis::OperatorBlocks<double> two_body_matrices;
  basis::SetOperatorToZero(two_body_sectors, two_body_matrices);

  // read xpn file
  // ReadXPNFile(two_body_space, two_body_sectors, two_body_matrices, run_parameters.input_filename);
  
  // write h2 file
  shell::OutH2Stream output_stream(
      run_parameters.h2_filename,
      orbital_space, two_body_space, two_body_sectors,
      run_parameters.output_h2_format
    );
  std::cout << output_stream.DiagnosticStr();

  for (std::size_t sector_index=0; sector_index<two_body_sectors.size(); ++sector_index)
    {
      // make reference to target sector
      const basis::TwoBodySectorsJJJPN::SectorType& two_body_sector
        = two_body_sectors.GetSector(sector_index);
      const auto& matrix = two_body_matrices[sector_index];
        
      output_stream.WriteSector(
          sector_index,
          matrix,
          basis::NormalizationConversion::kNone
        );
      std::cout << "." << std::flush;
    }
  output_stream.Close();
  std::cout << std::endl;

  std::exit(EXIT_SUCCESS);
}
