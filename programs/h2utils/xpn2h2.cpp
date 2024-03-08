/******************************************************************************

  xpn2h2.cpp -- convert BIGSTICK xpn to h2

  Syntax:
    + xpn2h2 orbital_filename input_filename output_filename

    programs/h2utils/xpn2h2 --truncation 6 6 ncci-tb-6.dat JISP16-tb-6-20.int JISP16-tb-6-20.dat

  Assumed file format:

    + comment lines beginning with hash ('#')

    + header (possibly wrapped over multiple lines):

        num_me spe1 spe2 ...

        num_me (int): number of matrix elements to follow
        spe1, ... (float): single particle energies

    + lines of form

        a b c d J T ME

        a, b, c, d (int): 1-based orbital indices; for pn matrix elements (ab)
        and (cd) are in order proton-neutron

        J (int): J

        T (int): nominally isospin, but ignored

        ME (float): matrix element

  Limitations:
    + Input single particle energies or one-body matrix elements are ignored.
    + Input XPN files with extra numbers after the SPEs (specifying the scaling) are *not* supported.
    + Output h2 format is fixed as Version 15099.

  See Sec 4.3.2 "Proton-neutron and other isospin-breaking formats" in
  C. W. Johnson et al., "BIGSTICK: A flexible configuration-interaction
  shell-model code", arXiv:1801.08432.

  Mark A. Caprio
  University of Notre Dame

  + 03/01/24 (mac): Created.
  + 03/06/24 (mac): Complete basic support for TBME conversion
  + 03/07/24 (mac):
    - Provide --truncation option to specify output weight truncation.
    - Factor out canonicalization routines to basis/jjjpn_operator.
      
******************************************************************************/

#include <fstream>

#include "basis/jjjpn_operator.h"
#include "fmt/format.h"
#include "mcutils/eigen.h"
#include "mcutils/parsing.h"
#include "tbme/h2_io.h"

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

struct RunParameters
// Stores simple parameters for run
{
  // filenames
  std::string orbital_filename;
  std::string input_filename;
  std::string output_filename;
  std::string obme_filename;
  // format
  shell::H2Format output_format;
  // output truncation
  basis::WeightMax weight_max;

  // default constructor
  RunParameters()
    : orbital_filename(""), input_filename(""), output_filename(""), obme_filename(""),
      output_format(shell::kVersion15099), weight_max(0,0)
  {}
};

void PrintUsage(const char **argv) {
  std::cout << "Usage: " << argv[0]
            << " [--truncation w1max w2max] [--obme-filename obme_file]" << std::endl
            << "       orbital_file input_file output_file"
            << std::endl;
}

void ProcessArguments(int argc, const char *argv[], RunParameters& run_parameters)
{
  // usage message
  if (argc-1 == 0)
    {
      PrintUsage(argv);
      std::exit(EXIT_SUCCESS);
    }

  int arg = 1;
  while (arg < argc && argv[arg][0] == '-')
    {
      std::istringstream parameter_stream(argv[arg++]);

      if (parameter_stream.str() == "--help" || parameter_stream.str() == "-h")
        {
          PrintUsage(argv);
          std::exit(EXIT_SUCCESS);
        }
      else if (parameter_stream.str() == "--truncation")
        {
          if (argc-arg < 2)
            {
              PrintUsage(argv);
              std::cerr << "Insufficient arguments for --truncation" << std::endl;
              std::exit(EXIT_FAILURE);
            }

          float w1max, w2max;
          std::istringstream w1max_stream(argv[arg++]);
          w1max_stream >> w1max;
          if (!w1max_stream) {
            PrintUsage(argv);
            std::cerr << "Invalid w1max" << std::endl;
            std::exit(EXIT_FAILURE);
          }
          std::istringstream w2max_stream(argv[arg++]);
          w2max_stream >> w2max;
          if (!w2max_stream) {
            PrintUsage(argv);
            std::cerr << "Invalid w2max" << std::endl;
            std::exit(EXIT_FAILURE);
          }
          run_parameters.weight_max = basis::WeightMax(w1max, w1max, w2max, w2max, w2max);
        }
      else if (parameter_stream.str() == "--obme-filename")
        {
          if (argc-arg < 1)
            {
              PrintUsage(argv);
              std::cerr << "Insufficient arguments for --obme-filename" << std::endl;
              std::exit(EXIT_FAILURE);
            }

          run_parameters.obme_filename = argv[arg++];
        }
      else
        {
          PrintUsage(argv);
          std::cerr << "Unrecognized option '" << parameter_stream.str() << "'" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        
    }

  if (argc-arg < 3)
    {
      PrintUsage(argv);
      std::cerr << "Insufficient arguments" << std::endl;
      std::exit(EXIT_FAILURE);
    }

  // orbital filename
  run_parameters.orbital_filename = argv[arg++];
  mcutils::FileExistCheck(run_parameters.orbital_filename, true, false);

  // input filename
  run_parameters.input_filename = argv[arg++];
  mcutils::FileExistCheck(run_parameters.input_filename, true, false);

  // output filename
  run_parameters.output_filename = argv[arg++];

}

struct XPNTBMEDatum
// Stores raw data from XPN file me data line
{
  int a, b, c, d, J, T;
  double me;
};

void ReadXPNFile(
    const std::string& filename, std::size_t num_orbitals,
    std::vector<double>& spe_data, std::vector<XPNTBMEDatum>& tbme_data
  )
// Read raw data from BIGSTICK XPN file.
//
// Arguments:
//   filename (std::string): XPN filename
//   num_orbitals (std::size_t): number of orbitals (per species)
//   spe_data (std::vector<double>): vector to store SPEs
//   tbme_data (std::vector<XPNTBMEDatum>): vector to store TBMEs
{

  // open xpn file
  std::ifstream is(filename);
  if (!is)
    {
      std::cout << "ERROR: Failure opening sps file" << std::endl;
      std::exit(EXIT_SUCCESS);
    }

  // set up line counter for use in error messages
  std::string line;
  int line_count = 0;

  // variables for tracking read progress
  std::size_t num_mes;
  bool have_read_num_mes = false;
  std::size_t num_spes_read = 0;
  bool have_read_spes = false;
  std::size_t num_mes_read = 0;
  bool have_read_mes = false;

  while (!have_read_mes)
    {

      // get line
      mcutils::GetLine(is,line,line_count);
      std::istringstream line_stream(line);
      // std::cout << line_count << ": " << line << std::endl;
      
      // read header
      if (!have_read_spes)
        {
          if (!have_read_num_mes)
            {
              line_stream >> num_mes;
              mcutils::ParsingCheck(line_stream,line_count,line);
              have_read_num_mes = true;
            }
          while (!have_read_spes && !line_stream.eof())
            {
              double spe;
              line_stream >> spe;
              mcutils::ParsingCheck(line_stream,line_count,line);
              spe_data.push_back(spe);
              ++num_spes_read;
              have_read_spes = num_spes_read == num_orbitals;
            }
        }
      //read mes
      else
        {
          XPNTBMEDatum datum;
          line_stream >> datum.a >> datum.b >> datum.c >> datum.d
                      >> datum.J >> datum.T
                      >> datum.me;
          mcutils::ParsingCheck(line_stream,line_count,line);
          tbme_data.push_back(datum);
          ++num_mes_read;
          have_read_mes = num_mes_read == num_mes;
          // std::cout << "  read me: " << num_mes_read << " of " << num_mes << std::endl;
        }
    }

  // diagnostic data dump
  // for (std::size_t spe_index=0; spe_index<spe_data.size(); ++spe_index)
  //   {
  //     double spe = spe_data[spe_index];
  //     std::cout << fmt::format(
  //         "{:6d} {:e}",
  //         spe_index,
  //         spe
  //       ) << std::endl;
  //   }
  // for (std::size_t datum_index=0; datum_index<tbme_data.size(); ++datum_index)
  //   {
  //     const auto& datum = tbme_data[datum_index];
  //     std::cout << fmt::format(
  //         "{:6d} {:2d} {:2d} {:2d} {:2d} {:2d} {:2d} {:e}",
  //         datum_index,
  //         datum.a, datum.b, datum.c, datum.d, datum.J, datum.T, datum.me
  //       ) << std::endl;
  //   }
}

std::tuple<std::size_t,std::size_t,double>
LookUpStateFromXPNLabels(
    const basis::OrbitalSpacePN& orbital_space,
    const basis::TwoBodySpaceJJJPN& two_body_space,
    int a, int b, int J
  )
// Look up state within two body space matching given XPN labels.
//
// Arguments:
//   orbital_space (basis::OrbitalSpacePN): orbital space, for orbital qn lookups
//   two_body_space (basis::TwoBodySpaceJJJPN): two-body space, for state lookup
//   a, b (int): raw XPN indices for orbitals
//   J (int): two-particle state J, to calculate canonicalization factor
//
// Returns:
//   subspace_index (std::size_t): subspace index for two-body state
//   state_index (std::size_t): state index for two-body state
//   canonicalization_factor (double): canonicalization factor from canonicalizing orbitals
{

  std::size_t num_orbitals = orbital_space.dimension()/2;  // orbitals for single species
  
  // identify orbitals
  std::size_t orbital_state_index1 = (a-1) % num_orbitals;
  std::size_t orbital_subspace_index1 = (a-1) / num_orbitals;
  const basis::OrbitalStatePN& orbital_state1 = orbital_space.GetSubspace(orbital_subspace_index1).GetState(orbital_state_index1);
  std::size_t orbital_state_index2 = (b-1) % num_orbitals;
  std::size_t orbital_subspace_index2 = (b-1) / num_orbitals;
  const basis::OrbitalStatePN& orbital_state2 = orbital_space.GetSubspace(orbital_subspace_index2).GetState(orbital_state_index2);

  // canonicalize orbital indices
  if (orbital_subspace_index1 > orbital_subspace_index2)
    {
      std::cout << "ERROR: XPN file format specification does not support orbital indices in order np" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  double canonicalization_factor;
  std::tie(
      orbital_subspace_index1, orbital_subspace_index2,
      orbital_state_index1, orbital_state_index2,
      canonicalization_factor
    ) = CanonicalizeTwoBodyOrbitalIndicesJJJPN(
        orbital_space, J,
        orbital_subspace_index1, orbital_subspace_index2,
        orbital_state_index1, orbital_state_index2
      );

  // deduce subspace labels
  int g = (orbital_state1.g() + orbital_state2.g()) % 2;
  basis::TwoBodySpeciesPN two_body_species;
  if (orbital_subspace_index1 == 0 && orbital_subspace_index2 == 0)
    two_body_species = basis::TwoBodySpeciesPN::kPP;
  else if (orbital_subspace_index1 == 1 & orbital_subspace_index2 == 1)
    two_body_species = basis::TwoBodySpeciesPN::kNN;
  else if (orbital_subspace_index1 == 0 && orbital_subspace_index2 == 1)
    two_body_species = basis::TwoBodySpeciesPN::kPN;

  // look up indices
  std::size_t two_body_subspace_index =
    two_body_space.LookUpSubspaceIndex(basis::TwoBodySubspaceJJJPN::LabelsType(two_body_species, J, g));
  const basis::TwoBodySubspaceJJJPN& two_body_subspace =
    two_body_space.GetSubspace(two_body_subspace_index);
  std::size_t two_body_state_index =
    two_body_subspace.LookUpStateIndex(basis::TwoBodyStateJJJPN::LabelsType(orbital_state_index1, orbital_state_index2));

  // diagnostics
  // const basis::TwoBodyStateJJJPN& two_body_state =
  //   two_body_subspace.GetState(two_body_state_index);
  // std::cout
  //   << fmt::format("indices {} {} {} {} {}",
  //                  orbital_subspace_index1, orbital_subspace_index2,
  //                  orbital_state_index1, orbital_state_index2,
  //                  canonicalization_factor
  //     )
  //   << std::endl;
  // std::cout
  //   << fmt::format(
  //     "{} {} => subspace {} state {} {}",
  //     orbital_state1.LabelStr(), orbital_state2.LabelStr(),
  //     two_body_subspace.LabelStr(), two_body_state.LabelStr(),
  //     canonicalization_factor
  //   )
  //   << std::endl;

  return std::tuple<std::size_t,std::size_t,double>(
      two_body_subspace_index,
      two_body_state_index,
      canonicalization_factor
    );
}

void StoreTBMEs(
    const basis::OrbitalSpacePN& orbital_space,
    const std::vector<XPNTBMEDatum>& tbme_data,
    const basis::TwoBodySpaceJJJPN& two_body_space,
    const basis::TwoBodySectorsJJJPN& two_body_sectors,
    basis::OperatorBlocks<double>& two_body_matrices
  )
// Store raw XPN tbmes into standard JJJPN storage structures.
//
// Arguments:
//    orbital_space (basis::OrbitalSpacePN, input): orbitals
//    tbme_data (std::vector<XPNTBMEDatum>, input): raw data from input xpn file
//    two_body_space (basis::TwoBodySpaceJJJPN, input): space for storage of tbmes
//    two_body_sectors (basis::TwoBodySectorsJJJPN, input): sectors for storage of tbmes
//    two_body_matrices (basis::OperatorBlocks<double>, output): matrices for storage of tbmes
{

  for (std::size_t datum_index=0; datum_index<tbme_data.size(); ++datum_index)
    {
      // retrieve XPN TMBE datum
      const auto& datum = tbme_data[datum_index];

      // std::cout << fmt::format(
      //     "{:6d} {:2d} {:2d} {:2d} {:2d} {:2d} {:2d} {:e}",
      //     datum_index,
      //     datum.a, datum.b, datum.c, datum.d, datum.J, datum.T, datum.me
      //   ) << std::endl;

      // look up bra and ket states
      std::size_t subspace_index_bra, subspace_index_ket;
      std::size_t state_index_bra, state_index_ket;
      double canonicalization_factor_bra, canonicalization_factor_ket;
      std::tie(subspace_index_bra, state_index_bra, canonicalization_factor_bra)
        = LookUpStateFromXPNLabels(orbital_space, two_body_space, datum.a, datum.b, datum.J);
      std::tie(subspace_index_ket, state_index_ket, canonicalization_factor_ket)
        = LookUpStateFromXPNLabels(orbital_space, two_body_space, datum.c, datum.d, datum.J);
      // std::cout << fmt::format("before canonicalization: subspace indices {} {} state indices {} {}", subspace_index_bra, subspace_index_ket, state_index_bra, state_index_ket) << std::endl;

      // canonicalize matrix element labels
      int J0 = 0, g0 = 0;
      double canonicalization_factor;
      std::tie(
          subspace_index_bra, subspace_index_ket,
          state_index_bra, state_index_ket,
          canonicalization_factor
        ) =
        CanonicalizeIndicesJJJPN(
            two_body_space, J0, g0,
            subspace_index_bra, subspace_index_ket,
            state_index_bra, state_index_ket
          );
      // std::cout << fmt::format("after canonicalization: subspace indices {} {} state indices {} {}", subspace_index_bra, subspace_index_ket, state_index_bra, state_index_ket) << std::endl;
      
      // store ME
      double value = canonicalization_factor * canonicalization_factor_bra * canonicalization_factor_ket * datum.me;
      std::size_t sector_index = two_body_sectors.LookUpSectorIndex(subspace_index_bra, subspace_index_ket);
      if (sector_index == basis::kNone)
        {
          std::cout << fmt::format("ERROR: matrix element in nonexistent sector (subspace indices {} {} => sector index {}", subspace_index_bra, subspace_index_ket, sector_index) << std::endl;
        }
      // std::cout
      //   << fmt::format(
      //       "subspace indices {} {} sector index {} state indices {} {} value {}",
      //       subspace_index_bra, subspace_index_ket, sector_index, state_index_bra, state_index_ket, value
      //     )
      //   << std::endl;      
      two_body_matrices[sector_index](state_index_bra, state_index_ket) = value;
      
    }

  
}
int main(int argc, const char *argv[])
{
  // header
  std::cout << std::endl;
  std::cout << "xpn2h2  -- xpn to h2 conversion" << std::endl;
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  // read parameters
  RunParameters run_parameters;
  ProcessArguments(argc,argv,run_parameters);

  // read orbitals
  std::cout << fmt::format("Reading orbital file {}...", run_parameters.orbital_filename) << std::endl;
  std::ifstream orbital_stream(run_parameters.orbital_filename);
  basis::OrbitalPNList orbitals =
    basis::ParseOrbitalPNStream(orbital_stream, true);
  basis::OrbitalSpacePN orbital_space(orbitals);

  // set up operator storage
  // Note: Both xpn and h2 use NAS storage, so use NAS internally as well.
  int J0 = 0;
  int g0 = 0;
  int Tz0 = 0;
  const basis::TwoBodySpaceJJJPN two_body_space(
      orbital_space,
      run_parameters.weight_max,
      shell::kH2SpaceOrdering.at(run_parameters.output_format)
    );
  const basis::TwoBodySectorsJJJPN two_body_sectors(
      two_body_space,
      J0, g0, Tz0
    );
  basis::OperatorBlocks<double> two_body_matrices;
  basis::SetOperatorToZero(two_body_sectors, two_body_matrices);

  // read xpn file
  // Note: If memory limitations become severe, can combine the read and storage steps.
  std::cout << fmt::format("Reading XPN file {}...", run_parameters.input_filename) << std::endl;
  std::vector<double> spe_data;
  std::vector<XPNTBMEDatum> tbme_data;
  ReadXPNFile(run_parameters.input_filename, orbital_space.dimension(), spe_data, tbme_data);

  // accumulate XPN matrix elements into operator matrices
  StoreTBMEs(
    orbital_space,
    tbme_data,
    two_body_space,
    two_body_sectors,
    two_body_matrices
    );

  // process OBMEs
  if (run_parameters.obme_filename != "")
    {
      std::cout << fmt::format("Writing OBME file {}...", run_parameters.obme_filename) << std::endl
                << "Implementation in progress!" << std::endl
                << std::endl;
      std::exit(EXIT_FAILURE);
      //WIP
    }
  
  // write h2 file
  std::cout << fmt::format("Writing h2 file {}...", run_parameters.output_filename) << std::endl
            << std::endl;
  shell::OutH2Stream output_stream(
      run_parameters.output_filename,
      orbital_space, two_body_space, two_body_sectors,
      run_parameters.output_format
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
