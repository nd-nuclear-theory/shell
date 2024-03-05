/******************************************************************************

  xpn2h2.cpp -- convert BIGSTICK xpn to h2

  Syntax:
    + xpn2h2 orbital_filename input_filename output_filename

    programs/h2utils/xpn2h2 ncci-tb-6.dat JISP16-tb-6-20.int JISP16-tb-6-20.dat

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
    + Output h2 format is fixed as Version 15099.

  See Sec 4.3.2 "Proton-neutron and other isospin-breaking formats" in
  C. W. Johnson et al., "BIGSTICK: A flexible configuration-interaction
  shell-model code", arXiv:1801.08432.

  Mark A. Caprio
  University of Notre Dame

  + 03/01/23 (mac): Created, based on me2j2h2.

******************************************************************************/

#include <fstream>

#include "basis/operator.h"
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
  // format
  shell::H2Format output_format;
  
  // default constructor
  RunParameters()
    : orbital_filename(""), output_filename(""), input_filename(""),
      output_format(shell::kVersion15099)
  {}
};

void PrintUsage(const char **argv) {
  std::cout << "Usage: " << argv[0]
            << " orbital_file input_file output_file"
            << std::endl;
}

void ProcessArguments(int argc, const char *argv[], RunParameters& run_parameters)
{
  // usage message
  if (argc-1 < 3)
    {
      PrintUsage(argv);
      std::exit(EXIT_SUCCESS);
    }

  // orbital filename
  run_parameters.orbital_filename = argv[1];
  mcutils::FileExistCheck(run_parameters.orbital_filename, true, false);

  // input filename
  run_parameters.input_filename = argv[2];
  mcutils::FileExistCheck(run_parameters.input_filename, true, false);

  // output filename
  run_parameters.output_filename = argv[3];
  
}

struct XPNTBMEDatum
// Stores raw data from XPN file me data line
{
  int a, b, c, d, J, T;
  double me;
};

std::vector<XPNTBMEDatum> ReadXPNFile(const std::string& filename, std::size_t num_orbitals)
// Read raw BIGSTICK xpn file tbme data.
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

  // data
  std::size_t num_mes;
  std::vector<XPNTBMEDatum> tbme_data;
  
  // progress flags
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

  return tbme_data;
}

//basis::TwoBodyStateJJJPN
void LookUpStateFromXPNLabels(
    const basis::OrbitalSpacePN& orbital_space,
    const basis::TwoBodySpaceJJJPN& two_body_space,
    int a, int b, int J
  )
// Look up state within two body space matching given XPN labels.
//
// CAVEAT: Need to implement canonicalization.
{

  std::size_t num_orbitals = orbital_space.dimension()/2;  // orbitals for single species
  
  // convert state indices
  std::size_t index1 = (a-1) % num_orbitals;
  const basis::OrbitalSpeciesPN& species1 = basis::OrbitalSpeciesPN((a-1) / num_orbitals);
  const basis::OrbitalStatePN& orbital_state1 = orbital_space.LookUpSubspace(basis::OrbitalSubspacePNLabels(species1)).GetState(index1);
  std::size_t index2 = (b-1) % num_orbitals;
  const basis::OrbitalSpeciesPN& species2 = basis::OrbitalSpeciesPN((b-1) / num_orbitals);
  const basis::OrbitalStatePN& orbital_state2 = orbital_space.LookUpSubspace(basis::OrbitalSubspacePNLabels(species2)).GetState(index2);

  // look up subspace
  int g = (orbital_state1.g() + orbital_state2.g()) % 2;
  std::cout << fmt::format("{} {} {}", orbital_state1.LabelStr(), orbital_state2.LabelStr(), g)
            << std::endl;

  //return basis::TwoBodyStateJJJPN();
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
      const auto& datum = tbme_data[datum_index];

      // find bra state
      
      std::cout << fmt::format(
          "{:6d} {:2d} {:2d} {:2d} {:2d} {:2d} {:2d} {:e}",
          datum_index,
          datum.a, datum.b, datum.c, datum.d, datum.J, datum.T, datum.me
        ) << std::endl;
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
  std::ifstream orbital_stream(run_parameters.orbital_filename);
  basis::OrbitalPNList orbitals =
    basis::ParseOrbitalPNStream(orbital_stream, true);
  basis::OrbitalSpacePN orbital_space(orbitals);

  // read xpn file
  std::vector<XPNTBMEDatum> tbme_data = ReadXPNFile(run_parameters.input_filename, orbital_space.dimension());
  for (std::size_t datum_index=0; datum_index<tbme_data.size(); ++datum_index)
    {
      const auto& datum = tbme_data[datum_index];
      std::cout << fmt::format(
          "{:6d} {:2d} {:2d} {:2d} {:2d} {:2d} {:2d} {:e}",
          datum_index,
          datum.a, datum.b, datum.c, datum.d, datum.J, datum.T, datum.me
        ) << std::endl;
    }
  
  // set truncation
  // basis::WeightMax weight_max(2,2,2,2,2);  // WIP
  basis::WeightMax weight_max(6, 6);  // WIP

  // set up operator storage
  // Note: Both xpn and h2 use NAS storage, so use NAS internally as well.
  int J0 = 0;
  int g0 = 0;
  int Tz0 = 0;
  const basis::TwoBodySpaceJJJPN two_body_space(
      orbital_space,
      weight_max,
      shell::kH2SpaceOrdering.at(run_parameters.output_format)
    );
  const basis::TwoBodySectorsJJJPN two_body_sectors(
      two_body_space,
      J0, g0, Tz0
    );
  basis::OperatorBlocks<double> two_body_matrices;
  basis::SetOperatorToZero(two_body_sectors, two_body_matrices);

  LookUpStateFromXPNLabels(orbital_space, two_body_space, 1, 1, 0);
  std::exit(EXIT_SUCCESS);
  
  // write h2 file
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
