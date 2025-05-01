/******************************************************************************

  h22xpn.cpp -- convert H2 to BIGSTICK XPN

  Syntax:

    h22xpn input_filename output_filename

  Limitations:

  This code currently imposes the (overly strict) requirement that the
  single-particle orbitals be in an oscillator-like truncation.  For xpn format
  (and Bigstick), it suffices to require that the proton and neutron orbital
  sets are identical.  The check should be relaxed to support, e.g., traditional
  shell model TBME files.  But then we also would need to worry about specifying
  nonzero single-particle energies.

  This code dumps the orbital list into the xpn file header comment, in the
  format expected for a Bigstick sps file (in "iso" mode).  But it does not
  actually write on an sps file.  This functionality could be added.  Or a
  separate workflow h2->orbitals->sps could be provided.

  Mark A. Caprio
  University of Notre Dame

  + 04/30/25 (mac): Created, drawing on xpn2h2 and h22me2j.

******************************************************************************/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>

#include "basis/jjjpn_scheme.h"
#include "basis/jjjpn_operator.h"
#include "fmt/format.h"
#include "mcutils/parsing.h"
#include "mcutils/profiling.h"

#include "tbme/h2_io.h"

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

struct RunParameters
// Stores simple parameters for run
{
  // filenames
  std::string input_filename;
  std::string output_filename;

  // file format
  std::size_t float_size;

  // default constructor
  RunParameters()
    : input_filename(""), output_filename(""), float_size(4)
  {}
  
};

void PrintUsage(const char **argv) {
  std::cout << "Usage: " << argv[0]
            << " input_filename output_filename"
            << std::endl;
}

void ProcessArguments(int argc, const char *argv[], RunParameters& run_parameters)
{
  
  int arg = 1;

  // process options
  while (arg < argc && argv[arg][0] == '-')
    {
      std::istringstream parameter_stream(argv[arg++]);

      if (parameter_stream.str() == "--help" || parameter_stream.str() == "-h")
        {
          PrintUsage(argv);
          std::exit(EXIT_SUCCESS);
        }
      else
        {
          PrintUsage(argv);
          std::cerr << "Unrecognized option '" << parameter_stream.str() << "'" << std::endl;
          std::exit(EXIT_FAILURE);
        }
    }
  
  // process fixed arguments
  if (argc-arg < 2)
    {
      PrintUsage(argv);
      std::cerr << "Insufficient arguments" << std::endl;
      std::exit(EXIT_FAILURE);
    }

  // input filename
  run_parameters.input_filename = argv[arg++];
  mcutils::FileExistCheck(run_parameters.input_filename, true, false);

  // output filename
  run_parameters.output_filename = argv[arg++];
  
  // orbital filename
  // run_parameters.orbital_filename = argv[arg++];
  // mcutils::FileExistCheck(run_parameters.orbital_filename, true, false);
}


////////////////////////////////////////////////////////////////
// h2 input
/////////////////////////////////////////////////////////////////

void ReadH2File(
    const std::string& filename,
    basis::OrbitalSpacePN& orbital_space, basis::TwoBodySpaceJJJPN& two_body_space,
    basis::TwoBodySectorsJJJPN& two_body_sectors, basis::OperatorBlocks<double>& two_body_matrices,
    basis::NormalizationConversion conversion_mode
  )
// Read all data from h2 file.
//
// Arguments:
//   filename (std::string): filename
//   orbital_space (basis::OrbitalSpacePN, output): orbitals
//   two_body_space (basis::TwoBodySpaceJJJPN, output): two-body space
//   two_body_sectors (basis::TwoBodySectorsJJJPN, output): two-body sectors
//   two_body_matrices (basis::OperatorBlocks<double>, output): TBME matrices
//   conversion_mode (basis::NormalizationConversion, optional): selects AS/NAS conversion mode
//
// Returns:
//   num_tbmes (std::size_t): Number of TBMEs in h2 file
{

  // initialize stream

  std::cout << "Input stream" << std::endl;
  shell::InH2Stream input_stream(filename);
  std::cout << input_stream.DiagnosticStr();
  std::cout << std::endl;

  // initialize data structures
  orbital_space = input_stream.orbital_space();
  two_body_space = input_stream.space();
  two_body_sectors = input_stream.sectors();
  const basis::TwoBodySpaceJJJPNOrdering space_ordering =
    shell::kH2SpaceOrdering.at(input_stream.h2_format());
  two_body_space = basis::TwoBodySpaceJJJPN(
      orbital_space,
      two_body_space.weight_max(),
      space_ordering
    );
  two_body_sectors = basis::TwoBodySectorsJJJPN(
      two_body_space, two_body_sectors.J0(), two_body_sectors.g0(), two_body_sectors.Tz0()
    );
  two_body_matrices.resize(two_body_sectors.size());
  
  // read sectors
  for (std::size_t sector_index = 0; sector_index < input_stream.num_sectors(); ++sector_index)
    {
      auto& matrix = two_body_matrices[sector_index];
      input_stream.ReadSector(sector_index, matrix, conversion_mode);
    }

  // close stream
  input_stream.Close();

}


////////////////////////////////////////////////////////////////
// xpn output
/////////////////////////////////////////////////////////////////

std::tuple<int,int,int,int>
GenerateXPNLabelsFromState(
    const basis::TwoBodyStateJJJPN& two_body_state,
    std::size_t num_orbitals
  )
// Generate XPN labels given state within two body space.
//
// Arguments:
//   two_body_state (basis::TwoBodyStateJJJPN): two-body state
//   num_orbitals (std::size_t): number of orbitals for single species
//
// Returns:
//   a, b (int): raw XPN indices for orbitals
//   Tz, J (int): two-particle state Tz and J
{

  // extract subspace labels
  basis::TwoBodySpeciesPN two_body_species = two_body_state.two_body_species();
  int Tz = two_body_state.Tz();
  int J = two_body_state.J();

  // convert orbitals
  //
  // Convert to 1-based with offset for neutrons.
  int orbital_offset1 = 1;
  if (two_body_species == basis::TwoBodySpeciesPN::kNN)
    orbital_offset1 += num_orbitals;
  int orbital_offset2 = 1;
  if (two_body_species == basis::TwoBodySpeciesPN::kPN || (two_body_species == basis::TwoBodySpeciesPN::kNN))
    orbital_offset2 += num_orbitals;
  int a = orbital_offset1 + two_body_state.index1();
  int b = orbital_offset2 + two_body_state.index2();

  return std::tuple<int, int, int, int>(a, b, Tz, J);
}

void WriteXPNFile(
    const std::string& filename,
    const basis::OrbitalSpacePN& orbital_space, const basis::TwoBodySpaceJJJPN& two_body_space,
    const basis::TwoBodySectorsJJJPN& two_body_sectors, const basis::OperatorBlocks<double>& two_body_matrices
  )
// Write all data to xpn file.
//
// Arguments:
//   filename (std::string): filename
//   orbital_space (basis::OrbitalSpacePN, output): orbitals
//   two_body_space (basis::TwoBodySpaceJJJPN, output): two-body space
//   two_body_sectors (basis::TwoBodySectorsJJJPN, output): two-body sectors
//   two_body_matrices (basis::OperatorBlocks<double>, output): TBME matrices
//   num_tbmes (std::size_t): number of TBMEs 
{

  std::size_t num_orbitals = orbital_space.dimension()/2;  // orbitals for single species
  
  // open xpn file
  std::ofstream os(filename);
  if (!os)
    {
      std::cout << "ERROR: Failure opening xpn file" << std::endl;
      std::exit(EXIT_SUCCESS);
    }
  // shell::InH2Stream input_stream(filename);
  // std::cout << input_stream.DiagnosticStr();
  // std::cout << std::endl;

  // write header comment
  os << "# XPN file written by h22xpn" << std::endl;
  os << "#" << std::endl;
  os << "# Entries are of the form:" << std::endl;
  os << "#   a b c d J T ME" << std::endl;
  os << "# where T is a dummy field (set here to |Tz|)." << std::endl;
  os << "#" << std::endl;
  os << "# Orbitals:" << std::endl;
  os << fmt::format("#   {}", "iso") << std::endl;
  os << fmt::format("#   {}", num_orbitals) << std::endl;
  basis::OrbitalPNList proton_orbitals = orbital_space.GetSubspace(0).OrbitalInfo();
  for (const basis::OrbitalPNInfo& orbital : proton_orbitals)
    {
      os << fmt::format("#   {} {} {} {}", orbital.n, orbital.l, float(orbital.j), orbital.weight) << std::endl;
    }

  // write header
  // write number of TBMEs
  std::size_t num_tbmes = basis::UpperTriangularEntries(two_body_sectors);
  os << fmt::format("{}", num_tbmes) << std::endl;
  // write dummy SPEs
  for (std::size_t count=0; count < 2*num_orbitals; ++count)
    {
      os << "0.0" << std::endl;
    }

  // iterate over sectors
  //
  // Modeled on control loop for basis::WriteTwoBodyOperatorJJJPN().
  
  for (std::size_t sector_index = 0; sector_index < two_body_sectors.size(); ++sector_index)
    {

      // extract sector
      const typename basis::TwoBodySectorsJJJPN::SectorType& sector = two_body_sectors.GetSector(sector_index);
      const typename basis::TwoBodySectorsJJJPN::SubspaceType& bra_subspace = sector.bra_subspace();
      const typename basis::TwoBodySectorsJJJPN::SubspaceType& ket_subspace = sector.ket_subspace();

      // verify that sector is canonical
      //
      // This is a check that the caller's sector construction
      // followed the specification that only "upper triangle"
      // sectors are stored.
      assert(sector.IsUpperTriangle());

      // iterate over matrix elements
      for (std::size_t bra_index=0; bra_index<bra_subspace.size(); ++bra_index)
        for (std::size_t ket_index=0; ket_index<ket_subspace.size(); ++ket_index)
          {

            // diagonal sector: restrict to upper triangle
            if (sector.IsDiagonal())
              if (!(bra_index<=ket_index))
                continue;

            // retrieve states
            const basis::TwoBodyStateJJJPN bra(bra_subspace,bra_index);
            const basis::TwoBodyStateJJJPN ket(ket_subspace,ket_index);

            // extract state labels
            int a, b, Tz_bra, J_bra;
            std::tie(a, b, Tz_bra, J_bra) = GenerateXPNLabelsFromState(bra, num_orbitals);
            int c, d, Tz_ket, J_ket;
            std::tie(c, d, Tz_ket, J_ket) = GenerateXPNLabelsFromState(ket, num_orbitals);
            assert((Tz_bra==Tz_ket) && (J_bra==J_ket));
            
            // extract matrix element
            const double matrix_element = two_body_matrices[sector_index](bra_index,ket_index);

            // generate output line
            int T = abs(Tz_bra);  // dummy pseudo-T value
            int J = J_bra;
            os << fmt::format(" {:4d} {:4d} {:4d} {:4d} {:4d} {:4d} {:13.8f}", a, b, c, d, J, T, matrix_element)
               << std::endl;

          }

    }

  // close stream
  os.close();

}


////////////////////////////////////////////////////////////////
// main program
/////////////////////////////////////////////////////////////////

int main(int argc, const char **argv)
{

  ////////////////////////////////////////////////////////////////
  // initialization
  ////////////////////////////////////////////////////////////////

  // header
  std::cout << std::endl;
  std::cout << "h22xpn -- convert H2 to BIGSTICK XPN" << std::endl;
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  // read parameters
  RunParameters run_parameters;
  ProcessArguments(argc, argv, run_parameters);

  // start timing
  mcutils::SteadyTimer total_time;
  total_time.Start();

  ////////////////////////////////////////////////////////////////
  // conversion
  ////////////////////////////////////////////////////////////////
  
  // read h2

  basis::OrbitalSpacePN orbital_space;
  basis::TwoBodySpaceJJJPN two_body_jjjpn_space;
  basis::TwoBodySectorsJJJPN two_body_jjjpn_sectors;
  basis::OperatorBlocks<double> two_body_jjjpn_matrices;
  ReadH2File(
      run_parameters.input_filename,
      orbital_space, two_body_jjjpn_space,
      two_body_jjjpn_sectors, two_body_jjjpn_matrices,
      basis::NormalizationConversion::kNone
    );

  // extract source operator information
  
  // validate and extract orbital truncation
  // TODO: relax constraint to simple mirror-symmetric orbital set
  if (!orbital_space.is_oscillator_like())
    {
      std::cerr << "ERROR: Input h2 file defines orbital set which is not oscillator-like" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  int Nmax_orb = int(orbital_space.weight_max());
  
  // extract and validate operator quantum numbers
  int J0 = two_body_jjjpn_sectors.J0();
  int g0 = two_body_jjjpn_sectors.g0();
  int Tz0 = two_body_jjjpn_sectors.Tz0();
  if (J0!=0 || g0!=0 || Tz0!=0)
    {
      std::cerr << "ERROR: Input h2 file must specify a Hamiltonian-like (J0,g0,Tz0)=(0,0,0) operator" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  
  // write output
  std::cout << "Output stream" << std::endl;
  std::cout << fmt::format("  File: {}", run_parameters.output_filename) << std::endl;
  WriteXPNFile(
      run_parameters.output_filename,
      orbital_space, two_body_jjjpn_space,
      two_body_jjjpn_sectors, two_body_jjjpn_matrices
    );
  std::cout << std::endl;
  
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
