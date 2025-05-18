/******************************************************************************

  obme2dens.cpp -- convert obme to density tabulation

  Densities follow conventions of equation (5.1) of C. W. Johnson, "BIGSTICK: A
  flexible configuration-interaction shell-model code", arxiv:1801.08432.  These
  differ from the MFDn ROBDMEs by an angular momentum factor 1/sqrt(2*J0+1),
  where J0 is the multipolarity.  The time-reversal ("tilde") convention on the
  annilation operator in this convention is such that diagonal scalar densities,
  proportional to orbital occupations, are positive.

  An MFDn obdme file should first be converted to shell obme format (with
  obdme-conv), to then be input to obme2dens.

    obme2dens input_filename output_filename

  Mark A. Caprio
  University of Notre Dame

  + 05/14/25 (mac): Created.

******************************************************************************/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>

#include "am/am.h"
#include "basis/nlj_orbital.h"
#include "basis/operator.h"
#include "fmt/format.h"
#include "mcutils/parsing.h"
#include "mcutils/profiling.h"
#include "obme/obme.h"
#include "obme/obme_io.h"
#include "obme/radial.h"

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

struct RunParameters
// Stores simple parameters for run
{
  // filenames
  std::string input_filename;
  std::string output_filename;

  // default constructor
  RunParameters()
    : input_filename(""), output_filename("")
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
  mcutils::FileExistCheck(run_parameters.output_filename, false, true);

}


////////////////////////////////////////////////////////////////
// h2 input
/////////////////////////////////////////////////////////////////

void ReadOBMEFile(
    const std::string& filename,
    basis::OrbitalSpaceLJPN& orbital_space,
    basis::OrbitalSectorsLJPN& sectors, basis::OperatorBlocks<double>& matrices
  )
// Read all data from h2 file.
//
// Arguments:
//   filename (std::string): filename
//   orbital_space (basis::OrbitalSpaceLJPN, output): orbitals
//   sectors (basis::OrbitalSectorsLJPN, output): sectors
//   matrices (basis::OperatorBlocks<double>, output): OBME matrices
{

  // initialize stream
  std::cout << "Reading OBME file..." << std::endl;
  shell::InOBMEStream operator_stream(filename);
  basis::OrbitalSpaceLJPN ket_orbital_space;
  operator_stream.SetToIndexing(
      orbital_space, ket_orbital_space, sectors
    );
  assert(orbital_space.OrbitalInfo() == ket_orbital_space.OrbitalInfo());
  std::cout << fmt::format("  Filename: {}", filename) << std::endl;
  int J0 = sectors.J0();
  int g0 = sectors.g0();
  int Tz0 = sectors.Tz0();
  std::cout << fmt::format("  Operator properties: J0 {} g0 {} Tz0 {}", J0, g0, Tz0) << std::endl;
  
  // read matrices
  operator_stream.Read(matrices);

  // diagnostic
  std::cout
    << fmt::format(
        "  Diagnostics: sectors {} matrix elements {}",
        sectors.size(),
        basis::AllocatedEntries(matrices)
      )
    << std::endl;
  std::cout << std::endl;
  
  // close file
  operator_stream.Close();

}


////////////////////////////////////////////////////////////////
// densities table output
/////////////////////////////////////////////////////////////////

void WriteDensitiesFile(
    const std::string& filename,
    const basis::OrbitalSpaceLJPN& orbital_space,
    const basis::OrbitalSectorsLJPN& sectors, const basis::OperatorBlocks<double>& matrices
  )
// Write all data to densities table file.
//
// Arguments:
//   filename (std::string): filename
//   orbital_space (basis::OrbitalSpaceLJPN, input): orbitals
//   sectors (basis::OrbitalSectorsLJPN, input): sectors
//   matrices (basis::OperatorBlocks<double>, input): OBME matrices
{

  // open densities file
  std::ofstream os(filename);
  if (!os)
    {
      std::cout << "ERROR: Failure opening densities file" << std::endl;
      std::exit(EXIT_SUCCESS);
    }

  // write header comment
  os << "# density tabulation written by obme2dens" << std::endl
     << "#" << std::endl
     << "# Densities follow conventions of equation (5.1) of C. W. Johnson, \"BIGSTICK: A" << std::endl
     << "# flexible configuration-interaction shell-model code\", arxiv:1801.08432." << std::endl
     << "#" << std::endl
     << "# Isospin convention: Tz=+1/2 (proton); Tz=-1/2 (neutron)" << std::endl
     << "#" << std::endl
     << "# Entries are of the form:" << std::endl
     << "#  na   la 2*ja 2*Tza    nb   lb 2*jb 2*Tzb   J0   g0  Tz0              rho" << std::endl;

  for (std::size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
    {

      // extract sector
      const auto& sector = sectors.GetSector(sector_index);
      const auto& bra_subspace = sector.bra_subspace();
      const auto& ket_subspace = sector.ket_subspace();

      // iterate over matrix elements
      for (std::size_t bra_index=0; bra_index<bra_subspace.size(); ++bra_index)
        for (std::size_t ket_index=0; ket_index<ket_subspace.size(); ++ket_index)
          {

            // retrieve states
            const auto bra = bra_subspace.GetState(bra_index);
            const auto ket = ket_subspace.GetState(ket_index);

            // extract matrix element
            const double matrix_element = matrices[sector_index](bra_index, ket_index);

            // generate output line
            os << fmt::format(
                " {:4d} {:4d} {:4d} {:+4d}   {:4d} {:4d} {:4d} {:+4d}   {:4d} {:4d} {:+4d}   {:+13.8f}",
                bra.n(), bra.l(), bra.j().TwiceValue(), bra.Tz().TwiceValue(),
                ket.n(), ket.l(), ket.j().TwiceValue(), ket.Tz().TwiceValue(),
                sectors.J0(), sectors.g0(), sectors.Tz0(),
                matrix_element
              )
               << std::endl;

          }

    }

  // close output file
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
  std::cout << "obme2dens -- convert obme to density tabulation" << std::endl;
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
  
  // read obme
  basis::OrbitalSpaceLJPN orbital_space;
  basis::OrbitalSectorsLJPN sectors;
  basis::OperatorBlocks<double> matrices;
  ReadOBMEFile(
      run_parameters.input_filename,
      orbital_space, sectors, matrices
    );

  
  // write output
  std::cout << "Output stream" << std::endl;
  std::cout << fmt::format("  File: {}", run_parameters.output_filename) << std::endl;
  WriteDensitiesFile(
      run_parameters.output_filename,
      orbital_space, sectors, matrices
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
