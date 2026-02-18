/***************************************

  number-op-gen.cpp -- create obme file representing orbital number operator

  Syntax:

    + number-op-gen orbital_filename n l 2*j particle_species N_max output_filename

        orbital_filename (str): filename for orbitals.dat file used to make h2 file

        n, l, 2*j (int): quantum numbers for the orbital of desired pair-counting operator

        particle_species (int or str): particle species associated with desired pair-counting operator; 1 or 'p' if proton, 2, -1, or 'n' if neutron

        N_max (int): N_max beyond which to truncate (one-body truncation)
        
        output_filename (str): filename for output file

  Assumed file format for orbital file (orbitals.dat):

    + comment lines beginning with hash ('#') or bang ('!')

    + header (possibly wrapped over multiple lines):

        format norb_p norb_n

        format (str): h2 format of file
        
        norb_p, norb_n (int): number of proton (neutron) orbitals contained in file

    + lines of form

        index n l 2*j species weight

        index (int): 1-based orbital index

        n, l, 2*j (int): quantum numbers of orbital associated with that index

        species (int): 1 if proton, 2 if neutron if version 15099; 2T_z (HEP convention) if version 15200

        weight (float): N = 2 * n + l for harmonic oscillator orbitals
      
  Kayla E. O'Donnell
  Massachusetts Institute of Technology

***************************************/
#include <cmath>    // sqrt function

#include "basis/nlj_operator.h"
#include "mcutils/parsing.h"
#include "mcutils/profiling.h"
#include "obme/obme_io.h"

struct RunParameters
// store basic parameters for run
{
  // command-line inputs
  std::string orbital_filename;
  int n, l;
  HalfInt j;
  basis::OrbitalSpeciesPN particle_species;
  std::string output_filename;
  basis::WeightMax weight_max;
  
  // default constructor
  RunParameters()
  {
    orbital_filename = "";
    n = 0;
    l = 0;
    j = HalfInt(1, 2);
    particle_species = basis::OrbitalSpeciesPN::kP;
    output_filename = "";
    weight_max = basis::WeightMax(0, 0);
  }
};

void ProcessArguments(int argc, char **argv, RunParameters& run_parameters)
// process command-line arguments
{
  // usage message
  if (argc - 1 != 7)
  {
    std::cout << "7 arguments expected, " << argc - 1 << " arguments found." << std::endl;
    std::cout << "Usage: number-op-gen orbital_filename n l 2*j particle_species N_max output_filename" << std::endl;
    std::exit(EXIT_SUCCESS);
  }
  
  // input 1 (orbitals.dat filename)
  run_parameters.orbital_filename = argv[1];
  mcutils::FileExistCheck(run_parameters.orbital_filename, true, false);
  
  // inputs 2-4 (n, l, 2*j)
  std::istringstream parameter_stream_2(argv[2]);
  parameter_stream_2 >> run_parameters.n;
  if (!parameter_stream_2)
  {
    std::cerr << "ERROR: Invalid n \"" << argv[2] << "\"; expected integer" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::istringstream parameter_stream_3(argv[3]);
  parameter_stream_3 >> run_parameters.l;
  if (!parameter_stream_3)
  {
    std::cerr << "ERROR: Invalid l \"" << argv[3] << "\"; expected integer" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  int twice_j;
  std::istringstream parameter_stream_4(argv[4]);
  parameter_stream_4 >> twice_j;
  if (!parameter_stream_4)
  {
    std::cerr << "ERROR: Invalid 2*j \"" << argv[4] << "\"; expected integer" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  run_parameters.j = HalfInt(twice_j, 2);
  
  // input 5 (particle species of pair-counting operator)
  std::string s5 = std::string(argv[5]);
  if (s5 == "1" || s5 == "p")
  {
    run_parameters.particle_species = basis::OrbitalSpeciesPN::kP;
  }
  else if (s5 == "2" || s5 == "-1" || s5 == "n")
  {
    run_parameters.particle_species = basis::OrbitalSpeciesPN::kN;
  }
  else
  {
    std::cerr << "ERROR: Invalid particle type \"" << argv[5] << "\"" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  
  // input 6 (truncation cutoff)
  int truncation_cutoff;
  std::istringstream parameter_stream_6(argv[6]);
  parameter_stream_6 >> truncation_cutoff;
  if (!parameter_stream_6)
  {
    std::cerr << "ERROR: Invalid truncation cutoff \"" << argv[6] << "\"; expected integer" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  run_parameters.weight_max = basis::WeightMax(basis::Rank::kOneBody, truncation_cutoff);
  
  // input 7 (output filename)
  run_parameters.output_filename = argv[7];
}

int main(int argc, char** argv)
{
  // header
  std::cout << std::endl;
  std::cout << "number-op-gen -- generating obme files for orbital number operators" << std::endl;
  std::cout << "version: " << VCS_REVISION << std::endl;
  std::cout << std::endl;

  // read parameters
  RunParameters run_parameters;
  ProcessArguments(argc, argv, run_parameters);
  
  // start timing
  mcutils::SteadyTimer total_time;
  total_time.Start();
  
  // define quantum numbers for N_k
  const int J0 = 0;
  const int g0 = 0;
  const int Tz0 = 0;
  
  // read and process orbitals.dat file
  std::cout << "Processing " << run_parameters.orbital_filename << " ... ";
  std::ifstream is(run_parameters.orbital_filename);
  basis::OrbitalPNList orbitals_list_raw = basis::ParseOrbitalPNStream(is, true);
  is.close();
  if (run_parameters.weight_max.one_body[0] > basis::OrbitalSpacePN(orbitals_list_raw).GetSubspace(0).weight_max())
  {
    std::cerr << "\nERROR: Desired N_max is unsupported by provided orbitals.dat file" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  basis::OrbitalPNList orbitals_list = basis::TruncateOrbitalList(run_parameters.weight_max, orbitals_list_raw);
  std::cout << "Done" << std::endl;
  
  // set up one-body space
  std::cout << "Setting up one-body space ... ";
  basis::OrbitalSpaceLJPN one_body_space(orbitals_list);
  basis::OrbitalSubspaceLJPNLabels one_body_subspace_labels(run_parameters.particle_species, run_parameters.l, run_parameters.j);
  basis::OrbitalSubspaceLJPN one_body_subspace = one_body_space.LookUpSubspace(one_body_subspace_labels);
  basis::OrbitalStateLJPNLabels orbital_state_labels(run_parameters.n);
  if (!(one_body_subspace.ContainsState(orbital_state_labels)))
  {
    std::cerr << "\nERROR: Orbital index entered is either unsupported by orbitals.dat file entered or truncated away by N_max entered" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  basis::OrbitalSectorsLJPN one_body_sectors(one_body_space, one_body_space, J0, g0, Tz0);
  basis::OperatorBlocks<double> one_body_matrices;
  basis::SetOperatorToZero(one_body_sectors, one_body_matrices);
  std::cout << "Done" << std::endl;
  
  // locate and modify appropriate matrix element
  std::cout << "Setting matrix element for " << (run_parameters.particle_species == basis::OrbitalSpeciesPN::kP ? "proton" : "neutron") 
    << " orbital with n = " << run_parameters.n << ", l = " << run_parameters.l << ", j = " << run_parameters.j << " ... ";
  std::size_t one_body_subspace_index = one_body_space.LookUpSubspaceIndex(one_body_subspace_labels);
  std::size_t one_body_state_index = one_body_subspace.LookUpStateIndex(orbital_state_labels);
  if (one_body_state_index == basis::kNone)
  {
    std::cout << "\nERROR: Orbital index entered has been truncated away by N_max entered" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::size_t sector_index = one_body_sectors.LookUpSectorIndex(one_body_subspace_index, one_body_subspace_index);
  one_body_matrices[sector_index](one_body_state_index, one_body_state_index) = 1;
  std::cout << "Done" << std::endl;
  
  // write obme file
  std::cout << "Writing to " << run_parameters.output_filename << " ... ";
  shell::OutOBMEStream output_stream(run_parameters.output_filename, one_body_space, one_body_space, one_body_sectors, basis::OneBodyOperatorType::kSpherical);
  output_stream.Write(one_body_matrices);
  std::cout << "Done" << std::endl;
  
  // end timing
  total_time.Stop();
  std::cout << "(Total time: " << total_time.ElapsedTime() << ")" << std::endl;
  std::cout << std::endl;

  // exit
  return EXIT_SUCCESS;
}
