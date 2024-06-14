/***************************************

  h2pair.cpp -- create h2 file representing NAS pair-counting operator

  Syntax:

    + h2pair orbital_filename n l 2*j particle_species truncation_rank N_max output_filename output_format

				orbital_filename (str): filename for orbitals.dat file used to make h2 file

				n, l, 2*j (int): quantum numbers for the orbital of desired pair-counting operator

				particle_species (int): particle species associated with desired pair-counting operator; 1 or 'p' if proton, 2, -1, or 'n' if neutron

				truncation_rank (int): 1 if one-body truncation, 2 if two-body truncation

				N_max (int): N_max beyond which to truncate according to truncation_rank
				
				output_filename (str): filename for output file
				
				output_format (str): h2 format of the file to be created (15099 or 15200)

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

#include "basis/jjjpn_operator.h"
#include "mcutils/parsing.h"
#include "mcutils/profiling.h"
#include "tbme/h2_io.h"

struct RunParameters
// store basic parameters for run
{
	// command-line inputs
	std::string orbital_filename;
	int n, l;
	HalfInt j;
	basis::OrbitalSpeciesPN particle_species;
	std::string output_filename;
	shell::H2Format output_format;
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
		output_format = shell::kVersion15099;
		weight_max = basis::WeightMax(0, 0);
	}
};

void ProcessArguments(int argc, char **argv, RunParameters& run_parameters)
// process command-line arguments
{
	// usage message
	if (argc - 1 != 9)
	{
		std::cout << "9 arguments expected, " << argc - 1 << " arguments found." << std::endl;
		std::cout << "Usage: h2pair orbital_filename n l 2*j particle_species truncation_rank N_max output_filename output_format" << std::endl;
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
	
	// input 6 (truncation rank)
	basis::Rank truncation_rank;
	std::string s6 = std::string(argv[6]);
	if (s6 == "1")
	{
		truncation_rank = basis::Rank::kOneBody;
	}
	else if (s6 == "2")
	{
		truncation_rank = basis::Rank::kTwoBody;
	}
	else
	{
		std::cerr << "ERROR: Invalid truncation rank \"" << argv[6] << "\"; expected 1 or 2" << std::endl;
		std::exit(EXIT_FAILURE);
	}
	
	// input 7 (truncation cutoff)
	int truncation_cutoff;
	std::istringstream parameter_stream_7(argv[7]);
  parameter_stream_7 >> truncation_cutoff;
  if (!parameter_stream_7)
	{
		std::cerr << "ERROR: Invalid truncation cutoff \"" << argv[7] << "\"; expected integer" << std::endl;
		std::exit(EXIT_FAILURE);
	}
	run_parameters.weight_max = basis::WeightMax(truncation_rank, truncation_cutoff);
	
	// input 8 (output filename)
	run_parameters.output_filename = argv[8];
	
	// input 9 (output format)
	std::istringstream parameter_stream_9(argv[9]);
	parameter_stream_9 >> run_parameters.output_format;
	if (!parameter_stream_9)
	{
		std::cerr << "ERROR: Invalid output format \"" << argv[9] << "\"; expected \"15099\" or \"15200\"" << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

int main(int argc, char** argv)
{
	// header
  std::cout << std::endl;
  std::cout << "h2pair -- generating H2 files for pair-counting operators" << std::endl;
  std::cout << "version: " << VCS_REVISION << std::endl;
  std::cout << std::endl;

  // read parameters
  RunParameters run_parameters;
  ProcessArguments(argc, argv, run_parameters);
	
	// start timing
  mcutils::SteadyTimer total_time;
  total_time.Start();
	
	// read and process orbitals.dat file and set up one-body orbital space
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
	basis::OrbitalSpacePN orbital_space(orbitals_list);
	basis::OrbitalSubspacePN orbital_subspace = 
		orbital_space.LookUpSubspace(basis::OrbitalSubspacePN::SubspaceLabelsType(run_parameters.particle_species));
	basis::OrbitalStatePN::StateLabelsType orbital_state_labels(run_parameters.n, run_parameters.l, run_parameters.j);
	if (!(orbital_subspace.ContainsState(orbital_state_labels)))
	{
		std::cerr << "\nERROR: Orbital index entered is either unsupported by orbitals.dat file entered or truncated away by N_max entered" << std::endl;
		std::exit(EXIT_FAILURE);
	}
	std::size_t orbital_index = orbital_subspace.LookUpStateIndex(orbital_state_labels);
	std::cout << "Done" << std::endl;
	
	// define quantum numbers for pair-counting operator
  const int J0 = 0;
  const int g0 = 0;
  const int Tz0 = 0;
	
	// set up two-body space
	std::cout << "Setting up two-body space ... ";
	const basis::TwoBodySpaceJJJPN two_body_space(orbital_space, run_parameters.weight_max, shell::kH2SpaceOrdering.at(run_parameters.output_format));
  const basis::TwoBodySectorsJJJPN two_body_sectors(two_body_space, J0, g0, Tz0);
	basis::OperatorBlocks<double> two_body_matrices;
  basis::SetOperatorToZero(two_body_sectors, two_body_matrices);
	std::cout << "Done" << std::endl;
	
	// locate and modify appropriate matrix element
	std::cout << "Setting matrix element for " << (run_parameters.particle_species == basis::OrbitalSpeciesPN::kP ? "proton" : "neutron") 
		<< " orbital with n = " << run_parameters.n << ", l = " << run_parameters.l << ", j = " << run_parameters.j << " ... ";
	basis::TwoBodySpeciesPN two_body_species;
	two_body_species = (run_parameters.particle_species == basis::OrbitalSpeciesPN::kP) ? basis::TwoBodySpeciesPN::kPP : basis::TwoBodySpeciesPN::kNN;
	basis::TwoBodySubspaceJJJPN::LabelsType two_body_subspace_labels(two_body_species, 0, 0); // operator nonzero only for J=0 and g=0
	std::size_t two_body_subspace_index = two_body_space.LookUpSubspaceIndex(two_body_subspace_labels);
	basis::TwoBodySubspaceJJJPN two_body_subspace = two_body_space.GetSubspace(two_body_subspace_index);
	std::size_t sector_index = two_body_sectors.LookUpSectorIndex(two_body_subspace_index, two_body_subspace_index);
	basis::TwoBodyStateJJJPN::StateLabelsType two_body_state_labels(orbital_index, orbital_index);
	std::size_t two_body_state_index = two_body_subspace.LookUpStateIndex(two_body_state_labels);
	if (two_body_state_index == basis::kNone)
	{
		std::cout << "\nERROR: Orbital index entered has been truncated away by N_max entered" << std::endl;
		std::exit(EXIT_FAILURE);
	}
	two_body_matrices[sector_index](two_body_state_index, two_body_state_index) = -2; // nonzero element for NAS pair-counting operator is always -2
	std::cout << "Done" << std::endl;
	
	// write h2 file
	std::cout << "Writing to " << run_parameters.output_filename << " ... ";
  shell::OutH2Stream output_stream(run_parameters.output_filename, orbital_space, two_body_space, two_body_sectors, run_parameters.output_format);
  for (std::size_t sector_index = 0; sector_index < two_body_sectors.size(); ++sector_index)
	{
		const basis::TwoBodySectorsJJJPN::SectorType& two_body_sector = two_body_sectors.GetSector(sector_index);
		const auto& matrix = two_body_matrices[sector_index];
		output_stream.WriteSector(sector_index, matrix, basis::NormalizationConversion::kNone);
	}
	std::cout << "Done" << std::endl;
	
	// end timing
  total_time.Stop();
  std::cout << "(Total time: " << total_time.ElapsedTime() << ")" << std::endl;
  std::cout << std::endl;

  // exit
  return EXIT_SUCCESS;
}