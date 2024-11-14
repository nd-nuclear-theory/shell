/**********************************************************************************
  nomixer.cpp --Normal ordered one body matrix element generation
  Beta version to improve nomixer_V0.cpp
  
  Shwetha Vittal
  University of Notre Dame

*********************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "am/am.h"
#include "am/halfint.h"
#include "am/wigner_gsl.h"
#include "basis/nlj_orbital.h"
#include "basis/operator.h"
#include "basis/proton_neutron.h"
#include "density/obdme_io.h"
#include "fmt/format.h"
#include "mcutils/parsing.h"
#include "obme/obme_operator.h"
#include "obme/obme_io.h"
#include "tbme/h2_io.h"

// Stores parameters for run
struct RunParameters {
  // filenames
  std::string output_filename; //(slv) Needed for saving matrix elements
  std::string interaction_filename;
  
  std::ios_base::openmode file_mode; 
  // one body space
  basis::OrbitalSpaceLJPN space;
  std::vector<std::unique_ptr<shell::InOBDMEStream>> density_streams;

};

// Temporary function to read parameters
void ReadParameters(RunParameters& run_parameters, std::string input_filename) {
  std::string line;
  int line_count = 0;
  std::ifstream inFile;
  inFile.open(input_filename);

    if(inFile.is_open()){
      while (getline(inFile,line,'\n')) {
	line_count = line_count +1;
	std::istringstream line_stream(line);
	std::string keyword;
	line_stream >> keyword;
	
	// select action based on keyword
	if (keyword == "set-output-file") {
	  line_stream >> run_parameters.output_filename;
	  run_parameters.file_mode = std::ios_base::trunc;
	  
	  if (!line_stream.eof()) {
	    std::string mode;
	    line_stream >> mode;
	    if (mode == "append")
	      run_parameters.file_mode = std::ios_base::app;
	  }
	  mcutils::ParsingCheck(line_stream, line_count, line);
	  mcutils::FileExistCheck(
				  run_parameters.output_filename, false,
				  run_parameters.file_mode==std::ios_base::trunc
				  );
	}
	else if (keyword == "set-indexing") {
	  std::string orbital_filename;
	  line_stream >> orbital_filename;
	  mcutils::ParsingCheck(line_stream, line_count, line);
	  mcutils::FileExistCheck(orbital_filename, true, false);

	  std::ifstream orbital_stream(orbital_filename);
	  std::vector<basis::OrbitalPNInfo> input_orbitals =
	    basis::ParseOrbitalPNStream(orbital_stream, true);
	  run_parameters.space = basis::OrbitalSpaceLJPN(input_orbitals);
	}
	// Not reading the interaction file here. Just doing preliminary check whether it exists.
	else if (keyword == "interaction-file") {
	  std::string interaction_filename;
	  line_stream >> interaction_filename;
	  mcutils::ParsingCheck(line_stream, line_count, line);
	  mcutils::FileExistCheck(interaction_filename, true, false);
	  run_parameters.interaction_filename = interaction_filename;
  	}
		
	else if (keyword == "define-densities") {
	  //std::unique_ptr<shell::InOBDMEStream> density_stream;
	  float Jf, Ji;
	  int gf, nf, gi, ni;
	  std::string robdme_filename, robdme_info_filename="";
	  line_stream >> Jf >> gf >> nf >> Ji >> gi >> ni >> robdme_filename;
	  mcutils::ParsingCheck(line_stream, line_count, line);
	  mcutils::FileExistCheck(robdme_filename, true, false);
	  // std::cout << "robdme_filename : " << robdme_filename <<std::endl;
	  if (!line_stream.eof()) {
	    line_stream >> robdme_info_filename;
	    // std::cout << "robdme_info_filename : " << robdme_info_filename <<std::endl;
	    mcutils::ParsingCheck(line_stream, line_count, line);
	    mcutils::FileExistCheck(robdme_info_filename, true, false);
	    // std::cout << "Calling INOBDMEStreamMulti" << std::endl;
	    // construct multi-file stream
	    run_parameters.density_streams.emplace_back(new shell::InOBDMEStreamMulti(
            robdme_info_filename, robdme_filename, run_parameters.space,
            HalfInt(2*Jf,2), gf, nf, HalfInt(2*Ji, 2), gi, ni
	    ));
	  } else {
	    // construct single-file stream
	    run_parameters.density_streams.emplace_back(new shell::InOBDMEStreamSingle(
            robdme_filename, run_parameters.space,
            HalfInt(2*Jf,2), gf, nf, HalfInt(2*Ji, 2), gi, ni
	    ));
	  }
	}
      }
    } else {
      std::cout<< "Unable to open the file : " << input_filename << std::endl;
    }
  
}

double GetInteractionMatrixElement(shell::InH2Stream& input_stream,// (TODO (slv) should not pass an entire stream)
				   // But the reference to the object is not same as the object itself
				   std::size_t state_index_a,
				   std::size_t state_index_b,
				   std::size_t state_index_c,
				   std::size_t state_index_d,
				   const basis::TwoBodySpeciesPN two_body_species,
				   HalfInt J,HalfInt JPrime,  int g
				   ){
  // references
  const basis::OrbitalSpacePN& orbital_space = input_stream.orbital_space();
  const basis::TwoBodySpaceJJJPN& input_space = input_stream.space();
  const basis::TwoBodySectorsJJJPN& input_sectors = input_stream.sectors();

  // Hard coded right now but output_h2_format  can be passed later
  shell::H2Format output_h2_format = shell::kVersion15099;
  
  basis::TwoBodySpaceJJJPNOrdering space_ordering =
    shell::kH2SpaceOrdering.at(output_h2_format);

  const basis::TwoBodySpaceJJJPN tb_space = basis::TwoBodySpaceJJJPN(
      orbital_space,
      input_space.weight_max(),
      space_ordering
    );

  // This gives all the subspaces in a two body space.
  // std::cout<< tb_space.DebugStr() << std::endl;
  
  const basis::TwoBodySectorsJJJPN tb_sectors = basis::TwoBodySectorsJJJPN(
      tb_space, input_sectors.J0(), input_sectors.g0(), input_sectors.Tz0()
    );

  // look up indices
  std::size_t tb_subspace_index_bra =
    tb_space.LookUpSubspaceIndex(basis::TwoBodySubspaceJJJPN::LabelsType(two_body_species, J , g));
  std::size_t tb_subspace_index_ket =
    tb_space.LookUpSubspaceIndex(basis::TwoBodySubspaceJJJPN::LabelsType(two_body_species, JPrime , g));
  
  const basis::TwoBodySubspaceJJJPN& tb_subspace_bra = tb_space.GetSubspace(tb_subspace_index_bra);
  const basis::TwoBodySubspaceJJJPN& tb_subspace_ket = tb_space.GetSubspace(tb_subspace_index_ket);

  std::size_t sector_index = tb_sectors.LookUpSectorIndex(tb_subspace_index_bra, tb_subspace_index_ket);
  
  std::size_t state_index_bra =
    tb_subspace_bra.LookUpStateIndex(basis::TwoBodyStateJJJPN::StateLabelsType(state_index_a, state_index_b));
  std::size_t state_index_ket =
    tb_subspace_ket.LookUpStateIndex(basis::TwoBodyStateJJJPN::StateLabelsType(state_index_c, state_index_d));

  double matrix_element = 0.0;

  if (state_index_bra != basis::kNone && state_index_ket != basis::kNone && sector_index != basis::kNone ){

    // Access relevant matrix and get the matrix element
    Eigen::MatrixXd matrix;
    std::cout<< "tb sector index " << sector_index << std::endl;
    
    input_stream.ReadSector(sector_index, matrix);
  
    matrix_element = matrix(state_index_bra, state_index_ket);
    /*  
    std::cout<< "found matrix element with state indices : " << state_index_bra << "," << state_index_ket << " )"
	     << "sector_index : " << sector_index
	     << "state indices a, b " << state_index_a << " , " << state_index_b
	     << "state indices c, d " << state_index_c << " , " << state_index_d <<std::endl;
	     }
	     else{
	     std::cout<< "no matrix element with state indices : ( " << state_index_bra << "," << state_index_ket << " )"
	     << "sector_index : " << sector_index
	     << "state indices a, b " << state_index_a << " , " << state_index_b
	     << "state indices c, d " << state_index_c << " , " << state_index_d <<std::endl;      
    */
  }
  return matrix_element;
}



double CalculateMatrixElement(const basis::OrbitalSpaceLJPN& space,
			      shell::InH2Stream& input_stream,
			      basis::OrbitalSectorsLJPN& density_sectors,
			      basis::OperatorBlocks<double>& density_blocks,
			      std::size_t state_index_b,
			      std::size_t state_index_d,
			      basis::OrbitalSubspaceLJPN subspace_b,// 
			      basis::OrbitalSubspaceLJPN subspace_d,
			      int species_b, int species_d, int J0, int g0
			      ){

  double matrix_element_output = 0.0;

  for (std::size_t subspace_index_a=0; subspace_index_a<space.size(); ++subspace_index_a)
  //  for (basis::OrbitalSpaceLJPN::iterator subspace_iterator_a = space.begin(); subspace_iterator_a != space.end(); subspace_iterator_a++)
    {
     
      for (std::size_t subspace_index_c=0; subspace_index_c<space.size(); ++subspace_index_c)
	{
		  
	  const auto& subspace_a = space.GetSubspace(subspace_index_a);
	  const auto& subspace_c = space.GetSubspace(subspace_index_c);

	  int species_a =int(subspace_a.orbital_species());
	  int species_c =int(subspace_c.orbital_species());
	  
	  basis::TwoBodySpeciesPN two_body_species;
			      
	  // (slv) Double check the following conditions
	  if(species_a==0 && species_b==0 && species_c==0 && species_d==0){
	    two_body_species= basis::TwoBodySpeciesPN::kPP;
	  }
	  else if(species_a==1 && species_b==1 && species_c==1 && species_d==1){
	two_body_species = basis::TwoBodySpeciesPN::kNN;
	  }
	  else if(species_a != species_b && species_c != species_d){
	    two_body_species = basis::TwoBodySpeciesPN::kPN;
	  }
	  else continue;
	   

	  //The following condition ensures the coupled angular momentum =0.
	  if(subspace_a.j() != subspace_c.j()) continue;
	  
	  //HalfInt::vector valid_J = am::ProductAngularMomenta(subspace_a.j(),subspace_b.j());
	  //HalfInt::vector valid_JPrime = am::ProductAngularMomenta(subspace_c.j(),subspace_d.j());

	  // The following is temporarily  set to test whether RestrictedCalculateMatrixElement() and 
	  // CalculateMatrixElement() gives the same result when J = 0 and J' = 0
	  std::vector<int> valid_J = {0};
	  std::vector<int> valid_JPrime = {0};

	  // Here I need to be decent to use iterators
	  for(int x=0; x<valid_J.size(); x++)
	    {
	      for(int y=0; y<valid_JPrime.size();y++)
		{
		  // This is where triangle condition is checked
		  // Here J is not required to be same as J'
		  // The following is zero condition for the 6J symbol
		  // (J,   j_b,   j_a)
		  // (j_d, J' , lambda)
				  
		  if (!am::AllowedTriangle(valid_J[x], subspace_b.j(), subspace_a.j())||
		      !am::AllowedTriangle(subspace_d.j(), valid_JPrime[y], subspace_a.j())||
		      !am::AllowedTriangle(subspace_d.j(), subspace_b.j(), J0 )||
		      !am::AllowedTriangle(valid_J[y], valid_JPrime[y], J0)){
		    continue;
		  }

		  // check for 6J symbol to be zero; skip the rest of the loop if accidental zero
		  double cg_coeff_6j = am::Wigner6J(valid_J[x], subspace_b.j(), subspace_a.j(),
						    subspace_d.j(), valid_JPrime[y], J0);
		  
		  if (std::abs(cg_coeff_6j) < 1e-8) continue;
		  
		  double phase_factor = double(ParitySign( subspace_a.j() +
							   subspace_b.j()+
							   valid_J[x]+ J0));
		  phase_factor /= Hat(subspace_a.j());
		  phase_factor *= Hat(valid_J[x]);
		  phase_factor *= Hat(valid_JPrime[y]);

		  auto sector_index = density_sectors.LookUpSectorIndex(subspace_index_a, subspace_index_c);
					
		  if (sector_index == basis::kNone) continue;
				  
		  for (std::size_t state_index_a = 0; state_index_a < subspace_a.size(); ++state_index_a)
		    {
		      for (std::size_t state_index_c = 0; state_index_c < subspace_c.size(); ++state_index_c)
			{
			  
			  double interaction_matrix_element =
			    GetInteractionMatrixElement(input_stream, state_index_a, state_index_b, state_index_c,
							state_index_d, two_body_species, valid_J[x],valid_JPrime[y], g0 );
			  if(interaction_matrix_element != 0.0)
			    {
			      
			      matrix_element_output +=  phase_factor * cg_coeff_6j *
				density_blocks[sector_index](state_index_a, state_index_c) *
				interaction_matrix_element;
			    }
			}
		    }
		}
	    }
	}
    } 
  
  return matrix_element_output;
}

double RestrictedCalculateMatrixElement(const basis::OrbitalSpaceLJPN& space,
			      shell::InH2Stream& input_stream,
			      basis::OrbitalSectorsLJPN& density_sectors,
			      basis::OperatorBlocks<double>& density_blocks,
			      std::size_t state_index_b,
			      std::size_t state_index_d,
			      basis::OrbitalSubspaceLJPN subspace_b,// 
			      basis::OrbitalSubspaceLJPN subspace_d,
			      int species_b, int species_d, int J0, int g0
			      ){
  double matrix_element_output = 0.0;

  for (std::size_t subspace_index_a=0; subspace_index_a<space.size(); ++subspace_index_a)
    {
      for (std::size_t subspace_index_c=0; subspace_index_c<space.size(); ++subspace_index_c)
	{
		  
	  const auto& subspace_a = space.GetSubspace(subspace_index_a);
	  const auto& subspace_c = space.GetSubspace(subspace_index_c);

	  int species_a =int(subspace_a.orbital_species());
	  int species_c =int(subspace_c.orbital_species());
	  
	  basis::TwoBodySpeciesPN two_body_species;
			      
	  // (slv) Double check the following conditions
	  if(species_a==0 && species_b==0 && species_c==0 && species_d==0){
	    two_body_species= basis::TwoBodySpeciesPN::kPP;
	  }
	  else if(species_a==1 && species_b==1 && species_c==1 && species_d==1){
	    two_body_species = basis::TwoBodySpeciesPN::kNN;
	  }
	  else if(species_a != species_b && species_c != species_d){
	    two_body_species = basis::TwoBodySpeciesPN::kPN;
	  }
	  else continue;

	  //The following condition ensures the coupled angular momentum =0.
	  if(subspace_a.j() != subspace_c.j()) continue;

	  //(slv): The following initializations makes no sense because the expressions were derived
	  // explicitly for the case when J and J' are equal to zero
	  // HalfInt::vector valid_J = am::ProductAngularMomenta(subspace_a.j(),subspace_b.j());
	  // HalfInt::vector valid_JPrime = am::ProductAngularMomenta(subspace_c.j(),subspace_d.j());

	  // The following is set to test whether RestrictedCalculateMatrixElement() and CalculateMatrixElement
	  // gives the same result when J = 0 and J' = 0
	  
	  std::vector<int> valid_J = {0};
	  std::vector<int> valid_JPrime = {0};

	  // Here I need to be decent to use iterators
	  for(int x=0; x<valid_J.size(); x++)
	    {
	      for(int y=0; y<valid_JPrime.size();y++)
		{				  
		  if (!am::AllowedTriangle(valid_J[x], subspace_b.j(), subspace_a.j())||
		      !am::AllowedTriangle(subspace_d.j(), valid_JPrime[y], subspace_a.j())||
		      !am::AllowedTriangle(subspace_d.j(), subspace_b.j(), J0 )||
		      !am::AllowedTriangle(valid_J[y], valid_JPrime[y], J0)){
		    continue;
		  }
		  
		  		  
		  double phase_factor = 1.0;
		  phase_factor /= Hat(subspace_a.j());
		  phase_factor /= Hat(subspace_c.j());
		  phase_factor /= Hat(valid_J[x]);


		  auto sector_index = density_sectors.LookUpSectorIndex(subspace_index_a, subspace_index_c);
					
		  if (sector_index == basis::kNone) continue;
				  
		  for (std::size_t state_index_a = 0; state_index_a < subspace_a.size(); ++state_index_a)
		    {
		      for (std::size_t state_index_c = 0; state_index_c < subspace_c.size(); ++state_index_c)
			{	   
			  double interaction_matrix_element =
			    GetInteractionMatrixElement(input_stream, state_index_a, state_index_b, state_index_c,
							state_index_d, two_body_species, valid_J[x],valid_JPrime[y], g0 );
			  
			  
			  if(interaction_matrix_element != 0.0)
			    {
			      
			      matrix_element_output +=  phase_factor * 
				density_blocks[sector_index](state_index_a, state_index_c) *
				interaction_matrix_element;
			    }
			}
		    }
		}
	    }
	}
    } 
  
  return matrix_element_output;
}

int main(){
  
  RunParameters run_parameters;
  run_parameters.output_filename = "nomixer-test.out";
  std::string input_filename= "test.in";
  ReadParameters(run_parameters, input_filename);
  const basis::OrbitalSpaceLJPN& space = run_parameters.space;// ROBDME stored in LJPN scheme 

  std::cout << "Input stream" << std::endl;
  shell::InH2Stream input_stream(run_parameters.interaction_filename);
  std::cout << input_stream.DiagnosticStr();
  std::cout << std::endl;

  const int J0= 0; // This is the rank of the interaction tensor
  const int g0 = 0; // setting positive parity

  // open output
  std::ofstream out_stream(run_parameters.output_filename, run_parameters.file_mode);
  mcutils::StreamCheck( bool(out_stream), run_parameters.output_filename,
		       "Failure opening file for output"
		       );
  std::ostringstream section_stream;
  section_stream << "[One-body observable]" << std::endl;
  section_stream << fmt::format(
				"# {:>3} {:>3}  {:s}",
				"J0", "g0", "name"
				) << std::endl;
  section_stream << fmt::format(
				"  {:>3d} {:>3d}  {:s}",
				J0, g0 , "Hamiltonian-Like" // (slv) :I need to input this
				) << std::endl;
  section_stream << fmt::format(
				"# {:>4} {:>3} {:>3}  {:>4} {:>3} {:>3}  {:>15s}",
				"Jf", "gf", "nf",
				"Ji", "gi", "ni",
				"rme"
				) << std::endl;

  // Access the density matrix elements (there is just one statrobdme file now)
  for (const auto& density_stream : run_parameters.density_streams)
    {
      // Print the bra and ket states of the static ROBDME
      std::cout << fmt::format(
        "  {:>4.1f} {:>3d} {:>3d}  {:>4.1f} {:>3d} {:>3d}  ",
        float(density_stream->J_bra()), density_stream->g_bra(), density_stream->n_bra(),
	float(density_stream->J_ket()), density_stream->g_ket(), density_stream->n_ket())
	<< std::endl;
      
      const shell::InOBDMEStream& obdme_s = *density_stream;
      
      // get necessary density sectors
      basis::OrbitalSectorsLJPN density_sectors;
      basis::OperatorBlocks<double> density_blocks;
      
      
      //This stores density_sectors and  and density_blocks of multipole J0
      obdme_s.GetMultipole(J0, density_sectors, density_blocks);
      double matrix_element_output =0.0;
      
      for (std::size_t subspace_index_b=0; subspace_index_b<space.size(); ++subspace_index_b)
	{
	  for (std::size_t subspace_index_d=0; subspace_index_d<space.size(); ++subspace_index_d)
	    {

	      const auto& subspace_b = space.GetSubspace(subspace_index_b);
	      const auto& subspace_d = space.GetSubspace(subspace_index_d);
	      int species_b =int(subspace_b.orbital_species());
	      int species_d =int(subspace_d.orbital_species());

	      // if(species_b != species_d) continue;

	      // The below condition is not applicable if transition OBDMEs are available
	      
	      if(subspace_b.j() != subspace_d.j()) continue;

	      for (std::size_t state_index_b = 0; state_index_b < subspace_b.size(); ++state_index_b)
			{
			for (std::size_t state_index_d = 0; state_index_d < subspace_d.size(); ++state_index_d)
				{
				// The following should be updated if transition OBDMEs are allowed
				// Here the density_sectors and density_blocks are common for both [c_dag_a c_til_c]_0 and [c_dag_b c_til_d]_lambda
				auto sector_index = density_sectors.LookUpSectorIndex(subspace_index_b, subspace_index_d);
				std::cout<< "Sector_index : " << sector_index << std::endl;

				if(sector_index == basis::kNone) continue;

				matrix_element_output = density_blocks[sector_index](state_index_b, state_index_d);
				
				/*
				matrix_element_output *= CalculateMatrixElement(space, input_stream, density_sectors,
										density_blocks, state_index_b,
										state_index_d, subspace_b, subspace_d,
										species_b,species_d, J0, g0 );
				*/
				matrix_element_output *= RestrictedCalculateMatrixElement(space, input_stream, density_sectors,
										density_blocks, state_index_b,
										state_index_d, subspace_b, subspace_d,
										species_b,species_d, J0, g0 );      
				/*
				section_stream << fmt::format(
									"  {:>4.1f} {:>3d} {:>3d}  {:>4.1f} {:>3d} {:>3d}  {:15.8e}",
									float(subspace_b.j()), subspace_b.g(), state_index_b,
									float(subspace_d.j()), subspace_d.g(), state_index_d,
								matrix_element_output * pow(Hat(J0),2)
								) << std::endl;
				std::cout<< "< "<< state_index_b << " , " << subspace_index_b << "| " << " V " << " |"
					<< state_index_d << " , " << subspace_index_d  << "> = "
					<< std::pow(Hat(J0),2) *matrix_element_output << std::endl;

				*/
				section_stream << fmt::format(
									"  {:>4.1f} {:>3d} {:>3d}  {:>4.1f} {:>3d} {:>3d}  {:15.8e}",
									float(subspace_b.j()), subspace_b.g(), state_index_b,
									float(subspace_d.j()), subspace_d.g(), state_index_d,
								matrix_element_output
								) << std::endl;
				std::cout<< "< "<< state_index_b << " , " << subspace_index_b << "| " << " V " << " |"
					<< state_index_d << " , " << subspace_index_d  << "> = "
					<< matrix_element_output << std::endl;

				}
			}

	    }
	}
      out_stream << section_stream.str() << std::flush;
    }

  return EXIT_SUCCESS;
}
