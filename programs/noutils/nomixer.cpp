/***************************************************************************
  nomixer.cpp --Normal ordered one body matrix elements generation
  

  Shwetha Vittal
  University of Notre Dame


*************************************************************************/
#include <sys/stat.h>
#include <cmath>

#include <omp.h>
#include <cstdlib>
#include <iostream>
#include <fstream> // addition
#include <iomanip>
#include <limits>
#include <set>
#include <string>
#include <vector>

#include "am/am.h"
#include "basis/nlj_orbital.h"
#include "basis/operator.h"
// #include "basis/jjjpn_operator.h"
#include "am/halfint.h"
#include "am/wigner_gsl.h"
#include "basis/proton_neutron.h"
#include "fmt/format.h"
#include "density/obdme_io.h"
#include "mcutils/parsing.h"
#include "obme/obme_operator.h"
#include "obme/obme_io.h"
#include "tbme/h2_io.h"

// (slv): adapted from obutils/obscalc-ob.cpp
////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////
// Store one-body operators
/* // (slv): I will need this to store the matrix elements
struct OneBodyOperator {
-  std::string name;
  basis::OrbitalSpaceLJPN space;
  basis::OrbitalSectorsLJPN sectors;
  basis::OperatorBlocks<double> blocks;

  explicit OneBodyOperator(const std::string& name__, const std::string& filename)
      : name(name__)
  {
    // read operator
    shell::InOBMEStream is(filename);
    is.Read(blocks);

    // get indexing
    basis::OrbitalSpaceLJPN ket_space;
    is.SetToIndexing(space, ket_space, sectors);
    assert(space.OrbitalInfo() == ket_space.OrbitalInfo());
    is.Close();
  }
};
*/

// Stores parameters for run
struct RunParameters {
  // filenames
  std::string output_filename; //(slv) Needed for saving matrix elements later
  std::string interaction_filename;
  std::ios_base::openmode file_mode; //(slv)TODO: Check what this does
  basis::OrbitalSpaceLJPN space;
  // std::vector<OneBodyOperator> operators;  // (slv) Not needed for storing ROBDMEs
  std::vector<std::unique_ptr<shell::InOBDMEStream>> density_streams;
  // shell::H2Format output_h2_format;
};

void PrintUsage(char** argv) { std::cout << "Usage: " << argv[0] << std::endl; }


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
    // std::cout<< "line count : " << line_count << std::endl;
    // std::cout<< line << std::endl;
  
}

double GetInteractionMatrixElement(shell::InH2Stream& input_stream,// (TODO (slv) should not pass an entire)
				   std::size_t state_index_a,
				   std::size_t state_index_b,
				   std::size_t state_index_c,
				   std::size_t state_index_d,
				   const basis::TwoBodySpeciesPN two_body_species,
				   HalfInt J, int g
				   ){
  // references
  const basis::OrbitalSpacePN& orbital_space = input_stream.orbital_space();
  const basis::TwoBodySpaceJJJPN& input_space = input_stream.space();
  const basis::TwoBodySectorsJJJPN& input_sectors = input_stream.sectors();

   // TODO (slv) This should be selected according to what species a,b,c and d belong
  //  basis::TwoBodySpeciesPN two_body_species = basis::TwoBodySpeciesPN::kPN; // selecting pp interaction

  // objects
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
  // tb_subspace_index = subspace_index_bra = subspace_index_ket
  std::size_t tb_subspace_index =
    tb_space.LookUpSubspaceIndex(basis::TwoBodySubspaceJJJPN::LabelsType(two_body_species, J , g));
  // std::cout<< "tb_subspace index " << tb_subspace_index << std::endl;

  const basis::TwoBodySubspaceJJJPN& tb_subspace = tb_space.GetSubspace(tb_subspace_index);

  // std::cout<< tb_subspace.DebugStr() << std::endl;

  // std::cout<< "tb_subspace size " << tb_subspace.size() << std::endl;
  // Here I am printing exactly what I have passed in the arguments but now
  // I am using an accessor to access
  // std::cout<< "tb_subspace J " << tb_subspace.J() << std::endl;
  // std::cout<< "tb_subspace g " << tb_subspace.g() << std::endl;
  
  std::size_t sector_index = tb_sectors.LookUpSectorIndex(tb_subspace_index, tb_subspace_index);
  // std::cout<< "sector_index " << sector_index << std::endl;

  std::size_t state_index_bra =
    tb_subspace.LookUpStateIndex(basis::TwoBodyStateJJJPN::StateLabelsType(state_index_a, state_index_b));
  std::size_t state_index_ket =
    tb_subspace.LookUpStateIndex(basis::TwoBodyStateJJJPN::StateLabelsType(state_index_c, state_index_d));

  double matrix_element = 0.0;

  if (state_index_bra != basis::kNone && state_index_ket != basis::kNone ){
    // std::cout<< "state_index_bra  " << state_index_bra << std::endl;
    // std::cout<< "state_index_ket " << state_index_ket << std::endl;
    
    // Access relevant matrix and get the matrix element
    Eigen::MatrixXd matrix;
    input_stream.ReadSector(sector_index, matrix);
  
    matrix_element = matrix(state_index_bra, state_index_ket);
  }
  
  return matrix_element;
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

  /*
  
  // testing with hardcoded state_indices
  std::size_t state_index_a = 0;
  std::size_t state_index_b = 0;
  std::size_t state_index_c = 0;
  std::size_t state_index_d = 0;
  int J = 1;
  int g = 0;
  double matrix_element = GetInteractionMatrixElement(input_stream, state_index_a, state_index_b, state_index_c, state_index_d, J, g );
  std::cout << "Matrix element " << matrix_element << std::endl;

  */


  /*
  
  // Accessing the density matrix elements
  for (const auto& density_stream : run_parameters.density_streams)
    {
      
      std::cout << fmt::format(
        "  {:>4.1f} {:>3d} {:>3d}  {:>4.1f} {:>3d} {:>3d}  ",
        float(density_stream->J_bra()), density_stream->g_bra(), density_stream->n_bra(),
	float(density_stream->J_ket()), density_stream->g_ket(), density_stream->n_ket())
	<< std::endl;
      const shell::InOBDMEStream& obdme_s = *density_stream;
      
      // get necessary density sectors
      basis::OrbitalSectorsLJPN density_sectors;
      basis::OperatorBlocks<double> density_blocks;
      //(slv) Here I should have used an accessor like op.sectors.J0() but density_sectors.J0()
      // throws an error; Assertion `(J0 >= J0_min()) && (J0 <= J0_max())' failed.
      int J0=0;
      obdme_s.GetMultipole(J0, density_sectors, density_blocks);
      for (std::size_t subspace_index_a=0; subspace_index_a<space.size(); ++subspace_index_a)
	{
	  for (std::size_t subspace_index_b=0; subspace_index_b<space.size(); ++subspace_index_b)
	    {
	      const auto& subspace_a = space.GetSubspace(subspace_index_a);
	      const auto& subspace_b = space.GetSubspace(subspace_index_b);
	      auto sector_index = density_sectors.LookUpSectorIndex(subspace_index_a, subspace_index_b);
	      if (sector_index == basis::kNone) continue;

	      for (std::size_t state_index_a = 0; state_index_a < subspace_a.size(); ++state_index_a) {
		for (std::size_t state_index_b = 0; state_index_b < subspace_b.size(); ++state_index_b) {
		  std::cout<< density_blocks[sector_index](state_index_a, state_index_b) << " ";
		}
		std::cout<< "\n";
	      }
	      std::cout<< "End of subspace index {" << subspace_index_a << "," << subspace_index_b << "}\n";
	    }
	}
    }


    */

  
  
  //========================================================================
  //
  //Implementing Normal Ordering using the exercises above to access the
  //matrix elements from density and interaction files.
  //
  //========================================================================
  

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
      
      //(slv) Here I should perhaps use an accessor like op.sectors.J0() but density_sectors.J0()
      // throws an error; Assertion `(J0 >= J0_min()) && (J0 <= J0_max())' failed.
      // I presume this is equivalent to multipolarity K of the density operator, ref mfdn.rppobdme.info
      // this is my lambda
      int J0=0;
      int g = 0; // setting positive parity
      
      //This stores density_sectors and  and density_blocks of multipole J0
      obdme_s.GetMultipole(J0, density_sectors, density_blocks);

      

      std::cout<< "Size of the LJPN space : " << space.size() << std::endl;
      
      // the outer for loop is to compute each matrix element
      // j_b and j_d have to be same because we are taking lambda = 0
      // but |b> and |d> need not be
      for (std::size_t subspace_index_b=0; subspace_index_b<space.size(); ++subspace_index_b)
	{
	  for (std::size_t subspace_index_d=0; subspace_index_d<space.size(); ++subspace_index_d)
	    {

	      const auto& subspace_b = space.GetSubspace(subspace_index_b);
	      const auto& subspace_d = space.GetSubspace(subspace_index_d);
	      if(subspace_b.j() != subspace_d.j()){
			continue;
		      }
	      auto sector_index_output = density_sectors.LookUpSectorIndex(subspace_index_b, subspace_index_d);
	      std::cout<< "Size of subspace b with index : " << subspace_index_b << " in LJPN space of size : " << subspace_b.size() << std::endl;
	      std::cout<< "Size of subspace d with index : " << subspace_index_d << " in LJPN space of size : " << subspace_d.size() << std::endl;

	      for (std::size_t state_index_b = 0; state_index_b < subspace_b.size(); ++state_index_b)
		{
		  for (std::size_t state_index_d = 0; state_index_d < subspace_d.size(); ++state_index_d)
		    {
		      
		      // Set the output matrix element to be zero
		      double matrix_element_output = 0.0;

		      //the following is the summation

		      for (std::size_t subspace_index_a=0; subspace_index_a<space.size(); ++subspace_index_a)
			{
			  for (std::size_t subspace_index_c=0; subspace_index_c<space.size(); ++subspace_index_c)
			    {
		  
			      const auto& subspace_a = space.GetSubspace(subspace_index_a);
			      const auto& subspace_c = space.GetSubspace(subspace_index_c);

			      int species_a =int(subspace_a.orbital_species());
			      int species_b =int(subspace_b.orbital_species());
			      int species_c =int(subspace_c.orbital_species());
			      int species_d =int(subspace_d.orbital_species());
			      basis::TwoBodySpeciesPN two_body_species;
			      
			      if(species_a != species_c || species_b != species_d) continue;

			      if(species_a==0 && species_b==0){
				two_body_species= basis::TwoBodySpeciesPN::kPP;
			      }
			      else if(species_a==1 && species_b==0 ||species_a==0 && species_b==1  ){
				two_body_species = basis::TwoBodySpeciesPN::kPN;
			      }
			      else if(species_a==1 && species_b==1){
				two_body_species = basis::TwoBodySpeciesPN::kNN;
			      }
				
			      if(subspace_a.j() != subspace_c.j()) continue;

			      HalfInt::vector valid_J = am::ProductAngularMomenta(subspace_a.j(),subspace_b.j());
			      // Here I need to be decent to use iterators
				for(int i=0; i<valid_J.size(); i++)
				{
				  // This is where I have to check the triangle condition
				  // Here I am taking J = J'
				  // The following is zero condition for the 6J symbol
				  // (J,   j_b,   j_a)
				  // (j_d, J , lambda)
				  
				  if (!am::AllowedTriangle(valid_J[i], subspace_b.j(), subspace_a.j())||
				      !am::AllowedTriangle(subspace_d.j(), valid_J[i], subspace_a.j())||
				      !am::AllowedTriangle(subspace_d.j(), subspace_b.j(), J0 )||
				      !am::AllowedTriangle(valid_J[i], valid_J[i], J0)){
				    continue;
				  }

				  // check for 6J symbol to be zero; skip the rest of the loop if accidental zero
				  double cg_coeff_6j = am::Wigner6J(valid_J[i], subspace_b.j(), subspace_a.j(),
								    subspace_d.j(), valid_J[i], J0);

				  if (std::abs(cg_coeff_6j) < 1e-8) continue;
			      
			      
				  // TO DO : Must figure out how to include this phase factor
				  double phase_factor = 1; // - 2*((3 * subspace_a.j() + 3*J0 + 2*subspace_b.j()+ 2* J + J +subspace_d.j()) % 2);

				  phase_factor *= Hat(subspace_a.j());
				  phase_factor *= Hat(J0);
				
				
				  auto sector_index = density_sectors.LookUpSectorIndex(subspace_index_a, subspace_index_c);
		  
				  if (sector_index == basis::kNone) continue;

				  for (std::size_t state_index_a = 0; state_index_a < subspace_a.size(); ++state_index_a) {
				    for (std::size_t state_index_c = 0; state_index_c < subspace_c.size(); ++state_index_c) {

				      double interaction_matrix_element =
					GetInteractionMatrixElement(input_stream, state_index_a, state_index_b, state_index_c, state_index_d, two_body_species,valid_J[i], g );
				      if(interaction_matrix_element != 0.0){

					// std::cout<< "density matrix element : " << std::endl;
					// std::cout<< "< " << state_index_a <<"| ca_dag cc |" << state_index_c << ">  : " 
					//	     << density_blocks[sector_index](state_index_a, state_index_c) << " \n";

					// std::cout<< "interaction_matrix_element : " << "< " << state_index_a
					//	     << ", " << state_index_b << "| " << interaction_matrix_element
					//	     << " |" << state_index_c << ", " << state_index_d
					//	     << "> " << std::endl;
					matrix_element_output += phase_factor * cg_coeff_6j * density_blocks[sector_index](state_index_a, state_index_c) *
					  interaction_matrix_element;
				      }
				    }
				    // std::cout<< "\n";
				  }
		      
				}
			      
			    }
			}
		      std::cout<< "< "<< state_index_b << " , " << subspace_index_b << "| " << " V " << " |"
			       << state_index_d << " , " << subspace_index_d  << "> = " << matrix_element_output << std::endl;
		    }
		}
	    }
	}
    }
  
  return EXIT_SUCCESS;
}


