/****************************************************************
  relative-filter.cpp

  Copy relative matrix elements, with filtering operations on their values.

  These filtering operations may effectively accomplish "truncation" of the
  operator, by zeroing out RMEs, but the actual truncation of the storage is unchanged between source and target operators.

  See lsjt_operator.h for documentation of operator storage and the relative
  operator file format.

  Standard input:
    source_filename target_filename
    filter_name [filter_parameter_1 ...]

  The supported filter names and accompanying parameters are:

    identity -- no filter (simply copy)
    Nrelmax <cutoff> -- cutoff on bra and ket relative quanta
    N0max <cutoff> -- cutoff on quanta carried by operator (i.e., |delta Nrel|) 

  Example (relative-filter.in):
    coulomb_Nmax20_p_rel.dat coulomb_Nmax20_p_Nrelmax02_rel.dat
    Nrelmax 2

  Language: C++14

  Mark A. Caprio
  University of Notre Dame

  05/15/19 (mac/jh): Created, absorbing code in relative-gen.cpp.

****************************************************************/

#include "basis/lsjt_operator.h"
#include "basis/proton_neutron.h"
#include "fmt/format.h"
#include "mcutils/parsing.h"
#include "relative/relative_me.h"

////////////////////////////////////////////////////////////////
// parameter input
////////////////////////////////////////////////////////////////

struct Parameters
// Container for run input parameters.
{
  std::string source_filename;
  std::string target_filename;
  std::string filter_name;

  // optional parameters for specific operators
  int cutoff;
};

void ReadParameters(Parameters& parameters)
// Read run parameters from stdin.
//
// Arguments:
//   parameters (Parameters, output) :
//     container for input parameters
{

  // set up line counter for use in error messages
  std::string line;
  int line_count = 0;
  
  // line 1: operator filenames
  {
    ++line_count;
    std::getline(std::cin,line); // commented out by J.H.
    std::istringstream line_stream(line);
    line_stream >> parameters.source_filename
                >> parameters.target_filename;
    mcutils::ParsingCheck(line_stream,line_count,line);
  }

  // line 2: filtering parameters
  {
    ++line_count;
    std::getline(std::cin,line);
    std::istringstream line_stream(line);
    line_stream >> parameters.filter_name
                >> parameters.cutoff;
    mcutils::ParsingCheck(line_stream,line_count,line);
  }

}

////////////////////////////////////////////////////////////////
// operator work
////////////////////////////////////////////////////////////////

void FilterOperator(
    const Parameters& parameters,
    const basis::RelativeSpaceLSJT& space,
    const basis::RelativeOperatorParametersLSJT& operator_parameters,
    const std::array<basis::RelativeSectorsLSJT,3>& component_sectors,
    const std::array<basis::OperatorBlocks<double>,3>& source_component_blocks,
    std::array<basis::OperatorBlocks<double>,3>& target_component_blocks,
    bool verbose
  )
// Copy and filter operator.
//
// Filtering operations:
//   If the filter name is "identity", no filtering is performed (the file is just copied).
//   If the filter name is "Nrelmax", the RMEs for which the number of oscillator quanta of relative motion (Nrel)
//     carried by the bra or ket state is greater than the cutoff parameter are zeroed out.
//   If the filter name is "N0max", the RMEs for which the number of oscillator quanta carried by the
//     operator, i.e. |Nrel_bra - Nrel_ket|, is greater than the cutoff parameter are zeroed out.
//
// Arguments:
//   parameters (input): includes tensorial properties of operator
//      choice of operator to use
//   space (input): space for source and target operators
//   component_sectors (input): sectors for source and target operators
//   source_component_blocks (input): RMEs for source operator
//   target_component_blocks (output): RMEs for target operator
//   verbose (input): verbosity flag
{

  // iterate over operator components
  for (int T0=operator_parameters.T0_min; T0<=operator_parameters.T0_max; ++T0)
      // for each isospin component

    for (int sector_index = 0; sector_index < component_sectors[T0].size(); ++sector_index)
      // for each sector
      {
        // extract sector
        const typename basis::RelativeSectorsLSJT::SectorType& sector = component_sectors[T0].GetSector(sector_index);
	const typename basis::RelativeSectorsLSJT::SubspaceType& bra_subspace = sector.bra_subspace();
	const typename basis::RelativeSectorsLSJT::SubspaceType& ket_subspace = sector.ket_subspace();

        // retrieve source block
        const basis::OperatorBlock<double>& source_block = source_component_blocks[T0][sector_index];

        // allocate target block
        basis::OperatorBlock<double> target_block = basis::OperatorBlock<double>::Zero(bra_subspace.size(),ket_subspace.size());
        
        // copy matrix elements
        for (int bra_index = 0; bra_index < bra_subspace.size(); ++bra_index)
          for (int ket_index = 0; ket_index < ket_subspace.size(); ++ket_index)
            {
              // diagonal sector: restrict to upper triangle
              if (sector.IsDiagonal())
                if (!(bra_index<=ket_index))
                  continue;

              // define states (for easy access to quantum numbers)
	      const basis::RelativeStateLSJT bra(bra_subspace,bra_index);
	      const basis::RelativeStateLSJT ket(ket_subspace,ket_index);

              // copy matrix element

              if(parameters.filter_name=="identity")
	        target_block(bra_index,ket_index) = source_block(bra_index,ket_index);
	      else if(parameters.filter_name=="Nrelmax")
                if((bra.N()<=parameters.cutoff)&&(ket.N()<=parameters.cutoff))
                  target_block(bra_index,ket_index) = source_block(bra_index,ket_index);
	        else
                  target_block(bra_index,ket_index) = 0.0;
	      else if(parameters.filter_name=="N0max")
		if(abs(bra.N()-ket.N())<=parameters.cutoff)
                  target_block(bra_index,ket_index) = source_block(bra_index,ket_index);
	        else
	          target_block(bra_index,ket_index) = 0.0;
            }

        // diagnostics
        if (verbose)
          {
            // 
            std::cout 
              << " sector " 
              << sector_index
              << std::endl;
          }

        // store block
        target_component_blocks[T0].push_back(target_block);
      }
}

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  // header
  std::cout << std::endl;
  std::cout << "relative-filter -- copy and filter relative operator matrix elements" << std::endl;
  std::cout << "version: " VCS_REVISION << std::endl;
  std::cout << std::endl;

  // parameter input
  Parameters parameters;
  ReadParameters(parameters);

  // read operator
  basis::RelativeSpaceLSJT space;
  basis::RelativeOperatorParametersLSJT operator_parameters;
  std::array<basis::RelativeSectorsLSJT,3> component_sectors;
  std::array<basis::OperatorBlocks<double>,3> source_component_blocks;
  basis::ReadRelativeOperatorLSJT(
      parameters.source_filename,
      space,operator_parameters,
      component_sectors,source_component_blocks,
      true // verbose
    );

  // filter operator
  std::array<basis::OperatorBlocks<double>,3> target_component_blocks;
  FilterOperator(
      parameters,space,operator_parameters,component_sectors,
      source_component_blocks,target_component_blocks,
      false  // verbose
    );

  // write operator
  basis::WriteRelativeOperatorLSJT(
      parameters.target_filename,
      space,operator_parameters,
      component_sectors,target_component_blocks,
      true  // verbose
    );

  // termination
  return 0;
}
