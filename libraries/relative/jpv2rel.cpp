/****************************************************************
  jpv2rel.cpp

  Populate Hamiltonian-like operator from Iowa State ("JPV") format
  relative operator file.

  Initial support is for isoscalar operators, though the format has
  also been extended to non-isoscalar operators.

  See lsjt_operator.h for documentation of operator storage and the
  relative operator file format.

  Standard input:
    Nmax Jmax
    source_filename
    target_filename

  Language: C++11
                                 
  Mark A. Caprio
  University of Notre Dame

  7/25/16 (mac): Created, based upon writerel.cpp.

****************************************************************/

#include <fstream>

#include "basis/lsjt_operator.h"
#include "cppformat/format.h"
#include "mcpp/parsing.h"
#include "relative/construct_relative.h"

////////////////////////////////////////////////////////////////
// parameter input
////////////////////////////////////////////////////////////////

struct Parameters
// Container for run input parameters.
{
  int Nmax, Jmax;
  std::string source_filename;
  std::string target_filename;
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


  // line 1: truncation properties
  {
    ++line_count;
    std::getline(std::cin,line);
    std::istringstream line_stream(line);
    line_stream >> parameters.Nmax
                >> parameters.Jmax;
    ParsingCheck(line_stream,line_count,line);
  }

  // line 2: source filename
  {
    ++line_count;
    std::getline(std::cin,line);
    std::istringstream line_stream(line);
    line_stream >> parameters.source_filename;
    ParsingCheck(line_stream,line_count,line);
  }

  // line 3: target filename
  {
    ++line_count;
    std::getline(std::cin,line);
    std::istringstream line_stream(line);
    line_stream >> parameters.target_filename;
    ParsingCheck(line_stream,line_count,line);
  }

}

////////////////////////////////////////////////////////////////
// operator work
////////////////////////////////////////////////////////////////

void ReadJPVOperator(
    const std::string& source_filename,
    const basis::RelativeSpaceLSJT& relative_space,
    const basis::OperatorLabelsJT& operator_labels,
    const std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
    std::array<basis::MatrixVector,3>& relative_component_matrices
  )
// Define operator.
//
// Arguments:
//   source_filename (std::string) : input filename
//   relative_space (...) : target space
//   relative_component_sectors (...) : target sectors
//   relative_component_matrices (..., output) : target matrices to be populated
{

  // open stream for reading
  std::cout
    << "Reading relative operator file (JPV format)..." << std::endl
    << "  Filename: " << source_filename << std::endl;
  std::ifstream is(source_filename.c_str());
  OpenCheck(bool(is),source_filename);

  // set up references for convenience
  const basis::RelativeSectorsLSJT& sectors = relative_component_sectors[0];
  basis::MatrixVector& matrices = relative_component_matrices[0];

  // read source file
  std::string line;
  int line_count = 0;
  bool done = false;
  while (!done)
    {

      // parse sector header line
      //
      // J, S, L, Lp, Nmax, dimension, mn, hw, identifier
      int J, S, L, Lp, Nmax, dimension;
      double mn, hw;
      int identifier;
      {
        ++line_count;
        std::getline(is,line);
        std::istringstream line_stream(line);

        // check first for terminating line
        line_stream >> J;
        ParsingCheck(line_stream,line_count,line);
        if (J==99999)
          {
            done = true;
            std::cout << "  Reached termination marker." << std::endl;
            continue;
          }

        // read rest of header line
        line_stream >> S >> L >> Lp >> Nmax >> dimension >> mn >> hw >> identifier;
        ParsingCheck(line_stream,line_count,line);
      }
      std::cout << fmt::format("  Input sector (raw labels): J {} S {} L {} Lp {} ipcut {} dimension {} mn {} hw {} ident {}",J,S,L,Lp,Nmax,dimension,mn,hw,identifier)
                << std::endl;

      // deduce sector cutoffs
      int nmax = dimension-1;
      int Nmax_bra = 2*nmax+Lp;
      int Nmax_ket = 2*nmax+L;
      int expected_matrix_elements = dimension*(dimension+1)/2;  // JPV always stores just one triangle
      std::cout << fmt::format("  ==> expected {} nmax {} Nmax_bra {} Nmax_ket {}",expected_matrix_elements,nmax,Nmax_bra,Nmax_ket)
                << std::endl;
      assert(std::max(Nmax_bra,Nmax_ket)<=relative_space.Nmax());

      // deduce implied sectors labels
      int T = (L+S+1)%2; // isospin forced by L+S+T~1
      int Tp = (Lp+S+1)%2;
      assert(Tp==T);  // assumed delta T = 0 in JPV format
      int g = L%2;  // parity forced by L~g
      int gp = Lp%2;

      // look up corresponding sector in our internal representation
      int subspace_index_bra = relative_space.LookUpSubspaceIndex(
          basis::RelativeSubspaceLSJTLabels(Lp,S,J,Tp,gp)
        );
      int subspace_index_ket = relative_space.LookUpSubspaceIndex(
          basis::RelativeSubspaceLSJTLabels(L,S,J,T,g)
        );
      assert(subspace_index_bra<=subspace_index_ket);
      int sector_index = sectors.LookUpSectorIndex(subspace_index_bra,subspace_index_ket);
      const basis::RelativeSectorsLSJT::SectorType& sector = sectors.GetSector(sector_index);

      // print sector diagnostics
      std::cout
        << fmt::format("  Sector storage: sector_index {} bra_subspace {} ket_subspace {} diagonal {} dimensions {}x{}",
                       sector_index,
                       sector.bra_subspace().LabelStr(),
                       sector.ket_subspace().LabelStr(),
                       int(sector.IsDiagonal()),
                       sector.bra_subspace().size(),
                       sector.ket_subspace().size()
          )
        << std::endl;

      // read lines until satiated for matrix elements
      int matrix_element_count = 0;
      while (matrix_element_count < expected_matrix_elements)
        {

          // read one line of matrix elements
          ++line_count;
          std::getline(is,line);
          std::istringstream line_stream(line);
          assert(is);  // can fail if read past EOF

          std::cout << fmt::format("matrix_element_count {}",matrix_element_count) << std::endl;
          std::cout << line_count << " : " << line << std::endl;

          while (
              (matrix_element_count < expected_matrix_elements)
              && line_stream
            )
            {
              // attempt to extract matrix element
              double matrix_element;
              line_stream >> matrix_element;
              
              // save matrix element
              if (line_stream)
                {
                  // deduced matrix element indices
                  int np = matrix_element_count / nmax;  // bra index less rapidly varying
                  int n = matrix_element_count % nmax;  // ket index more rapidly varying

                  // save matrix element
                  matrices[sector_index](np,n) = matrix_element;

                  // advance counter
                  ++matrix_element_count;
                }
            }


//            for (int np=0; np<=nmax; np++)
//                for(int n=0; n<=np; n++)
//                {
//                  sector(np,n)=matrix_elements[pos];
//                  sector(n,np)=matrix_elements[pos];
//                  ++pos;
//                }
        }

    }
      
}

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  // parameter input
  Parameters parameters;
  ReadParameters(parameters);

  // set up zero operator
  std::cout << "Operator setup..." << std::endl;
  basis::RelativeSpaceLSJT relative_space(parameters.Nmax,parameters.Jmax);
  basis::OperatorLabelsJT operator_labels;
  operator_labels.J0 = 0;
  operator_labels.g0 = 0;
  operator_labels.T0_min = 0;
  operator_labels.T0_max = 0;
  operator_labels.symmetry_phase_mode = basis::SymmetryPhaseMode::kHermitian;
  std::array<basis::RelativeSectorsLSJT,3> relative_component_sectors;
  std::array<basis::MatrixVector,3> relative_component_matrices;
  relative::ConstructDiagonalConstantOperator(
      basis::RelativeOperatorParametersLSJT(operator_labels,parameters.Nmax,parameters.Jmax),
      relative_space,relative_component_sectors,relative_component_matrices,
      0.
    );

  // operator diagnostics
  std::cout << "  Truncation:"
            << " Nmax " << parameters.Nmax
            << " Jmax " << parameters.Jmax
            << std::endl;
  std::cout << "  Matrix elements:";
  for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
    std::cout << " " << basis::UpperTriangularEntries(relative_component_sectors[T0]);
  std::cout << std::endl;
  std::cout << "  Allocated:";
  for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
    std::cout << " " << basis::AllocatedEntries(relative_component_matrices[T0]);
  std::cout << std::endl;
        
  // populate matrix elements
  ReadJPVOperator(
      parameters.source_filename,
      relative_space,operator_labels,relative_component_sectors,relative_component_matrices
    );

  // write operator
  basis::WriteRelativeOperatorLSJT(
      parameters.target_filename,
      relative_space,
      operator_labels,  // only need operator labels
      relative_component_sectors,
      relative_component_matrices,
      true  // verbose
    );

  // termination
  return 0;
}
