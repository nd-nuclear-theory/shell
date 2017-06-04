/****************************************************************
  jpv_io.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "relative/jpv_io.h"

#include <fstream>

#include "basis/jjjpn_scheme.h"  // for TwoBodySpecies enum typedef
#include "cppformat/format.h"
#include "mcutils/eigen.h"  // for debugging output
#include "mcutils/parsing.h"

namespace relative {

  void ReadJPVOperator(
      const std::string& source_filename,
      const basis::RelativeSpaceLSJT& relative_space,
      const basis::OperatorLabelsJT& operator_labels,
      const std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::MatrixVector,3>& relative_component_matrices,
      bool verbose
    )
  {

    // validate operator labels

    // const basis::OperatorLabelsJT operator_labels_isoscalar(0,0,0,0,basis::SymmetryPhaseMode::kHermitian);
    // assert(operator_labels==operator_labels_isoscalar);  // no == defined on OperatorLabels
    assert(
        (operator_labels.J0==0) && (operator_labels.g0==0)
        && (operator_labels.T0_min==0) && (operator_labels.T0_max==0)
        && (operator_labels.symmetry_phase_mode==basis::SymmetryPhaseMode::kHermitian)
      );

    // open stream for reading
    std::cout
      << "Reading relative operator file (JPV format)..." << std::endl
      << "  Filename: " << source_filename << std::endl;
    std::ifstream is(source_filename);
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
        if (verbose)
          std::cout << fmt::format("  Input sector (raw labels): J {} S {} L {} Lp {} ipcut {} dimension {} mn {} hw {} ident {}",J,S,L,Lp,Nmax,dimension,mn,hw,identifier)
                    << std::endl;

        // canonicalize "lower triangle" sector
        //
        // We must flag the need to transpose (np,n) labels for
        // individual matrix elements below as well.
        bool flip = (Lp>L);
        if (flip)
          std::swap(Lp,L);
        if (verbose)
          std::cout << fmt::format("    After transposition to upper triangle: (Lp,L)=({},{}) flip {}",Lp,L,int(flip))
                    << std::endl;

        // deduce sector cutoffs
        int nmax = dimension-1;
        int Nmax_bra = 2*nmax+Lp;
        int Nmax_ket = 2*nmax+L;
        int expected_matrix_elements = dimension*(dimension+1)/2;  // JPV always stores just one triangle
        if (verbose)
          std::cout << fmt::format("    Input sector properties: expected m.e. {} nmax {} Nmax_bra {} Nmax_ket {}",expected_matrix_elements,nmax,Nmax_bra,Nmax_ket)
                    << std::endl;

        // check viability of sector
        //
        // * make sure target space includes the given J for this sector
        //   -- otherwise we have to be more careful with sector lookups
        //
        // * make sure target sector can hold all matrix elements --
        //   OMIT since we are okay with having a smaller target space
        //   than source space
        assert(J<=relative_space.Jmax());
        //    assert(std::max(Nmax_bra,Nmax_ket)<=relative_space.Nmax());


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
        // short circuit if subspace falls outside our target truncation
        if ((subspace_index_bra==basis::kNone)||(subspace_index_ket==basis::kNone))
          {
            std::cout << "ERROR: Input sector contains LSJT subspace not present in target truncation" << std::endl;
            std::exit(EXIT_FAILURE);
          }
        assert(subspace_index_bra<=subspace_index_ket);  // subspaces should be canonical after our L swap
        int sector_index = sectors.LookUpSectorIndex(subspace_index_bra,subspace_index_ket);
        const basis::RelativeSectorsLSJT::SectorType& sector = sectors.GetSector(sector_index);

        // look up target matrix dimensions
        int dimension_bra = sector.bra_subspace().size();
        int dimension_ket = sector.ket_subspace().size();
        // print sector diagnostics
        int sector_size;
        if (sector.IsDiagonal())
          sector_size = dimension_ket * (dimension_ket + 1);
        else
          sector_size = dimension_bra * dimension_ket;
        if (verbose)
          {
            std::cout
              << fmt::format("    Subspace labels: bra {} (dimension {}) ket {}  (dimension {})",
                             sector.bra_subspace().LabelStr(),
                             dimension_bra,
                             sector.ket_subspace().LabelStr(),
                             dimension_ket
                )
              << std::endl
              << fmt::format("    Sector storage: sector_index {} subspace_index_bra {} subspace_index_ket {} entries {}",
                             sector_index,
                             subspace_index_bra,
                             subspace_index_ket,
                             sector_size
                )
              << std::endl;
          }

        // reading matrix elements for sector:
        //
        //   repeat:
        //     read line
        //     repeat:
        //       try to extract matrix element from line
        //     until (read enough matrix elements) or (fail due to end of line)
        //   until (read enough matrix elements)

        int matrix_element_count = 0;
        int row_index = 0;
        int column_index = 0;
        bool done = (matrix_element_count == expected_matrix_elements);
        while (!done)
          {
            // read one line
            ++line_count;
            std::getline(is,line);
            std::istringstream line_stream(line);
            StreamCheck(bool(is),source_filename,"Failure reading matrix elements");  // can fail if there are not enough matrix elements and we read past EOF
            // std::cout << fmt::format("matrix_element_count {}",matrix_element_count) << std::endl;
            // std::cout << line_count << " : " << line << std::endl;

            // extract matrix elements from line
            while (!done)
              {
                // attempt to extract matrix element
                double matrix_element;
                line_stream >> matrix_element;

                // quit if end of line
                if (!line_stream)
                  break;

                // deduce matrix element indices
                //
                // Recall we have adopted the interpretation that matrix
                // elements (np,n) are column major in upper triangle.
                int np = row_index;
                int n = column_index;
                if (flip)
                  std::swap(np,n);
              
                // save matrix element
                if ((np<dimension_bra)&&(n<dimension_ket))
                  matrices[sector_index](np,n) = matrix_element;

                // advance counter and check for completion
                ++matrix_element_count;
                done = (matrix_element_count == expected_matrix_elements);

                // advance to next matrix element (row,column) indices
                //
                // column major ordering in upper triangle
                ++row_index;
                if (row_index>column_index)
                  {
                    row_index = 0;
                    ++column_index;
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

  void ReadJPVOperatorPN(
      const std::array<std::string,3>& source_filenames,
      const basis::RelativeSpaceLSJT& relative_space,
      const basis::RelativeOperatorParametersLSJT& operator_parameters,
      const std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::MatrixVector,3>& relative_component_matrices,
      bool verbose
    )
  {

    // validate operator labels
    assert(
        (operator_parameters.J0==0) && (operator_parameters.g0==0)
        && (operator_parameters.T0_min==0) && (operator_parameters.T0_max==2)
        && (operator_parameters.symmetry_phase_mode==basis::SymmetryPhaseMode::kHermitian)
      );

    // Relation of JT-reduced matrix elements to pp/nn/pn matrix elements
    //
    //   < T || A^{T0} || T >  vs.  < TTz | A | TTz>
    //
    // For T=0 sectors:
    //
    //   <0||A0||0> = <00|A|00>
    //
    // For T=1 sectors:
    // 
    // {<1||A0||1>, <1||A1||1>, <1||A2||1>}
    // = 1/3. * {
    //           {1,1,1},
    //           {sqrt(9./2.),-sqrt(9./2.),0},
    //           {sqrt(5./2.),sqrt(5./2.),-sqrt(10.)}
    //         }
    //   * {<1+1|A|1+1>, <1-1|A|1-1>, <10|A|10>}
    //
    // We have listed Tz sectors in the order pp/nn/pn to match
    // the TwoBodySpecies enum ordering.
    
    // transformation matrix
    //
    // indexed by (T,int(two_body_species))
    static Eigen::Matrix3d kIsospinCoefficientMatrixTzToTForT1;
    kIsospinCoefficientMatrixTzToTForT1
      << 1, 1, 1,
      std::sqrt(9./2.), -std::sqrt(9./2.), 0,
      std::sqrt(5./2.), std::sqrt(5./2.), -std::sqrt(10.);
    kIsospinCoefficientMatrixTzToTForT1 *= 1/3.;

    // read matrix elements by two-body species
    for (basis::TwoBodySpeciesPN two_body_species : {basis::TwoBodySpeciesPN::kPP,basis::TwoBodySpeciesPN::kNN,basis::TwoBodySpeciesPN::kPN})
      {

        // convert two body species to Tz
        int Tz = basis::kTwoBodySpeciesPNCodeTz[int(two_body_species)];

        // set up storage for input matrix elements
        //
        // Recall that they are stored in the JPV files as if they
        // were "isoscalar" MEs, so only T0=0 component of each will
        // be populated

        std::array<basis::RelativeSectorsLSJT,3> relative_component_sectors_input;
        std::array<basis::MatrixVector,3> relative_component_matrices_input;
        basis::ConstructZeroOperatorRelativeLSJT(
            basis::RelativeOperatorParametersLSJT(operator_parameters,operator_parameters.Nmax,operator_parameters.Jmax),
            relative_space,relative_component_sectors_input,relative_component_matrices_input
          );

        // read matrix elements
        const std::string& source_filename = source_filenames[int(two_body_species)];
        const basis::OperatorLabelsJT operator_labels_isoscalar(0,0,0,0,basis::SymmetryPhaseMode::kHermitian);
        ReadJPVOperator(
            source_filename,
            relative_space,
            operator_labels_isoscalar,
            relative_component_sectors_input,
            relative_component_matrices_input,
            verbose
          );

        // accumulate matrix elements
        const int num_sectors = relative_component_sectors_input[0].size();
        for (int sector_index=0; sector_index < num_sectors; ++sector_index)
          // for each source "isoscalar operator" sector
          {
            // set up aliases for convenience
            const typename basis::RelativeSectorsLSJT::SectorType& input_sector
              = relative_component_sectors_input[0].GetSector(sector_index);
            const Eigen::MatrixXd& input_matrix
              = relative_component_matrices_input[0][sector_index];
          
            // extract sector isospin labels
            int bra_T = input_sector.bra_subspace().T();
            int ket_T = input_sector.ket_subspace().T();
            if (bra_T!=ket_T)
              continue;  // short circuit known vanishing T-changing sectors
            int T = ket_T;
          
            if (T==0)
              // sector with (T'T)=(0,0)
              {
                // T=0 sectors only relevant in Tz=0 file; they are zeroed out in Tz!=0 files
                if(Tz==0)
                  {
                    // Simple copy to corresponding target sector.  Though we
                    // can write this as an accumulation for consistency.
                    relative_component_matrices[0][sector_index] += input_matrix;
                  }
              }
            else if (T==1)
              // sector with (T'T)=(1,1)
              {
                for (int T0=0; T0<=2; ++T0)
                  {
                    // look up target sector
                    int target_sector_index = relative_component_sectors[T0].LookUpSectorIndex(input_sector.bra_subspace_index(),input_sector.ket_subspace_index());
                    assert(sector_index!=basis::kNone);
                    const typename basis::RelativeSectorsLSJT::SectorType& target_sector
                      = relative_component_sectors[T0].GetSector(target_sector_index);

                    // accumulate matrix for sector
                    Eigen::MatrixXd& matrix
                      = relative_component_matrices[T0][target_sector_index];
                    double isospin_coefficient = kIsospinCoefficientMatrixTzToTForT1(T0,int(two_body_species));
                    matrix += isospin_coefficient * input_matrix;

                    // debugging output
                    // std::cout
                    //   << fmt::format(
                    //       "input sector Tz={:+d} {:s}x{:s} -> target sector T0={:d} {:s}x{:s}   coefficient {:.4f}",
                    //       Tz,
                    //       input_sector.bra_subspace().LabelStr(),input_sector.ket_subspace().LabelStr(),
                    //       T0,
                    //       target_sector.bra_subspace().LabelStr(),target_sector.ket_subspace().LabelStr(),
                    //       isospin_coefficient
                    //     )
                    //   << std::endl;
                    // std::cout
                    //   << mcutils::FormatMatrix(input_matrix,".7f","   ")
                    //   << std::endl;
                    // std::cout << "  ... accumulating target ... " << std::endl;
                    // std::cout
                    //   << mcutils::FormatMatrix(matrix,".7f","   ")
                    //   << std::endl;




                  }

          
              }
          }
      }
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
