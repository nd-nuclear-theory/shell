/****************************************************************
  jpv_io.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "relative/jpv_io.h"

#include <fstream>

#include "cppformat/format.h"
#include "mcutils/parsing.h"

namespace relative {

  void ReadJPVOperator(
      const std::string& source_filename,
      const basis::RelativeSpaceLSJT& relative_space,
      const basis::OperatorLabelsJT& operator_labels,
      const std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::MatrixVector,3>& relative_component_matrices
    )
  {

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
        std::cout << fmt::format("  Input sector (raw labels): J {} S {} L {} Lp {} ipcut {} dimension {} mn {} hw {} ident {}",J,S,L,Lp,Nmax,dimension,mn,hw,identifier)
                  << std::endl;

        // canonicalize "lower triangle" sector
        //
        // We must flag the need to transpose (np,n) labels for
        // individual matrix elements below as well.
        bool flip = (Lp>L);
        if (flip)
          std::swap(Lp,L);
        std::cout << fmt::format("    After transposition to upper triangle: (Lp,L)=({},{}) flip {}",Lp,L,int(flip))
                  << std::endl;

        // deduce sector cutoffs
        int nmax = dimension-1;
        int Nmax_bra = 2*nmax+Lp;
        int Nmax_ket = 2*nmax+L;
        int expected_matrix_elements = dimension*(dimension+1)/2;  // JPV always stores just one triangle
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

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
