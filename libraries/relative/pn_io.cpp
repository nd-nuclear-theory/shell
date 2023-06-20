/****************************************************************
  pn_io.cpp

  Anna E. McCoy
  TRIUMF

****************************************************************/

#include "relative/pn_io.h"

#include <fstream>

#include "basis/jjjpn_scheme.h"  // for TwoBodySpecies enum typedef
#include "fmt/format.h"
#include "mcutils/eigen.h"  // for debugging output
#include "mcutils/parsing.h"

namespace relative
{


  void ReadPNOperator(
      const std::string& source_filename,
      const basis::RelativeSpaceLSJT& relative_space,
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSectorsLSJT& sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
      bool verbose
    )
  {
    // validate operator labels
    assert(
        (operator_labels.J0==0) && (operator_labels.g0==0)
        && (operator_labels.symmetry_phase_mode==basis::SymmetryPhaseMode::kHermitian)
      );

    // open stream for reading
    std::cout
      << "Reading relative operator file (Petr Navratil format)..." << std::endl
      << "  Filename: " << source_filename << std::endl;
    std::ifstream is(source_filename);
    mcutils::StreamCheck(bool(is),source_filename,"Failed to open relative operator file");

    // read source file
    std::string line;
    int line_count = 0;
    int np,Lp,n,L,S,Tz,J;
    double me;
    // int identifier;

    while (std::getline(is,line))
      {
        ++line_count;
        std::istringstream line_stream(line);
        
        // read rest of header line
        line_stream >> np >> Lp >> n >> L >> S >> J >> Tz >> me;
        mcutils::ParsingCheck(line_stream,line_count,line);

        // deduce implied sectors labels
        int T = (L+S+1)%2; // isospin forced by L+S+T~1
        int Tp = (Lp+S+1)%2;
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

        // set up references for convenience
        basis::OperatorBlocks<double>& matrices = relative_component_matrices[Tz+1];


        int sector_index;
        // swap bra and ket in cases where the supplied matrix element is the adjoint 
        if(subspace_index_bra>subspace_index_ket)
          {

            sector_index = sectors.LookUpSectorIndex(subspace_index_ket,subspace_index_bra);
            std::swap(np,n);
          }
        else
          sector_index = sectors.LookUpSectorIndex(subspace_index_bra,subspace_index_ket);

      
        const basis::RelativeSectorsLSJT::SectorType& sector = sectors.GetSector(sector_index);

        // look up target matrix dimensions
        int dimension_bra = sector.bra_subspace().size();
        int dimension_ket = sector.ket_subspace().size();

        // save matrix element (only keep if within chosen Nmax truncation?)
        if ((np<dimension_bra)&&(n<dimension_ket))
        {
          matrices[sector_index](np,n) = me;
          // If diagonal subspace
          if(subspace_index_bra==subspace_index_ket)
            matrices[sector_index](n,np) = me; 
        } 
      }
  }
  ////////////////////////////////////////////////////////////////

  void ReadPNOperatorPN(
        const std::string& source_filename,
        const basis::RelativeSpaceLSJT& relative_space,
        const basis::RelativeOperatorParametersLSJT& operator_parameters,
        const std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
        std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
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
    // We have listed Tz sectors in the order pp/pn/nn to match
    // Tz=-1,0,1 ordering.
    
    // transformation matrix
    //
    // indexed by (T,Tz)
    static Eigen::Matrix3d kIsospinCoefficientMatrixTzToTForT1;
    kIsospinCoefficientMatrixTzToTForT1
      << 1, 1, 1,
      std::sqrt(9./2.),  0, -std::sqrt(9./2.),
      std::sqrt(5./2.), -std::sqrt(10.), std::sqrt(5./2.) ;
    kIsospinCoefficientMatrixTzToTForT1 *= 1/3.;

    // set up storage for input matrix elements
    std::array<basis::RelativeSectorsLSJT,3> relative_component_sectors_input;
    std::array<basis::OperatorBlocks<double>,3> relative_component_matrices_input;
    basis::ConstructZeroOperatorRelativeLSJT(
        basis::RelativeOperatorParametersLSJT(operator_parameters,operator_parameters.Nmax,operator_parameters.Jmax),
        relative_space,relative_component_sectors_input,relative_component_matrices_input
      );

    //Hack to get correct zero initialized matrices and sectors 
    // PN matrix elements all treated as T0=0 sectors 
    const basis::RelativeSectorsLSJT& input_sectors = relative_component_sectors[0];
    relative_component_matrices_input[1]=relative_component_matrices_input[0];
    relative_component_matrices_input[2]=relative_component_matrices_input[0];
  
    //For proton-neutron RMEs 'pretend' operator is isoscalar
    //The three components of relative_component_matrices_input will be 0: Tz=-1, 1: Tz=0 and 2: Tz=1
    const basis::OperatorLabelsJT operator_labels_isoscalar(0,0,0,0,basis::SymmetryPhaseMode::kHermitian);
    relative::ReadPNOperator(
        source_filename,
        relative_space,
        operator_labels_isoscalar,
        input_sectors,
        relative_component_matrices_input,
        verbose
      );

    // accumulate matrix elements
    const int num_sectors = relative_component_sectors_input[0].size();
    for (int sector_index=0; sector_index < num_sectors; ++sector_index)
      {
        // set up aliases for convenience
        const typename basis::RelativeSectorsLSJT::SectorType& 
          input_sector = input_sectors.GetSector(sector_index);

        // extract sector isospin labels
        int bra_T = input_sector.bra_subspace().T();
        int ket_T = input_sector.ket_subspace().T();
        if (bra_T!=ket_T)
          continue;  // short circuit known vanishing T-changing sectors
        int T = ket_T;

        for(int Tz=-1; Tz<=1; ++Tz)
          {
            const Eigen::MatrixXd& 
              input_matrix = relative_component_matrices_input[Tz+1][sector_index];

            if (T==0)
              // sector with (T'T)=(0,0)
              {
                // T=0 sectors only relevant in Tz=0 file; they are zeroed out in Tz!=0 files
                if(Tz==0)
                  {
                    // Simple copy to corresponding target sector.  Though we
                    // can write this as an accumulation for consistency.
                    int T0=0;
                    relative_component_matrices[T0][sector_index] += input_matrix;
                  }
              }
            else if (T==1)
              // sector with (T'T)=(1,1)
              {
                for (int T0=0; T0<=2; ++T0)
                  {
                    // look up target sector
                    int target_sector_index 
                      = relative_component_sectors[T0].LookUpSectorIndex(input_sector.bra_subspace_index(),input_sector.ket_subspace_index());
                    assert(sector_index!=basis::kNone);
                    const typename basis::RelativeSectorsLSJT::SectorType& target_sector
                      = relative_component_sectors[T0].GetSector(target_sector_index);

                    // accumulate matrix for sector
                    Eigen::MatrixXd& matrix
                      = relative_component_matrices[T0][target_sector_index];
                    double isospin_coefficient = kIsospinCoefficientMatrixTzToTForT1(T0,Tz+1);
                    matrix += isospin_coefficient * input_matrix;
                    
                  }
              }
          }
      }
  }

  ////////////////////////////////////////////////////////////////
} // namespace
