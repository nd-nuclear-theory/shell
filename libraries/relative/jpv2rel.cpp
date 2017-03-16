/****************************************************************
  jpv2rel.cpp

  Populate Hamiltonian-like operator from Iowa State ("JPV") format
  relative operator file.

  See lsjt_operator.h for documentation of operator storage and the
  relative operator file format.

  Notes on the interpretation of JPV relative files are provided in
  comments below.

  The basic format as considered here provides support only for
  isoscalar operators, though the format has also been extended to
  non-isoscalar operators through the use of separate "pp", "nn", and
  "pn" files.

  Standard input:
    Nmax Jmax
    source_filename
    target_filename

  Language: C++11
                                 
  Mark A. Caprio
  University of Notre Dame

  7/25/16 (mac): Created, based upon writerel.cpp.
  8/10/16 (mac): Fix input indexing.
  10/9/16 (pjf): Rename mcpp -> mcutils.
  10/19/16 (mac): Remove superflous debugging options.

****************************************************************/

// JPV relative interaction readme file from NERSC m94 project space
// (2014):
//
// c**    For the relative HO matrix element files stored in
// c**    this directory, the following parameters and 
// c**    conventions are employed.
// 
//       idmaxp = 120
// 
// c**   For uncoupled channels and subchannels of coupled channels:
//       ipsiz = (2*idmaxp-min0(LR,LL))/2
// 
// c**   Thus, the full coupled-channel matrix is dimensioned at (2*ipsiz) X (2*ipsiz)
// 
// c**   Conventions for storing the coupled-channels matrices
// c**   and writing them out to the file.
// c**
// c**   Commands to write the file and their formats are provided
// c**   below.  Each square of the coupled channels matrices and 
// c**   the uncoupled channel matrices is stored in upper-triangle
// c**   form.  The ordering of the orbital angular momenta in the
// c**   bras and the kets is shown in the following figure.
// c**
// c**     -------------------------------
// c**     |               |             |
// c**     |   IPART= 1    |      2      |
// c**     |     {L,L}     |   {L,L+2}   |
// c**     -------------------------------
// c**     |               |             |
// c**     |       3       |      4      |
// c**     |   {L+2,L}     |  {L+2,L+2}  |
// c**     -------------------------------
// c**
// c**     Within a square subchannel, the ordering of the stored
// c**     matrix elements is:
// c**
// C**         1 -> 2   4   7  11
// c**              3   5   8  etc
// c**                  6   9
// c**                     10
// c**
// c**   OMEGF = nucleon mass in MeV
// c**   HOP   = HO energy
// c**   ipcut = ignorable
// c**   IDENT = ignorable
// c**
// c**   Stop reading when you encounter JT = 99999 on input.
// c**
// c**   Output the IPART = 1 subchannel                                                                                 
// c**   Note that ordering does not matter as it is symmetric                                                           
//       WRITE(12,1250) JT,IS,LR-2,LL-2,ipcut, ipsiz, OMEGF,HOP,IDENT
//  1250 FORMAT(6I5,2(1PE16.6),I5)
//       write(12,1300)((VECP(i,j),i=1,j),j=1,ipsiz)
//  1300  format(7(1PE15.6))
// 
// c**   Note that if IS = 1 and JT = "LR-2" + 1, this is a coupled-channel.
// c**   In this case all 4 "IPART" subsections follow in the file in the
// c**   sequence 1, 2, 3 & 4.
// 
//       ipsizc = 2*ipsiz - 1
// 
// c**   Output the IPART = 2 subchannel (as a square submatrix with                                                     
// c**   zeroes padding it)                                                                                              
//       WRITE(12,1250) JT,IS,LR-2,LL,ipcut, ipsiz, OMEGF,HOP,IDENT
//       write(12,1300)((VECP(i,j),j=ipsiz+1,i+ipsiz),i=1,ipsiz)
// 
// c**   Output the IPART = 3 subchannel (as a square submatrix with                                                     
// c**   zeroes padding it)                                                                                              
//       WRITE(12,1250) JT,IS,LR,LL-2,ipcut, ipsiz, OMEGF,HOP,IDENT
//       write(12,1300)((VECP(i,j),j=1,i-ipsiz),i=ipsiz+1,ipsizc+1)
// 
// c**   Output the IPART = 4 subchannel (as a square submatrix with                                                     
// c**   zeroes padding it)                                                                                              
// c**   Note that ordering does not matter as it is symmetric                                                           
//       WRITE(12,1250) JT,IS,LR,LL,ipcut, ipsiz, OMEGF,HOP,IDENT
//       write(12,1300)((VECP(i,j),i=ipsiz+1,j),j=ipsiz+1,ipsizc+1)
// 
// c**   Uncoupled channels are stored similar to IPART = 1 of the
// c**   coupled channels.
//       WRITE(12,12) JT,IS,LR,LL,ipcut, ipsiz, OMEGF,HOP,IDENT
//       write(12,1300)((VECP(i,j),i=1,j),j=1,ipsiz)
// 
// * * * * *  END OF THE README FILE * * * * * 

// JPV relative interaction readme file README_Heff_notes_v2.txt
// (2015):
//
// c**    Important updates added November 27, 2015 for Version 2
// c**    
// 
// c**    For the relative HO matrix element files stored in
// c**    this directory, the following parameters and 
// c**    conventions are employed.
// 
//       idmaxp = 120
// 
// c**   For uncoupled channels and subchannels of coupled channels:
//       ipsiz = (2*idmaxp-min0(LR,LL))/2
// 
// c**   Thus, the full coupled-channel matrix is dimensioned at (2*ipsiz) X (2*ipsiz)
// 
// c**   Conventions for storing the coupled-channels matrices
// c**   and writing them out to the file.
// c**
// c**   Commands to write the file and their formats are provided
// c**   below.  Each square of the coupled channels matrices and 
// c**   the uncoupled channel matrices is stored in a triangle
// c**   form (be careful about the row-col ordering - see below). 
// c** 
// c**   The ordering of the orbital angular momenta in the
// c**   bras and the kets is shown in the following figure.
// c**
// c**     -------------------------------
// c**     |               |             |
// c**     |   IPART= 1    |      2      |
// c**     |     {L,L}     |   {L,L+2}   |
// c**     -------------------------------
// c**     |               |             |
// c**     |       3       |      4      |
// c**     |   {L+2,L}     |  {L+2,L+2}  |
// c**     -------------------------------
// c**
// c**     For IPART = 1 and 4 
// c**     Within these square subchannels, the ordering of the stored
// c**     matrix elements is:
// c**
// c**         1 -> 2   4   7  11
// c**              3   5   8  etc
// c**                  6   9
// c**                     10
// c**
// c**	However, one may transpose each of these sub channels 
// c**     due to symmetry - i.e. interchange the rows and columns.
// c**
// c**     For the "off-diagonal subchannels" i.e. for IPART = 2 & 3
// c**     the row-col ordering is as follows (see the write statements
// c**     below for explicit indexing):
// c**
// c**         1
// c**         2    3
// c**         4    5   6
// c**         7    8   9  10
// c**         etc
// c**
// c**    From this stored triangle, the symmetric part of the matrix
// c**    in the complementary IPART area is obtained.
// c**
// c**   Variables defined:
// c**   OMEGF = nucleon mass in MeV
// c**   HOP   = HO energy
// c**   ipcut = ignorable
// c**   IDENT = ignorable
// c**
// c**   Stop reading when you encounter JT = 99999 on input.
// c**
// c**   Output the IPART = 1 subchannel                                                                                 
// c**   Note that ordering does not matter as it is symmetric                                                           
//       WRITE(12,1250) JT,IS,LR-2,LL-2,ipcut, ipsiz, OMEGF,HOP,IDENT
//  1250 FORMAT(6I5,2(1PE16.6),I5)
//       write(12,1300)((VECP(i,j),i=1,j),j=1,ipsiz)
//  1300  format(7(1PE15.6))
// 
// c**   Note that if IS = 1 and JT = "LR-2" + 1, this is a coupled-channel.
// c**   In this case all 4 "IPART" subsections follow in the file in the
// c**   sequence 1, 2, 3 & 4.
// 
//       ipsizc = 2*ipsiz - 1
// 
// c**   Output the IPART = 2 subchannel (as a square submatrix with                                                     
// c**   zeroes padding it)                                                                                              
//       WRITE(12,1250) JT,IS,LR-2,LL,ipcut, ipsiz, OMEGF,HOP,IDENT
//       write(12,1300)((VECP(i,j),j=ipsiz+1,i+ipsiz),i=1,ipsiz)
// 
// c**   Output the IPART = 3 subchannel (as a square submatrix with                                                     
// c**   zeroes padding it)                                                                                              
//       WRITE(12,1250) JT,IS,LR,LL-2,ipcut, ipsiz, OMEGF,HOP,IDENT
//       write(12,1300)((VECP(i,j),j=1,i-ipsiz),i=ipsiz+1,ipsizc+1)
// 
// c**   Output the IPART = 4 subchannel (as a square submatrix with                                                     
// c**   zeroes padding it)                                                                                              
// c**   Note that ordering does not matter as it is symmetric                                                           
//       WRITE(12,1250) JT,IS,LR,LL,ipcut, ipsiz, OMEGF,HOP,IDENT
//       write(12,1300)((VECP(i,j),i=ipsiz+1,j),j=ipsiz+1,ipsizc+1)
// 
// c**   Uncoupled channels are stored similar to IPART = 1 of the
// c**   coupled channels.
//       WRITE(12,12) JT,IS,LR,LL,ipcut, ipsiz, OMEGF,HOP,IDENT
//       write(12,1300)((VECP(i,j),i=1,j),j=1,ipsiz)
// 
// * * * * *  END OF THE README FILE * * * * * 


// Further comments in convert_MFDn_to_Hagen.f from NERSC m94 project
// space (2014):
//
// c**    Begin processing the first channel, or a new quadrant of a coupled channel
// c**    or a new channel depending on values of ICOUP and IPART
// c**       JT  -   total angular momentum = J
// c**       IS  -   total spin
// c**       LR  -   orbital angular momentum of the bra
// c**       LL  -   orbital angular momentum of the ket
// c**       ipcut - reserved for model space cutoff signifier
// c**       idmax - dimension of input (bare) matrix assumed symmetric
// c**               so number of matrix elements: kxmax = idmax*(idmax+1)/2
// c**               For uncoupled channels, idmax is full dimension while
// c**               for coupled channels, it is the dimension of the subchannel
// c**               assumed to be a square subchannel.
// c**       OMEGF - Fermion mass (MeV) - Formerly, the starting energy
// c**       HO    - Harmonic oscillator energy parameter = hbar*omega (MeV)
// c**       IDENT - Interaction identification number
// c**

// Notes on JPV relative file format:
//
// * Note the apparent counterintuitive definition of LR (="right"?)
//   as "bra" and LL (="left"?) as ket, in the comment in
//   convert_MFDn_to_Hagen.f.  This assignment is substantiated by the
//   README.  Consider, e.g., sector 2, which apparently has {bra,ket}
//   orbital angular momenta {L,L+2}, but these are written in the
//   header line as (LR-2,LL), suggesting that the entries in the
//   header are indeed written as {Lp,L} but confusingly denoted as
//   {"right","left"}.  However, in handwritten notes from discussion
//   at Iowa State 12/26/14, the contradictory note is made that the
//   header line contains "..., LR=l, LL=l', ...".  That would attach
//   the natural meanings to "left" and "right".  We shall explore
//   this issue further below.
//
// * All diagonal sectors, i.e., for which Lp=L, are stored in an
//   identical scheme.  For purposes of reading a sector, we don't have
//   to worry about whether a *diagonal* sector is part of a "coupled
//   channel" or not, provided we take the "dimension"=ipsiz=idmax from
//   the header line.
//
//   The only way in which the nature of the subchannel is of concern
//   apparently is if we tried to calculate this dimension directly from
//   the (Lp,L) values and the two-body cutoff, since the dimension
//   *might* be calculated based on the "LL" and "LR" values used to
//   label the four coupled channels overall (this must be verified).
//
//   Why define ipsizc=2*ipsizp-1 then use ipsizc+1 as the dimension?
//   Why not just define ipsizc=2*ipsizp?  That *is* the size of the
//   coupled channel matrix, as noted earlier in the README, and
//   presumably "siz" is meant to indicate size...
//
//   With the original 1-based FORTRAN indexing:
//
//     for each j in [1,dimension]:
//       for each i in [1,j]:
//         A(i,j)
//
//   Translating to 0-based numbering based on the radial quantum
//   numbers (i->np and j->n):
//
//     for each n in [0,dimension):
//       for each np in [0,n):
//         A(np,n)
//
//   This iteration is column major (vertical stripes) in the upper
//   triangle of the sector, as indicated by the diagram in the README
//   (pay attention to the sequence of numbers, not the arrow, which
//   could give a false impression as to the direction of traversal).
//
//   When comparing with the off-diagonal sectors in the discussion
//   below, it is worth noting that the diagonal sectors are
//   symmetric, so we could just as easily interpret the matrix
//   elements as being row major (horizontal stripes) in the lower
//   triangle of the sector.
//
//   I am not quite sure of the implications of the claim that the
//   upper diagonal subchannel is "as a square submatrix with zeroes
//   padding it".  Of course it's square, and why are there zeros
//   padding it?  At higher L, with the same two-body cutoff, radial
//   quantum numbers stop one step lower, so maybe that's the issue.
//
// * For the off-diagonal sectors, we have:
//
//   Lp<L (IPART=2): This sector is in the upper triangle in JPV's
//   readme and also in our canonical sense of (bra subspace) < (ket
//   subspace).
//
//   With the original 1-based FORTRAN indexing, but subtracting off
//   the dimension as needed to restart indexing from 1 in this
//   "subchannel":
//
//     for each i in [1,dimension]:
//       for each j in [1,i]:
//         A(i,j)
//
//   Translating to 0-based numbering based on the radial quantum
//   numbers (i->np and j->n):
//
//     for each np in [0,dimension):
//       for each n in [0,np):
//         A(np,n)
//
//   This iteration is row major (horizontal stripes) in the lower
//   triangle of the sector, in apparent contradiction to the diagram
//   in the README.
//
//   Lp>L (IPART=3): This sector is in the lower triangle in JPV's
//   readme and also in our canonical sense of (bra subspace) > (ket
//   subspace).
//
//   The iteration is exactly as for the Lp<L sector, i.e., row major
//   (horizontal stripes) in the lower triangle of the sector, in
//   apparent contradiction to the diagram in the README.
//
// * Relation of the off-diagonal ambiguities to symmetry:
//
//   We cannot simply flip (np,n) within an
//   L-changing sector.  The individual sector matrix is not
//   symmetric.  So, no, we cannot trivially flip each off-diagonal sector to be
//   column major in the upper triangle.
//
//   We *can* flip (np,n) and (Lp,L) simultaneously, since the overall
//   Hamiltonian matrix is symmetric (and thus the 2x2 block "coupled
//   channel" matrix sketched in the README is also symmetric).
//
//   So, we are free to reverse the meaning of the L values in the
//   header line from the meaning in the readmes
//
//     ... Lp L ...
//
//   to match the *handwritten* notes
//
//     ... L Lp ...
//
//   and simultaneously reinterpret the matrix elements in these
//   sectors as row-major upper triangle (thus contradicting the WRITE
//   statements in the README but matching the diagram in the README).
//
//   Let's give this a 50% chance of being right...
//
// * The code convert_MFDn_to_Hagen.f reads the header line (see READ
//   FORMAT 12) as
//
//     ... LR LP ...
//
//   identifies these (according to the "backwards" naming convention)
//   as
//
//     ... Lp L ...
//
//   then labels the matrix elements for *all* sectors on input as
//
//     for each np in [0,dimension):
//       for each n in [0,np):
//         A(np,n)
//
//   i.e., row major, lower triangle.  The output (see WRITE FORMAT
//   1600), indeed, seems to pair the label np with Lp and n with L.
//   Both symmetric partner matrix elements are written out.
//
//   If we believe convert_MFDn_to_Hagen.f, there are thus two
//   equivalent interpretations, by overall symmetry of the
//   Hamiltonian:
//
//   (A) Header gives (Lp,L).  Matrix elements are stored row major,
//   lower triangle for *all* sectors (for each np: for each n<np).  This
//   is the nominal reading of the code.
//
//   (B) Header gives (L,Lp).  Matrix elements are stored column
//   major, upper triangle for *all* sectors (for each np: for each
//   np<n).  This is the transpose of the nominal reading of the code.
//
//   We observe that the written code entries in the README file are
//   completely consistent with (A), once we take into account the
//   symmetry of the diagonal sectors.  It is only the diagram on the
//   ordering of the stored matrix elements which contradicts this
//   interpretation, i.e., is inconsistent with either the iteration
//   order or the L labeling convention (take your pick) in the
//   off-diagonal channels!
//
// * Our choice of interpretation: Interpretation (B) matches our
//   preference for "upper triangle" storage.  For the Lp>L sector, we
//   will canonicalize the sector labels by interchanging (Lp,L) and
//   (np,n) simultaneously.  We will thus use these entries to fill in
//   the *lower* triangle of the *upper* triangle (Lp<L) sector.  The
//   diagonal matrix elements should be redundant across these sectors
//   (we could try to check for equality, but this assumes we know the
//   order in which the Lp>L an dLp<L sectors are read from the input
//   file).
//
// * The calculation
//
//     ipsiz = (2*idmaxp-min0(LR,LL))/2
//
//   indicates that idmaxp (=120 for the example files) is almost but
//   not quite the relative oscillator *radial* quantum number cutoff.
//   Since 2*idmaxp is guaranteed even and the numerator is
//   (presumably) guaranteed nonnegative, floor division gives
//
//      dimension = idmaxp - [ min(Lp,L) / 2 with upward rounding]
//
//   If, indeed, the same dimension is used for all subchannels, and
//   (LL,LR) indicate the angular momenta for the higher-L diagonal
//   subchannel (IPART=4), as suggested by the readme, then LR=LP, so
//   it is not clear why a "min" is needed in this formula.
//
//   Taking idmaxp=120 as claimed means the highest nmax for a sector
//   is 119, for the L=0 sector, and the N for this sector is thus
//   Nmax=2*nmax=238.  Naive interpretation of the subchannel matrices
//   as "square" causes some to be interpreted as going up to
//   Nmax=239, but this is presumably just counting the last row or
//   column of "zero padding".

// Example:
//
//   For a full listing of the sectors in a JISP16 JPV input file and
//   the corresponding output sectors, see the results of this test
//   run:
//
//     % cd libraries/relative/test
//     % ../jpv2rel < jpv2rel_jisp16_Nmax238_hw20.0_rel.in > jpv2rel_jisp16_Nmax238_hw20.0_rel.out
//
//   The input file contents are:
//
//     238 4
//     Vrel_JISP16_bare_Jmax4.hw20
//     jisp16_Nmax238_hw20.0_rel.dat


#include <fstream>

#include "basis/lsjt_operator.h"
#include "cppformat/format.h"
#include "mcutils/parsing.h"
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
