/****************************************************************
  jpv_io.h

  Read JPV format relative matrix element files.

  Language: C++11
                                 
  Mark A. Caprio
  University of Notre Dame

  3/26/17 (mac): Extracted from jpv2rel (created 7/25/16).
  3/28/17 (mac):
    - Implement input of non-isoscalar operators.
    - Add verbose option to input functions.
  4/9/17 (mac): Finish implementing conversion of non-isoscalar operators.

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

// JPV operators are always represented within the storage appropriate
// to an isoscalar operator.  A non-isoscalar operator is stored
// (indexing-wise) as if it were three isoscalar operators.

#ifndef JPV_IO_H_
#define JPV_IO_H_

#include "basis/lsjt_operator.h"

namespace relative {

void ReadJPVOperator(
    const std::string& source_filename,
    const basis::RelativeSpaceLSJT& relative_space,
    const basis::OperatorLabelsJT& operator_labels,
    const std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
    std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
    bool verbose = false
  );
// Read a single isoscalar Hamiltonian-like relative operator in JPV
// format.
//
// The input operator is represented in terms of non-reduced matrix
// elements between good-JT (Tz=0) states, which, for an isoscalar
// operator, are identical to group-theory convention JT-reduced
// matrix elments.
//
// The input process carried out by this function also suffices to
// read the data for one of the pp, pn, or nn components of a
// non-isoscalar operator in JPV format, though further processing
// (see ReadJPVOperatorPN) is required to extract the actual isospin
// components.
//
// The operator_labels argument is only included for validation
// purposes.  It is checked with an assertion that the operator being
// read is a Hamiltonian-like operator with T0_max=0.  There is no
// further useful information in these labels.
//
// The returned relative_component_sectors and
// relative_component_matrices from this function will only have a
// T0=0 component.
//
// Arguments:
//   source_filename (input): input filename
//   relative_space (input): target space
//   operator_labels (input): operator labels
//   relative_component_sectors (input): target sectors
//   relative_component_matrices (output): target matrices to be populated
//   verbose (input, optional): verbosity

void ReadJPVOperatorPN(
    const std::array<std::string,3>& source_filenames,
    const basis::RelativeSpaceLSJT& relative_space,
    const basis::RelativeOperatorParametersLSJT& operator_parameters,
    const std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
    std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
    bool verbose = false
  );
// Read a generic (non-isoscalar) Hamiltonian-like relative operator in JPV
// format.
//
// The input operator is represented in terms of non-reduced matrix
// elements between good-JT (and Tz) states, for Tz = +1 (pp), 0 (pn),
// and -1 (nn).  The T=0 sectors for the pp and nn inputs are
// explicitly stored but contain all zero entries.
//
// The operator_parameters argument is used to allocate intermediate
// storage to hold the pp/nn/pn matrix elements.  The
// operator_parameters argument is also validated.  It is checked with
// an assertion that the operator being read is a Hamiltonian-like
// operator with T0_max=2.
//
// Arguments:
//   source_filenames (input): input filenames (in order pp, nn, pn)
//   relative_space (input): target space
//   operator_parameters (input): operator labels and truncation
//   relative_component_sectors (input): target sectors
//   relative_component_matrices (output): target matrices to be populated
//   verbose (input, optional): verbosity

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
