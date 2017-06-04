/****************************************************************
  moshinsky.cpp -- Performs Moshinsky transformation of 
  general relative operator in LSJT scheme.

  Mark A. Caprio, University of Notre Dame.
  Language: C++

  11/13/15 (mac): Created.

  11/21/15 (mac): Extract interaction I/O and manipulation functions
    to interaction_lsjt.  Allow input of generic relative interaction.

  11/26/15 (mac): Initial running version.

  01/06/16 (mac): Documentation updates.

  LEGACY version (note added 7/4/16, mac): Based on algorithm as
  sketched Fall 2015.  Works by summing over sectors of fixed Ncm.
  Supplanted by version working in blocks of fixed N.  Note that
  library calls are based on the old directory structure and before
  major overhauls to the am and basis modules.

****************************************************************/

#include <iomanip>

#include <am/wigner_gsl.h>

#include <shell/moshinsky_bracket.h>
#include <shell/indexing_lsjt.h>
#include <shell/interaction_lsjt.h>
#include <shell/indexing_jt.h>

////////////////////////////////////////////////////////////////
// transformation
////////////////////////////////////////////////////////////////

Eigen::MatrixXd
MoshinskyMatrix(
		const SpuriousRelativeSubspaceLSJT& spurious_relative_subspace,
		const TwoBodySubspaceLSJT& two_body_subspace
		)
// Generates Moshinsky transformation matrix in given L sector
// (rather, in a given LSJT sector, but the non-L labels are spectators), for
// fixed spuriosity (Ncm,lcm).
{
  Eigen::MatrixXd M(spurious_relative_subspace.Dimension(),two_body_subspace.Dimension());


  for (int i=0; i<spurious_relative_subspace.Dimension(); ++i)
    for (int j=0; j<two_body_subspace.Dimension(); ++j)
      {

	// retrieve spurious relative state
	SpuriousRelativeStateLSJT spurious_relative_state(spurious_relative_subspace,i);
	int Nr = spurious_relative_state.Nr();
	int lr = spurious_relative_state.lr();
	int nr = (Nr - lr)/2;
	int Ncm = spurious_relative_state.Ncm();
	int lcm = spurious_relative_state.lcm();
	int ncm = (Ncm - lcm)/2;
	int L = spurious_relative_state.L();

	// retrieve two-body state
	TwoBodyStateLSJT two_body_state(two_body_subspace,j);
	int N1 = two_body_state.N1();
	int l1 = two_body_state.l1();
	int n1 = (N1 - l1)/2;
	int N2 = two_body_state.N2();
	int l2 = two_body_state.l2();
	int n2 = (N2 - l2)/2;
	int L_check = two_body_state.L();
	assert(L==L_check);

	// evaluate bracket
	M(i,j) = MoshinskyBracket(nr,lr,ncm,lcm,n1,l1,n2,l2,L);
				  
      }

  //  std::cout << "Moshinsky " << M << std::endl;

  return M;
}


Eigen::MatrixXd
SpuriousRelativeMatrix(
		       const SpuriousRelativeSubspaceLSJT& spurious_relative_subspace2,
		       const SpuriousRelativeSubspaceLSJT& spurious_relative_subspace1,
		       const RelativeSpaceLSJT& relative_space,
		       const RelativeSectorsLSJT& relative_sectors,
		       const SectorMatrices& relative_matrices,
		       int J0
		       )
// Upgrades relative matrix to spurious relative matrix, for fixed
// spuriosity (Ncm,lcm).
//
// As(is2,is1) = sum of 
//              Unitary6JZ(lr2,lcm,L2,J2,S2,Jr2) * Unitary6JZ(lr1,lcm,L1,J1,S1,Jr1)
//		* RacahReductionFactorFirstSystem(Jr2,lcm,J2,Jr1,lcm,J1,J0)
//		* Ar(ir2,ir1);
{

  // initial matrix for accumulation
  Eigen::MatrixXd As = Eigen::MatrixXd::Zero(spurious_relative_subspace2.Dimension(),spurious_relative_subspace1.Dimension());

  // extract spurious relative (target) subspace labels
  int L2 = spurious_relative_subspace2.L();
  int S2 = spurious_relative_subspace2.S();
  int J2 = spurious_relative_subspace2.J();
  int T2 = spurious_relative_subspace2.T();
  int g2 = spurious_relative_subspace2.g();
  int L1 = spurious_relative_subspace1.L();
  int S1 = spurious_relative_subspace1.S();
  int J1 = spurious_relative_subspace1.J();
  int T1 = spurious_relative_subspace1.T();
  int g1 = spurious_relative_subspace1.g();

  // check delta condition on CM state
  int Ncm = spurious_relative_subspace2.Ncm();
  int lcm = spurious_relative_subspace2.lcm();
  int Ncm_check = spurious_relative_subspace1.Ncm();
  int lcm_check = spurious_relative_subspace1.lcm();
  assert(Ncm==Ncm_check);
  assert(lcm==lcm_check);

  // validate that relative (source) subspaces support required Nr
  int Nmax2 = spurious_relative_subspace2.Nmax();
  int Nmax1 = spurious_relative_subspace1.Nmax();
  int Nr_max = std::max(Nmax2,Nmax1) - Ncm;
  assert(relative_space.Nmax()>=Nr_max);

  // parities of relative (source) subspace determined by parity of spurious motion
  int gr2 = (g2+Ncm)%2;  
  int gr1 = (g1+Ncm)%2;

  // iterate over spurious relative (target) matrix elements
  for (int is2=0; is2<spurious_relative_subspace2.Dimension(); ++is2)
    for (int is1=0; is1<spurious_relative_subspace1.Dimension(); ++is1)
      {

	// retrieve spurious relative states
	SpuriousRelativeStateLSJT spurious_relative_state2(spurious_relative_subspace2,is2);
	int Nr2 = spurious_relative_state2.Nr();
	int lr2 = spurious_relative_state2.lr();
	SpuriousRelativeStateLSJT spurious_relative_state1(spurious_relative_subspace1,is1);
	int Nr1 = spurious_relative_state1.Nr();
	int lr1 = spurious_relative_state1.lr();

	// iterate over relevant relative sectors
	//   by triangle selection between relative and spurious motion
	for (int Jr2=abs(lcm-J2); Jr2<=(lcm+J2); ++Jr2)
	  for (int Jr1=abs(lcm-J1); Jr1<=(lcm+J1); ++Jr1)
	    {

	      // impose triangle selection on relative motion
	      if (!(AllowedTriangle(lr2,S2,Jr2)&&AllowedTriangle(lr1,S1,Jr1)))
		continue;

	      // impose triangle selection from relative operator
	      if (!AllowedTriangle(Jr2,J0,Jr1))
		continue;

	      // look up relative sector
	      int relative_subspace_index2 = relative_space.LookUpSubspaceIndex(RelativeSubspaceLSJT::SubspaceLabelsType(lr2,S2,Jr2,T2,gr2));
	      int relative_subspace_index1 = relative_space.LookUpSubspaceIndex(RelativeSubspaceLSJT::SubspaceLabelsType(lr1,S1,Jr1,T1,gr1));
	      int relative_sector_index = relative_sectors.LookUpSectorIndex(Sector(relative_subspace_index2,relative_subspace_index1));
	      //Eigen::MatrixXd& Ar = relative_matrices[relative_sector_index]; // fails as invalid initialization
	      Eigen::MatrixXd Ar = relative_matrices[relative_sector_index];

	      // locate matrix element
	      const RelativeSubspaceLSJT& relative_subspace2 = relative_space.GetSubspace(relative_subspace_index2);
	      const RelativeStateLSJT& relative_state2 = RelativeStateLSJT(relative_subspace2,RelativeSubspaceLSJT::StateLabelsType(Nr2));
	      int ir2 = relative_state2.index();
	      const RelativeSubspaceLSJT& relative_subspace1 = relative_space.GetSubspace(relative_subspace_index1);
	      const RelativeStateLSJT& relative_state1 = RelativeStateLSJT(relative_subspace1,RelativeSubspaceLSJT::StateLabelsType(Nr1));
	      int ir1 = relative_state1.index();

	      // accumulate matrix element
	      double contribution 
		= Unitary6JZ(lr2,lcm,L2,J2,S2,Jr2) * Unitary6JZ(lr1,lcm,L1,J1,S1,Jr1)
		* RacahReductionFactorFirstSystem(Jr2,lcm,J2,Jr1,lcm,J1,J0)
		* Ar(ir2,ir1);
	      As(is2,is1) += contribution;
	    }			  
      }

  return As;
}

Eigen::MatrixXd 
TransformedSector(
		  const TwoBodySubspaceLSJT& two_body_subspace2,
		  const TwoBodySubspaceLSJT& two_body_subspace1,
		  const RelativeSpaceLSJT& relative_space,
		  const RelativeSectorsLSJT& relative_sectors,
		  const SectorMatrices& relative_matrices,
		  int J0
		  )
// Carries out Moshinsky transform to generate the reduced matrix
// elements on a LSTJ-coupled two-body sector.  
// 
// Important: The input relative matrices must contain *reduced*
// matrix elements, not plain matrix elements.  The present code is
// oriented towards the tranformation of general spherical tensor
// relative operators, for which it is appropriate to work with
// reduced metrix elements.  In contrast, for simple scalar
// (interaction-type) relative operators, common convention is to
// specify the relative interaction in terms of its plain, unreduced
// matrix elements (which are M-independent) in the relative subspace.
//
// The output matrix element are antisymmetrized (AS) matrix elements,
// rather than normalized antisymmetrized (NAS) matrix elements (on
// the one hand) or distinguishable-particle matrix elements (on the
// other).
//
// Note inclusion of factor of 2 relative two "naive" application of
// Moshinsky transform on distinguishable-particle states, arising
// from the bookkeeping since our two-body states are antisymmetrized,
// rather than distinguishable-particle states.

{
  // define target matrix
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(two_body_subspace2.Dimension(),two_body_subspace1.Dimension());

  // extract subspace labels
  int L2 = two_body_subspace2.L();
  int S2 = two_body_subspace2.S();
  int J2 = two_body_subspace2.J();
  int T2 = two_body_subspace2.T();
  int g2 = two_body_subspace2.g();

  int L1 = two_body_subspace1.L();
  int S1 = two_body_subspace1.S();
  int J1 = two_body_subspace1.J();
  int T1 = two_body_subspace1.T();
  int g1 = two_body_subspace1.g();


  // obtain N truncations
  //
  // spurious spaces contribute for Nmax only up to lesser of the bra
  // and ket two-body space Nmax values
  int Nmax2 = two_body_subspace2.Nmax();
  int Nmax1 = two_body_subspace1.Nmax();
  int Ncm_max = std::min(Nmax2,Nmax1); 

  // accumulate cm contributions
  for (int Ncm=0; Ncm <= Ncm_max; ++Ncm)
    for (int lcm = Ncm%2; lcm <= Ncm; lcm+=2)
    {
      // define intermediate spurious relative subspaces
      const SpuriousRelativeSubspaceLSJT spurious_relative_subspace2(Ncm,lcm,L2,S2,J2,T2,g2,Nmax2);
      const SpuriousRelativeSubspaceLSJT spurious_relative_subspace1(Ncm,lcm,L1,S1,J1,T1,g1,Nmax1);

      // check that spurious relative subspace is nonempty
      if (spurious_relative_subspace2.Dimension()==0)
	continue;
      if (spurious_relative_subspace1.Dimension()==0)
	continue;

      // construct Ar_spurious
      // PLACEHOLDER: Eigen::MatrixXd As = Eigen::MatrixXd::Zero(spurious_relative_subspace2.Dimension(),spurious_relative_subspace1.Dimension());
      Eigen::MatrixXd As = SpuriousRelativeMatrix(
						  spurious_relative_subspace2,spurious_relative_subspace1,
						  relative_space,relative_sectors,relative_matrices,J0
						  );
      // std::cout << "As " << As << std::endl;

      // construct M2
      // PLACEHOLDER: Eigen::MatrixXd M2 = Eigen::MatrixXd::Zero(spurious_relative_subspace2.Dimension(),two_body_subspace2.Dimension());
      Eigen::MatrixXd M2 = MoshinskyMatrix(spurious_relative_subspace2,two_body_subspace2);
      // std::cout << "M2 " << M2 << std::endl;

      // construct M1
      // PLACEHOLDER: Eigen::MatrixXd M1 = Eigen::MatrixXd::Zero(spurious_relative_subspace1.Dimension(),two_body_subspace1.Dimension());
      Eigen::MatrixXd M1 = MoshinskyMatrix(spurious_relative_subspace1,two_body_subspace1);
      // std::cout << "M1 " << M1 << std::endl;
      
      // generate contribution to A
      Eigen::MatrixXd A_contribution = 2 * M2.transpose() * As * M1;

      A += A_contribution;
    }


  return A;

}



////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  // configuration parameters
  const int Nmax_relative = 2;
  const int J0 = 0;
  const int g0 = 0;
  const int Nmax_two_body = 2;

  // set up relative space
  std::cout << std::endl;
  std::cout << "Setting up relative space..." << std::endl;
  std::cout << "Nmax " << Nmax_relative << std::endl;
  RelativeSpaceLSJT relative_space(Nmax_relative);
  // WriteRelativeSubspaces(std::cout,relative_space);

  // set up relative operator
  std::cout << std::endl;
  std::cout << "Setting up relative sectors..." << std::endl;
  std::cout << "J0 " << J0 << " g0 " << g0 << std::endl;
  RelativeSectorsLSJT relative_sectors(relative_space,J0,g0);  // for scalar operator
  SectorMatrices Ar_matrices;
  std::cout << std::endl;

  // input relative operator
  std::cout << "Reading relative operator..." << std::endl;
  // ReadRelativeOperator(std::cin,relative_space,relative_sectors,Ar_matrices);
  SetOperatorToIdentityReduced(relative_space,relative_sectors,Ar_matrices);

  // output relative operator
  std::cout << std::endl;
  std::cout << "Relative operator" << std::endl;
  // caution: affects future output precision
  WriteRelativeOperatorMatrices(std::cout,relative_space,relative_sectors,Ar_matrices,10,7);

  // set up two-body space
  std::cout << std::endl;
  std::cout << "Setting up two-body space..." << std::endl;
  TwoBodySpaceLSJT two_body_space(Nmax_two_body);
  std::cout << "Nmax " << Nmax_two_body << std::endl;
  // WriteRelativeSubspaces(std::cout,relative_space);

  // construct two-body operator
  std::cout << std::endl;
  std::cout << "Transforming..." << std::endl;
  TwoBodySectorsLSJT two_body_sectors(two_body_space,J0,g0);  // for scalar operator
  SectorMatrices A_matrices;
  for (int two_body_sector_index = 0; two_body_sector_index < two_body_sectors.size(); ++two_body_sector_index)
    {
      Sector two_body_sector = two_body_sectors.GetSector(two_body_sector_index);
      const TwoBodySubspaceLSJT& two_body_subspace2 = two_body_space.GetSubspace(two_body_sector.index2());
      const TwoBodySubspaceLSJT& two_body_subspace1 = two_body_space.GetSubspace(two_body_sector.index1());
      Eigen::MatrixXd two_body_sector_matrix = TransformedSector(
								 two_body_subspace2,two_body_subspace1,
								 relative_space,relative_sectors,Ar_matrices,
								 J0 // for scalar operator
								 );

      A_matrices.push_back(two_body_sector_matrix);
    }

  // dump two_body operator
  std::cout << std::endl;
  std::cout << "Two-body operator" << std::endl;
  WriteTwoBodyOperator(std::cout,two_body_space,two_body_sectors,A_matrices);
  WriteTwoBodyOperatorMatrices(std::cout,two_body_space,two_body_sectors,A_matrices,10,7);

  

  // termination
  return 0;
}
