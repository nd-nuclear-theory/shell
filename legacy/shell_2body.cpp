/****************************************************************
  shell_2body.cpp

  Created by Mark A. Caprio, University of Notre Dame.
  Last modified 4/25/15.

****************************************************************/

#include <shell/shell_2body.h>
#include <cmath>

namespace shell {

  ////////////////////////////////////////////////////////////////
  // arithmetic constants
  ////////////////////////////////////////////////////////////////

  const double kSqrtTwo = std::sqrt(static_cast<double>(2));
  const double kSqrtTwoInverse = 1. / std::sqrt(static_cast<double>(2));


  ////////////////////////////////////////////////////////////////
  // (Nlj)^2 sector iterator implementation
  ////////////////////////////////////////////////////////////////

  SectorNljTzJP::SectorNljTzJP () {
    SetTruncation(0,0);
  }

  SectorNljTzJP::SectorNljTzJP (const TwoBodyBasisNljTzJP& Basis) {

    // initialize to zeroth sector
    k_ = 0;

    // save max J values for each state type
    for (TwoSpeciesStateType state_type = kTwoSpeciesStateTypeBegin; state_type <= kTwoSpeciesStateTypeEnd; ++state_type)
      J_max_.push_back(Basis.JMax(state_type));

    // store total num sectors
    num_sectors_ = 0;
    for (TwoSpeciesStateType state_type = kTwoSpeciesStateTypeBegin; state_type <= kTwoSpeciesStateTypeEnd; ++state_type)
      num_sectors_ += 2 * (J_max_[state_type]+1);

  }

  SectorNljTzJP::SectorNljTzJP (int N1b_max, int N2b_max) {
    SetTruncation(N1b_max, N2b_max);
  }


  void SectorNljTzJP::SetTruncation (int N1b_max, int N2b_max) {

    // initialize to zeroth sector
    k_ = 0;

    // initialize two-body maximal J parameters
    int J_max_like = std::min( (2 * N1b_max), N2b_max + 1);
    int J_max_unlike = std::min( (2 * N1b_max + 1), N2b_max + 1);
    J_max_.resize(3);
    J_max_[kPP] = J_max_like;
    J_max_[kPN] = J_max_unlike;
    J_max_[kNN] = J_max_like;

    // store total num sectors
    num_sectors_ = 0;
    for (TwoSpeciesStateType state_type = kTwoSpeciesStateTypeBegin; state_type <= kTwoSpeciesStateTypeEnd; ++state_type)
      num_sectors_ += 2 * (J_max_[state_type]+1);

  }

  int SectorNljTzJP::NumSectors () const {
    return num_sectors_;
  };

  int SectorNljTzJP::GetIndex () const {
    return k_;
  }

  int SectorNljTzJP::GetStateType () const {
    int sector_count = 0;
    for (TwoSpeciesStateType state_type = kTwoSpeciesStateTypeBegin; state_type <= kTwoSpeciesStateTypeEnd; ++state_type)
      {
	sector_count += 2 * (J_max_[state_type]+1);
	if (k_ < sector_count)
	  return state_type;
      }
  }

  int SectorNljTzJP::GetJ () const {
    int sector_count = 0;
    // for each state class
    for (TwoSpeciesStateType state_type = kTwoSpeciesStateTypeBegin; state_type <= kTwoSpeciesStateTypeEnd; ++state_type)
      {
	// for each J
	int J_max = J_max_[state_type];
	for (int J = 0; J <= J_max; ++J)
	  {
	    sector_count += 2;
	    if (k_ < sector_count)
	      return J;
	  }
      }
  }

  int SectorNljTzJP::GetGrade () const {
    return (k_ % 2);
  }

  SectorNljTzJP& SectorNljTzJP::operator ++ ()
  {
    ++k_;

    return *this;
  }

  SectorNljTzJP SectorNljTzJP::operator ++ (int)
  {
    SectorNljTzJP x = *this;
    ++k_;
    return x;
  }

  bool SectorNljTzJP::InRange () const {
    return k_ < NumSectors ();
  }

  bool SectorNljTzJP::IsFirstOfType () const {
    // iterate over state types, accumulating total number of sectors
    int sector_count = 0;
    for (TwoSpeciesStateType state_type = kTwoSpeciesStateTypeBegin; state_type <= kTwoSpeciesStateTypeEnd; ++state_type)
      {
	// flag if current sector lies at beginning of given state type
	if (k_ == sector_count)
	  return true;
	sector_count += 2 * (J_max_[state_type]+1);
      }
    // otherwise current sector does not lie at beginning of any state type
    return false;
  }

  bool SectorNljTzJP::IsLastOfType () const {
    // iterate over state types, accumulating total number of sectors
    int sector_count = 0;
    for (TwoSpeciesStateType state_type = kTwoSpeciesStateTypeBegin; state_type <= kTwoSpeciesStateTypeEnd; ++state_type)
      {
	sector_count += 2 * (J_max_[state_type]+1);
	// flag if current sector lies at end of given state type
	if (k_ == sector_count -1)
	  return true;
      }
    // otherwise current sector does not lie at end of any state type
    return false;
  }

  ////////////////////////////////////////////////////////////////
  // (Nlj)^2 basis
  ////////////////////////////////////////////////////////////////

  // InitializeN1bN2b -- following combined 1-body and 2-body phonon truncation scheme
  //   N1b_max -- maximum N (=2*n+l) for one-body space (shell cutoff)
  //   N2b_max -- maximum N (=N1+N2) for two-body space (phonon cutoff)
  // Note: For pure one-body cutoff (FCI-like), use N2b=2*N1b.

  void TwoBodyBasisNljTzJP::InitializeN1bN2b(int N1b_max, int N2b_max)
  {

    // save definition

    N1b_max_ = N1b_max;
    N2b_max_ = N2b_max;

    // initialize one-body space parameters

    int sp_dimension = LevelCountNlj(N1b_max);
    SPOrbitalNlj a_min=SPOrbitalNlj(0);
    SPOrbitalNlj a_max=SPOrbitalNlj(sp_dimension-1);

    // initialize two-body maximal J parameters
    int J_max_like = std::min( (2 * N1b_max), N2b_max + 1);
    int J_max_unlike = std::min( (2 * N1b_max + 1), N2b_max + 1);
    J_max_.resize(3);
    J_max_[kPP] = J_max_like;
    J_max_[kPN] = J_max_unlike;
    J_max_[kNN] = J_max_like;

    // clear any old values
    state_enumeration_.clear();
    state_lookup_.clear();

    // set up enumeration and lookup array framework

    // for each state class
    // TODO: write cleanly using sector iterator
    state_enumeration_.resize(3);
    state_lookup_.resize(3);
    for (TwoSpeciesStateType state_type = kTwoSpeciesStateTypeBegin; state_type <= kTwoSpeciesStateTypeEnd; ++state_type)
      {
	// set up prototype indexing array
	IndexingType indexing_type = (state_type == kPN) ? kSquare : kUpperTriangular;

	// for each J
	state_enumeration_[state_type].resize(J_max_[state_type]+1);
	state_lookup_[state_type].resize(J_max_[state_type]+1);
	for (int J = 0; J <= J_max_[state_type]; ++J)
	  {
	    // for each parity
	    //   Note: enumeration vector is initialized to null length by default constructor
	    state_enumeration_[state_type][J].resize(2);
	    state_lookup_[state_type][J].resize(2);
	    for (int g = 0; g <= 1; ++g)
	      state_lookup_[state_type][J][g].Initialize(indexing_type, sp_dimension, kNullIndex);
	  }
      }
	
    // enumerate and index all |a1,a2,J> states

    //   Note: lexicographical is retained, and only canonically
    //   ordered pairs (a1 <= a2) are retained for like particles

    for (SPOrbitalNlj a1 = a_min; a1 <= a_max; ++a1)
      for (SPOrbitalNlj a2 = a_min; a2 <= a_max; ++a2)
	for (int J = IValue(abs( a1.Getj() - a2.Getj() )); J <= IValue( a1.Getj() + a2.Getj() ); ++J) 
	  {
	    // triage state by parity
	    int g = (a1.Getl() + a2.Getl()) % 2;
				
	    for (TwoSpeciesStateType state_type = kTwoSpeciesStateTypeBegin; state_type <= kTwoSpeciesStateTypeEnd; ++state_type)
	      {

		// skip states exceeding N2b limit
		if ( (a1.GetN() + a2.GetN()) > N2b_max )
		  continue;

		// for like species: skip noncanonical-ordered or antisymmetry-forbidden states 
		bool like = ( (state_type == kPP) || (state_type == kNN) );
		if ( (like) &&  (a1 > a2) )
		  continue;
		if ( (like) &&  (a1 == a2) && ((J%2) == 1) )
		  continue;

		// store state and index
		int k = state_enumeration_[state_type][J][g].size();
		state_enumeration_[state_type][J][g].push_back(TwoBodyStateNlj(a1, a2, J));
		state_lookup_[state_type][J][g](a1.GetIndex(), a2.GetIndex()) = k;
	      }
	  }
  }

  ////////////////////////////////////////////////////////////////
  // (Nlj)^2 matrix dimensions
  ////////////////////////////////////////////////////////////////

  int TwoBodyMatrixNljTzJPDimension (const TwoBodyBasisNljTzJP& two_body_basis, TwoSpeciesStateType state_type, int J, int g)
  {
    int dimension = two_body_basis.GetDimension(state_type,J,g);
    return IndexingDimension(kUpperTriangular, dimension);
  }

  int TwoBodyMatrixNljTzJPDimension (const TwoBodyBasisNljTzJP& two_body_basis, TwoSpeciesStateType state_type)
  {
    int size = 0;

    // for each J
    int J_max = two_body_basis.JMax(state_type);
    for (int J = 0; J <= J_max; ++J)
      // for each parity
      for (int g = 0; g <= 1; ++g)
	size += TwoBodyMatrixNljTzJPDimension(two_body_basis, state_type, J, g);

    return size;
  }

  int TwoBodyMatrixNljTzJPDimension (const TwoBodyBasisNljTzJP& two_body_basis)
  {
    int size = 0;
    // allocate matrices

    // for each state class
    for (TwoSpeciesStateType state_type = kTwoSpeciesStateTypeBegin; state_type <=kTwoSpeciesStateTypeEnd; ++state_type)
      size += TwoBodyMatrixNljTzJPDimension(two_body_basis, state_type);
	
    return size;
  }

  ////////////////////////////////////////////////////////////////
  // (Nlj)^2 matrix member functions
  ////////////////////////////////////////////////////////////////

  TwoBodyMatrixNljTzJP::TwoBodyMatrixNljTzJP () 
  {
  }

  TwoBodyMatrixNljTzJP::TwoBodyMatrixNljTzJP (const TwoBodyBasisNljTzJP& two_body_basis)
  {
    SetBasis(two_body_basis);
  }

  void TwoBodyMatrixNljTzJP::SetBasis (const TwoBodyBasisNljTzJP& two_body_basis)
  {

    // save copy of basis
    two_body_basis_ = two_body_basis;

    // reset entire matrix
    matrix_.clear();

    // for each state class
    matrix_.resize(3);
    for (TwoSpeciesStateType state_type = kTwoSpeciesStateTypeBegin; state_type <= kTwoSpeciesStateTypeEnd; ++state_type)
      {
	// for each J
	int J_max = two_body_basis_.JMax(state_type);
	matrix_[state_type].resize(J_max + 1);
	for (int J = 0; J <= J_max; ++J)
	  {
	    // for each parity
	    matrix_[state_type][J].resize(2);
	  }
      }

  }

  void TwoBodyMatrixNljTzJP::Initialize (TwoSpeciesStateType state_type, int J, int g) {
    matrix_[state_type][J][g].Initialize(kUpperTriangular,two_body_basis_.GetDimension(state_type,J,g),0.);
  }

  void TwoBodyMatrixNljTzJP::Initialize (const SectorNljTzJP& sector) {
    Initialize(sector.GetStateType(),sector.GetJ(),sector.GetGrade());
  }


  void TwoBodyMatrixNljTzJP::Initialize () {
    // TODO: write cleanly using sector iterator
    for (TwoSpeciesStateType state_type = kTwoSpeciesStateTypeBegin; state_type <= kTwoSpeciesStateTypeEnd; ++state_type)
      {
	// for each J
	int J_max = two_body_basis_.JMax(state_type);
	for (int J = 0; J <= J_max; ++J)
	  // for each parity
	  for (int g = 0; g <= 1; ++g)
	    Initialize(state_type,J,g);
      }
  }

  void TwoBodyMatrixNljTzJP::Free (TwoSpeciesStateType state_type, int J, int g) {
    matrix_[state_type][J][g].Free();
  }

  void TwoBodyMatrixNljTzJP::Free (const SectorNljTzJP& sector){
    Free(sector.GetStateType(),sector.GetJ(),sector.GetGrade());
  }


  void TwoBodyMatrixNljTzJP::Free () {
    // TODO: write cleanly using sector iterator
    for (TwoSpeciesStateType state_type = kTwoSpeciesStateTypeBegin; state_type <= kTwoSpeciesStateTypeEnd; ++state_type)
      {
	// for each J
	int J_max = two_body_basis_.JMax(state_type);
	for (int J = 0; J <= J_max; ++J)
	  // for each parity
	  for (int g = 0; g <= 1; ++g)
	    Free(state_type,J,g);
      }
  }


  ////////////////////////////////////////////////////////////////
  // (Nlj)^2 matrix sector operations
  ////////////////////////////////////////////////////////////////

  // TODO: retrofit to use sector argument

  void TwoBodyMatrixSectorCopy (const TwoBodyMatrixNljTzJP& source_matrix, TwoBodyMatrixNljTzJP& destination_matrix, TwoSpeciesStateType state_type, int J, int g)
  {
    const int dimension = destination_matrix.GetDimension(state_type,J,g);
				
    // for canonical pairs of states in destination two-body space (in lexicographical order)
    for (int k1 = 0; k1 < dimension; ++k1)
      for (int k2 = k1; k2 < dimension; ++k2)
	{
	  TwoBodyStateNlj s1 = destination_matrix.GetTwoBodyBasis().GetState(state_type,J,g,k1);
	  TwoBodyStateNlj s2 = destination_matrix.GetTwoBodyBasis().GetState(state_type,J,g,k2);
			
	  double matrix_element;
	  matrix_element = source_matrix.GetMatrixElementNAS(state_type, s1, s2);
	  destination_matrix.SetMatrixElementNAS(state_type, s1, s2, matrix_element); 
	}
  }

  void TwoBodyMatrixSectorCopy (const TwoBodyMatrixNljTzJP& source_matrix, TwoBodyMatrixNljTzJP& destination_matrix, const SectorNljTzJP& sector)
  {
    TwoBodyMatrixSectorCopy (source_matrix, destination_matrix, sector.GetStateType(), sector.GetJ(), sector.GetGrade());
  }


  void TwoBodyMatrixSectorAdd (const TwoBodyMatrixNljTzJP& source_matrix, double scale, TwoBodyMatrixNljTzJP& destination_matrix, TwoSpeciesStateType state_type, int J, int g)
  {
    const int dimension = destination_matrix.GetDimension(state_type,J,g);
				
    // for canonical pairs of states in destination two-body space (in lexicographical order)
    for (int k1 = 0; k1 < dimension; ++k1)
      for (int k2 = k1; k2 < dimension; ++k2)
	{
	  TwoBodyStateNlj s1 = destination_matrix.GetTwoBodyBasis().GetState(state_type,J,g,k1);
	  TwoBodyStateNlj s2 = destination_matrix.GetTwoBodyBasis().GetState(state_type,J,g,k2);
			
	  double matrix_element;
	  matrix_element = destination_matrix.GetMatrixElementNAS(state_type, s1, s2);
			
	  if (source_matrix.HasMatrixElement(state_type, s1, s2))
	    matrix_element += scale * source_matrix.GetMatrixElementNAS(state_type, s1, s2);
	  destination_matrix.SetMatrixElementNAS(state_type, s1, s2, matrix_element); 
	}
  }

  void TwoBodyMatrixSectorAdd (const TwoBodyMatrixNljTzJP& source_matrix, double scale, TwoBodyMatrixNljTzJP& destination_matrix, const SectorNljTzJP& sector)
  {
    TwoBodyMatrixSectorAdd (source_matrix, scale, destination_matrix, sector.GetStateType(), sector.GetJ(), sector.GetGrade());
  }


  void TwoBodyMatrixSectorAddIdentity (double scale, TwoBodyMatrixNljTzJP& destination_matrix, TwoSpeciesStateType state_type, int J, int g)
  {
    const int dimension = destination_matrix.GetDimension(state_type,J,g);
				
    // for each diagonal entry
    for (int k1 = 0; k1 < dimension; ++k1)
      {
	TwoBodyStateNlj s1 = destination_matrix.GetTwoBodyBasis().GetState(state_type,J,g,k1);
	TwoBodyStateNlj s2 = s1;

	double matrix_element;
	matrix_element = destination_matrix.GetMatrixElementNAS(state_type, s1, s2);
	matrix_element += scale;
	destination_matrix.SetMatrixElementNAS(state_type, s1, s2, matrix_element); 
      }
  }

  void TwoBodyMatrixSectorAddIdentity (double scale, TwoBodyMatrixNljTzJP& destination_matrix, const SectorNljTzJP& sector)
  {
    TwoBodyMatrixSectorAddIdentity (scale, destination_matrix, sector.GetStateType(), sector.GetJ(), sector.GetGrade());
  }



  ////////////////////////////////////////////////////////////////
  // EXAMPLE: loop over subspaces
  ////////////////////////////////////////////////////////////////

  // OLD STYLE:
  //
  //	// for each state class
  //	for (TwoSpeciesStateType state_type = kTwoSpeciesStateTypeBegin; state_type <= kTwoSpeciesStateTypeEnd; ++state_type)
  //	{
  //		// for each J
  //		int J_max = two_body_basis_.JMax(state_type);
  //		for (int J = 0; J <= J_max; ++J)
  //		{
  //			// for each parity
  //			for (int g = 0; g <= 1; ++g)
  //				matrix_[state_type][J][g].Initialize(kUpperTriangular,two_body_basis_.GetDimension(state_type,J,g),0.);
  //		}
  //	}


  // NEW STYLE:
  //
  //	for (SectorNljTzJP sector(output_basis); sector.InRange(); ++sector)
  //	{
  //		// recover sector properties
  //		const TwoSpeciesStateType state_type = sector.GetStateType();
  //		const int J = sector.GetJ();
  //		const int g = sector.GetGrade();
  //		
  //		cout << "sector " << state_type << " " << J << " " << g << endl;
  //	}


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
