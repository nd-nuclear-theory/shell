// HEADER

////////////////////////////////////////////////////////////////
// (Nlj)^2 basis decomposed by Tz-J-P
////////////////////////////////////////////////////////////////

class TwoBodyBasisNljTzJP{

	public:

	////////////////////////////////////////////////////////////////
        // constructors
	////////////////////////////////////////////////////////////////
	
	// default constructor -- synthesized
	// Note: subsequent call to an initializer required for well-defined behavior

	// copy constructor -- synthesized

	////////////////////////////////////////////////////////////////
        // initializers
	////////////////////////////////////////////////////////////////

	void InitializeN1bN2b (int, int);

	////////////////////////////////////////////////////////////////
        // accessors
	////////////////////////////////////////////////////////////////

	// access of state by index
	const TwoBodyStateNlj& GetState (TwoSpeciesStateType state_type, int J, int g, int k) const;

	// lookup of index for state
	int GetIndex (TwoSpeciesStateType state_type, const TwoBodyStateNlj& state) const;

	////////////////////////////////////////////////////////////////
        // two-body basis information retrieval
	////////////////////////////////////////////////////////////////

	// subspace dimension retrieval
	int GetDimension (TwoSpeciesStateType state_type, int J, int g) const {
		return state_enumeration_[state_type][J][g].size();
	};

	// maximal J retrieval
	int JMax(TwoSpeciesStateType state_type) const {return J_max_[state_type];};

	////////////////////////////////////////////////////////////////
        // internal storage
	////////////////////////////////////////////////////////////////

	private:
	
	// index lookup initialization value
	static const int kNullIndex = -1;

	// maximal J values
	std::vector< int > J_max_;

	// state enumeration and index lookup data
	//    Note: data arranged by [state_type=kPP,kPN,kNN][J=0,1,...J_max][g=0,1]
 	typedef std::vector< std::vector< std::vector< TwoBodyStateSetNlj > > > StateEnumerationContainerType;
	typedef PairLookupArray< int > IndexLookupArrayType;
 	typedef std::vector< std::vector< std::vector< IndexLookupArrayType > > > StateLookupContainerType;
	StateEnumerationContainerType state_enumeration_;
	StateLookupContainerType state_lookup_;

};

inline
const TwoBodyStateNlj& TwoBodyBasisNljTzJP::GetState (TwoSpeciesStateType state_type, int J, int g, int k) const
{
	// check that state enumeration exists and index is in range for enumeration
	if ( !( ((0 <= k) && (k < state_enumeration_[state_type][J][g].size()) ) ) )
	{
		std::cerr << __FILE__ << ":" << __LINE__ << ": "
			  << "TwoBodyBasisNljTzJP::GetState: requested state index out of range for subspace"
			  << " (" << state_type << "," << J << "," << g <<")" 
			  << k
			  << std::endl;
		std::exit(EXIT_FAILURE);
	}

	return state_enumeration_[state_type][J][g][k];
};

inline 
int TwoBodyBasisNljTzJP::GetIndex (TwoSpeciesStateType state_type, const TwoBodyStateNlj& state) const
{

	// look up index by SP orbitals
	int J = state.J;
	int g = (state.a1.Getl() + state.a2.Getl()) % 2;
	int i1 = state.a1.GetIndex();
	int i2 = state.a2.GetIndex();

	// make sure orbitals are in range and in canonical order
	if ( (i1 >=state_lookup_[state_type][J][g].MinorDimension()) 
	     || (i2 >=state_lookup_[state_type][J][g].MinorDimension()) 
		)
	{
		std::cerr << __FILE__ << ":" << __LINE__ << ": "
			  << "TwoBodyBasisNljTzJP::GetIndex: requested orbitals above cutoff for subspace"
			  << " (" << state_type << "," << J << "," << g <<")" 
			  << " |" << state.a1 << state.a2 << ";" << J << ">"
			  << std::endl;
		std::exit(EXIT_FAILURE);
	}
	if ( (i1>i2) && (state_type != kPN) )
	{
		std::cerr << __FILE__ << ":" << __LINE__ << ": "
			  << "TwoBodyBasisNljTzJP::GetIndex: requested orbitals not in canonical order"
			  << " (" << state_type << "," << J << "," << g <<")" 
			  << " |" << state.a1 << state.a2 << ";" << J << ">"
			  << std::endl;
		std::exit(EXIT_FAILURE);
	}

	// look up index
	int index = state_lookup_[state_type][J][g](i1,i2);

	// make sure state was found in lookup
	if (index == kNullIndex)
	{
		std::cerr << __FILE__ << ":" << __LINE__ << ": "
			  << "TwoBodyBasisNljTzJP::GetIndex: requested orbital pair not found in subspace"
			  << " (" << state_type << "," << J << "," << g <<")" 
			  << " |" << state.a1 << state.a2 << ";" << J << ">"
			  << std::endl;
		std::exit(EXIT_FAILURE);
	}

	// return result
	return index;
};

////////////////////////////////////////////////////////////////
// (Nlj)^2 matrix dimensions decomposed by Tz-J-P
////////////////////////////////////////////////////////////////

// Notation:  NljTzJP refers to pn-scheme, with Nlj orbitals and conserved (Tz,J,P)

// expected size (number of matrix elements)
//   overloaded: (1) total, (2) then by state_type, (3) then by specific (J,g) space

int TwoBodyMatrixNljTzJPDimension(const TwoBodyBasisNljTzJP& two_body_basis, TwoSpeciesStateType state_type, int J, int g);
int TwoBodyMatrixNljTzJPDimension(const TwoBodyBasisNljTzJP& two_body_basis, TwoSpeciesStateType state_type);
int TwoBodyMatrixNljTzJPDimension(const TwoBodyBasisNljTzJP& two_body_basis);


////////////////////////////////////////////////////////////////
// (Nlj)^2 matrix elements decomposed by Tz-J-P
////////////////////////////////////////////////////////////////

// class for holding matrix elements
// Note: requires continued existence of TwoBodyBasisNljTzJP data

class TwoBodyMatrixNljTzJP{

	public:

	////////////////////////////////////////////////////////////////
	// type definitions
	////////////////////////////////////////////////////////////////

	typedef double MatrixElementType;

	////////////////////////////////////////////////////////////////
        // constructors
	////////////////////////////////////////////////////////////////
	
	// default constructor -- disabled
	TwoBodyMatrixNljTzJP ();


	// set up storage structure 

	// Note: We must specify two body basis *reference* in the
	// *constructor*, since the initializer list is the only legal
	// way to set the reference variable two_body_basis_.  This is
	// okay, since since do not expect need to use more than one
	// basis with the same matrix element container.

	// Note: This does *not* allocate memory for the TBMEs
	// themselves, which must be done by Initialize().

	TwoBodyMatrixNljTzJP (const TwoBodyBasisNljTzJP&);

	// copy constructor -- synthesized

	////////////////////////////////////////////////////////////////
        // allocation and deallocation
	////////////////////////////////////////////////////////////////

	// allocate and initialize single (Tz,J,g) subspace
	void Initialize (TwoSpeciesStateType state_type, int J, int g);

	// allocate and initialize all subspaces
	void Initialize ();

	// deallocate single (Tz,J,g) subspace
	void Free (TwoSpeciesStateType state_type, int J, int g);

	// deallocate all subspaces
	void Free ();

	////////////////////////////////////////////////////////////////
        // accessors
	////////////////////////////////////////////////////////////////

	// Storage is of "canonically ordered" matrix elements <s1|O|s2>
        // for "normalized antisymmetrized" (NAS) states.

	// Canonical ordering: 
	//    -- for like particles: s1.a1 <= s1.a2 by orbital index
	//    -- for like particles: s2.a1 <= s2.a2 by orbital index
	//    -- s1 <= s2 lexicographically by orbital indices

	// Normalization:
	//    NAS = normalized antisymmetrized
	//    UNAS = unnormalized antisymmetrized
	// These both refer to angular momentum coupling of the
	// antisymmetrized product states, or 
	//    (a^\dagger_{j_1} \times a^\dagger_{j_2})^J
	// acting on the vacuum in second quantized notation --
	// but then either *including* a normalization
	// factor 1/sqrt(2) for like particles (NAS) or *not*
	// including this factor (UNAS).

	// GetMatrixElement -- returns NAS ME, arguments must be canonical
	MatrixElementType GetMatrixElement (TwoSpeciesStateType state_type, const TwoBodyStateNlj& s1, const TwoBodyStateNlj& s2) const 
	{
		return InternalReference(state_type, s1, s2);
	};

	// GetMatrixElement -- returns UNAS ME, arguments must be canonical
	MatrixElementType GetMatrixElementUNAS (TwoSpeciesStateType state_type, const TwoBodyStateNlj& s1, const TwoBodyStateNlj& s2) const 
	{
		double factor = 1.;
		if (state_type != kPN)
		{ 
			if (s1.a1 == s1.a2)
				factor *= kSqrtTwo;
			if (s2.a1 == s2.a2)
				factor *= kSqrtTwo;
		}
		return factor * InternalReference(state_type, s1, s2);
	};

	// SetMatrixElement -- accepts NAS ME, arguments must be canonical
	void SetMatrixElement (TwoSpeciesStateType state_type, const TwoBodyStateNlj& s1, const TwoBodyStateNlj& s2, MatrixElementType value) {
		InternalReference(state_type, s1, s2) = value;
	};

	// SetMatrixElement -- accepts UNAS ME, arguments must be canonical
	void SetMatrixElementUNAS (TwoSpeciesStateType state_type, const TwoBodyStateNlj& s1, const TwoBodyStateNlj& s2, MatrixElementType value) {
		double factor = 1.;
		if (state_type != kPN)
		{ 
			if (s1.a1 == s1.a2)
				factor *= kSqrtTwoInverse;
			if (s2.a1 == s2.a2)
				factor *= kSqrtTwoInverse;
		}
		InternalReference(state_type, s1, s2) = factor * value;

	};


	////////////////////////////////////////////////////////////////
	// matrix and basis information retrieval
	////////////////////////////////////////////////////////////////

	// reference to underlying two-body basis
	const TwoBodyBasisNljTzJP& GetTwoBodyBasis () const {return two_body_basis_;};

	// Note: this is essentially a pass-through of information
	// which can be calculated from the underlying two_body basis

	// dimension information -- of (state_type,J,g) sector
	int GetDimension (TwoSpeciesStateType state_type, int J, int g) const {
		return two_body_basis_.GetDimension(state_type, J, g);
	};

	// size information -- for full matrix, of (state_type) sector, or of (state_type,J,g) sector
	int GetSize (TwoSpeciesStateType state_type, int J, int g) const {
		return TwoBodyMatrixNljTzJPDimension(two_body_basis_, state_type, J, g);
	};
	int GetSize (TwoSpeciesStateType state_type) const {
		return TwoBodyMatrixNljTzJPDimension(two_body_basis_, state_type);
	};
	int GetSize () const {
		return TwoBodyMatrixNljTzJPDimension(two_body_basis_);
	};

	private:

	////////////////////////////////////////////////////////////////
        // internal indexing
	////////////////////////////////////////////////////////////////

	// lookup of entry corresponding to given bra-ket pair
	MatrixElementType& InternalReference (TwoSpeciesStateType state_type, const TwoBodyStateNlj& s1, const TwoBodyStateNlj& s2);

	const MatrixElementType& InternalReference (TwoSpeciesStateType state_type, const TwoBodyStateNlj& s1, const TwoBodyStateNlj& s2) const;

	////////////////////////////////////////////////////////////////
        // internal storage
	////////////////////////////////////////////////////////////////

	// reference to underlying basis
	// CAVEAT: may be dangerous if two body basis goes out of scope and is destroyed
	const TwoBodyBasisNljTzJP& two_body_basis_;

	// matrix element data
	//    Note: data arranged by [state_type=kPP,kPN,kNN][J=0,1,...J_max][g=0,1]
	typedef PairLookupArray<MatrixElementType> MatrixType;
 	typedef std::vector< std::vector< std::vector< MatrixType > > > MatrixContainerType;
	MatrixContainerType matrix_;

};

////////////////////////////////////////////////////////////////
// (Nlj)^2 matrix elements (inline implementation)
////////////////////////////////////////////////////////////////

inline
TwoBodyMatrixNljTzJP::MatrixElementType& TwoBodyMatrixNljTzJP::InternalReference (TwoSpeciesStateType state_type, const TwoBodyStateNlj& s1, const TwoBodyStateNlj& s2) 
{

	// extract J for subspace
	int J = s1.J;

	// extract g for subspace
	int g = ( s1.a1.Getl() + s1.a2.Getl() ) % 2;

	// recall state indices
	int k1 = two_body_basis_.GetIndex(state_type,s1);
	int k2 = two_body_basis_.GetIndex(state_type,s2);

	// make sure states are in canonical order
	if ( !(k1 <= k2) ) 
	{
		std::cerr << __FILE__ << ":" << __LINE__ << ": "
			  << "TwoBodyMatrixNljTzJP::InternalReference: requested bra and ket not in canonical order"
			  << " (" << state_type << "," << J << "," << g <<")" 
			  << " <" << k1 << "|O|" << k2 << ">"
			  << std::endl;
		std::exit(EXIT_FAILURE);
	}
		
	// access matrix element
	return matrix_[state_type][J][g](k1,k2);
};

inline
const TwoBodyMatrixNljTzJP::MatrixElementType& TwoBodyMatrixNljTzJP::InternalReference (TwoSpeciesStateType state_type, const TwoBodyStateNlj& s1, const TwoBodyStateNlj& s2) const {


	// extract J for subspace
	int J = s1.J;

	// extract g for subspace
	int g = ( s1.a1.Getl() + s1.a2.Getl() ) % 2;

	// recall state indices
	int k1 = two_body_basis_.GetIndex(state_type,s1);
	int k2 = two_body_basis_.GetIndex(state_type,s2);

	// make sure states are in canonical order
	if ( !(k1 <= k2) ) 
	{
		std::cerr << __FILE__ << ":" << __LINE__ << ": "
			  << "TwoBodyMatrixNljTzJP::InternalReference: requested bra and ket not in canonical order"
			  << " (" << state_type << "," << J << "," << g <<")" 
			  << " <" << k1 << "|O|" << k2 << ">"
			  << std::endl;
		std::exit(EXIT_FAILURE);
	}
		
	// access matrix element
	return matrix_[state_type][J][g](k1,k2);
};

// BODY

////////////////////////////////////////////////////////////////
// (Nlj)^2 matrix
////////////////////////////////////////////////////////////////

// matrix size

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
	for (TwoSpeciesStateType state_type = kPP; state_type <= kNN; ++state_type)
		size += TwoBodyMatrixNljTzJPDimension(two_body_basis, state_type);
	
	return size;
};

// constructor

TwoBodyMatrixNljTzJP::TwoBodyMatrixNljTzJP (const TwoBodyBasisNljTzJP& two_body_basis)
	: two_body_basis_(two_body_basis)
{

	// for each state class
	matrix_.resize(3);
	for (TwoSpeciesStateType state_type = kPP; state_type <= kNN; ++state_type)
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

};

void TwoBodyMatrixNljTzJP::Initialize (TwoSpeciesStateType state_type, int J, int g) {
	matrix_[state_type][J][g].Initialize(kUpperTriangular,two_body_basis_.GetDimension(state_type,J,g),0.);
};


void TwoBodyMatrixNljTzJP::Initialize () {
	for (TwoSpeciesStateType state_type = kPP; state_type <= kNN; ++state_type)
	{
		// for each J
		int J_max = two_body_basis_.JMax(state_type);
		for (int J = 0; J <= J_max; ++J)
			// for each parity
			for (int g = 0; g <= 1; ++g)
				Initialize(state_type,J,g);
	}
};

void TwoBodyMatrixNljTzJP::Free (TwoSpeciesStateType state_type, int J, int g) {
	matrix_[state_type][J][g].Free();
};


void TwoBodyMatrixNljTzJP::Free () {
	for (TwoSpeciesStateType state_type = kPP; state_type <= kNN; ++state_type)
	{
		// for each J
		int J_max = two_body_basis_.JMax(state_type);
		for (int J = 0; J <= J_max; ++J)
			// for each parity
			for (int g = 0; g <= 1; ++g)
				Free(state_type,J,g);
	}
};

////////////////////////////////////////////////////////////////
// EXAMPLE: loop over subspaces
////////////////////////////////////////////////////////////////

//	// for each state class
//	for (TwoSpeciesStateType state_type = kPP; state_type <= kNN; ++state_type)
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






	// set up storage structure 

	// Note: We must specify two body basis *reference* in the
	// *constructor*, since the initializer list is the only legal
	// way to set the reference variable two_body_basis_.  This is
	// okay, since since do not expect need to use more than one
	// basis with the same matrix element container.

	// Note: This does *not* allocate memory for the TBMEs
	// themselves, which must be done by Initialize().

	TwoBodyMatrixNljTzJP (const TwoBodyBasisNljTzJP&);

	// reference to underlying basis
	// CAVEAT: may be dangerous if two body basis goes out of scope and is destroyed
	const TwoBodyBasisNljTzJP& two_body_basis_;

// constructor

TwoBodyMatrixNljTzJP::TwoBodyMatrixNljTzJP (const TwoBodyBasisNljTzJP& two_body_basis)
	: two_body_basis_(two_body_basis)
{

	// for each state class
	matrix_.resize(3);
	for (TwoSpeciesStateType state_type = kPP; state_type <= kNN; ++state_type)
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

};
