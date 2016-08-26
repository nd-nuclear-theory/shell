struct SPOrbitalNj {
	explicit SPLabelsNj(int k);
	SPLabelsNj(int N0, const HalfInt& j0);
	int N;
	int l;
	HalfInt j;
};

	// iterate over all orbitals for first index, in ascending order
	for (int N1 = 0; N1 <= Ni_max; ++N1)
		for (HalfInt j1 = HalfInt(1,2); j1 <= HalfInt(1,2) + N1; j1 += 2)

	// constructor: allocates state sets for J <= J_max
	explicit TwoBodyStateComplexNlj (size_type J_max);


	TwoBodyStateSetNlj TwoBodySpaceFCI(int Ni_max, int J, ParityGrade g, bool unlike)
{
	TwoBodyStateNlj state;
	TwoBodyStateSetNlj state_set;
	SPOrbitalNlj a_min=SPOrbital(0, HalfInt(1,2));
	SPOrbitalNlj a_max=SPOrbital(Ni_max, Ni_max + HalfInt(1,2));
	

// compare: enum TwoSpeciesStateType {kPP = 0, kPN = 1, kNN = 2};


	////////////////////////////////////////////////////////////////
	// data type definitions
	////////////////////////////////////////////////////////////////

	// Note: Data types are needed publically *only* to enable use
	// of size_type rather than int in member function arguments,
	// to conform pedantically to C++ norms.

 	typedef std::vector< std::vector< std::vector< TwoBodyStateSetNlj > > > TwoBodyStateComplexNljData;

	typedef TwoBodyStateComplexNljData::size_type size_type;

// single-likeness loop

	for (SPOrbitalNlj a1 = a_min; a1 <= a_max; ++a1)
		for (SPOrbitalNlj a2 = (like ? a1 : a_min); a2 <= a_max; ++a2)
			for (int J = IValue(abs( a1.Getj() - a2.Getj() )); J <= IValue( a1.Getj() + a2.Getj() ); ++J) 
			{
				// bypass symmetry-forbidden case |aa;J> (J odd) for like particles
				if ( (like) && (a1 == a2) && ((J%2) == 1) )
					continue;

				// triage state by parity
				int g = (a1.Getl() + a2.Getl()) % 2;
				
				TwoBodyStateNlj state;
				state.a1 = a1;
				state.a2 = a2;
				state.J = J;
				space(state_type,J,g).push_back(state);
			}


////////////////////////////////////////////////////////////////
// generic JP complex -- INC
////////////////////////////////////////////////////////////////

template <typename T>
class TwoBodyComplexJP{

	public:

	////////////////////////////////////////////////////////////////
        // type definitions
	////////////////////////////////////////////////////////////////

	typedef T value_type;
	typedef typename std::vector<T> vector_type;

	////////////////////////////////////////////////////////////////
        // constructors
	////////////////////////////////////////////////////////////////
	
	// default constructor -- synthesized
	// Note: subsequent call to Initialize() required for well-defined behavior

	// copy constructor -- synthesized

	////////////////////////////////////////////////////////////////
        // initializors
	////////////////////////////////////////////////////////////////

	void Initialize (int, int, int);

	////////////////////////////////////////////////////////////////
        // accessors
	////////////////////////////////////////////////////////////////

	// (J,pi) set access
	vector_type& operator() (TwoSpeciesStateType state_type, int J, int g) {
		return state_data_[state_type][J][g];
	};
	const vector_type& operator() (TwoSpeciesStateType state_type, int J, int g) const {
		return state_data_[state_type][J][g];
	};

	// maximal J retrieval
	int JMax(TwoSpeciesStateType state_type) const {return J_max_[state_type];};

	// entry index lookup

	int GetIndex(TwoSpeciesStateType state_type, int J, int g, const value_type& entry) const 
	{
		typename vector_type::iterator it;
		it = lower_bound(state_data_[state_type][J][g].begin(), state_data_[state_type][J][g].end(), key);
		return it;
	};

	////////////////////////////////////////////////////////////////
        // internal storage
	////////////////////////////////////////////////////////////////

	// Note: data arranged by [state_type][J=0,1,...J_max][g=0,1]
	 
	private:
	std::vector< int > J_max_;
 	typedef std::vector< std::vector< std::vector< vector_type > > > InternalDataComplexType;
	InternalDataComplexType state_data_;
};

// basis

typedef TwoBodyComplexJP<TwoBodyStateNlj> TwoBodyNljBasisComplex;
void ConstructTwoBodySpaceFCI(int, TwoBodyNljBasisComplex&);


////////////////////////////////////////////////////////////////
// triangle indexing type -- INC
////////////////////////////////////////////////////////////////



template <typename Tk, typename Ti>
class TriangleIndex {
public:
	////////////////////////////////////////////////////////////////
        // types
	////////////////////////////////////////////////////////////////

	typedef Ti minor_index_type;
	typedef Tk major_index_type;
	typedef std::pair<minor_index_type,minor_index_type> pair_type;

	////////////////////////////////////////////////////////////////
        // constructors
	////////////////////////////////////////////////////////////////

	// default -- disabled
	TriangleIndex();

	// by major index
        explicit TriangleIndex(major_index_type k) {k_ = k;};

	// by minor indices
	TriangleIndex(minor_index_type i, minor_index_type j) {k_ = UpperTriangleIndex(i,j);};	

	////////////////////////////////////////////////////////////////
        // accessors
	////////////////////////////////////////////////////////////////

	// label extraction
	minor_index_type Geti() const; //TODO
	minor_index_type Getj() const; //TODO
	pair_type Getij() const; //TODO

	// accessor for index
	major_index_type GetIndex() const {return k_;};

	////////////////////////////////////////////////////////////////
        // arithmetic operators
	////////////////////////////////////////////////////////////////

	// incrementor
	TriangleIndex& operator ++ (); // prefix ++, increments then returns reference
	TriangleIndex operator ++ (int); // postfix ++, returns value then increments

private:

	// // static label cache
	// static std::vector< pair > label_cache_;
	// static int max_cached_N_;
	// void ExtendCache() const;

	// index for this instance
        //   internal index is always 0-based
	major_index_type k_;
};

////////////////////////////////////////////////////////////////
// triangle indexing type implementation -- INC
////////////////////////////////////////////////////////////////

template <typename Tk, typename Ti>
	inline TriangleIndex<Tk,Ti>& TriangleIndex<Tk,Ti>::operator ++ ()
{
	++k_;
	//ExtendCache();

	return *this;
}

template <typename Tk, typename Ti>
inline TriangleIndex<Tk,Ti> TriangleIndex<Tk,Ti>::operator ++ (int)
{
	TriangleIndex x = *this;
	++k_;
	//ExtendCache();
	return x;
}




// variable base

	// static indexing convention
	static const int base_ = 1;

inline
SPOrbitalNlj::SPOrbitalNlj(int k)
{
	// validate argument
	if (k < base_)
	{
		std::cerr << "SPOrbitalNlj given index below base" << std::endl;
		std::exit(EXIT_FAILURE);
	}

	// store index
	k_ = k - base_;

	// extend cache
	ExtendCache();
}

inline
int SPOrbitalNlj::GetIndex() const
{
	return k_ + base_;
}


void ConstructTwoBodySpaceFCI(int Ni_max, TwoBodyStateComplexNlj& space)
{

	// initialize space parameters
	SPOrbitalNlj a_min=SPOrbitalNlj(0, HalfInt(1,2));
	SPOrbitalNlj a_max=SPOrbitalNlj(Ni_max, Ni_max + HalfInt(1,2));
	int J_max_like = (2 * Ni_max);
	int J_max_unlike = (2 * Ni_max + 1);

	// initialize space container
	space.Initialize(J_max_like, J_max_unlike, J_max_like); 
	
	// iterate over all (i1,i2,J) states

	//   Note: lexicographical is retained, and only canonically
	//   ordered pairs (i1 <= i2) are retained for like particles

	for (SPOrbitalNlj a1 = a_min; a1 <= a_max; ++a1)
		for (SPOrbitalNlj a2 = a_min; a2 <= a_max; ++a2)
			for (int J = IValue(abs( a1.Getj() - a2.Getj() )); J <= IValue( a1.Getj() + a2.Getj() ); ++J) 
			{
				// set up state
				TwoBodyStateNlj state(a1, a2, J);

				// triage state by parity
				int g = (a1.Getl() + a2.Getl()) % 2;
				
				for (TwoSpeciesStateType state_type = kPP; state_type <= kNN; ++state_type)
				{
					// skip noncanonical or symmetry-forbidden states for like species
					bool like = ( (state_type == kPP) || (state_type == kNN) );
					if ( (like) &&  (a1 > a2) )
 						continue;
					if ( (like) &&  (a1 == a2) && ((J%2) == 1) )
						continue;

					space(state_type,J,g).push_back(state);
				}
			}
};

class TwoBodyStateComplexNlj{

	public:

	////////////////////////////////////////////////////////////////
        // constructors
	////////////////////////////////////////////////////////////////
	
	// default constructor -- synthesized
	// Note: subsequent call to Initialize() required for well-defined behavior

	// copy constructor -- synthesized

	////////////////////////////////////////////////////////////////
        // initializors
	////////////////////////////////////////////////////////////////

	void Initialize (int, int, int);

	////////////////////////////////////////////////////////////////
        // accessors
	////////////////////////////////////////////////////////////////

	// (J,pi) set access
	TwoBodyStateSetNlj& operator() (TwoSpeciesStateType state_type, int J, int g) {
		return state_data_[state_type][J][g];
	};
	const TwoBodyStateSetNlj& operator() (TwoSpeciesStateType state_type, int J, int g) const {
		return state_data_[state_type][J][g];
	};

	// maximal J retrieval
	int JMax(TwoSpeciesStateType state_type) const {return J_max_[state_type];};

	////////////////////////////////////////////////////////////////
        // internal storage
	////////////////////////////////////////////////////////////////

	// Note: data arranged by [state_type][J=0,1,...J_max][g=0,1]
	 
	private:
 	typedef std::vector< std::vector< std::vector< TwoBodyStateSetNlj > > > TwoBodyStateComplexNljData;
	std::vector< int > J_max_;
	TwoBodyStateComplexNljData state_data_;
};


	SPOrbitalNlj a_min=SPOrbitalNlj(0, HalfInt(1,2));
	SPOrbitalNlj a_max=SPOrbitalNlj(N1b_max, N1b_max + HalfInt(1,2));


	// (J,pi) set access
	TwoBodyStateSetNlj& operator() (TwoSpeciesStateType state_type, int J, int g) {
		return state_enumeration_[state_type][J][g];
	};
	const TwoBodyStateSetNlj& operator() (TwoSpeciesStateType state_type, int J, int g) const {
		return state_enumeration_[state_type][J][g];
	};



	cout << "****" << endl;
	cout << "(Nlj)^2 FCI indexing tests" << endl;

	int Ni_max = 2;
	TwoBodyStateComplexNlj basis;

	ConstructTwoBodySpaceFCI(Ni_max, basis);
	for (TwoSpeciesStateType state_type = kPP; state_type <= kNN; ++state_type)
	{
		cout << "state_type " << state_type << " J_max " << basis.JMax(state_type) << endl;

		for (int J = 0; J <= basis.JMax(state_type); ++J)
			for (int g = 0; g <= 1; ++g)
			{
				cout << "**** J " << " " << J << " g " << g << " dimension " << basis(state_type,J,g).size() << endl;
				
				for (int i = 0; i < basis(state_type,J,g).size(); ++i)
					cout << basis(state_type,J,g)[i].a1.GetIndex1() << " "
					     << basis(state_type,J,g)[i].a1 << " "
					     << basis(state_type,J,g)[i].a2.GetIndex1() << " "
					     << basis(state_type,J,g)[i].a2 << " "
					     << endl;
			}
	}


// INC -- REPLACED

void TwoBodyBasisNljComplex::Initialize (int J_max_PP, int J_max_PN, int J_max_NN)
{
	// save maximal J values
	J_max_.resize(3);
	J_max_[kPP] = J_max_PP;
	J_max_[kPN] = J_max_PN;
	J_max_[kNN] = J_max_NN;

	// clear old arrays to allow for reuse
	state_enumeration_.clear();
	state_lookup_.clear();

	// set up array framework

	// for each state class
	state_enumeration_.resize(3);
	state_lookup_.resize(3);
	for (TwoSpeciesStateType state_type = kPP; state_type <= kNN; ++state_type)
	{
		// for each J
		state_enumeration_[state_type].resize(J_max_[state_type]+1);
		state_lookup_[state_type].resize(J_max_[state_type]+1);
		for (int J = 0; J <= J_max_[state_type]; ++J)
		{
			// for each parity
			state_enumeration_[state_type][J].resize(2);
			IndexingType indexing_type = (state_type == kPN) ? kSquare : kUpperTriangular;
			state_lookup_[state_type][J].Initialize(indexing_type,-99999);
	}
}



 	typedef std::vector< std::vector< std::vector< IndexLookupArray > > > StateLookupContainerType;
					state_lookup_[state_type][J][g].SetIndex(a1.GetIndex(), a2.GetIndex(), k);


////////////////////////////////////////////////////////////////
// two-index lookup
//   DEPRECATED --to absorb into PairLookupArray
////////////////////////////////////////////////////////////////

class IndexLookupArray {

public:

	////////////////////////////////////////////////////////////////
	// constants
	////////////////////////////////////////////////////////////////
	
	// flag for empty value
	static const int kNullIndex = -1;

	////////////////////////////////////////////////////////////////
	// constructors
	////////////////////////////////////////////////////////////////

	// default -- disabled
	IndexLookupArray();

	// copy -- synthesized

	// by minor index dimension
	IndexLookupArray (IndexingType indexing_type, int minor_dimension){
		indexing_type_ = indexing_type;
		minor_dimension_ = minor_dimension;
		major_dimension_ = IndexingDimension(indexing_type,minor_dimension);
		index_vector_.resize(major_dimension_,kNullIndex);
	};

	////////////////////////////////////////////////////////////////
	// accessors
	////////////////////////////////////////////////////////////////

	int GetIndex (int i, int j) const {
		int k = index_vector_[InternalIndex(i,j)];
		if (k == kNullIndex) {
			std::cerr << "IndexLookupArray::GetIndex: no index value defined for (" << i << "," << j << ")" << std::endl;
			std::exit(EXIT_FAILURE);
		}
		return k;
	};

	void SetIndex (int i, int j, int k) {
		index_vector_[InternalIndex(i,j)]=k;
	};

	////////////////////////////////////////////////////////////////
	// diagnostics
	////////////////////////////////////////////////////////////////

	int MinorDimension () const {return minor_dimension_;};
	int MajorDimension () const {return major_dimension_;};
	int InternalIndex (int i, int j)  const {return LexicographicIndex(indexing_type_,minor_dimension_,i,j);};

private:

	////////////////////////////////////////////////////////////////
	// data
	////////////////////////////////////////////////////////////////

	IndexingType indexing_type_;	
	int minor_dimension_, major_dimension_;
	std::vector<int> index_vector_;

};


 	IndexLookupArrayType index_lookup_array(indexing_type, sp_dimension, kNullIndex);



		if (s1.J != s2.J)
		{
			std::cerr << "TwoBodyMatrixNljTzJP::MatrixElement: inconsistent J values " 
				  << s1.J << " " << s2.J << std::endl;
			std::exit(EXIT_FAILURE);
		}

		if ( (( s1.a1.Getl() + s1.a2.Getl() ) + ( s2.a1.Getl() + s2.a2.Getl() ) ) % 2 != 0)
		{
			std::cerr << "TwoBodyMatrixNljTzJP::MatrixElement: inconsistent parities "
				  << std::endl;
			std::exit(EXIT_FAILURE);
		}


cout << "    " << " state_type " << state_type << " J " << J << " g " << g << endl;


							if (abs(coefficient) > 0.1)
							{
								cout << "<" << s1p.a1 << s1p.a2 << "|" << s2p.a1 << s2p.a2 << ">"
								     << "  <--  "
								     << "<" << s1.a1 << s1.a2 << "|" << s2.a1 << s2.a2 << ">"
								     << " " << coefficient
								     << endl;
								cout << " " << n11p << " " << n12p << " " << n21p << " " << n22p << endl
								     << " " << l11p << " " << l12p << " " << l21p << " " << l22p << endl
								     << " " << n11 << " " << n12 << " " << n21 << " " << n22 << endl
								     << " " << N11 << " " << N12 << " " << N21 << " " << N22 << endl;
								cout << "a22(N22,s2p.a2.Getj()) : " << N22 << " " << s2p.a2.Getj() << " " << a22 << SPOrbitalNlj(N22,s2p.a2.Getj()) << endl;
							}


	header.matrix_size_PP = matrix.MatrixSize(kPP);
	header.matrix_size_PN = matrix.MatrixSize(kPN);
	header.matrix_size_NN = matrix.MatrixSize(kNN);

	// dimension information -- of (state_type,J,g) sector
	int GetDimension (TwoSpeciesStateType state_type, int J, int g) const {
		return matrix_[state_type][J][g].MinorDimension();
	};

	//   Note: these are REDUNDANT to the TwoBodyMatrixNljTzJPDimension functions, which calculate
	//   the size directly by formula, rather than relying on the dimension information returned by
	//   the underlying matrix container.
	int GetSize (TwoSpeciesStateType state_type, int J, int g) const {
		return matrix_[state_type][J][g].MajorDimension();
	};
	int GetSize (TwoSpeciesStateType state_type) const {
		int size = 0;
		for (int J = 0; J <= two_body_basis_.JMax(state_type); ++J)
			for (int g = 0; g <= 1; ++g)
				size += MatrixSize(state_type, J, g);
		return size;
	};
	int GetSize () const {
		int size = 0;
		for (TwoSpeciesStateType state_type = kPP; state_type <= kNN; ++state_type)
			size += MatrixSize(state_type);
		return size;
	};
