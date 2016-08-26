inline SectorNljTzJP::SectorNljTzJP (const TwoBodyBasisNljTzJP& Basis, int k) {
	for (TwoSpeciesStateType state_type = kTwoSpeciesStateTypeBegin; state_type <= kTwoSpeciesStateTypeEnd; ++state_type)
		J_max_.push_back(Basis.JMax[state_type]);

	k_ = k;
}


	// TODO: convert TwoBodyStateNlj from struct to class
	// SPOrbitalNlj a = s1.Geta1();
	// SPOrbitalNlj b = s1.Geta2();
	// SPOrbitalNlj c = s2.Geta1();
	// SPOrbitalNlj d = s2.Geta2();

void ReadTwoBodyMatrixSectorMFDnH2 (std::istream& is, TwoBodyMatrixNljTzJP& matrix, TwoSpeciesStateType state_type, int J, int g)

void ReadTwoBodyMatrixSectorMFDnH2 (std::istream& is, TwoBodyMatrixNljTzJP& matrix, const SectorNljTzJP& sector)
{
	ReadTwoBodyMatrixSectorMFDnH2(is, matrix, sector.GetStateType(), sector.GetJ(), sector.GetGrade());
}

void ReadTwoBodyMatrixBodyMFDnH2 (std::istream& is, TwoBodyMatrixNljTzJP& matrix)
{
	for (SectorNljTzJP sector(output_basis); sector.InRange(); ++sector)
		ReadTwoBodyMatrixSectorMFDnH2(is, matrix, sector);
}

inline double OscillatorScale(double hw) 
{
	return sqrt(kMFDn_mc2 * hw) / kMFDn_hc;
}


// All bad...

					in_matrix_streams.resize[in_matrix_streams.size()+1];
					in_matrix_streams[in_matrix_streams.size()-1].filename = basename + "-r2.dat";
					in_matrix_streams.push_back(InMatrixStream(basename + "-r2.dat"));
					in_matrix_streams.push_back(basename + "-r1r2.dat");

//void ReadTwoBodyMatrixSectorMFDnH2 (std::istream&, TwoBodyMatrixNljTzJP&, TwoSpeciesStateType, int, int);
//void WriteTwoBodyMatrixSectorMFDnH2 (std::ostream&, const TwoBodyMatrixNljTzJP&, TwoSpeciesStateType, int, int);


typedef vector< vector<MixingAmplitude> > MixingAmplitudes;




		// store truncation information
		int J_max_like = std::min( (2 * header.N1b_max), header.N2b_max + 1);
		int J_max_unlike = std::min( (2 * header.N1b_max + 1), header.N2b_max + 1);
		it->J_max_input[kPP] = J_max_like;
		it->J_max_input[kPN] = J_max_unlike;
		it->J_max_input[kNN] = J_max_like;

	int N1b_max, N2b_max;

	enum IOMode { kIn, kOut };
	IOMode io_mode_; 

		stream_->write(static_cast<char*> &header_.format,int_size);
		stream_->write(static_cast<char*> &header_.num_types,int_size); 
		stream_->write(static_cast<char*> &header_.num_P,int_size);    
		stream_->write(static_cast<char*> &header_.num_N,int_size);
		stream_->write(static_cast<char*> &header_.N1b_max,int_size);
		stream_->write(static_cast<char*> &header_.N2b_max,int_size);
		stream_->write(static_cast<char*> &header_.matrix_size_PP,int_size);
		stream_->write(static_cast<char*> &header_.matrix_size_NN,int_size);
		stream_->write(static_cast<char*> &header_.matrix_size_PN,int_size);

		STREAMWRITE(*stream_,header_.format,int); 
		STREAMWRITE(*stream_,header_.num_types,int); 
		STREAMWRITE(*stream_,header_.num_P,int);    
		STREAMWRITE(*stream_,header_.num_N,int);
		STREAMWRITE(*stream_,header_.N1b_max,int);
		STREAMWRITE(*stream_,header_.N2b_max,int);
		STREAMWRITE(*stream_,header_.matrix_PP,int);
		STREAMWRITE(*stream_,header_.matrix_NN,int);
		STREAMWRITE(*stream_,header_.matrix_size_PN,int);

// STREAMWRITE(var,type) casts var to given type in temporary variable and passes to write method
#define STREAMWRITE(stream,var,t) ({t stream_write_tmp = var; (stream).write (static_cast<char*> &stream_write_tmp,sizeof(t));})

	template <typename T>
	void StreamWrite  (const std::basic_ofstream& os, const T value)
	{
		os.write(reinterpret_cast<const char*> (&value),sizeof(T));
	}


		StreamWrite<int>(*stream_,header_.format); 
		StreamWrite<int>(*stream_,header_.num_types); 
		StreamWrite<int>(*stream_,header_.num_P);    
		StreamWrite<int>(*stream_,header_.num_N);
		StreamWrite<int>(*stream_,header_.N1b_max);
		StreamWrite<int>(*stream_,header_.N2b_max);
		StreamWrite<int>(*stream_,header_.matrix_size_PP);
		StreamWrite<int>(*stream_,header_.matrix_size_NN);
		StreamWrite<int>(*stream_,header_.matrix_size_PN);
error: no matching function for call to Å‚Å†è∏StreamWrite(std::basic_ofstream<char, std::char_traits<char> >&, const int&)Å‚Å†èπ


 od -j 40 -t fF



		// for each lp
		// start at lp or lp+1, depending on parity of dl
		int l1 = lp + (dl_ % 2);  
		// end at lp+dl, but max at l_max_-1 or l_max_, depending on parity of dl 
		int l2 = std::min<int>(l_max_,lp+dl_);
		l2 = l2 - ((l2 - l1) - dl) %2;
		// deduce number of sectors to expect for this lp, noting l goes in step 2
		int dim = (l2 - l1)/2 + 1;
		radial_matrices_[lp].resize(dim);
		std::cout << "Reading radial matrices " << lp << " " << l1 << " " << l2 << " " << dim << std::endl;
		for (int i = 0; i < dim; ++i)
			{
				// for each l
				radial_matrices_[lp][i].Initialize(kSquare,n_max_+1,0.);
				for (int np = 0; np <= n_max_; ++np)
					for (int n = 0; n <= n_max_; ++n)
						is >> radial_matrices_[lp][i](np,n);
			}

	// temporary diagnostics
	std::cout << r2k2_radial_matrices.r2_p.GetMatrixElement (0,0,0,0) << " " << std::endl;
	std::cout << r2k2_radial_matrices.r2_p.GetMatrixElement (1,1,0,0) << " " << std::endl;






  class RelativeStateLSTJ
  {
    
  public:

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    // default constructor -- unsupported

    // copy constructor -- synthesized

    // value initialization constructors

    RelativeStateLSTJ (const RelativeSpaceLSTJ& space)
      // Construct state, defaulting to 0th state in space.
    {
      space_ = space;
      index_ = 0;
    };

    RelativeStateLSTJ(const RelativeSpaceLSTJ& space, int index)
      // Construct state, setting to index-th state in space.
    {
      space_ = space;
      index_ = index;
    };


    ////////////////////////////////////////////////////////////////
    // pass-through accessors
    ////////////////////////////////////////////////////////////////
    int L() const {return space_.L();};
    int S() const {return space_.S();};
    int T() const {return space_.T();};
    int J() const {return space_.J();};
    int grade() const {return space_.grade();};

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////
    int index() const {return index_;};

    ////////////////////////////////////////////////////////////////
    // computed indexing
    ////////////////////////////////////////////////////////////////

    int N() const 
    {
      int value = grade() + 2*index();
      return  value;
    };

    int ValidIndex() const
    {
      return index() < space_.Dimension();
    };

    ////////////////////////////////////////////////////////////////
    // iteration
    ////////////////////////////////////////////////////////////////

    RelativeStateLSTJ& operator ++ ()
       // prefix increment operator
      {
	++index_;
	return *this;
      }

    ////////////////////////////////////////////////////////////////
    // private storage
    ////////////////////////////////////////////////////////////////

  private:

    RelativeSpaceLSTJ space_;  // LSTJ space in which state lies
    int index_;   // 0-based index within space
  };



  const RelativeSpaceContainerLSTJ RelativeSpaceEnumerationLSTJ(int Nmax)
  {
    RelativeSpaceContainerLSTJ spaces;

    // iterate over L
    for (int L=0; L<=Nmax; ++L)
      {
	int grade = L%2;

	// iterate over S
	for (int S=0; S<=1; ++S)
	  {
	    int T = (L+S+1)%2;

	    // iterate over J
	    for (int J=abs(L-S); J<=L+S; ++J)
	      {
		RelativeSpaceLSTJ space(L,S,T,J,grade,Nmax);
		assert(space.Dimension()!=0);
		spaces.push_back(space);
	      }
	  }
      }
   
    return spaces;
  }



  const RelativeSpaceContainerLSTJ RelativeSpaceEnumerationLSTJ(int Nmax)
  {
    RelativeSpaceContainerLSTJ spaces;

    // iterate over L
    for (int L=0; L<=Nmax; ++L)
      {
	int g = L%2;

	// iterate over S
	for (int S=0; S<=1; ++S)
	  {
	    int T = (L+S+1)%2;

	    // iterate over J
	    for (int J=abs(L-S); J<=L+S; ++J)
	      {
		// downshift Nmax if necessary to match parity of subspace
		int Nmax_for_space = Nmax - (Nmax-g)%2;

		// std::cout 
		//        << std::setw(3) << L 
		// 	  << std::setw(3) << S 
		// 	  << std::setw(3) << T 
		// 	  << std::setw(3) << J 
		// 	  << std::setw(3) << g 
		// 	  << std::setw(3) << Nmax_for_space 
		// 	  << std::endl;

		RelativeSpaceLSTJ space(L,S,T,J,g,Nmax_for_space);
		assert(space.Dimension()!=0);
		spaces.push_back(space);
	      }
	  }
      }
   
    return spaces;
  }

  const SectorContainer
  RelativeSectorEnumerationLSTJ(RelativeSpaceContainerLSTJ& relative_spaces, int J0, int g0)
  // Enumerates sector pairs connected by an operator of given tensorial and parity character.
  //
  // Sectors are enumerated in "both directions", i.e., no asumption of hermiticity is made.
  // TODO: genericize to any space with LSTJ labels??  
  {
    SectorContainer relative_sectors;
    for (int s2=0; s2<relative_spaces.size(); ++s2)
      for (int s1=0; s1<relative_spaces.size(); ++s1)
	{
	  // verify triangularity and p 
	  int J2 = relative_spaces[s2].J();
	  int J1 = relative_spaces[s1].J();
	  int g2 = relative_spaces[s2].g();
	  int g1 = relative_spaces[s1].g();

	  if ( AllowedTriangle(J2,J0,J1) && ((g2+g0+g1)%2==0))
	    relative_sectors.push_back(std::make_pair(s2,s1));
	}

    return relative_sectors;
  }


  ////////////////////////////////////////////////////////////////
  // enumeration of subspaces
  ////////////////////////////////////////////////////////////////


  // container type for enumerating spaces
  typedef std::vector<RelativeSpaceLSTJ> RelativeSpaceContainerLSTJ;

  const RelativeSpaceContainerLSTJ RelativeSpaceEnumerationLSTJ(int Nmax);
  // Enumerates all relative LSTJ spaces of given dimension up to a
  // given Nmax cutoff.
  //
  // All parities, isospins, etc., are included.
  // The maximal L is determined by Nmax.
  //
  // Iteration order:
  //    -- increasing L (L=0,1,...,Nmax)
  //    -- [this determines P]
  //    -- increasing S (S=0,1)
  //    -- [this determines T]
  //    -- increasing J
  // subject to:
  //    -- triangularity of L,S,T
  //    -- [assertion: sanity check on nonzero dimension for space]
  // Note that ordering of spaces is therefore lexicographic by (L,S,J).




// Moshinsky bracket -- seed 
double SeedMoshinskyBracket (int, int, int, int, int, int, int);

// Moshinsky bracket -- special case N_CM = 0
//   subsumed under Moshinsky Bracket below
double CMMoshinskyBracket (int, int, int, int, int);
double CMMoshinskyBracket (const TwoBodyStateNl&);


////////////////////////////////////////////////////////////////
// CMMoshinskyBracket
////////////////////////////////////////////////////////////////

// now superfluous special case of MoshinskyBracket

double CMMoshinskyBracket (
	int n1, int l1, int n2, int l2,
	int Lambda
	)
{
	// bra labels fixed from ket
	int n1_dot = n1 + n2 + (l1+l2-Lambda)/2;
	int l1_dot = Lambda; 

	// angular-momentum forbidden case
	if ((l1 + l2 + Lambda)%2 != 0)
		return 0.;

	// recurse to value
	double value;

	if ( (n1 == 0) && (n2 == 0) )
		// case: at seed value on RHS
		value = SeedMoshinskyBracket(n1_dot,l1_dot,0,0,l1,l2,Lambda);
	else if (n1 == 0)
		// case: cannot recurse n1 further but can recurse n2
		value = ParitySign(Lambda) 
			* CMMoshinskyBracket(n2,l2,0,l1,Lambda);
	else
		// otherwise: recurse n1
		value = 1/2.
			* sqrt(
				( n1_dot * (n1_dot + l1_dot + 1/2.) )
				/ ( n1 * (n1 + l1 + 1/2.) )
				)
			* CMMoshinskyBracket(n1-1,l1,n2,l2,Lambda);
	
	return value;
}

double CMMoshinskyBracket (const TwoBodyStateNl& state)
{
	const int n1 = state.a1.Getn();
	const int n2 = state.a2.Getn();
	const int l1 = state.a1.Getl();
	const int l2 = state.a2.Getl();
	const int L = state.L;

	return CMMoshinskyBracket (n1,l1,n2,l2,L);
}



    typedef std::vector<TwoBodySpaceLSJT> SpaceContainer;


    // typedef std::pair<int,int> KeyType;
    // KeyType Key() const
    // {
    //   return KeyType(index2,index1);
    // }
    //bool operator < (Sector& s2) const
    //{
    //  return std::make_tuple(index2,index1) < std::make_tuple(s2.index2,s2.index1);
    //}

    ////////////////////////////////////////////////////////////////
    //  composite label structures
    ////////////////////////////////////////////////////////////////

    struct SubspaceLabels
    {
      SubspaceLabels (){};  // not sure why we need this (gcc 4.5.3), since should be synthesized
      SubspaceLabels(int L_, int S_, int J_, int T_, int g_)
      {
	L = L_; S = S_; J = J_; T = T_; g = g_;
      };
      typedef std::tuple<int,int,int,int> KeyType;
      KeyType Key() const
      {
	return KeyType(L,S,J,g);
      }

      //bool operator < (SubspaceLabels& s2) const
      //{
      //	return std::make_tuple(L,S,T,g) < std::make_tuple(s2.L,s2.S,s2.T,s2.g);
      //};

      int L, S, J, T, g;
    };

    struct StateLabels
    {
      StateLabels(int N_)
      {
	N = N_;
      };
      int N;
    };

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////
    int L() const {return labels_.L;}
    int S() const {return labels_.S;}
    int J() const {return labels_.J;}
    int T() const {return labels_.T;}
    int g() const {return labels_.g;}
    int Nmax() const {return Nmax_;}


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  // logging helper function for debugging -- TEMPORARY
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  // // tuple output
  // //  template <typename T>
  // //    void WriteMultiplet(std::ostream& os, typename std::initializer_list<T>& vals)
  // //    inline void WriteMultipletInt(std::ostream& os, std::initializer_list<int> vals)
  // template <typename T>
  //   void WriteMultiplet(std::ostream& os, T& vals)
  //   // Generalization of "print" function on Josuttis, "The C++
  //   // Standard Library", 2ed, page 16.  Probably to change into
  //   // string conversion function.
  //   {
  //     char delimiter_left[] = "(";
  //     char delimiter_middle[] = ",";
  //     char delimiter_right[] = ")";
  // 
  //     os << delimiter_left;
  // 
  //     bool at_start = true;
  //     for(auto p=vals.begin(); p != vals.end(); ++p) // make into range-based?
  // 	{
  // 	  if (!at_start)
  // 	    os << delimiter_middle;
  // 	  at_start = false;
  // 	  os << (*p);
  // 	}
  //     os << delimiter_right;
  // 
  //   }
  // //  WriteMultiplet<int>(std::cout,{1,2});
  // //	  WriteMultiplet< std::tuple<int,int> >(std::cout,std::make_tuple(1,2));
  // // error: no matching function for call to 'WriteMultiplet(std::ostream&, std::tuple<int, int>)'


    // special constructor -- quantum number lookup
    //
    // Note: may be supplanted by generic lookup from BaseState

    RelativeStateLSJT(const RelativeSubspaceLSJT& space, const StateLabels& labels)
      // Construct state, by reverse-lookup on labels (N).
      //
      // Example usage:
      //   RelativeStateLSJT relative_state = RelativeStateLSJT(RelativeStateLSJT::StateLabels(Nr));

      {
	space_ptr_ = &space;
	int index = (labels.N - g())/2;
	index_ = index;
      }



    // manual version
    // int N() const 
    // {
    //   int value = g() + 2*index();
    //   return  value;
    // }


  // logging output -- TEMPORARY

  inline
    std::ostream& operator<< (std::ostream& os, const RelativeSubspaceLSJT& subspace)
    {
      // FUTURE: replace with some form of standardized tuple output
      // on initializer list
      os << "(" 
	 << subspace.L() << "," << subspace.S() << "," << subspace.T() << "," << subspace.J() << "," << subspace.g() 
	 << ")";
	
      return os;
    }


  // relative operator construction

  void SetRelativeOperatorToZero(
				 const RelativeSpaceLSJT& space,
				 const RelativeSectorsLSJT& sectors,
				 SectorMatrices& matrices
				 );
  
  void SetRelativeOperatorToIdentity(
				     const RelativeSpaceLSJT& space,
				     const RelativeSectorsLSJT& sectors,
				     SectorMatrices& matrices
				     );






  void SetRelativeOperatorToZero(
				 const RelativeSpaceLSJT& space,
				 const RelativeSectorsLSJT& sectors,
				 SectorMatrices& matrices
				 )
  {
    matrices.clear();
    for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
      {
	const Sector& sector = sectors.GetSector(sector_index);
	const RelativeSubspaceLSJT& subspace2 = space.GetSubspace(sector.index2());
	const RelativeSubspaceLSJT& subspace1 = space.GetSubspace(sector.index1());

	Eigen::MatrixXd sector_matrix;
	sector_matrix = Eigen::MatrixXd::Zero(subspace2.Dimension(),subspace1.Dimension());
	matrices.push_back(sector_matrix);
      }
  }

  void SetRelativeOperatorToIdentity(
				     const RelativeSpaceLSJT& space,
				     const RelativeSectorsLSJT& sectors,
				     SectorMatrices& matrices
				     )
  {
    matrices.clear();
    for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
      {
	const Sector& sector = sectors.GetSector(sector_index);
	const RelativeSubspaceLSJT& subspace2 = space.GetSubspace(sector.index2());
	const RelativeSubspaceLSJT& subspace1 = space.GetSubspace(sector.index1());

	Eigen::MatrixXd sector_matrix;
	if (sector.index2() == sector.index1())
	  {
	    int J = subspace1.J();
	    sector_matrix = Hat(J) * Eigen::MatrixXd::Identity(subspace2.Dimension(),subspace1.Dimension());
	  }
	else
	  {
	    sector_matrix = Eigen::MatrixXd::Zero(subspace2.Dimension(),subspace1.Dimension());
	  }
	matrices.push_back(sector_matrix);
      }
  }


    // antisymmetry
    valid &= ((L()+S()+T())%2 == 1);

