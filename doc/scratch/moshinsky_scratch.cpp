
	
	// cerr << "<" <<  n1_dot << " " << l1_dot  << " -- -- | "
	//      << n1  << " " << l1 << " " << n2 << " " << l2 << " ; " << Lambda << ">" 
	//      << endl;
	

	// else if (n1_dot == 0)
	// 	// case: at seed value on LHS
	// {
	// 	cerr << "case 2" << endl;
	// 	value = ParitySign(l2+Lambda)
	// 		* SeedMoshinskyBracket(n1,l1,n2,l2,l1_dot,0,Lambda);
	// }


////////////////////////////////////////////////////////////////

// as raw labels

struct TwoBodyLabelsNl {
	int N1, l1, N2, l2, L;
};

typedef std::vector<TwoBodyLabelsNl> TwoBodyLabelSetNl;

TwoBodyLabelSetNl TwoBodySpaceNL(int, int);


				current_labels.N1 = N1;
				current_labels.N2 = N2;
				current_labels.l1 = l1;
				current_labels.l2 = l2;


////////////////////////////////////////////////////////////////
// version with TwoBodyStateNl indexing -- DEPRECATED/DISABLED
////////////////////////////////////////////////////////////////

//  double MoshinskyBracket (const TwoBodyStateNl& state_dot, const TwoBodyStateNl& state)
//  {
//    // extract parameters in Moshinsky HO book notation
//    const int n1_dot = state_dot.a1.Getn();
//    const int l1_dot = state_dot.a1.Getl();
//    const int n2_dot = state_dot.a2.Getn();
//    const int l2_dot = state_dot.a2.Getl();
//    const int L_dot = state_dot.L;
//
//    const int n1 = state.a1.Getn();
//    const int l1 = state.a1.Getl();
//    const int n2 = state.a2.Getn();
//    const int l2 = state.a2.Getl();
//    const int L = state.L;
//
//    // check equality of L values
//    if ( L_dot != L )
//      return 0.;
//
//    // invoke basic version of function
//    return MoshinskyBracket (n1_dot,l1_dot,n2_dot,l2_dot,n1,l1,n2,l2,L);
//  }

