
	
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


