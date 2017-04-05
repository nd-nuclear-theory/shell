
        ////////////////////////////////
	// simple loop by shell cutoff
	////////////////////////////////

	const int N_smax = 30;
	
	double dummy_sum = 0;
	int bracket_count = 0;
	for(int N1 = 0; N1 <= N_smax; ++N1) 
		for(int N2 = 0; N2 <= N1; ++N2) 
			for (int n1 = 0; 2*n1 <= N1; ++n1)
				for (int n2 = 0; 2*n2 <= N2; ++n2)
				{
					int l1 =  N1 - 2*n1;
					int l2 =  N2 - 2*n2;
					int Lambda_min = abs(l1-l2);
					int Lambda_max = abs(l1+l2);				
					for (int Lambda = Lambda_min; Lambda <= Lambda_max; Lambda += 2)
					{
						++bracket_count;
						double bracket = CMMoshinskyBracket(n1,l1,n2,l2,Lambda);
						if (abs(bracket) > 1)
						{
							cout << bracket << endl;
							exit(0);
						}
						dummy_sum += abs(bracket);
					}
				}
	cout << bracket_count << " " << dummy_sum << endl;
					
