/****************************************************************
  halfint.cpp                       
                                    
  Created by Ke Cai on 6/18/10.    
                                    
  Last changes by mac 2/22/11       
                                    
****************************************************************/

#include "halfint.h"

using namespace std;

////////////////////////////////////////////////////////////////
// stream output
////////////////////////////////////////////////////////////////

ostream& operator<< (ostream& os, const HalfInt& h)
{
	int num = h.TwiceValue();
	if (num % 2 == 0)
		os << num/2;
	else
		os << num << "/" << "2";

	return os;
}

