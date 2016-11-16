/****************************************************************
  halfint_bound.cpp                       
                                    
  Created by Ke Cai on 6/25/10.    
                                    
  Last changes by mac 3/4/11.
                                    
****************************************************************/

#include "halfint_bound.h"

using namespace std;

HalfIntBound :: HalfIntBound(){
	upper_ = HalfInt(-998, 2);
	lower_ = HalfInt(-999, 2);
}

HalfIntBound :: HalfIntBound(const HalfInt& h1, const HalfInt& h2){
	if (h1 >= h2) {
		upper_ = h1;
		lower_ = h2;
	}
	else {
		upper_ = h2;
		lower_ = h1;
	}
}


void HalfIntBound:: operator = (const HalfIntBound& hb){
	this->lower_ = hb.lower_;
	this->upper_ = hb.upper_;
}

HalfIntBound HalfIntBound::operator * (const HalfIntBound& hb){
	HalfIntBound result;
	//HalfIntBound nullBound(zero,zero);
	if ((this->lower_)>= hb.lower_) {
		result.lower_ = this->lower_;
	}
	else {
		result.lower_ = hb.lower_;
	}

	if ((this->upper_)>= hb.upper_) {
		result.upper_ = hb.upper_;
	}
	else {
		result.upper_ = this->upper_;
	}

	return result;
}


////////////////////////////////
// Stream output
////////////////////////////////

ostream& operator<< (ostream& os, const HalfIntBound& hb)
{
	cout<< "(" << hb.lower() << "," << hb.upper() << ")";
	return os;
}

