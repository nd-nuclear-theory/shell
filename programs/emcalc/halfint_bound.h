/****************************************************************
  halfint_bound.h                 
                                  
  Created by Ke Cai on 6/25/10.   
                                  
  10/04/10 stylistic changes (vh): 
    (1) arguments passed by constant reference instead of by value for
        constructor, overloaded operator, and triangleBound()
    (2) copy constructor doesn't need to be defined :
        the default copy constructor & assignment operator provide memberwise
        copy & assignment
    (3) "get" member functions and print member function made const member function
  2/16/11 (mac): 
     - file renamed to halfint_bound
     - member ShowBound() replaced with stream output overload
     - move out function TriangleBound() since specific to angular momentum
     - make accessors const
  3/4/11 (mac):
     - fix missing return value for operator<<

****************************************************************/


#ifndef halfint_bound_h
#define halfint_bound_h

#include <iostream>

#include "halfint.h"

//HalfIntBound objects contain an upper bound and a lower bound, both of which are halfints
class HalfIntBound{
public:
	//CONSTRUCTORS
	//Default constructor: returns (-999/2, -998/2).
	HalfIntBound();

	//Constructs a halfBound object with the smaller argument as the lower bound and the greater argument as the upper bound
	HalfIntBound(const HalfInt&, const HalfInt&);


	//Assignment operator: assigns the upper bound and lower bound of the argument to the upper bound and lower bound of the object, respectively
	void operator = (const HalfIntBound&);

	//Intersection: take the intersection of two bounds; lower bound of the return value is the larger one between the lower bounds of object and argument;
	//upper bound of the return value is the smaller one between the upper bounds of object and argument.
	//If the lower bound of the returned value is greater than its upper bound, an error message will be printed and the program will be terminated
	HalfIntBound operator * (const HalfIntBound&);

	//Returns the lower bound of a HalfIntBound object
	HalfInt lower() const {return lower_;};

	//Returns the upper bound of a HalfIntBound object
	HalfInt upper() const {return upper_;};

private:
	//Two members.
	HalfInt lower_;
	HalfInt upper_;

};

////////////////////////////////
// Stream output
////////////////////////////////

// Textual output to stream
// EX: HalfIntBound(HalfInt(1,2),HalfInt(3,2)) -> "(1/2,3/2)"
std::ostream& operator<< (std::ostream&, const HalfIntBound&);

#endif
