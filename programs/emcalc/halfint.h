/****************************************************************
  halfint.h                       

  Defines arithmetic type HalfInt which stores integer or half-integer values,
  as needed, e.g., for angular momentum quantum numbers.
                                  
  Created by Ke Cai on 6/18/10.   

  10/04/10 stylistic changes (vh): 
    (1) removed "using namespace std" in Halfint.h and added to Halfint.cpp
    (2) copy constructor & copy assignment operator don't need to be defined : 
        the default copy constructor & assignment operator provide memberwise
        copy & assignment
    (3) changed arguments passed by value to arguments passed by const reference
        for various functions
    (4) getnum() and display() made const member functions
  11/27/10+ (mac) full overhaul of functionality plus stylistic changes:
    Goals:
      - Make HalfInt behave as much like a built-in arithmetic type as possible
      - Eliminate redundant features
    - Remove function hat() in favor of global function hat()
    - Remove member function absolute() in favor of global function abs()
    - Remove set() member function as redundant to assignment operator
    - Remove maximum() member function as redundant to standard algorithm max()
    - Implement full set of unary operators +, -, ++, --
    - Implement full set of arithmetic assignment operators += and -=
    - Implement stream output operator << for HalfInt
    - Change binary arithmetic operators + and - to global functions
        for aesthetics and to allow use with automatic integer->HalfInt 
        conversion, and reimplement in terms of += and -=
    - Remove static constant definitions zero and one, since are redundant to
      integers 0 and 1 (due to automatic conversions from
      int to HalfInt)
    - Rename data member numerator to twice_value_, with accessor twice_value().
      The new choice of name is since "numerator" is ambiguous.  It assumes
      implicitly that we always mean numeratorrelative to
      a "denominator" of 2 (not 1!)
    - Inline all member functions
    - Pull out angular momentum related functions to angular_momentum.h
  12/8/10 (mac):
    - Rename value accessors to TwiceValue() and DValue()
  12/15/10 (mac):
    - Add value accessor IValue() which gives *truncated* integer value
  12/24/10 (mac):
    - Eliminate member function display() as redundant to stream output 
    - Fix sign error just introduced to abs(HalfInt)
  2/16/11 (mac):
    - Add nonmember versions of accessor functions
  2/21/11 (mac):
    - Replace exception generation with exit()
  2/22/11 (mac):
    - Fix operator<< to return output stream

****************************************************************/

#ifndef halfint_h
#define halfint_h

#include <cstdlib>
#include <iostream>

class HalfInt{

public:
	////////////////////////////////////////////////////////////////
	// constructors
	////////////////////////////////////////////////////////////////

	// note: Copy constructor is synthesized copy constructor 
	// since only data member needs copying

	// default constructor: initializes value to zero
	HalfInt () : twice_value_(0) {}; 

	// conversion from integer
	// EX: HalfInt(1) constructs 2/2, HalfInt(2) constructs 4/2
	HalfInt (int);

	// construct from numerator and denominator
	// EX: HalfInt(1,2) constructs 1/2, HalfInt(2, 2) constructs 1,
        //     HalfInt(2,1) constructs 4/2
	HalfInt (int, int);

	////////////////////////////////////////////////////////////////
	// accessors
	////////////////////////////////////////////////////////////////

	// return twice value as type int
	int TwiceValue() const {return twice_value_;};
	
	// return value as type int, by truncation
	int IValue() const {
		if ((twice_value_ % 2) != 0)
		{
			std::cerr << "HalfInt::IValue called for nonintegral HalfInt value" << std::endl;
			std::exit(EXIT_FAILURE);
		}
		return twice_value_/2;
	};

	// return value as type double
	double DValue() const {return static_cast<double>(twice_value_)/2;};

	// note on DValue: Defining an operator double(), for
	// automatic conversion to double, naively seems desirable.
	// However, this would lead to ambiguous overloads, in
	// combination of automatic conversion from int to HalfInt,
	// e.g., 1 + 1/2 -> 2/2 + 1/2 -> 3/2 vs. 1 + 1/2 -> 1.0 + 0.5
	// -> 1.5.

	////////////////////////////////////////////////////////////////
        // arithmetic assignment operators
	////////////////////////////////////////////////////////////////

	// note: operator = is synthesized assignment operator 

        // operators return reference to allow for chained (a+=b)+=c structures
	HalfInt& operator += (const HalfInt&);
	HalfInt& operator -= (const HalfInt&);

	////////////////////////////////////////////////////////////////
	// unary arithmetic operators
	////////////////////////////////////////////////////////////////

	HalfInt operator + () const;
	HalfInt operator - () const;

	HalfInt& operator ++ (); // prefix ++, increments then returns reference
	HalfInt operator ++ (int); // postfix ++, returns value then increments
	HalfInt& operator -- (); // prefix --, increments then returns reference
	HalfInt operator -- (int); // postfix --, returns value then increments

private:
	// contains one member, which stores twice the actual value, since this is
        // guaranteed to be an integer.
	int twice_value_;

};

////////////////////////////////////////////////////////////////
// constructors
////////////////////////////////////////////////////////////////

inline HalfInt::HalfInt (int value)
{
	twice_value_ = 2*value;
}

inline HalfInt::HalfInt (int numerator, int denominator)
{
	if (!((denominator==1)||(denominator==2)))
	{
		std::cerr << "HalfInt constructed with denominator not 1 or 2" << std::endl;
		std::exit(EXIT_FAILURE);
	}
	twice_value_ = (2/denominator)*numerator;
}

////////////////////////////////////////////////////////////////
// arithmetic assignment operators
////////////////////////////////////////////////////////////////

inline HalfInt& HalfInt::operator += (const HalfInt& b)
{
	twice_value_ += b.twice_value_;
	return *this;
}

inline HalfInt& HalfInt::operator -= (const HalfInt& b)
{
	twice_value_ -= b.twice_value_;
	return *this;
}

////////////////////////////////////////////////////////////////
// nonmember accessor functions
////////////////////////////////////////////////////////////////


inline int TwiceValue(const HalfInt& h) 
{
	return (h.TwiceValue());
}

inline int IValue(const HalfInt& h) 
{
	return (h.IValue());
}

inline double DValue(const HalfInt& h) 
{
	return (h.DValue());
}

////////////////////////////////////////////////////////////////
// unary arithmetic operators
////////////////////////////////////////////////////////////////

inline HalfInt HalfInt::operator + () const 
{
	return *this;
};

inline HalfInt HalfInt::operator - () const 
{
	HalfInt negative;
	negative.twice_value_ = -twice_value_;
	return negative;
};

inline HalfInt& HalfInt::operator ++ ()
{
	twice_value_ += 2;	
	return *this;
}

inline HalfInt HalfInt::operator ++ (int)
{
	HalfInt h = *this;
	twice_value_ += 2;
	return h;
}

inline HalfInt& HalfInt::operator -- ()
{
	twice_value_ -= 2;	
	return *this;
}

inline HalfInt HalfInt::operator -- (int)
{
	HalfInt h = *this;
	twice_value_ -= 2;
	return h;
}


////////////////////////////////////////////////////////////////
// binary arithmetic operators
////////////////////////////////////////////////////////////////

inline HalfInt operator + (const HalfInt& a, const HalfInt& b) 
{
	HalfInt sum(a);
	sum += b;
	return sum;
}

inline HalfInt operator - (const HalfInt& a, const HalfInt& b) 
{
	HalfInt sum(a);
	sum -=b;
	return sum;
}

////////////////////////////////////////////////////////////////
// relational operators
////////////////////////////////////////////////////////////////

inline bool operator == (const HalfInt& h1, const HalfInt& h2)
{
	return h1.TwiceValue() == h2.TwiceValue();
}

inline bool operator < (const HalfInt& h1, const HalfInt& h2)
{
	return h1.TwiceValue() < h2.TwiceValue();
}

inline bool operator > (const HalfInt& h1, const HalfInt& h2)
{
	return !((h1==h2)||(h1<h2));
}

inline bool operator >= (const HalfInt& h1, const HalfInt& h2)
{
	return !(h1 < h2);
}

inline bool operator <= (const HalfInt& h1, const HalfInt& h2)
{
	return ((h1 < h2)||(h1 == h2));
}

inline bool operator != (const HalfInt& h1, const HalfInt& h2)
{
	return !(h1 == h2);
}

////////////////////////////////////////////////////////////////
// arithmetic functions
////////////////////////////////////////////////////////////////

inline HalfInt abs(const HalfInt& h) 
{
	if (h.TwiceValue() < 0)
		return -h;
	else 
		return h;
}

////////////////////////////////////////////////////////////////
// stream output
////////////////////////////////////////////////////////////////

// textual output to stream
// EX: HalfInt(3) or HalfInt(6,2) -> "3", HalfInt(3,2) -> "3/2"
std::ostream& operator<< (std::ostream&, const HalfInt&);

#endif
