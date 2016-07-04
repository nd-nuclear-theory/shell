/****************************************************************

  vector_tuple.h

  Template library for fixed-length tuples of arithmetic types, 
  based on storage in a vector.  Supported arithmetic operations:
    -- unary + -
    -- entrywise + - * /
    -- lexicographical comparison

  Created by Mark A. Caprio 11/27/10
  Inspired by Tomas Dytrych CTuple.h, but using vector rather than array.
                                  
  12/15/10 (mac):
  -- add option template parameter B=0 as indexing base for VectorTuple.
        VectorTuple<T,N,0> for C-style 0-based indexing (default)
        VectorTuple<T,N,1> for mathematical 1-based indexing
  12/24/10 (mac):
  -- add ostream output of VectorTuple, with configurable delimiters
  2/17/11 (mac):
  -- fix size()
  2/23/11 (mac):
  -- fix namespace issues
                                  
****************************************************************/


#ifndef VECTOR_TUPLE_H_
#define VECTOR_TUPLE_H_

#include <cstddef>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

template<class T, size_t N, size_t B=0>
class VectorTuple{

public:

	////////////////////////////////
	// type definitions
	////////////////////////////////

        // standard vector definitions (subset)
        typedef typename std::vector<T>::iterator iterator;
	typedef typename std::vector<T>::const_iterator const_iterator;
	typedef typename std::vector<T>::size_type size_type;
	typedef typename std::vector<T>::difference_type difference_type;
	typedef typename std::vector<T>::value_type value_type;
	typedef typename std::vector<T>::reverse_iterator reverse_iterator;
        typedef typename std::vector<T>::const_reverse_iterator const_reverse_iterator;

	// related vector type
	typedef typename std::vector<T> vector_type;

	////////////////////////////////
	// constructors
	////////////////////////////////

	// default -- default-initialize vector of correct length
	VectorTuple() : values_(N){};

        // copy -- synthesized constructor copies value

	// initialization by fill with constant value
	explicit VectorTuple(const T& a) : values_(N,a) {};
	
	// conversion from vector 
	VectorTuple(const vector_type& v) : values_(v) {};

	////////////////////////////////
        // accessors
	////////////////////////////////

	// entry access by indexing (nonconst and const references)
	T& operator[](size_type i) {return values_[i-B];};
	const T& operator[](size_type i) const {return values_[i-B];};

	// copy of entries as a vector
	// std::vector<T> GetVector() const {return values_;};

	// iterators
	iterator begin() {return values_.begin();};
	iterator end() {return values_.end();};
	const_iterator begin() const {return values_.begin();};
	const_iterator end() const {return values_.end();};
	reverse_iterator rbegin() {return values_.rbegin();};
	reverse_iterator rend() {return values_.rend();};
	const_reverse_iterator rbegin() const {return values_.rbegin();};
	const_reverse_iterator rend() const {return values_.rend();};

	// size
        size_type size() const {return values_.size();};

	////////////////////////////////
        // arithmetic assignment operators
	////////////////////////////////

	// assignment --  synthesized assignment copies value

	VectorTuple& operator += (const VectorTuple&);
	VectorTuple& operator -= (const VectorTuple&);
	VectorTuple& operator *= (const VectorTuple&);
	VectorTuple& operator /= (const VectorTuple&);
	
	////////////////////////////////
	// unary arithmetic operators
	////////////////////////////////

	VectorTuple operator + () const;
	VectorTuple operator - () const;

	////////////////////////////////
	// ostream output
	////////////////////////////////

        // output operator -- friend declaration for access to delimiters
        template <class TX, size_t NX, size_t BX>
        friend std::ostream& operator<< (std::ostream&, const VectorTuple<TX,NX,BX>&);

        // configuring delimiter strings
        //   static member function sets delimiters for *all* VectorTuple 
        //   instances with the given template parameters
        // EX: VectorTuple<int,3,1>::SetDelimiters("< "," "," >"); 
        static void SetDelimiters(const std::string&, const std::string&, const std::string&);

private:

        // ostream delimiters
        static std::string delimiter_left_;
        static std::string delimiter_middle_;
        static std::string delimiter_right_;

        // instance data storage
        vector_type values_;

};

////////////////////////////////
// arithmetic utilities
////////////////////////////////

// note: didn't get <functional> functors to work in this role

template<class T> T BinaryPlus(const T& a, const T& b)
{
	return a+b;
}
template<class T> T BinaryMinus(const T& a, const T& b)
{
	return a-b;
}
template<class T> T BinaryTimes(const T& a, const T& b)
{
	return a*b;
}
template<class T> T BinaryDivide(const T& a, const T& b)
{
	return a/b;
}
template<class T> T UnaryMinus(const T& a)
{
	return -a;
}

////////////////////////////////
// arithmetic assignment operators
////////////////////////////////

template <class T, size_t N, size_t B>
VectorTuple<T,N,B>& VectorTuple<T,N,B>::operator += (const VectorTuple& b)
{
	std::transform(values_.begin(),values_.end(),b.values_.begin(),values_.begin(),BinaryPlus<T>);
	return *this;
}

template <class T, size_t N, size_t B>
VectorTuple<T,N,B>& VectorTuple<T,N,B>::operator -= (const VectorTuple& b)
{
	std::transform(values_.begin(),values_.end(),b.values_.begin(),values_.begin(),BinaryMinus<T>);
	return *this;
}

template <class T, size_t N, size_t B>
VectorTuple<T,N,B>& VectorTuple<T,N,B>::operator *= (const VectorTuple& b)
{
	std::transform(values_.begin(),values_.end(),b.values_.begin(),values_.begin(),BinaryTimes<T>);
	return *this;
}

template <class T, size_t N, size_t B>
VectorTuple<T,N,B>& VectorTuple<T,N,B>::operator /= (const VectorTuple& b)
{
	std::transform(values_.begin(),values_.end(),b.values_.begin(),values_.begin(),BinaryDivide<T>);
	return *this;
}

////////////////////////////////
// unary arithmetic operators
////////////////////////////////

template <class T, size_t N, size_t B>
VectorTuple<T,N,B> VectorTuple<T,N,B>::operator + () const 
{
	return *this;
};

template <class T, size_t N, size_t B>
VectorTuple<T,N,B> VectorTuple<T,N,B>::operator - () const 
{
	VectorTuple<T,N,B> c;
	std::transform(values_.begin(),values_.end(),c.values_.begin(),UnaryMinus<T>);
	return c;
};


////////////////////////////////
// binary arithmetic operators
////////////////////////////////

template <class T, size_t N, size_t B>
VectorTuple<T,N,B> operator + (const VectorTuple<T,N,B>& a, const VectorTuple<T,N,B>& b) 
{
	VectorTuple<T,N,B> c(a);
	c += b;
	return c;
}

template <class T, size_t N, size_t B>
VectorTuple<T,N,B> operator - (const VectorTuple<T,N,B>& a, const VectorTuple<T,N,B>& b) 
{
	VectorTuple<T,N,B> c(a);
	c -=b;
	return c;
}

template <class T, size_t N, size_t B>
VectorTuple<T,N,B> operator * (const VectorTuple<T,N,B>& a, const VectorTuple<T,N,B>& b) 
{
	VectorTuple<T,N,B> c(a);
	c *= b;
	return c;
}

template <class T, size_t N, size_t B>
VectorTuple<T,N,B> operator / (const VectorTuple<T,N,B>& a, const VectorTuple<T,N,B>& b) 
{
	VectorTuple<T,N,B> c(a);
	c /= b;
	return c;
}

////////////////////////////////
// relational operators
////////////////////////////////

template <class T, size_t N, size_t B>
bool operator == (const VectorTuple<T,N,B>& v1, const VectorTuple<T,N,B>& v2)
{
	return equal(v1.begin(),v1.end(),v2.begin());
}

template <class T, size_t N, size_t B>
bool operator < (const VectorTuple<T,N,B>& v1, const VectorTuple<T,N,B>& v2)
{
	return lexicographical_compare(v1.begin(),v1.end(),v2.begin(),v2.end());
}

template <class T, size_t N, size_t B>
bool operator > (const VectorTuple<T,N,B>& v1, const VectorTuple<T,N,B>& v2)
{
	return !((v1==v2)||(v1<v2));
}

template <class T, size_t N, size_t B>
bool operator >= (const VectorTuple<T,N,B>& v1, const VectorTuple<T,N,B>& v2)
{
	return !(v1 < v2);
}

template <class T, size_t N, size_t B>
bool operator <= (const VectorTuple<T,N,B>& v1, const VectorTuple<T,N,B>& v2)
{
	return ((v1 < v2)||(v1 == v2));
}

template <class T, size_t N, size_t B>
bool operator != (const VectorTuple<T,N,B>& v1, const VectorTuple<T,N,B>& v2)
{
	return !(v1 == v2);
}

////////////////////////////////
// stream output
////////////////////////////////

// initialize delimiter variables
template <class T, size_t N, size_t B>
	std::string VectorTuple<T,N,B>::delimiter_left_("(");
template <class T, size_t N, size_t B>
	std::string VectorTuple<T,N,B>::delimiter_middle_(",");
template <class T, size_t N, size_t B>
	std::string VectorTuple<T,N,B>::delimiter_right_(")");

// delimiter configuration
template <class T, size_t N, size_t B>
void VectorTuple<T,N,B>::SetDelimiters(
	const std::string& left, 
	const std::string& middle, 
	const std::string& right
	)
{
	VectorTuple<T,N,B>::delimiter_left_ = left;
	VectorTuple<T,N,B>::delimiter_middle_ = middle;
	VectorTuple<T,N,B>::delimiter_right_ = right;
};


// output operator
template <class T, size_t N, size_t B>
std::ostream& operator<< (std::ostream& os, const VectorTuple<T,N,B>& v)
{
	os << VectorTuple<T,N,B>::delimiter_left_;
	for(typename VectorTuple<T,N,B>::const_iterator i = v.begin(); i != v.end(); ++i)
	{
		if (i != v.begin())
			os << VectorTuple<T,N,B>::delimiter_middle_;
		os << (*i);
	}
	os << VectorTuple<T,N,B>::delimiter_right_;

	return os;
}


#endif
