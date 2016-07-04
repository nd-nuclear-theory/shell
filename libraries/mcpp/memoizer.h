/****************************************************************

  memoizer.h

  Template library for generic "memoization"/caching structures.

  Created by Mark A. Caprio, 2/23/11.

****************************************************************/

#ifndef MEMOIZER_H_
#define MEMOIZER_H_

#include <cstddef>
#include <algorithm>
#include <iostream>
#include <map>
#include <vector>

// MEMOIZE(m,x,y) applies Memoizer m to key given by expression x
//   - if x already exists as key, a const reference to the stored "y" 
//     value for key x is returned
//   - else expression y is evaluated, stored as the "y" value for key x,
//     and the "y" value is returned
// Optimization: Note key expression x is evaluated twice by the macro.
// On the other hand, it cannot be saved within m, because of possible
// reentrance problems if m is invoked in evaluating y.

#define MEMOIZE(m,x,y) ( (m.Seek(x)) ? m.GetValue() : m.SetValue(x,y) )

template<typename Key, typename T, 
	typename Compare = std::less<Key> ,
	typename Alloc = std::allocator<std::pair<const Key, T> > >
class Memoizer{

public:

	////////////////////////////////
	// type definitions
	////////////////////////////////

	// underlying map type
	typedef typename std::map<Key,T,Compare,Alloc> map_type;

        // standard map type definitions
	typedef typename map_type::key_type key_type;
	typedef typename map_type::mapped_type mapped_type;
	typedef typename map_type::value_type value_type;
	typedef typename map_type::key_compare key_compare;
	typedef typename map_type::allocator_type allocator_type;
	typedef typename map_type::reference reference;
	typedef typename map_type::const_reference const_reference;
	typedef typename map_type::iterator iterator;
	typedef typename map_type::const_iterator const_iterator;
	typedef typename map_type::size_type size_type;
	typedef typename map_type::difference_type difference_type;
	typedef typename map_type::pointer pointer;
	typedef typename map_type::const_pointer const_pointer;
	typedef typename map_type::reverse_iterator reverse_iterator;
	typedef typename map_type::const_reverse_iterator const_reverse_iterator;
	typedef typename map_type::value_compare value_compare;

	// vector type for key-value pair dump
	typedef typename std::vector<value_type> vector_type;

	////////////////////////////////
	// constructors
	////////////////////////////////

	// default
	//   default-initialize empty map
	//   enable or disable caching
	//   provide current_reference_ with useless but syntactically-required
	//     initial value (to a temporary object)
	Memoizer() : cache_enabled_(true) {};
	Memoizer(bool b) {cache_enabled_ = b;};

        // copy -- synthesized constructor copies members

	////////////////////////////////
        // accessors
	////////////////////////////////

	// Seek(x) seeks entry for key x, returns true if found
	//   as side effect, saves reference to "y" value for subsequent rapid
 	//   access by GetValue()
	bool Seek(const Key& x);

	// GetValue() returns the "y" value for 
	//   the most recently sought key
	//   note: this requires that no recursive call be made 
	//   before GetValue() invoked
	T GetValue() const;

	// SetValue(y) stores y value 
        //   and returns a copy of this value
	T SetValue(const Key& x, const T& y);

	// Known(x) determines whether or not a value for x is already stored
	//   for debugging and diagnostic use -- not part of MEMOIZE call
	bool Known(const Key& x) const {return (values_.count(x) == 1);}

	////////////////////////////////
        // iterators
	////////////////////////////////

	// iterators for read access only are defined

	const_iterator begin() const {return values_.begin();};
	const_iterator end() const {return values_.end();};
	const_reverse_iterator rbegin() const {return values_.rbegin();};
	const_reverse_iterator rend() const {return values_.rend();};

	////////////////////////////////
	// bulk access
	////////////////////////////////

        size_type size() const {return values_.size();};
        void clear() {return values_.clear();};

	////////////////////////////////
	// ostream output
	////////////////////////////////

        // output operator -- friend declaration for access to delimiters
	template<typename KeyX, typename TX, typename CompareX, typename AllocX>
	friend std::ostream& operator<< (std::ostream&, const Memoizer<KeyX, TX, CompareX, AllocX>&);

	////////////////////////////////
	// configuration
	////////////////////////////////

	// mode flags
	void EnableCaching(bool b) {cache_enabled_ = b;};

        // configuring delimiter strings
        //   static member function sets delimiters for *all* Memoizer
        //   instances with the given template parameters
        // EX: Memoizer<...>::SetDelimiters(" ( ", " -> ", " )\n"); 
	static void SetDelimiters(const std::string&, const std::string&, const std::string&);

private:

	////////////////////////////////
	// configuration data
	////////////////////////////////

        // ostream delimiters (static)
        static std::string delimiter_left_;
        static std::string delimiter_middle_;
        static std::string delimiter_right_;

	// mode variables
	bool cache_enabled_;

	////////////////////////////////
	// caching data
	////////////////////////////////

        // cache map container
        map_type values_;

	// current entry access
	// Key current_key_;
	T current_result_;

};

template<typename Key, typename T, typename Compare, typename Alloc>
inline
bool Memoizer<Key,T,Compare,Alloc>::Seek(const Key& x)
{
	
	if (cache_enabled_)
		// cache enabled
	{
		// look for key in map
		iterator it = values_.find(x);
		
		// process key result
		if ( it != values_.end()) 
			// key found
		{
			// store associated "y" value for rapid access
			current_result_ = it->second;
			
			return true;
		}
		else
			// key not found
			return false;
	}
	else
		// cache not enabled
		return false;
}

template<typename Key, typename T, typename Compare, typename Alloc>
inline
T Memoizer<Key,T,Compare,Alloc>::GetValue() const
{
	// return "y" value
	return current_result_;
}

template<typename Key, typename T, typename Compare, typename Alloc>
inline
T Memoizer<Key,T,Compare,Alloc>::SetValue(const Key& x, const T& y)
{
	if (cache_enabled_)
	{
		values_.insert(value_type(x,y));
	}

	return y;
};


////////////////////////////////
// stream output
////////////////////////////////

// initialize delimiter variables
template<typename Key, typename T, typename Compare, typename Alloc>
	std::string Memoizer<Key,T,Compare,Alloc>::delimiter_left_("  ");
template<typename Key, typename T, typename Compare, typename Alloc>
	std::string Memoizer<Key,T,Compare,Alloc>::delimiter_middle_("->");
template<typename Key, typename T, typename Compare, typename Alloc>
	std::string Memoizer<Key,T,Compare,Alloc>::delimiter_right_("\n");

// delimiter configuration
template<typename Key, typename T, typename Compare, typename Alloc>
void Memoizer<Key,T,Compare,Alloc>::SetDelimiters(
	const std::string& left, 
	const std::string& middle, 
	const std::string& right
	)
{
	Memoizer<Key,T,Compare,Alloc>::delimiter_left_ = left;
	Memoizer<Key,T,Compare,Alloc>::delimiter_middle_ = middle;
	Memoizer<Key,T,Compare,Alloc>::delimiter_right_ = right;
};


// output operator
template<typename Key, typename T, typename Compare, typename Alloc>
std::ostream& operator<< (std::ostream& os, const Memoizer<Key,T,Compare,Alloc>& m)
{
	for(
		typename Memoizer<Key,T,Compare,Alloc>::const_iterator i = m.begin(); 
		i != m.end(); 
		++i
		)
	{
		os << Memoizer<Key,T,Compare,Alloc>::delimiter_left_;
		os << i->first;
		os << Memoizer<Key,T,Compare,Alloc>::delimiter_middle_;
		os << i->second;
		os << Memoizer<Key,T,Compare,Alloc>::delimiter_right_;
	}

	return os;
}


#endif
