/****************************************************************
  shell_indexing_nlj.cpp

  Mark A. Caprio, University of Notre Dame.
  Last modified 4/25/15.

****************************************************************/

#include <sstream>

#include <shell/shell_indexing_nlj.h>

namespace shell {

  ////////////////////////////////////////////////////////////////
  // SPOrbitalNlj cache static variable definitions
  ////////////////////////////////////////////////////////////////

  std::vector<SPOrbitalNlj::SPLabelsNlj> SPOrbitalNlj::label_cache_;
  int SPOrbitalNlj::max_cached_N_ = -1;

  ////////////////////////////////////////////////////////////////
  // SPOrbitalNlj implementation
  ////////////////////////////////////////////////////////////////

  std::string SPOrbitalNlj::String() const
  {
    std::ostringstream ss;

    ss << GetIndex() << "(" << Getn() << "," << Getl() << "," << Getj() << ")";

    return ss.str();
  }

  std::ostream& operator<< (std::ostream& os, const SPOrbitalNlj& orbital)
  {
    os << orbital.String();
	
    return os;
  }


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
