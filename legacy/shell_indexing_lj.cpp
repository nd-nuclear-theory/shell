/****************************************************************
  shell_indexing_lj.cpp

  Mark A. Caprio, University of Notre Dame.
  Last modified 4/25/15.

****************************************************************/

#include <sstream>

#include <shell/shell_indexing_lj.h>

namespace shell {

  ////////////////////////////////////////////////////////////////
  // SPSpacelj implementation
  ////////////////////////////////////////////////////////////////

  std::string SPSpacelj::String() const
  {
    std::ostringstream ss;

    ss << GetIndex() << "(" << Getl() << "," << Getj() << ")";

    return ss.str();
  }

  std::ostream& operator<< (std::ostream& os, const SPSpacelj& orbital)
  {
    os << orbital.String();
	
    return os;
  }


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
