/****************************************************************
  moshinsky_bracket.h 

  Defines Moshinsky bracket calculations.

  Reference: M. Moshinsky and Y. F. Smirnov, The harmonic oscillator
  in modern physics (Harwood Academic Publishers, Amsterdam, 1996).

  Language: C++11
                                 
  Mark A. Caprio, University of Notre Dame.

  2/15/11 (mac): Created.
  11/14/15 (mac): Updated header file.  Removed special Ncm=0 functions.
  11/26/15 (mac): Renamed from moshinsky to moshinsky_bracket.
  3/5/16 (mac): Eradicate use of shell_indexing_nl.
  7/4/16 (mac):
    - Update to current usage of am module (namespace conventions, use 
      of HalfInt::pair, use of int conversion).
    - Update #include guard and include files.
  10/9/16 (pjf): Rename mcpp -> mcutils.

FUTURE: Reimplement using C++11 STL tuples in place of home-grown
VectorTuple.  Perhaps eradicate wrappers for shell_indexing_nl state
arguments.  Expurgate "using" statements.

****************************************************************/

#ifndef MOSHINSKY_BRACKET_H_
#define MOSHINSKY_BRACKET_H_

#include <vector>

namespace moshinsky {


  ////////////////////////////////////////////////////////////////
  // bracket calculation
  ////////////////////////////////////////////////////////////////

  // debugging flag variable
  extern bool trace_moshinsky;

  // Moshinsky bracket -- general form
  double MoshinskyBracket(
      int n1_dot, int l1_dot, 
      int n2_dot, int l2_dot, 
      int n1, int l1, int n2, int l2,
      int Lambda
    );

  // Returns Moshinsky bracket between relative/cm (dotted) and
  // single-particle (undotted) oscillator product states.
  //
  // Args:
  //   n1_dot (int) : relative *radial* qn
  //   ...
  //   n2_dot (int) : center-of-mass *radial* qn
  //   ...
  //   n1 (int) : particle 1 *radial* qn
  //   ...
  //   n1 (int) : particle 1 *radial* qn
  //   ...
  //   Lambda (int) : coupled orbital angular momentum


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
