/****************************************************************
  moshinsky_bracket.h

  Implements recursive calculation of harmonic oscillator Moshinsky brackets.

  References:

  [Moshinsky1959] M. Moshinsky, Transformation brackets for harmonic oscillator
  functions, Nucl. Phys. 13, 104 (1959). DOI 10.1016/0029-5582(59)90143-9.

  [HO] M. Moshinsky and Y. F. Smirnov, The harmonic oscillator in modern physics
  (Harwood Academic Publishers, Amsterdam, 1996).

  [TTB] T. A. Brody and M. Moshinsky, Tables of transformation brackets for
  nuclear shell-model calculations, Monografias del Instituto de Fisica,
  Universidad Nacional Autonoma de Mexico, Mexico, 1960.

  Language: C++11

  Mark A. Caprio, University of Notre Dame.

  + 02/15/11 (mac): Created.
  + 11/14/15 (mac): Update header file.  Remove special Ncm=0 functions.
  + 11/26/15 (mac): Rename from moshinsky to moshinsky_bracket.
  + 03/05/16 (mac): Eradicate use of shell_indexing_nl.
  + 07/04/16 (mac):
    - Update to current usage of am module (namespace conventions, use
      of HalfInt::pair, use of int conversion).
    - Update #include guard and include files.
  + 10/09/16 (pjf): Rename mcpp -> mcutils.
  + 05/16/26 (mac): Change memoization key type from VectorTuple to std::tuple. 

****************************************************************/

#ifndef MOSHINSKY_BRACKET_H_
#define MOSHINSKY_BRACKET_H_

namespace moshinsky {

  ////////////////////////////////////////////////////////////////
  // bracket calculation
  ////////////////////////////////////////////////////////////////

  // debugging flag variable
  
  extern bool trace_moshinsky;

  // Moshinsky bracket
  
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
  //   n1_dot (int): Relative radial quantum number.
  //   l1_dot (int): Relative orbital angular momentum.
  //   n2_dot (int): Center-of-mass radial quantum number.
  //   l2_dot (int): Center-of-mass orbital angular momentum.
  //   n1 (int): Particle 1 radial quantum number.
  //   l1 (int): Particle 1 orbital angular momentum.
  //   n2 (int): Particle 2 radial quantum number.
  //   l2 (int): Particle 2 orbital angular momentum.
  //   Lambda (int): Coupled orbital angular momentum.
  //
  //  Returns:
  //    Moshinsky bracket.

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
