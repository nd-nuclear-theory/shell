/****************************************************************
  moshinsky_bracket.h

  Implements recursive calculation of harmonic oscillator Moshinsky brackets.

  References (Moshinsky bracket):

    [Moshinsky (1959)] M. Moshinsky, Transformation brackets for harmonic
    oscillator functions, Nucl. Phys. 13, 104 (1959). DOI
    10.1016/0029-5582(59)90143-9.

    [HOMP] M. Moshinsky and Y. F. Smirnov, The harmonic oscillator in modern
    physics (Harwood Academic Publishers, Amsterdam, 1996).

    [TTB] T. A. Brody and M. Moshinsky, Tables of transformation brackets for
    nuclear shell-model calculations, Monografias del Instituto de Fisica,
    Universidad Nacional Autonoma de Mexico, Mexico, 1960.

  References (generalized Moshinsky bracket):

    [Trlifaj (1972)] PRC 5, 5 (1972).  [In fact, PRC 5(5), 5...]

    [Bevelacqua (1978)] CJP 57, 1136 (1979).

    [Kamuntavicius (2001)] NPA 695, 191 (2001).

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
  // Return Moshinsky bracket between relative/cm (dotted) and single-particle
  // (undotted) oscillator product states.
  //
  // Args:
  //   n1_dot (input): Relative radial quantum number.
  //   l1_dot (input): Relative orbital angular momentum.
  //   n2_dot (input): Center-of-mass radial quantum number.
  //   l2_dot (input): Center-of-mass orbital angular momentum.
  //   n1 (input): Particle 1 radial quantum number.
  //   l1 (input): Particle 1 orbital angular momentum.
  //   n2 (input): Particle 2 radial quantum number.
  //   l2 (input): Particle 2 orbital angular momentum.
  //   Lambda (input): Coupled orbital angular momentum.
  //
  //  Returns:
  //    Moshinsky bracket.

  double GeneralizedMoshinskyBracket(
      int n1_dot, int l1_dot,
      int n2_dot, int l2_dot,
      int n1, int l1, int n2, int l2,
      int Lambda,
      double d
    );
  // Return generalized Moshinsky bracket between relative/cm (dotted) and
  // single-particle (undotted) oscillator product states, that is, generalized
  // to arbitrary mass ratio.
  //
  // Limitation: Initial implementation is restricted to case in which either
  // l1_dot=0 or l2_dot=0.
  //
  // Args:
  //   n1_dot (input): Relative radial quantum number.
  //   l1_dot (input): Relative orbital angular momentum.
  //   n2_dot (input): Center-of-mass radial quantum number.
  //   l2_dot (input): Center-of-mass orbital angular momentum.
  //   n1 (input): Particle 1 radial quantum number.
  //   l1 (input): Particle 1 orbital angular momentum.
  //   n2 (input): Particle 2 radial quantum number.
  //   l2 (input): Particle 2 orbital angular momentum.
  //   Lambda (input): Coupled orbital angular momentum.
  //   d (input): Mass ratio parameter.
  //
  //  Returns:
  //    Generalized Moshinsky bracket.
  
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
