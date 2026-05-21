/****************************************************************
  moshinsky_bracket.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/


#include "moshinsky_bracket.h"

#include <algorithm>
#include <cassert>
#include <iostream>  // for debugging
#include <tuple>

#include "gsl/gsl_sf_gamma.h"  // for factorial
#include "gsl/gsl_math.h" // for integer powers

#include "am/am.h"
#include "am/wigner_gsl.h"
#include "mcutils/arithmetic.h"  // for ONLYIF
#include "mcutils/memoizer.h"

namespace moshinsky {

  ////////////////////////////////////////////////////////////////
  // debugging flag
  ////////////////////////////////////////////////////////////////

  bool trace_moshinsky = false;

  ////////////////////////////////////////////////////////////////
  // Moshinsky bracket seed generation
  ////////////////////////////////////////////////////////////////

  double MoshinskyACoefficient(
      int l1, int l1_dot, int l2, int l2_dot, int kappa
    )
  // Calculate coefficient A(l1, l1_dot, l2, l2_dot, kappa) entering into
  // expression for seed Moshinsky coefficient with n1=n2=0.
  //
  // Defined in (64) of Moshinsky (1959) or (2.12) of TTB.
  {
    // compute prefactor
    double prefactor = sqrt(
        ( gsl_sf_fact(l1+l1_dot+kappa+1) * gsl_sf_fact(l1+l1_dot-kappa)
          * gsl_sf_fact(l1-l1_dot+kappa) )
        / ( gsl_sf_fact(l1_dot-l1+kappa) )
        * ( gsl_sf_fact(l2+l2_dot+kappa+1) * gsl_sf_fact(l2+l2_dot-kappa)
            * gsl_sf_fact(l2-l2_dot+kappa) )
        / ( gsl_sf_fact(l2_dot-l2+kappa) )
      );

    // determine summation range for q
    int q_min = kappa + ( (kappa + l1_dot + l1) % 2);
    int q_max = std::min(l1_dot + l1, l2_dot + l2);

    // accumulate sum
    double sum = 0.;
    for (int q = q_min; q <= q_max; q += 2)
      {
        // std::cerr << "  q: " << q << std::endl;
        sum += ParitySign((l1_dot-l1+q)/2)
          / ( gsl_sf_fact(q-kappa) * gsl_sf_fact(q+kappa+1) )
          * ( gsl_sf_fact(l1_dot-l1+q) )
          / ( gsl_sf_fact((l1_dot-l1+q)/2) * gsl_sf_fact((l1_dot+l1-q)/2) )
          * ( gsl_sf_fact(l2_dot-l2+q) )
          / ( gsl_sf_fact((l2_dot-l2+q)/2) * gsl_sf_fact((l2_dot+l2-q)/2) )
          ;
      }

    return prefactor * sum;

  }

  double SeedMoshinskyBracket(
      int n1_dot, int l1_dot, int n2_dot, int l2_dot,
      int l1, int l2,
      int Lambda
    )
  // Calculate seed Moshinsky coefficient with n1=n2=0.
  //
  // Expression is provided in (63) of Moshinsky (1959) or (2.11) of TTB.  This
  // is a simplification of (10.30) of HO.
  {
    // compute prefactor
    double prefactor = ParitySign(n1_dot+l1_dot+l2_dot-Lambda)
      * sqrt(
          ( gsl_sf_fact(l1) * gsl_sf_fact(l2) )
          / ( gsl_sf_fact(2*l1) * gsl_sf_fact(2*l2) )
          * ( (2*l1_dot+1) * (2*l2_dot+1) )
          / ( gsl_pow_int(2.,l1_dot+l2_dot) )
          * ( gsl_sf_fact(n1_dot+l1_dot) )
          / ( gsl_sf_fact(n1_dot) * gsl_sf_fact(2*n1_dot+2*l1_dot+1) )
          * ( gsl_sf_fact(n2_dot+l2_dot) )
          / ( gsl_sf_fact(n2_dot) * gsl_sf_fact(2*n2_dot+2*l2_dot+1) )
        );

    // determine summation range for kappa
    HalfInt::pair kappa_bound = am::AngularMomentumRangeIntersection(
        am::ProductAngularMomentumRange(l1,l1_dot),
        am::ProductAngularMomentumRange(l2,l2_dot)
      );
    int kappa_min = int(kappa_bound.first);
    int kappa_max = int(kappa_bound.second);

    // accumulate sum
    double sum = 0.;
    for (int kappa = kappa_min; kappa <= kappa_max; ++kappa)
      {
        // std::cerr << "kappa: " << kappa << std::endl;
        sum += (2*kappa+1)
          * MoshinskyACoefficient(l1,l1_dot,l2,l2_dot,kappa)
          / ( Hat(Lambda) * Hat(kappa) )
          * am::Unitary6J(l1_dot,l2_dot,Lambda,l2,l1,kappa);
      }

    return prefactor * sum;
  }

  
  ////////////////////////////////////////////////////////////////
  // MoshinskyBracket
  ////////////////////////////////////////////////////////////////

  double MoshinskyBracket(
      int n1_dot, int l1_dot,
      int n2_dot, int l2_dot,
      int n1, int l1, int n2, int l2,
      int Lambda
    )
  // Implements recurrence relation (10.36) of HOMP, for recurrence on n1.  When
  // n1 vanishes, the symmetry relation (10.38) of HOMP is used to swap n1 and
  // n2.  Thus, both n1 and n2 are reduced to zero, reducing the problem to
  // evaluation of seed Moshinsky coefficients with n1=n2=0.
  //
  // See also mac notes 11/09/16.
  //
  // A static (persistent) cache is used to memoize computed values.
  {

    // tracing output
    if (trace_moshinsky)
      {
        std::cerr << "<" <<  n1_dot << " " << l1_dot << " "  << n2_dot << " " << l2_dot << " ; " << Lambda
                  << " | "
                  << n1  << " " << l1 << " " << n2 << " " << l2 << " ; " << Lambda << ">"
                  << std::endl;
      }

    // validate bracket
    //   return 0. for negative argument cases
    if ( (n1_dot < 0) || (l1_dot < 0) || (n2_dot < 0) || (l2_dot < 0)
         || (n1 < 0) || (l1 < 0) || (n2 < 0) || (l2 < 0) )
      {
        if (trace_moshinsky)
          std::cerr << "Moshinsky negative arg???" << std::endl;
        return 0.;
      }

    // compute phonon parameters
    int rho_dot = 2*n1_dot + l1_dot + 2*n2_dot + l2_dot;
    int rho = 2*n1 + l1 + 2*n2 + l2;

    // set up caching key
    typedef std::tuple<int, int, int, int, int, int, int, int, int> KeyType;
    static mcutils::Memoizer<KeyType,double> m;
    KeyType memo_key(n1_dot, l1_dot, n2_dot, l2_dot, n1, l1, n2, l2, Lambda);
    
    // evaluate bracket
    double value;
    if (m.Seek(memo_key))
      // key found
      {
        // retrieve stored value
        value = m.GetValue();
      }
    else
      {
        // calculate new value
        if ( (rho_dot != rho)
             || !am::AllowedTriangle(l1_dot,l2_dot,Lambda)
             || !am::AllowedTriangle(l1,l2,Lambda)
          )
          // case: phonon or angular-momentum forbidden
          {
            if (trace_moshinsky) std::cerr << "   " <<  "case: forbidden" << std::endl;
            value = 0.;
          }
        else if ( (n1 == 0) && (n2 == 0) )
          // case: at seed value on RHS
          {
            if (trace_moshinsky) std::cerr << "   " <<  "case: seed" << std::endl;
            value = SeedMoshinskyBracket(n1_dot,l1_dot,n2_dot,l2_dot,l1,l2,Lambda);
          }
        else if (n1 == 0)
          // case: cannot recurse n1 further but can recurse n2
          {
            if (trace_moshinsky) std::cerr << "   " <<  "case: reverse ket" << std::endl;
            value = ParitySign(l2_dot-Lambda)
              * MoshinskyBracket(n1_dot,l1_dot,n2_dot,l2_dot,n2,l2,n1,l1,Lambda);
          }
        else
          // otherwise: recurse n1
          {
            if (trace_moshinsky) std::cerr << "   " <<  "case: recurse" << std::endl;

            double prefactor = 1 / sqrt( n1 * (n1 + l1 + 1/2.) ) ;
            double terms =
              ONLYIF( n1_dot != 0 ,
                      1/2.*sqrt(n1_dot*(n1_dot+l1_dot+1/2.))
                      * MoshinskyBracket(
                          n1_dot-1, l1_dot, n2_dot, l2_dot,
                          n1-1, l1, n2, l2,
                          Lambda
                        )
                )
              + ONLYIF( n2_dot != 0 ,
                        1/2.*sqrt(n2_dot*(n2_dot+l2_dot+1/2.))
                        * MoshinskyBracket(
                            n1_dot, l1_dot, n2_dot-1, l2_dot,
                            n1-1, l1, n2, l2,
                            Lambda
                          )
                )
              + ONLYIF( (n1_dot != 0) && (n2_dot != 0) ,
                        ParitySign(l1_dot+l2_dot+Lambda)*sqrt((n1_dot)*(n2_dot)*(l1_dot+1)*(l2_dot+1))
                        * am::Wigner6J(l1_dot, l1_dot+1, 1, l2_dot+1, l2_dot, Lambda)
                        * MoshinskyBracket(
                            n1_dot-1, l1_dot+1, n2_dot-1, l2_dot+1,
                            n1-1, l1, n2, l2,
                            Lambda
                          )
                )
              + ONLYIF( (n1_dot != 0) && (l2_dot != 0) ,
                        ParitySign(l1_dot+l2_dot+Lambda)*sqrt((n1_dot)*(n2_dot+l2_dot+1/2.)*(l1_dot+1)*(l2_dot))
                        * am::Wigner6J(l1_dot, l1_dot+1, 1, l2_dot-1, l2_dot, Lambda)
                        * MoshinskyBracket(
                            n1_dot-1, l1_dot+1, n2_dot, l2_dot-1,
                            n1-1, l1, n2, l2,
                            Lambda
                          )
                )
              + ONLYIF( (l1_dot != 0) && (n2_dot != 0) ,
                        ParitySign(l1_dot+l2_dot+Lambda)*sqrt((n1_dot+l1_dot+1/2.)*(n2_dot)*(l1_dot)*(l2_dot+1))
                        * am::Wigner6J(l1_dot, l1_dot-1, 1, l2_dot+1, l2_dot, Lambda)
                        * MoshinskyBracket(
                            n1_dot, l1_dot-1, n2_dot-1, l2_dot+1,
                            n1-1, l1, n2, l2,
                            Lambda
                          )
                )
              + ONLYIF( (l1_dot != 0) && (l2_dot != 0) ,
                        ParitySign(l1_dot+l2_dot+Lambda)*sqrt((n1_dot+l1_dot+1/2.)*(n2_dot+l2_dot+1/2.)*(l1_dot)*(l2_dot))
                        * am::Wigner6J(l1_dot, l1_dot-1, 1, l2_dot-1, l2_dot, Lambda)
                        * MoshinskyBracket(
                            n1_dot, l1_dot-1, n2_dot, l2_dot-1,
                            n1-1, l1, n2, l2,
                            Lambda
                          )
                );
            value = prefactor * terms;
          }

        // store new value
        m.SetValue(memo_key,value);
      }

    if (trace_moshinsky) std::cerr << "   " <<  value << std::endl;

    return value;

  }


  ////////////////////////////////////////////////////////////////
  // generalized Moshinsky bracket seed generation
  ////////////////////////////////////////////////////////////////

  double SeedGeneralizedMoshinskyBracket(
      int n1_dot, int l1_dot, int n2_dot, int l2_dot,
      int l1, int l2,
      int Lambda,
      double d
    )
  // Calculate seed generalized Moshinsky coefficient with n1=n2=0.
  //
  // Limitation: Initial implementation is only for cases in which either l1_dot=0, as in (11) of Trlifaj
  // (1972), or l2_dot=0, which can then be obtained by 1<->2 interchange, with the phase factor of (10.38) of HO,
  //  as verified for generalized Moshinsky brackets in
  // (15) of Kamuntavicius (2001).  Note that center-of-mass free (CMF)
  // Moshinsky brackets, for which n2_dot=l2_dot=0, are obtained as an important special case.
  //
  // For future general implementation, one would specialize either (23) of
  // Trlifaj (1972) or (27) of Kamuntavicius (2001) to n1=n2=0, and then apply
  // simplifications as in going from (60) to (63) of Moshinsky (1959).
  {
    double value = 0.0;

    // short circuit: perform 1<->2 swap if needed
    if ((l2_dot != 0) && (l1_dot ==0))
      {
        if (trace_moshinsky) std::cerr << "   " <<  "case: perform 1<->2 swap for l1_dot==0" << std::endl;
        value = ParitySign(l1-Lambda)
          * SeedGeneralizedMoshinskyBracket(n2_dot, l2_dot, n1_dot, l1_dot, l1, l2, Lambda, d);
        return value;
      }
    
    assert( l2_dot == 0 );

    // PLACEHOLDER: use seed for ordinary Moshinsky bracket
    if (d == 1.) {
      value = SeedMoshinskyBracket(
          n1_dot, l1_dot, n2_dot, l2_dot,
          l1, l2,
          Lambda
        );
      return value;
    }

    // restrict to case (l2_dot=0) of (11) of Trlifaj (1972)
    // with n1=n2=0, which corespond to the seed
    if ((2*n2_dot+2*n1_dot+l1_dot) == (l1+l2)) {
        value  = std::pow(-1, n1_dot) * 0.5 * std::sqrt(M_PI) * std::pow(1+d, -0.5*l1) * std::pow(d/(1+d), 0.5*l2)
                * am::ClebschGordan(l1, 0, l2, 0, l1_dot, 0)
                * std::sqrt((2*l1+1) * (2*l2+1) * gsl_sf_fact(n2_dot))
                / std::sqrt((2*l1_dot+1) * gsl_sf_fact(n1_dot) * gsl_sf_gamma(n2_dot+1.5) * gsl_sf_gamma(n1_dot+l1_dot+1.5))
                * gsl_sf_fact(0.5*(l1+l2-l1_dot)) * gsl_sf_gamma(0.5*(l1+l2+l1_dot) + 1.5)
                / gsl_sf_fact((l1+l2-l1_dot)/2 - n1_dot) 
                / std::sqrt(gsl_sf_gamma(l1+1.5) * gsl_sf_gamma(l2+1.5));
    } else {
        value = 0.;
    }
      
    return value;
  };

  
  ////////////////////////////////////////////////////////////////
  // GeneralizedMoshinskyBracket
  ////////////////////////////////////////////////////////////////

  double GeneralizedMoshinskyBracket(
      int n1_dot, int l1_dot,
      int n2_dot, int l2_dot,
      int n1, int l1, int n2, int l2,
      int Lambda,
      double d
    )
  // Implements recurrence relation (10.36) of HOMP, as modified for generalized Moshinsky brackets in
  // Bevelacqua (1979), for recurrence on n1.  When n1 vanishes, the symmetry
  // relation (10.38) of HOMP, as verified for generalized Moshinsky brackets in
  // (14) of Kamuntavicius (2001), is used to swap n1 and n2.  Thus, both n1 and n2
  // are reduced to zero, reducing the problem to evaluation of seed Moshinsky
  // coefficients with n1=n2=0.
  //
  // A persistent cache is used to memoize computed values.  This cache cannot
  // be shared between calculations with distinct d.  Beware that caching
  // retrieval relies on floating point comparison of d, which should be
  // successful if d is not independently recalculated between different calls
  // to GeneralizedMoshinskyBracket.
  {

    // tracing output
    if (trace_moshinsky)
      {
        std::cerr << "<" <<  n1_dot << " " << l1_dot << " "  << n2_dot << " " << l2_dot << " ; " << Lambda
                  << " | "
                  << n1  << " " << l1 << " " << n2 << " " << l2 << " ; " << Lambda << ">"
                  << std::endl;
      }

    // validate bracket
    //   return 0. for negative argument cases
    if ( (n1_dot < 0) || (l1_dot < 0) || (n2_dot < 0) || (l2_dot < 0)
         || (n1 < 0) || (l1 < 0) || (n2 < 0) || (l2 < 0) )
      {
        if (trace_moshinsky)
          std::cerr << "Moshinsky negative arg???" << std::endl;
        return 0.;
      }

    // compute phonon parameters
    int rho_dot = 2*n1_dot + l1_dot + 2*n2_dot + l2_dot;
    int rho = 2*n1 + l1 + 2*n2 + l2;

    // set up caching key
    typedef std::tuple<int, int, int, int, int, int, int, int, int, double> KeyType;
    static mcutils::Memoizer<KeyType,double> m;
    KeyType memo_key(n1_dot, l1_dot, n2_dot, l2_dot, n1, l1, n2, l2, Lambda, d);
    
    // evaluate bracket
    double value;
    if (m.Seek(memo_key))
      // key found
      {
        // retrieve stored value
        value = m.GetValue();
      }
    else
      {
        // calculate new value
        if ( (rho_dot != rho)
             || !am::AllowedTriangle(l1_dot,l2_dot,Lambda)
             || !am::AllowedTriangle(l1,l2,Lambda)
          )
          // case: phonon or angular-momentum forbidden
          {
            if (trace_moshinsky) std::cerr << "   " <<  "case: forbidden" << std::endl;
            value = 0.;
          }
        else if ( (n1 == 0) && (n2 == 0) )
          // case: at seed value on RHS
          {
            if (trace_moshinsky) std::cerr << "   " <<  "case: seed" << std::endl;
            value = SeedGeneralizedMoshinskyBracket(n1_dot,l1_dot,n2_dot,l2_dot,l1,l2,Lambda,d);
          }
        else if (n1 == 0)
          // case: cannot recurse n1 further but can recurse n2
          {
            if (trace_moshinsky) std::cerr << "   " <<  "case: reverse ket" << std::endl;
            value = ParitySign(l2_dot-Lambda)
              * GeneralizedMoshinskyBracket(n1_dot,l1_dot,n2_dot,l2_dot,n2,l2,n1,l1,Lambda,d);
          }
        else
          // otherwise: recurse n1
          {
            if (trace_moshinsky) std::cerr << "   " <<  "case: recurse" << std::endl;

            // TODO (mac): Add prefactors from Bevelacqua
            double prefactor = (
                2 / (1+d)  // modification factor for generalized Moshinsky bracket
                * 1 / sqrt( n1 * (n1 + l1 + 1/2.) )
              );
            double terms =
              ONLYIF( n1_dot != 0 ,
                      d *  // modification factor for generalized Moshinsky bracket
                      1/2.*sqrt(n1_dot*(n1_dot+l1_dot+1/2.))
                      * GeneralizedMoshinskyBracket(
                          n1_dot-1, l1_dot, n2_dot, l2_dot,
                          n1-1, l1, n2, l2,
                          Lambda, d
                        )
                )
              + ONLYIF( n2_dot != 0 ,
                        1. *  // modification factor for generalized Moshinsky bracket
                        1/2.*sqrt(n2_dot*(n2_dot+l2_dot+1/2.))
                        * GeneralizedMoshinskyBracket(
                            n1_dot, l1_dot, n2_dot-1, l2_dot,
                            n1-1, l1, n2, l2,
                            Lambda, d
                          )
                )
              + ONLYIF( (n1_dot != 0) && (n2_dot != 0) ,
                        sqrt(d) *  // modification factor for generalized Moshinsky bracket
                        ParitySign(l1_dot+l2_dot+Lambda)*sqrt((n1_dot)*(n2_dot)*(l1_dot+1)*(l2_dot+1))
                        * am::Wigner6J(l1_dot, l1_dot+1, 1, l2_dot+1, l2_dot, Lambda)
                        * GeneralizedMoshinskyBracket(
                            n1_dot-1, l1_dot+1, n2_dot-1, l2_dot+1,
                            n1-1, l1, n2, l2,
                            Lambda, d
                          )
                )
              + ONLYIF( (n1_dot != 0) && (l2_dot != 0) ,
                        sqrt(d) *  // modification factor for generalized Moshinsky bracket
                        ParitySign(l1_dot+l2_dot+Lambda)*sqrt((n1_dot)*(n2_dot+l2_dot+1/2.)*(l1_dot+1)*(l2_dot))
                        * am::Wigner6J(l1_dot, l1_dot+1, 1, l2_dot-1, l2_dot, Lambda)
                        * GeneralizedMoshinskyBracket(
                            n1_dot-1, l1_dot+1, n2_dot, l2_dot-1,
                            n1-1, l1, n2, l2,
                            Lambda, d
                          )
                )
              + ONLYIF( (l1_dot != 0) && (n2_dot != 0) ,
                        sqrt(d) *  // modification factor for generalized Moshinsky bracket
                        ParitySign(l1_dot+l2_dot+Lambda)*sqrt((n1_dot+l1_dot+1/2.)*(n2_dot)*(l1_dot)*(l2_dot+1))
                        * am::Wigner6J(l1_dot, l1_dot-1, 1, l2_dot+1, l2_dot, Lambda)
                        * GeneralizedMoshinskyBracket(
                            n1_dot, l1_dot-1, n2_dot-1, l2_dot+1,
                            n1-1, l1, n2, l2,
                            Lambda, d
                          )
                )
              + ONLYIF( (l1_dot != 0) && (l2_dot != 0) ,
                        sqrt(d) *  // modification factor for generalized Moshinsky bracket
                        ParitySign(l1_dot+l2_dot+Lambda)*sqrt((n1_dot+l1_dot+1/2.)*(n2_dot+l2_dot+1/2.)*(l1_dot)*(l2_dot))
                        * am::Wigner6J(l1_dot, l1_dot-1, 1, l2_dot-1, l2_dot, Lambda)
                        * GeneralizedMoshinskyBracket(
                            n1_dot, l1_dot-1, n2_dot, l2_dot-1,
                            n1-1, l1, n2, l2,
                            Lambda, d
                          )
                );
            value = prefactor * terms;
          }

        // store new value
        m.SetValue(memo_key,value);
      }

    if (trace_moshinsky) std::cerr << "   " <<  value << std::endl;

    return value;

  }
  
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
