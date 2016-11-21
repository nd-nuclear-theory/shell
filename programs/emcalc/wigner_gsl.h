/****************************************************************
  wigner_gsl.h                       

  Defines Wigner coupling and recoupling symbols as wrappers 
  for GSL angular momentum functions.  

  Includes support for HalfInt arguments.
                                  
  Created by Mark A. Caprio, University of Notre Dame, 2/16/10.

  2/16/10 (mac): initiated.

****************************************************************/

#ifndef gsl_wigner_h
#define gsl_wigner_h

// #include </usr/common/usg/gsl/1.16/pgi/gsl-1.16/gsl/gsl_sf_coupling.h>
// March 3rd 2016. Changing the path of the library..

#include </usr/common/usg/gsl/1.16/gnu/include/gsl/gsl_sf_coupling.h>

#include "halfint.h"

#include "angular_momentum.h"

// Naming convention: 
//   - Function names *not* ending in '2' accept HalfInt arguments J.
//   - Function names ending in '2' accept integer arguments 2*J.

// Wigner3J(ja,jb,jc,ma,mb,mc)
//   returns Wigner 3-J symbol
//   wrapper for gsl_sf_coupling_3j

inline 
double Wigner3J(
	const HalfInt& ja, const HalfInt& jb, const HalfInt& jc, 
	const HalfInt& ma, const HalfInt& mb, const HalfInt& mc
	)
{
	return gsl_sf_coupling_3j(
		TwiceValue(ja), TwiceValue(jb), TwiceValue(jc), 
		TwiceValue(ma), TwiceValue(mb), TwiceValue(mc)
		);
}

// ClebschGordan(ja,ma,jb,mb,jc,mc)
//   returns Clebsch-Gordan coefficient
//   wrapper for gsl_sf_coupling_3j

inline 
double ClebschGordan(
	const HalfInt& ja, const HalfInt& ma, 
	const HalfInt& jb, const HalfInt& mb, 
	const HalfInt& jc, const HalfInt& mc
	)
{
	return Hat(jc)*ParitySign(ja-jb+mc)
		*gsl_sf_coupling_3j(
		TwiceValue(ja), TwiceValue(jb), TwiceValue(jc), 
		TwiceValue(ma), TwiceValue(mb), -TwiceValue(mc)
		);
;
}

// Wigner6J(ja,jb,jc,jd,je,jf)
//   returns Wigner 6-J symbol
//   wrapper for gsl_sf_coupling_6j

inline 
double Wigner6J(
	const HalfInt& ja, const HalfInt& jb, const HalfInt& jc, 
	const HalfInt& jd, const HalfInt& je, const HalfInt& jf
	)
{
	return gsl_sf_coupling_6j(
		TwiceValue(ja), TwiceValue(jb), TwiceValue(jc), 
		TwiceValue(jd), TwiceValue(je), TwiceValue(jf)
		);
}

// Unitary6J(ja,jb,jc,jd,je,jf)
//   wrapper for gsl_sf_coupling_6j
//   returns unitary recoupling symbol for (12)3-1(23) recoupling

inline 
double Unitary6J(
	const HalfInt& ja, const HalfInt& jb, const HalfInt& jc, 
	const HalfInt& jd, const HalfInt& je, const HalfInt& jf
	)
{
	return ParitySign(ja+jb+jd+je)*Hat(jc)*Hat(jf)
		*gsl_sf_coupling_6j(
		TwiceValue(ja), TwiceValue(jb), TwiceValue(jc), 
		TwiceValue(jd), TwiceValue(je), TwiceValue(jf)
		);
}

// Wigner9J(ja,jb,jc,jd,je,jf,jg,jh,ji)
//   returns Wigner 9-J symbol
//   wrapper for gsl_sf_coupling_9j

inline 
double Wigner9J(
	const HalfInt& ja, const HalfInt& jb, const HalfInt& jc, 
	const HalfInt& jd, const HalfInt& je, const HalfInt& jf,
	const HalfInt& jg, const HalfInt& jh, const HalfInt& ji
	)
{
	return gsl_sf_coupling_9j(
		TwiceValue(ja), TwiceValue(jb), TwiceValue(jc), 
		TwiceValue(jd), TwiceValue(je), TwiceValue(jf),
		TwiceValue(jg), TwiceValue(jh), TwiceValue(ji)
		);
}

// Unitary9J(ja,jb,jc,jd,je,jf,jg,jh,ji)
//   returns unitary 9-J symbol
//   wrapper for gsl_sf_coupling_9j

inline 
double Unitary9J(
	const HalfInt& ja, const HalfInt& jb, const HalfInt& jc, 
	const HalfInt& jd, const HalfInt& je, const HalfInt& jf,
	const HalfInt& jg, const HalfInt& jh, const HalfInt& ji
	)
{
	return Hat(jc)*Hat(jf)*Hat(jg)*Hat(jh)
		*gsl_sf_coupling_9j(
			TwiceValue(ja), TwiceValue(jb), TwiceValue(jc), 
			TwiceValue(jd), TwiceValue(je), TwiceValue(jf),
			TwiceValue(jg), TwiceValue(jh), TwiceValue(ji)
			);
}

#endif
