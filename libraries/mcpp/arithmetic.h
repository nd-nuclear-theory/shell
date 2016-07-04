/****************************************************************
  arithmetic.h                       

  Defines arithmetic shorthands.
                                  
  Created by Mark A. Caprio on 2/17/11.   
  Drawing upon libmc/mcutils C code.

  2/23/11 (mac): Renamed from mc_arithmetic to arithmetic.

****************************************************************/

#ifndef ARITHMETIC_H
#define ARITHMETIC_H

// ONLYIF(cond,x) evaluates and returns x only if cond is true
#define ONLYIF(cond,x) ( (cond) ? (x) : 0)    

// sqr(x) returns the arithmetic square of x by self-multiplication
//   Note: Use of inline template avoids double evaluation of x which
//   would occur in a macro implementation.

template <typename T>
inline
T sqr(const T& x) 
{
	return x*x;
}


#endif
