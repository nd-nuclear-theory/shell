/****************************************************************
  mfdn_scaling.h                       

  Defines NCSM scaling parameters, using MFDn values for physical constants.
                                  
  Mark A. Caprio, University of Notre Dame.

  3/14/12 (mac): Initiated.
  4/25/15 (mac): Source file reformatted.

****************************************************************/

#ifndef mfdn_scaling_h
#define mfdn_scaling_h

#include <cmath>
namespace shell {

  extern const double kMFDn_mc2;  // nucleon mass
  extern const double kMFDn_hc;   // hbar-c

  // OscillatorLength returns the oscillator length parameter (in fm)
  //    which is inverse to the oscillator scale parameter

  inline double OscillatorLength(double hw) 
  {
    return kMFDn_hc / sqrt(kMFDn_mc2 * hw);
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
