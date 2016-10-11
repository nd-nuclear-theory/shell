/****************************************************************
  oscillator_parameters.h

  Defines NCSM scaling parameters, using MFDn values for physical constants.
                                  
  Mark A. Caprio
  University of Notre Dame

  3/14/12 (mac): Initiated (mfdn_scaling).
  4/25/15 (mac): Source file reformatted.
  10/11/16 (mac): Renamed oscillator_parameters.

****************************************************************/

#ifndef OSCILLATOR_PARAMETERS_H_
#define OSCILLATOR_PARAMETERS_H_

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
