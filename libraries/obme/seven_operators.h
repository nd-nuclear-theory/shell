/****************************************************************
  @file seven_operators.h

  Defines functions for computing single-particle reduced matrix
  elements of the "seven operators" arising in semileptonic
  electroweak interactions with nuclei.

  These operators appear in the multipole expansion of nuclear
  currents at order (v/c):
    M_J(qr)          -- Coulomb/charge operator          (normal parity)
    Δ_J(qr)          -- longitudinal operator            (abnormal parity)
    Δ'_J(qr)         -- transverse electric operator     (normal parity)
    Σ_J(qr)          -- spin operator                    (normal parity)
    Σ'_J(qr)         -- transverse magnetic operator     (abnormal parity)
    Σ''_J(qr)        -- spin-longitudinal operator       (abnormal parity)
    Ω_J(qr)          -- spin-velocity operator           (abnormal parity)
    Ω'_J(qr)         -- combined spin-velocity           (abnormal parity)


  References:
    [1] Donnelly & Haxton, ADNDT 23 (1979) 103.
        "Multipole operators in semileptonic weak and electromagnetic
         interactions with nuclei: HO single-particle matrix elements."
    [2] Haxton & Lunardini, CPC 179 (2008) 345.
        "SevenOperators, a Mathematica script for HO nuclear matrix
         elements arising in semileptonic electroweak interactions."
    [3] Serot, NPA 308 (1978) 457.
        "Semileptonic weak and electromagnetic interactions with
         nuclei: Nuclear current operators through order (v/c)^2."


  Victor Duménil
  University of Notre Dame & LPC Caen
  

  TO DO : 
     - Add the other 'seven operators' in ref [3]

****************************************************************/

#ifndef OBME_SEVEN_OPERATORS_H_
#define OBME_SEVEN_OPERATORS_H_

#include "basis/nlj_orbital.h"
#include "obme/obme_operator.h"

namespace shell {


////////////////////////////////////////////////////////////////
// Seven operator type enumeration
////////////////////////////////////////////////////////////////

enum class SevenOperatorType {
  kMJ,       // M_J(qr)    -- normal parity
  kDeltaJ,   // Δ_J(qr)    -- abnormal parity
  kDeltaJP,  // Δ'_J(qr)   -- normal parity
  kSigmaJ,   // Σ_J(qr)    -- normal parity
  kSigmaJP,  // Σ'_J(qr)   -- abnormal parity
  kSigmaJPP, // Σ''_J(qr)  -- abnormal parity
  kOmegaJ,   // Ω_J(qr)    -- abnormal parity
  kOmegaJP,  // Ω'_J(qr)   -- abnormal parity
};


// Return parity change g0 = (l_f + l_i + J) mod 2 expected for a
// given operator type and rank J:
//   0  =>  normal parity  (MJ, DeltaJP, SigmaJ)
//   1  =>  abnormal parity (DeltaJ, SigmaJP, SigmaJPP, OmegaJ, OmegaJP)
int SevenOperatorParityChange(SevenOperatorType operator_type, int J);


////////////////////////////////////////////////////////////////

double SevenOperator(
    SevenOperatorType operator_type,
    int ni, int li, double ji, double bi,
    int nf, int lf, double jf, double bf,
    int J, double q
  );
// Compute the reduced matrix element
//    <nf lf jf || T_J(q) || ni li ji>
// in the harmonic oscillator basis.
//
// Arguments:
//   operator_type : which of the eight SevenOperators to evaluate
//   ni, li, ji    : ket quantum numbers (n, l, j)
//   nf, lf, jf    : bra quantum numbers (n', l', j')
//   bi, bf        : harmonic oscillator length parameter
//   J             : operator rank
//   q             : momentum transfer (units: 1/b if b given explicitly)
//
// Returns 0 if the parity / triangle selection rules are not satisfied.


void SevenOperatorsOneBodyOperator(
    SevenOperatorType operator_type,
    int J,
    double q,
    double b,
    const basis::OrbitalSpaceLJPN& space,
    const basis::OrbitalSectorsLJPN& sectors,
    basis::OperatorBlocks<double>& matrices
  );
// Fill reduced matrix elements of one of the seven electroweak
// single-particle operators into the standard sector/matrix structure.
//
// The caller must provide sectors constructed with
//   J0  = J
//   g0  = SevenOperatorParityChange(operator_type, J)
//   Tz0 = 0
//
// For each sector (bra subspace lf,jf; ket subspace li,ji) the
// matrix element matrices[s](row, col) is set to
//    <nf lf jf || T_J(q) || ni li ji>
// where row = index of nf in bra subspace, col = index of ni in ket.
//
// RMEs are in Rose convention (consistent with the rest of obme/).
//
// Arguments:
//   operator_type : which SevenOperator to evaluate
//   J             : operator rank
//   q             : momentum transfer
//   b             : harmonic oscillator length parameter
//   space         : one-body orbital space
//   sectors       : pre-built sectors with matching J0/g0/Tz0
//   matrices      : output operator blocks (resized and filled here)
	
} // end namespace shell
#endif  // OBME_SEVEN_OPERATORS_H_
