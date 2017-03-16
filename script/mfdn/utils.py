"""utils.py -- helper functions for MFDn runs

    - 2/20/17 (pjf): Created, extracted from mfdn.py.
"""

import math

################################################################
# physical constants
################################################################

k_mN_csqr = 938.92  # (m_N c^2)~938.92 MeV
k_hbar_c = 197.327  # (hbar c)~197.327 Mev fm

################################################################
# utility calculations
################################################################

def weight_max_string(truncation):
    """ Convert (rank,cutoff) to "wp wn wpp wnn wpn" string.

    Valid truncations:
        ("ob",w1b)
        ("tb",w2b)
        (w1b,w2b)
        (wp,wn,wpp,wnn,wpn) -- TODO

    >>> weight_max_string(("ob",4))
        "4 4 8 8 8"
    >>> weight_max_string(("tb",4))
        "4 4 4 4 4"
    """

    if (truncation[0] == "ob"):
        (code, N) = truncation
        cutoffs = (N,N,2*N,2*N,2*N)
    elif (truncation[0] == "tb"):
        (code, N) = truncation
        cutoffs = (N,N,N,N,N)
    elif (len(truncation)==2):
        (w1,w2) = truncation
        cutoffs = (w1,w1,w2,w2,w2)
    elif (len(truncation)==5):
        cutoffs = truncation

    return "{0[0]} {0[1]} {0[2]} {0[3]} {0[4]}".format(cutoffs)

def oscillator_length(hw):
    """ Calculate oscillator length for given oscillator frequency.

    b(hw) = (hbar c)/[(m_N c^2) (hbar omega)]^(1/2)

    Arguments:
        hw (numeric): hbar omega in MeV

    Returns:
        (float): b in fm
    """

    return k_hbar_c/math.sqrt(k_mN_csqr*hw)

def natural_orbital_indicator(natural_orbital_iteration):
    """Construct natural orbital indicator string."""
    if (natural_orbital_iteration is None):
        indicator = ""
    else:
        indicator = "-no{:1d}".format(natural_orbital_iteration)
    return indicator
