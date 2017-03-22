"""operators.py -- define two-body operators for h2mixer input

    - 2/18/17 (pjf): Created.
"""
import mcscript.utils
from . import utils

################################################################
# identity operator
################################################################

def identity():
    return mcscript.utils.CoefficientDict(identity=1.)

################################################################
# radial kinematic operators
################################################################

kinematic_operator_set = {
    "identity", "Ursqr", "Vr1r2", "Uksqr", "Vk1k2"
}

def Ursqr():
    return mcscript.utils.CoefficientDict(Ursqr=1.)

def Vr1r2():
    return mcscript.utils.CoefficientDict(Vr1r2=1.)

def Uksqr():
    return mcscript.utils.CoefficientDict(Uksqr=1.)

def Vk1k2():
    return mcscript.utils.CoefficientDict(Vk1k2=1.)

################################################################
# angular momentum operators
################################################################

angular_momentum_operator_set = {
    "L", "Sp", "Sn", "S", "J"
}

def L():
    return mcscript.utils.CoefficientDict(L=1.)

def Sp():
    return mcscript.utils.CoefficientDict(Sp=1.)

def Sn():
    return mcscript.utils.CoefficientDict(Sn=1.)

def S():
    return mcscript.utils.CoefficientDict(S=1.)

def J():
    return mcscript.utils.CoefficientDict(J=1.)

################################################################
# interactions
################################################################

def VNN():
    return mcscript.utils.CoefficientDict(VNN=1.)

def VC_unscaled():
    return mcscript.utils.CoefficientDict(VC_unscaled=1.)

def VC(bsqr_coul=1.0):
    return VC_unscaled() / math.sqrt(bsqr_coul)

################################################################
# common observables
################################################################

def rrel2(A, hw, **kwargs):
    """Relative r^2 two-body operator.

    Arguments:
        A (int): mass number
        hw (float): length parameter

    Returns:
        CoefficientDict containing coefficients for rrel2 operator.
    """
    out = mcscript.utils.CoefficientDict()
    out += ((A-1)*(utils.oscillator_length(hw)/A)**2) * Ursqr()
    out += (-2*(utils.oscillator_length(hw)/A)**2) * Vr1r2()
    return out

def Ncm(A, bsqr, **kwargs):
    """Number of oscillator quanta in the center-of-mass.

    Arguments:
        A (int): mass number
        bsqr (float): beta squared

    Returns:
        CoefficientDict containing coefficients for Ncm operator.
    """
    out = mcscript.utils.CoefficientDict()
    out += (1/(2*A*bsqr)) * Ursqr()
    out += (1/(A*bsqr)) * Vr1r2()
    out += ((1/(2*A))*bsqr) * Uksqr()
    out += ((1/A)*bsqr) * Vk1k2()
    out += -3/2 * identity()
    return out

def Ntotal(A, bsqr, **kwargs):
    """Total number of oscillator quanta.

    Arguments:
        A (int): mass number
        bsqr (float): beta squared

    Returns:
        CoefficientDict containing coefficients for N operator.
    """
    out = mcscript.utils.CoefficientDict()
    out += (1/(2*bsqr)) * Ursqr()
    out += ((1/2)*bsqr) * Uksqr()
    out += (-3/2*A) * identity()

def Nintr(A, bsqr, **kwargs):
    """Number of oscillator quanta in the intrinsic frame.

    Arguments:
        A (int): mass number
        bsqr (float): beta squared

    Returns:
        CoefficientDict containing coefficients for Nintr operator.
    """
    return Ntotal(A, bsqr) - Ncm(A, bsqr)

def Trel(A, hw, **kwargs):
    """Two-body kinetic energy operator.

    Arguments:
        A (int): mass number
        bsqr (float): beta squared

    Returns:
        CoefficientDict containing coefficients for Trel operator.
    """
    out = mcscript.utils.CoefficientDict()
    out += ((A-1)/(2*A)*hw) * Uksqr()
    out += (-1/A*hw) * Vk1k2()
    return out

################################################################
# standard Hamiltonian
################################################################

def Hamiltonian(A, hw, a_cm, bsqr_intr=1.0, use_coulomb=True, bsqr_coul=1.0, **kwargs):
    """A standard Hamiltonian for shell-model runs.

    Arguments:
        A (int): mass number
        hw (float): oscillator basis parameter
        a_cm (float): Lawson term coefficient
        bsqr_intr (float, default 1.0): beta-squared for Lawson term
        use_coulomb (bool, default True): include Coulomb interaction
        bsqr_coul (float, optional): beta-squared for Coulomb scaling
    """
    kinetic_energy = Trel(A, hw)
    lawson_term = a_cm * Ncm(A, bsqr_intr)
    interaction = VNN()
    if use_coulomb:
        coulomb_interaction = VC(bsqr_coul)
    else:
        coulomb_interaction = mcscript.utils.CoefficientDict()
    return (kinetic_energy + interaction + coulomb_interaction + lawson_term)
