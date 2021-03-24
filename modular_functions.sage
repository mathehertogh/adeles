"""
This file contains the implementation of some modular functions.
"""


def weber_f(x):
    r"""
    Evaluate Weber's `f` function in ``x``
    """
    zeta48 = CF(exp(2*CF(pi)*CF(I)/48))
    return eta((x+1)/2) / (zeta48 * eta(x))

def weber_gamma_2(x):
    r"""
    Evaluate Weber's `\gamma_2` function in ``x``
    """
    fx = weber_f(x)
    return (fx^24 - 16) / fx^8