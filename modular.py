
from sage.rings.complex_mpfr import ComplexField
from sage.functions.log import exp
from sage.symbolic.constants import pi, I
from sage.misc.functional import eta

def weber_f(x, prec=53):
    r"""
    Evaluate Weber's `f` function in ``x`` with ``prec`` bits of precision
    """
    CC = ComplexField(prec)
    x = CC(x)
    zeta48 = CC(exp(2*CC(pi)*CC(I)/48))
    return eta((x+CC(1))/CC(2)) / (zeta48 * eta(x))

def weber_f1(x, prec=53):
    r"""
    Evaluate Weber's `f_1` function in ``x`` with ``prec`` bits of precision
    """
    CC = ComplexField(prec)
    x = CC(x)
    return eta(x/CC(2)) / eta(x)

def weber_f2(x, prec=53):
    r"""
    Evaluate Weber's `f_2` function in ``x`` with ``prec`` bits of precision
    """
    CC = ComplexField(prec)
    x = CC(x)
    sqrt2 = CC(sqrt(2))
    return sqrt2 * eta(CC(2) * x) / eta(x)

def gamma_2(x, prec=53):
    r"""
    Evaluate Weber's `\gamma_2` function in ``x`` with ``prec`` bits of precision
    """
    CC = ComplexField(prec)
    fx = weber_f(x, prec)
    return (fx**CC(24) - CC(16)) / fx**CC(8)

def print_action_on_gamma_2(STs):
    r"""
    Print the action of ST-product (as returned by ST_factor()) on gamma_2

    INPUT:

    - ``STs`` -- an element of the free group on {S, T}; usually obtained from
                 an `SL_2(\ZZ/3\ZZ)`-matrix passed into :func:`ST_factor`

    OUTPUT:

    Does not return anything.
    Prints a string describing the action.

    ALGORITHM:

    We use the fact that `S` acts trivially on `\gamma_2` and `T` acts as
    multiplication by `\zeta_3^{-1}`, where `\zeta_3 = exp(2i\pi/3)`.
    """
    from sage.rings.finite_rings.integer_mod_ring import Zmod
    from sage.functions.generalized import sign
    G = STs.parent()
    # S acts trivially on gamma_2, so we are only interested in the number of
    # T's in U.
    n_Ts = sum([sign(STs) for STs in STs.substitute(S=G.one()).Tietze()])
    # T acts as multiplication by zeta^-1 on gamma_2, so U acts as
    # multiplication by zeta^-n_Ts on gamma_2
    e = Zmod(3)(-n_Ts)
    if e == 0:
        print("  gamma_2 ]--> gamma_2")
    elif e == 1:
        print("  gamma_2 ]--> zeta_3 * gamma_2")
    else:
        print("  gamma_2 ]--> zeta_3^2 * gamma_2")

def print_action_on_f2(STs):
    r"""
    Print the action of ST-product (as returned by ST_factor()) on Weber's f_2
    modular function

    INPUT:

    - ``STs`` -- an element of the free group on {S, T}; usually obtained from
                 an `SL_2(\ZZ/48\ZZ)`-matrix passed into :func:`ST_factor`

    OUTPUT:

    Does not return anything.
    Prints a string describing the action.

    ALGORITHM:

    We use the fact that `S` and `T` act as follows on
    `\QQ(\zeta_48, f, f_1, f_2)`, with `\zeta_3 = exp(2i\pi/48)`:

    - `S: (f, f_1, f_2) \mapsto (f, f_2, f_1)`;
    - `T: (f, f_1, f_2) \mapsto (\zeta_48^{-1} f_1, \zeta_48^{-1} f, \zeta_48^2 f_2)`.
    """
    from sage.rings.finite_rings.integer_mod_ring import Zmod
    from sage.functions.generalized import sign
    G = STs.parent()
    F = FreeGroup(['zeta48', 'f', 'f1', 'f2'], abelian=True)
    F.inject_variables()
    phi_S = Hom(F, F)([zeta48, f, f2, f1])
    phi_T = Hom(F, F)([zeta48, zeta48**ZZ(-1)*f1, zeta48**ZZ(-1)*f, zeta48**ZZ(2)*f2])
    phi_T_inv = Hom(F, F)([zeta48, zeta48*f1, zeta48*f, zeta48**ZZ(-2)*f2])
    # Note that the inverse "phi_S_inv" is equal to phi_S itself

    f2_image = f2
    for R in STs.Tietze():
        if R in [1, -1]: # R == S^+/-1
            f2_image = phi_S(f2_image)
        elif R == 2: # R == T
            f2_image = phi_T(f2_image)
        else: # R == T^-1
            f2_image = phi_T_inv(f2_image)

    print("  f2 ]--> {}".format(f2_image))
