r"""
Some Modular Functions and their Actions

This file implements (the numerical evaluation of) the modular functions

- `\gamma_2` -- the cube root of the `j`-invariant satisfying
  `\gamma_2(i) = 12`, which is a modular function of level `3`;
- `\mathfrak{f}` -- Weber's modular `\mathfrak{f}` function of level `48`;
- `\mathfrak{f_1}` -- Weber's modular `\mathfrak{f_1}` function of level `48`;
- `\mathfrak{f_2}` -- Weber's modular `\mathfrak{f_2}` function of level `48`.

We also implement functions printing the action of `GL_2(\ZZ/3\ZZ)` on
`\gamma_2` and the action of `GL_2(\ZZ/48\ZZ)` on `\mathfrak{f}` and
`\mathfrak{f_2}`.

REFERENCES:

[Her2021] Mathé Hertogh, Computing with adèles and idèles, master's thesis,
Leiden University, 2021.

This file is part of the master's thesis [Her2021]. See Chapter 9 for context
on the utility of these functions, in particular Section 9.3.4.

AUTHORS:

- Mathé Hertogh (2021-07): initial version based on [Her2021]
"""

# ****************************************************************************
#       Copyright (C) 2021 Mathé Hertogh <m.c.hertogh@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.groups.free_group import FreeGroup
from sage.categories.homset import Hom
from sage.rings.complex_mpfr import ComplexField
from sage.misc.functional import eta
from sage.functions.other import sqrt
from sage.modular.arithgroup.congroup_sl2z import SL2Z
from sage.modular.arithgroup.arithgroup_perm import sl2z_word_problem
from sage.modular.local_comp.liftings import lift_matrix_to_sl2z
from sage.matrix.special import diagonal_matrix
from sage.misc.functional import det
from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.functions.generalized import sign


def gamma_2(x, prec=53):
    r"""
    Evaluate the modular function `\gamma_2` of level `3` in ``x``

    The function `\gamma_2` is the cube root of the `j`-invariant satisfying
    `\gamma_2(i) = 12`. It can be given in terms of Weber's `f` function (cf.
    :func:`weber_f`) as follows, for `z \in \CC` with `Im(z)>0`:

    .. MATH::

        \gamma_2(z) = (f(z)^{24} - 16) / f(z)^8.

    INPUT:
    
    - ``x`` -- complex point in the upper half plane
    - ``prec`` -- positive integer (default: 53); number of bits of precision to
      use

    OUTPUT:

    The value `\gamma_2(x)` as an element of ``ComplexField(prec)``.

    EXAMPLES::

        sage: gamma_2(I)
        12.0000000000000 + 6.72179444308389e-15*I
        sage: gamma_2(1+I, prec=20)
        -6.0000 - 10.392*I
        sage: gamma_2(-1.2345+10*I, prec=100)
        -1.0590687206704318964160589464e9 + 6.5818683636617232191130875508e8*I
    """
    CC = ComplexField(prec)
    x = CC(x)
    fx = weber_f(x, prec)
    return (fx**CC(24) - CC(16)) / fx**CC(8)

def weber_f(x, prec=53):
    r"""
    Evaluate Weber's `f` function in ``x``

    Weber's `f` function is a modular function of level 48 given for
    `z \in \CC` with `Im(z)>0` by

    .. MATH::

        f(z) = \zeta_{48}^{-1} \eta((z+1)/2) / \eta(z)

    where `\zeta_{48} = \exp(2 \pi i/48)` and `\eta` denotes the Dedekind
    `\eta`-function (cf. :func:`sage.misc.functional.eta`).

    INPUT:
    
    - ``x`` -- complex point in the upper half plane
    - ``prec`` -- positive integer (default: 53); number of bits of precision to
      use

    OUTPUT:

    The value `f(x)` as an element of ``ComplexField(prec)``.

    EXAMPLES::

        sage: weber_f(I)
        1.18920711500272 + 2.77555756156289e-17*I
        sage: weber_f(1+I, prec=20)
        1.0812 - 0.14234*I
        sage: weber_f(-1.2345+10*I, prec=100)
        3.6542217197057340054214271813 + 0.59570067343661179281710012847*I
    """
    CC = ComplexField(prec)
    x = CC(x)
    return eta((x+CC(1))/CC(2)) / (CC.zeta(48) * eta(x))

def weber_f1(x, prec=53):
    r"""
    Evaluate Weber's `f_1` function in ``x``

    Weber's `f_1` function is a modular function of level 48 given for
    `z \in \CC` with `Im(z)>0` by

    .. MATH::

        f_1(z) = \eta(x/2) / \eta(x)

    where `\eta` denotes the Dedekind `\eta`-function (cf.
    :func:`sage.misc.functional.eta`).

    INPUT:
    
    - ``x`` -- complex point in the upper half plane
    - ``prec`` -- positive integer (default: 53); number of bits of precision to
      use

    OUTPUT:

    The value `f_1(x)` as an element of ``ComplexField(prec)``.
    
    EXAMPLES::

        sage: weber_f1(I)
        1.09050773266526
        sage: weber_f1(1+I, prec=20)
        1.1790 - 0.15522*I
        sage: weber_f1(-1.2345+10*I, prec=100)
        3.6542217197058751251562566155 + 0.59570067343652031920472325046*I
    """
    CC = ComplexField(prec)
    x = CC(x)
    return eta(x/CC(2)) / eta(x)

def weber_f2(x, prec=53):
    r"""
    Evaluate Weber's `f_2` function in ``x``

    Weber's `f_2` function is a modular function of level 48 given for
    `z \in \CC` with `Im(z)>0` by

    .. MATH::

        f_2(x) = \sqrt{2} \eta(2x) / \eta(x)

    where `\eta` denotes the Dedekind `\eta`-function (cf.
    :func:`sage.misc.functional.eta`).

    INPUT:
    
    - ``x`` -- complex point in the upper half plane
    - ``prec`` -- positive integer (default: 53); number of bits of precision to
      use

    OUTPUT:

    The value `f_2(x)` as an element of ``ComplexField(prec)``.

    EXAMPLES::

        sage: weber_f2(I)
        1.09050773266526
        sage: weber_f2(1+I, prec=20)
        1.0533 + 0.28224*I
        sage: weber_f2(-1.2345+10*I, prec=100)
        0.097824329747547433152328660118 - 0.032764790050361566443041641916*I
    """
    CC = ComplexField(prec)
    x = CC(x)
    sqrt2 = CC(sqrt(2))
    return sqrt2 * eta(CC(2) * x) / eta(x)

def apply_fractional_linear_transformation(M, tau):
    r"""
    Apply the fractional linear transformation given by ``M`` to ``tau``

    INPUT:

    - ``M`` -- a matrix in `GL_2(\QQ)` with positive determinant
    - ``tau`` -- a non-zero field element

    OUTPUT:
    
    Writing ``M = (a, b; c, d)``, return the element ``(a*tau+b)/(c*tau+d)``.

    EXAMPLES::

        sage: M = matrix([[1, 0], [0, 2]])
        sage: apply_fractional_linear_transformation(M, I)
        1/2*I
        sage: N = matrix([[1, 2], [3, 4]])
        sage: apply_fractional_linear_transformation(N, I)
        -2/25*I + 11/25
    """
    a, b, c, d = M.list()
    return (a*tau + b) / (c*tau + d)

def ST_factor(A, return_homomorphism=False):
    r"""
    Factor `A \in SL_2(\ZZ/N\ZZ)` into a product of the standard generators `S`
    and `T`

    Here we have `S = (0, -1; 1, 0)` and `T = (1, 1; 0, 1)`. Together they
    generate `SL_2(\ZZ)` and hence also `SL_2(\ZZ/N\ZZ)` for any integer `N`.

    INPUT:

    - ``A`` -- a matrix in `SL_2(\ZZ/N\ZZ)` for some integer `N`
    - ``return_homomorphism`` -- boolean (default: ``False``); whether or not to
      return the homomorphism `f` described below as well

    OUPUT:

    An element of the free multiplicative group G generated by `S` and `T` which
    is mapped to ``A`` by the homomorphism `f: G \to SL_2(\ZZ/N\ZZ)` that maps S
    to (0, -1; 1, 0) and T to (1, 1; 0, 1).

    If ``return_homomorphism`` is ``True``, also returns `f`.

    EXAMPLES::

        sage: G = SL(2, Zmod(10))
        sage: A = G([7, 5, 8, 3]); A
        [7 5]
        [8 3]
        sage: STs, f = ST_factor(A, True); STs
        T^7*(T*(S^3*T^-1*S)^2)^2
        sage: STs.parent()
        Free Group on generators {S, T}
        sage: f
        Group morphism:
          From: Free Group on generators {S, T}
          To:   Special Linear Group of degree 2 over Ring of integers modulo 10
        sage: f(STs) == A
        True

    TESTS::

        sage: A = G.random_element()
        sage: STs, f = ST_factor(A, True)
        sage: f(STs) == A
        True
    """
    SL2ZmodN = A.parent()
    N = A.base_ring().order()
    G = FreeGroup(names=('S', 'T',)); (S, T,) = G._first_ngens(2)

    # First we convert A to a list of 4 elements and lift it to SL_2(ZZ)
    A = A.list()
    if len(A) != 4:
        A = sum([row for row in A], [])
    A_lift = SL2Z(lift_matrix_to_sl2z(A, N))

    # We want to factor A_lift into a product of S's and T's. The function
    # sl2z_word_problem() already factors it into L's and R's for us, where L=T
    # and R=(1, 0; 1, 1)=S^3*T^-1*S.
    L = T
    R = S**ZZ(3)  * T**ZZ(-1)  * S
    factorization = G.one()
    for is_R, e in sl2z_word_problem(A_lift):
        if is_R:
            factorization = factorization * R**e
        else:
            factorization = factorization * L**e

    if not return_homomorphism:
        return factorization

    # Construct the homomorphism f: G --> SL_2(Z/NZ) with f(S) = (0, -1; -1, 0)
    # and f(T) = (1, 1; 0, 1).
    Sm = SL2ZmodN([0, -1, 1, 0])
    Tm = SL2ZmodN([1, 1, 0, 1])
    f = Hom(G, SL2ZmodN)([Sm, Tm])

    return factorization, f

def print_action_on_gamma_2(B):
    r"""
    Print the action of `B \in GL_2(\ZZ/3\ZZ)` on the function `\gamma_2`

    See :func:`\gamma_2` for information on the modular function `\gamma_2` of
    level `3`.

    INPUT:

    - ``B`` -- a `GL_2(\ZZ/3\ZZ)`-matrix

    ALGORITHM:

    We factor `B` as `B = (1, 0; 0, d) \cdot U` where with `d = \det(B)` and
    `U \in SL_2(\ZZ)`. Next we write `U` on the standard generators
    `S = (0, -1; 1, 0)` and `T = (1, 1; 0, 1)` of `SL_2(\ZZ)`. Then we use the
    fact that `(\ZZ/3\ZZ)^*` and `S` act trivially on `\gamma_2`, while `T` acts
    as multiplication by `\zeta_3 = \exp(2 \pi i/3)`.

    EXAMPLES::

        sage: S = matrix(Zmod(3), [[0, -1], [-1, 0]])
        sage: print_action_on_gamma_2(S)
          gamma_2 ]--> gamma_2
        sage: T = matrix(Zmod(3), [[1, 1], [0, 1]])
        sage: print_action_on_gamma_2(T)
          gamma_2 ]--> zeta_3^2 * gamma_2
        sage: B = matrix(Zmod(3), [[0, 1], [1, 2]])
        sage: print_action_on_gamma_2(B)
          gamma_2 ]--> zeta_3 * gamma_2

    TESTS::

        sage: B = matrix(Zmod(3), [[1, 1], [2, 2]])
        sage: print_action_on_gamma_2(B)
        Traceback (most recent call last):
        ...
        ValueError: B does not lie in GL_2(ZZ/48ZZ)
    """
    d = det(B)
    if not (B.parent().base_ring() is Zmod(3) and d.is_unit()):
        raise ValueError("B does not lie in GL_2(ZZ/48ZZ)")
    U = diagonal_matrix([1, 1/d]) * B
    STs = ST_factor(U)

    G = STs.parent()
    # As S acts trivially on gamma_2, we are only interested in the number of
    # T's in U.
    n_Ts = sum([sign(STs) for STs in STs.substitute(S=G.one()).Tietze()])
    # T acts as multiplication by zeta_3^-1 on gamma_2, so U acts as
    # multiplication by zeta^-n_Ts on gamma_2.
    e = Zmod(3)(-n_Ts)
    if e == 0:
        print("  gamma_2 ]--> gamma_2")
    elif e == 1:
        print("  gamma_2 ]--> zeta_3 * gamma_2")
    else:
        print("  gamma_2 ]--> zeta_3^2 * gamma_2")

def print_action_on_weber_f(B):
    r"""
    Print the action of `B \in GL_2(\ZZ/48\ZZ)` on Weber's `f` function

    See :func:`weber_f` for information on Weber's `f` function.

    INPUT:

    - ``B`` -- a `GL_2(\ZZ/48\ZZ)`-matrix

    ALGORITHM:

    We factor `B` as `B = (1, 0; 0, d) \cdot U` where with `d = \det(B)` and
    `U \in SL_2(\ZZ)`. Next we write `U` on the standard generators
    `S = (0, -1; 1, 0)` and `T = (1, 1; 0, 1)` of `SL_2(\ZZ)`. Then we use the
    following (hard-coded) knowledge.

    Writing `\zeta_{48} = exp(2 \pi i/48)`, the generators `S` and `T` act as
    follows on `\QQ(\zeta_{48}, f, f_1, f_2)`:

    - `S: (f, f_1, f_2) \mapsto (f, f_2, f_1)`;
    - `T: (f, f_1, f_2) \mapsto (\zeta_{48}^{-1} f_1, \zeta_{48}^{-1} f, \zeta_{48}^2 f_2)`.
    
    and `(\ZZ/48\ZZ)^*` acts trivially on `f`.

    EXAMPLES::

        sage: S = matrix(Zmod(48), [[0, -1], [-1, 0]])
        sage: print_action_on_weber_f(S)
          f       ]--> f
        sage: T = matrix(Zmod(48), [[1, 1], [0, 1]])
        sage: print_action_on_weber_f(T)
          f       ]--> zeta48^-1*f1
        sage: B = matrix(Zmod(48), [[-1, 4], [3, 7]])
        sage: print_action_on_weber_f(B)
          f       ]--> zeta48^-23*f2

    TESTS::

        sage: B = matrix(Zmod(48), [[1, 0], [0, 2]])
        sage: print_action_on_weber_f(B)
        Traceback (most recent call last):
        ...
        ValueError: B does not lie in GL_2(ZZ/48ZZ)
    """
    d = det(B)
    if not (B.parent().base_ring() is Zmod(48) and d.is_unit()):
        raise ValueError("B does not lie in GL_2(ZZ/48ZZ)")
    U = diagonal_matrix([1, 1/d]) * B
    STs = ST_factor(U)

    G = STs.parent()
    F = FreeGroup(['zeta48', 'f', 'f1', 'f2'])
    zeta48, f, f1, f2 = F.gens()
    phi_S = Hom(F, F)([zeta48, f, f2, f1])
    phi_T = Hom(F, F)([zeta48, zeta48**ZZ(-1)*f1, zeta48**ZZ(-1)*f, zeta48**ZZ(2)*f2])
    phi_T_inv = Hom(F, F)([zeta48, zeta48*f1, zeta48*f, zeta48**ZZ(-2)*f2])
    # Note that the inverse "phi_S_inv" is equal to phi_S itself

    f_image = f

    for R in STs.Tietze():
        if R in [1, -1]: # R == S or R == S^-1
            f_image = phi_S(f_image)
        elif R == 2: # R == T
            f_image = phi_T(f_image)
        else: # R == T^-1
            f_image = phi_T_inv(f_image)

    print("  f       ]--> {}".format(f_image))

def print_action_on_weber_f2(B):
    r"""
    Print the action of `B \in GL_2(\ZZ/48\ZZ)` on Weber's `f_2` function

    See :func:`weber_f2` for information on Weber's `f_2` function.

    INPUT:

    - ``B`` -- a `GL_2(\ZZ/48\ZZ)`-matrix

    ALGORITHM:

    We factor `B` as `B = (1, 0; 0, d) \cdot U` where with `d = \det(B)` and
    `U \in SL_2(\ZZ)`. Next we write `U` on the standard generators
    `S = (0, -1; 1, 0)` and `T = (1, 1; 0, 1)` of `SL_2(\ZZ)`. Then we use the
    following (hard-coded) knowledge.

    Writing `\zeta_{48} = exp(2 \pi i/48)`, the generators `S` and `T` act as
    follows on `\QQ(\zeta_{48}, f, f_1, f_2)`:

    - `S: (f, f_1, f_2) \mapsto (f, f_2, f_1)`;
    - `T: (f, f_1, f_2) \mapsto (\zeta_{48}^{-1} f_1, \zeta_{48}^{-1} f, \zeta_{48}^2 f_2)`.
    
    and the action of `(\ZZ/48\ZZ)^*` is given by `f_2^d = -f_2` for `d \equiv
    3, 5 \mod 8` and `f_2` is invariant under all other `d \in (\ZZ/48\ZZ)^*`.

    EXAMPLES::

        sage: S = matrix(Zmod(48), [[0, -1], [-1, 0]])
        sage: print_action_on_weber_f2(S)
          f2    ]--> f1
        sage: T = matrix(Zmod(48), [[1, 1], [0, 1]])
        sage: print_action_on_weber_f2(T)
          f2    ]--> zeta48^2*f2
        sage: B = matrix(Zmod(48), [[8, 3], [-33, 5]])
        sage: print_action_on_weber_f2(B)
          f2    ]--> zeta48^21*f

    TESTS::

        sage: B = matrix(Zmod(48), [[10, 33], [-21, 9]])
        sage: print_action_on_weber_f2(B)
        Traceback (most recent call last):
        ...
        ValueError: B does not lie in GL_2(ZZ/48ZZ)
    """
    d = det(B)
    if not (B.parent().base_ring() is Zmod(48) and d.is_unit()):
        raise ValueError("B does not lie in GL_2(ZZ/48ZZ)")
    U = diagonal_matrix([1, 1/d]) * B
    STs = ST_factor(U)

    G = STs.parent()
    F = FreeGroup(['zeta48', 'f', 'f1', 'f2'])
    zeta48, f, f1, f2 = F.gens()
    phi_S = Hom(F, F)([zeta48, f, f2, f1])
    phi_T = Hom(F, F)([zeta48, zeta48**ZZ(-1)*f1, zeta48**ZZ(-1)*f, zeta48**ZZ(2)*f2])
    phi_T_inv = Hom(F, F)([zeta48, zeta48*f1, zeta48*f, zeta48**ZZ(-2)*f2])
    # Note that the inverse "phi_S_inv" is equal to phi_S itself

    f2_image = f2

    if d.lift() % 8 in [3, 5]:
        f2_image = zeta48**ZZ(24) * f2

    for R in STs.Tietze():
        if R in [1, -1]: # R == S or R == S^-1
            f2_image = phi_S(f2_image)
        elif R == 2: # R == T
            f2_image = phi_T(f2_image)
        else: # R == T^-1
            f2_image = phi_T_inv(f2_image)

    print("  f2    ]--> {}".format(f2_image))

