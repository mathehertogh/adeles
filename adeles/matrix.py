r"""
Adèlic Matrix Factorization

This file implements the adèlic matrix factorization algorithms described in
Chapter 8 of [Her2021].

The function :func:`factor_GLQhat` implements Algorithm 8.4 of [Her2021], which
factors an adèlic matrix `M \in GL_n(\hat{\QQ})` as `M = BA` with
`B \in GL_n(\hat{\ZZ})` and `A \in GL_n^+(\QQ)`.

The function :func:`factor_GSpQhat` implements Algorithm 8.9 of [Her201], which
factors a general symplectic matrix `M \in GSp_{2g}(\hat{\QQ})` as `M = BA` with
`B \in GSp_{2g}(\hat{\ZZ})` and `A \in GSp_{2g}^+(\QQ)`.

REFERENCES:

[Her2021] Mathé Hertogh, Computing with adèles and idèles, master's thesis,
Leiden University, 2021.

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

from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.arith.functions import lcm
from sage.arith.misc import gcd
from sage.matrix.constructor import matrix
from sage.matrix.special import identity_matrix, block_matrix, diagonal_matrix
from sage.matrix.symplectic_basis import symplectic_basis_over_ZZ
from sage.misc.functional import det

from .profinite_number import Qhat


def value_matrix(M):
    r"""
    Return the value matrix of the `\hat{\QQ}`-matrix ``M``

    This is the matrix consisting of the values of the entries of ``M``.

    EXAMPLES::

        sage: M = matrix(Qhat, [[Qhat(1, 7/6), Qhat(1/3, 6)],
        ....:                   [ Qhat(0, 12),  Qhat(5, 18)]])
        sage: value_matrix(M)
        [  1 1/3]
        [  0   5]
    """
    return matrix ([[entry.value() for entry in row] for row in M])

def denominator_matrix(M):
    r"""
    Return the denominator matrix of the `\hat{\QQ}`-matrix ``M``

    This is the diagonal matrix `D` with `D_{jj}` equal to to least common
    multiple of the denominators of the entries of the `j`-th column of ``M``.

    EXAMPLES::

        sage: M = matrix(Qhat, [[Qhat(1, 7/6), Qhat(1/3, 6)],
        ....:                   [ Qhat(0, 12),  Qhat(5, 18)]])
        sage: denominator_matrix(M)
        [6 0]
        [0 3]
    """
    diag = [lcm([entry.denominator() for entry in column]) for column in M.columns()]
    return diagonal_matrix(diag)

def denominator(M):
    r"""
    Return the denominator of the `\hat{\QQ}`-matrix ``M``

    This is the least common multiple of the denominators of the entries of
    ``M``.

    EXAMPLES::

        sage: M = matrix(Qhat, [[Qhat(1, 7/6), Qhat(1/3, 6)],
        ....:                   [ Qhat(0, 12),  Qhat(5, 18)]])
        sage: denominator(M)
        6
    """
    return lcm([entry.denominator() for entry in M.list()])

def modulus(M):
    r"""
    Return the modulus (i.e. precision) of the `\hat{\QQ}`-matrix ``M``

    This is the greatest common divisor of the moduli of the entries of ``M``.

    EXAMPLES::

        sage: M = matrix(Qhat, [[Qhat(1, 7/6), Qhat(1/3, 6)],
        ....:                   [ Qhat(0, 12),  Qhat(5, 18)]])
        sage: modulus(M)
        1/6
    """
    return gcd([entry.modulus() for entry in M.list()])

def matrix_modulo(A, N):
    r"""
    Project a matrix over `\hat{\ZZ}` to a matrix over `\ZZ/N\ZZ`

    INPUT:

    - ``A`` -- a matrix with entries *integral* elements of `\hat{\QQ}` and with
      modulus divisible by ``N``
    - ``N`` -- an integer
    """
    rows_modN = []
    for row in A:
        row_modN = []
        for a in row:
            if not N.divides(a.modulus()):
                raise ValueError("not every entry of A is defined modulo {}".format(N))
            row_modN.append(Zmod(N)(a.value()))
        rows_modN.append(row_modN)
    return matrix(Zmod(N), rows_modN)

def has_good_precision(M, Delta):
    r"""
    Return wether or not the sqaure matrix ``M`` over ``Qhat`` has good
    precision with respect to ``Delta``

    By `M` having *good precision* with respect to `\Delta \in \QQ_{>0}` we mean
    the following. Let `D` be the denominator matrix of `M` and let `m`
    be the precision (i.e. modulus) of `MD`. Then `M` has *good precision* with
    respect to `\Delta` if for all prime numbers `p` we have
    `ord_p(m) \geq ord_p(\Delta \det(D))`.

    INPUT:

    - ``M`` -- a square matrix over the profinite `\QQ`-integers
    - ``Delta`` -- a positive rational number

    EXAMPLES::

        sage: M = matrix([[Qhat(1/3, 10/7), Qhat(4, 20)],
        ....:             [Qhat(1/5, 20/7), Qhat(5/3, 30)]])
        sage: has_good_precision(M, 2/21)
        True
        sage: has_good_precision(M, 4/21)
        False
        sage: has_good_precision(M, 2/7)
        False
        sage: has_good_precision(M, 2/3)
        False
        sage: has_good_precision(M, 1/21)
        True
    """
    D = denominator_matrix(M)
    m = modulus(M * D)
    n = Delta * det(D)
    d = lcm(m.denominator(), n.denominator())
    m, n = ZZ(d*m), ZZ(d*n)
    return n.divides(m)

def factor_GLQhat(M, detM):
    r"""
    Factor the `GL_n(\hat{\QQ})`-matrix `M` as `M = BA` with
    `B \in GL_n(\hat{\ZZ})` and `A \in GL_n^+(\QQ)`

    Here `GL_n` denotes the general linear group, i.e. the group of invertible
    `n \times n`-matrices. `GL_n^+(\QQ)` denotes the subgroup of `GL_n(\QQ)`
    of matrices with positive determinant.

    INPUT:

    - ``M`` -- an `n \times n`-matrix over the ring of profinite `\QQ`-numbers
      having good precision with respect to ``detM``
      (cf. :func:`has_good_precision`) and representing at least one element of
      `GL_n(\hat{\QQ})` of determinant `d` satisfying
      `d\hat{\ZZ} = detM\hat{\ZZ}`.
    - ``detM`` -- a positive rational number

    OUTPUT:

    A matrix `A \in GL_n^+(\QQ)` such that for every matrix
    `N \in GL_n(\hat{\QQ})` that `M` represents and satisfying
    `\det(N)\hat{\ZZ} = detM \hat{\ZZ}` we have `N A^{-1} \in GL_n(\hat{\ZZ})`.

    ALGORITHM:

    We perform Algorithm 8.4 of [Her2021].

    EXAMPLES:

    Let `I` be the `2 \times 2`-identity matrix, i.e. `(1, 0; 0, 1)`. Let
    `N \in GL_2(\hat{\QQ})` be given by `N_2 = (1, 0; 0, 2)` and `N_p = I` for
    all odd `p` (where `N_p` denotes the projection of `N` to `GL_2(\QQ_p)`).
    Note that `\det(N)\hat{\ZZ} = 2\hat{\ZZ}`. The following matrix ``M``
    represents `N` and has good precision with respect to `2`. ::

            sage: M = matrix([[Qhat(1, 2), Qhat(0, 2)],
            ....:             [Qhat(0, 2), Qhat(0, 2)]]); M
            [1 mod 2 0 mod 2]
            [0 mod 2 0 mod 2]
            sage: has_good_precision(M, 2)
            True
    
    This function returns the following `A`::

        sage: A = factor_GLQhat(M, 2); A
        [1 0]
        [0 2]

    Now it is indeed the case that `N A^{-1} \in GL_2(\hat{\ZZ})`: at `2` this
    matrix equals `I` and at every odd prime `p` it equals the inverse of `A`
    and we have `A \in GL_2(\ZZ_p)`.

    Next we let `N \in GL_2(\hat{\QQ})` at `2` be given by ::

        sage: N_2 = matrix(QQ, [[1, 0], [0, 1/4]]); N_2
        [  1   0]
        [  0 1/4]
    
    and at `3` by ::

        sage: N_3 = matrix(QQ, [[0, 2], [3, 1/5]]); N_3
        [  0   2]
        [  3 1/5]

    and at the other primes by `N_p = I`. In this case we have
    `\det(N)\hat{\ZZ} = 3/4\hat{\ZZ}`. A matrix representing `N` and having good
    precision with respect to `3/4` is::

        sage: M = matrix([[Qhat(21, 60), Qhat(20, 60)],
        ....:             [Qhat(0, 60), Qhat(209/4, 60)]]); M
        [   21 mod 60    20 mod 60]
        [    0 mod 60 209/4 mod 60]
        sage: has_good_precision(M, 3/4)
        True

    We obtain the following `A \in GL_2^+(\QQ)`::

        sage: A = factor_GLQhat(M, 3/4); A
        [  3   0]
        [  0 1/4]
    
    We check that this is correct by checking that `N_p A^{-1} \in GL_2(\ZZ_p)`
    for all prime numbers `p`. This is clear for `p \geq 5` since we then have
    `I, A \in GL_2(\ZZ_p)`. At `2` and `3` we have ::

        sage: N_2 * ~A, N_3 * ~A
        (
        [1/3   0]  [  0   8]
        [  0   1], [  1 4/5]
        )

    which indeed lie in `GL_2(\ZZ_2)` and `GL_3(\ZZ_3)` respectively.

    Let's do an example with `n=3`; consider `N \in GL_3(\hat{\QQ})` given at
    `3` by ::

        sage: N_3 = matrix(QQ, [[14, 1/3, 1], [7, 2, 0], [5, 3, 2/3]]); N_3
        [ 14 1/3   1]
        [  7   2   0]
        [  5   3 2/3]

    and at `5` by ::

        sage: N_5 = matrix(QQ, [[15, 5, 0], [9, 11, -1], [5, -9, 1]]); N_5
        [15  5  0]
        [ 9 11 -1]
        [ 5 -9  1]
    
    and at all other primes `p` by the `3 \times 3`-identity matrix `I`. We have
    ::

        sage: det(N_3).valuation(3), det(N_5).valuation(5)
        (-2, 1)

    and so `det(N) \hat{\ZZ} = 5/9\hat{\ZZ}`. The following matrix represents
    `N` and has good precision with respect to `5/9`::

        sage: M = matrix([[Qhat([Qp(3)(N_3[i,j], 0), Qp(5)(N_5[i,j], 1)])
        ....:             for j in range(3)] for i in range(3)])
        sage: M
        [   0 mod 5 10/3 mod 5    0 mod 5]
        [   4 mod 5    1 mod 5    4 mod 5]
        [   0 mod 5    1 mod 5  8/3 mod 5]
        sage: has_good_precision(M, 5/9)
        True

    Hence we compute::

        sage: A = factor_GLQhat(M, 5/9); A
        [  1   0 1/3]
        [  0 1/3 1/3]
        [  0   0 5/3]

    and we check correctness by checking that `N_3 A^{-1} \in GL_3(\ZZ_3)` and
    `N_5 A^{-1} \in GL_5(\ZZ_5)`::

        sage: N_3 * ~A, N_5 * ~A
        (
        [   14     1 -12/5]  [ 15  15  -6]
        [    7     6 -13/5]  [  9  33  -9]
        [    5     9 -12/5], [  5 -27   5]
        )
        sage: det(N_3*~A).valuation(3) == 0
        True
        sage: det(N_5*~A).valuation(5) == 0
        True

    If the matrix does not have good precision, we throw an exception. ::

        sage: M = matrix([[Qhat(1,2), Qhat(0,2)], [Qhat(0, 2), Qhat(0,2)]])
        sage: factor_GLQhat(M, 4)
        Traceback (most recent call last):
        ...
        ValueError: matrix precision too low!
    """
    if not has_good_precision(M, detM):
        raise ValueError("matrix precision too low!")

    n = M.nrows()

    # Step 1: Compute the denominator matrix of M.
    D = denominator_matrix(M)

    # Step 2: Construct J
    MD = M * D
    detMD = ZZ(detM * det(D))
    detMD_I = detMD * identity_matrix(n)
    vMD = value_matrix(MD).change_ring(ZZ)
    J = block_matrix([[vMD], [detMD_I]])

    # Step 3: Compute the Hermite Normal Form A_0 of J.
    A_0 = J.hermite_form(include_zero_rows=False)

    # Step 4: Ensure that the determinant of A_0 is positive.
    if det(A_0) < 0:
        A_0 = diagonal_matrix([1 for i in range(n-1)] + [-1]) * A_0

    # Step 5: Compute and return A.
    A = A_0 * D.inverse()
    return A

def Omega(g):
    r"""
    Return the matrix `\Omega_{2g}`

    The matrix `\Omega_{2g}` is the `2g \times 2g`-matrix given in
    `g \times g`-blocks by `(0, 1; -1, 0)`.

    EXAMPLES::

        sage: Omega(3)
        [ 0  0  0  1  0  0]
        [ 0  0  0  0  1  0]
        [ 0  0  0  0  0  1]
        [-1  0  0  0  0  0]
        [ 0 -1  0  0  0  0]
        [ 0  0 -1  0  0  0]
    """
    I = identity_matrix(g)
    return block_matrix([[0*I, I], [-I, 0*I]], subdivide=False)

def factor_GSpQhat(M, detM):
    r"""
    Factor the `GSp_{2g}(\hat{\QQ})`-matrix `M` as `M = BA` with
    `B \in GSp_{2g}(\hat{\ZZ})` and `A \in GSp_{2g}^+(\QQ)`

    Here `GSp_{2g}(R)` denotes the general symplectic group over the ring `R`
    and `GSp_{2g}^+(\QQ)` denotes the subgroup of matrcies with positive
    mulitplier. See Chapter 8 of [Her2021] for details.

    INPUT:

    - ``M`` -- a matrix over the ring of profinite `\QQ`-numbers having good
      precision with respect to ``detM`` (cf. :func:`has_good_precision`) and
      representing at least one element in `GSp_{2g}(\hat{\QQ})` of determinant
      `d` satisfying `d\hat{\ZZ} = detM\hat{\ZZ}`
    - ``detM`` -- a positive rational number

    OUTPUT:

    A matrix `A \in GSp_{2g}^+(\QQ)` such that for every matrix
    `N \in GSp_{2g}(\hat{\QQ})` that `M` represents and satisfying
    `\det(N)\hat{\ZZ} = detM \hat{\ZZ}` we have
    `N A^{-1} \in GSp_{2g}(\hat{\ZZ})`.

    ALGORITHM:

    We perform Algorithm 8.9 of [Her2021].

    EXAMPLES:

    Consider the following symplectic matrix over `\QQ`::

        sage: I = identity_matrix(2)
        sage: S = matrix([[1, 2], [2, 3]])
        sage: N_2 = block_matrix([[I, S], [0*I, 1/2*I]], subdivide=False); N_2
        [  1   0   1   2]
        [  0   1   2   3]
        [  0   0 1/2   0]
        [  0   0   0 1/2]
        sage: N_2.transpose() * Omega(2) * N_2 == 1/2 * Omega(2)
        True

    It has multiplier `1/2` and hence determinant `1/4`. Now consider the matrix
    `N \in GSp_4(\hat{\QQ})` given at `2` by `N_2` and at all other primes by
    the `4 \times 4`-identity matrix. We have `\det(N)\hat{\ZZ} = 1/4\hat{\ZZ}`.

    The following matrix represents `N` and has good precision with respect to
    `1/4`::

        sage: z = Qhat(0, 2)
        sage: e = Qhat(1, 2)
        sage: h = Qhat(1/2, 2)
        sage: M = matrix([[e,z,e,z],[z,e,z,e],[z,z,h,z],[z,z,z,h]]); M
        [  1 mod 2   0 mod 2   1 mod 2   0 mod 2]
        [  0 mod 2   1 mod 2   0 mod 2   1 mod 2]
        [  0 mod 2   0 mod 2 1/2 mod 2   0 mod 2]
        [  0 mod 2   0 mod 2   0 mod 2 1/2 mod 2]
        sage: has_good_precision(M, 1/4)
        True

    Hence we can compute::

        sage: A = factor_GSpQhat(M, 1/4); A
        [  1   0   0   0]
        [  0   1   0   0]
        [  0   0 1/2   0]
        [  0   0   0 1/2]

    Then indeed `A` is symplectic::

        sage: A.transpose() * Omega(2) * A == 1/2 * Omega(2)
        True

    and for every odd prime number `p` we have `A \in GSp_4(\ZZ_p)`. Hence
    to check that `N A^{-1} \in GSp_4(\hat{\ZZ})` is suffices to check that
    `N_2 A^{-1} \in GSp_4(\ZZ_2)`::

        sage: N_2 * ~A
        [1 0 2 4]
        [0 1 4 6]
        [0 0 1 0]
        [0 0 0 1]
        sage: det(N_2*~A)
        1

    Let's now do a bigger example. We construct two rational symplectic matrices
    ::

        sage: I = identity_matrix(3)
        sage: S = matrix([[1/2, 1, 2], [1, 2, 3], [2, 3, 1/5]])
        sage: N_2 = block_matrix([[6*I, S], [0*I, 2*I]], subdivide=False); N_2
        [  6   0   0 1/2   1   2]
        [  0   6   0   1   2   3]
        [  0   0   6   2   3 1/5]
        [  0   0   0   2   0   0]
        [  0   0   0   0   2   0]
        [  0   0   0   0   0   2]
        sage: N_2.transpose() * Omega(3) * N_2 == 12 * Omega(3)
        True


    and ::

        sage: U = matrix([[3, 2, 1/2], [1, -3, -3], [4, -2, 5]])
        sage: X = block_matrix([[I, 0*I], [S, 1/3*I]], subdivide=False)
        sage: Y = block_matrix([[U, 0*I], [0*I, ~U.transpose()]], subdivide=False)
        sage: N_3 = X*Y; N_3
        [      3       2     1/2       0       0       0]
        [      1      -3      -3       0       0       0]
        [      4      -2       5       0       0       0]
        [   21/2      -6    29/4    7/92  17/276  -5/138]
        [     17     -10    19/2  11/276 -13/276  -7/138]
        [   49/5   -27/5      -7   3/184 -19/552  11/276]
        sage: N_3.transpose() * Omega(3) * N_3 == 1/3 * Omega(3)
        True
    
    We let `N \in GSp_6(\hat{\QQ})` at `2` be given by `N_2`, at `3` by `N_3`
    and by the identity matrix at all other primes. Then because we have ::

        sage: det(N_2).valuation(2)
        6
        sage: det(N_3).valuation(3)
        -3

    the determinant of `N` satisfies `det(N)\hat{\ZZ} = 64/27\hat{\ZZ}`.

    We construct a representation `M` of `N` of good precision::

        sage: M = matrix([[Qhat([Qp(2)(N_2[i,j], 7), Qp(3)(N_3[i,j], 1)])
        ....:             for j in range(6)] for i in range(6)])
        sage: M
        [    6 mod 384   128 mod 384   128 mod 384 513/2 mod 384   129 mod 384   258 mod 384]
        [  256 mod 384     6 mod 384     0 mod 384   129 mod 384   258 mod 384     3 mod 384]
        [  256 mod 384   256 mod 384   134 mod 384   258 mod 384     3 mod 384   333 mod 384]
        [    0 mod 384     0 mod 384   128 mod 384     2 mod 384 256/3 mod 384 256/3 mod 384]
        [  128 mod 384   128 mod 384   128 mod 384 640/3 mod 384 646/3 mod 384 128/3 mod 384]
        [  128 mod 384     0 mod 384   128 mod 384     0 mod 384 128/3 mod 384 262/3 mod 384]
        sage: has_good_precision(M, 64/27)
        True

    Now we can use this factorization function::

        sage: A = factor_GSpQhat(M, 64/27); A
        [  2   0   0 1/6 1/3   0]
        [  0   2   0 1/3   0 1/3]
        [  0   0   2   0 1/3 1/3]
        [  0   0   0 2/3   0   0]
        [  0   0   0   0 2/3   0]
        [  0   0   0   0   0 2/3]

    Let's check that `A` does indeed lie in `GSp_6^+(\QQ)`::
        
        sage: A.transpose() * Omega(3) * A == 4/3 * Omega(3)
        True

    From the entries of `A` and the fact that `A` has multiplier `4/3` it
    follows that `A \in GSp_6(\ZZ_p)` for all prime numbers bigger than `3`.
    At the primes `2` and `3` we have ::

        sage: min([entry.valuation(2) for entry in (N_2 * ~A).list()])
        0
        sage: det(N_2*~A).valuation(2)
        0
        sage: min([entry.valuation(3) for entry in (N_3 * ~A).list()])
        0
        sage: det(N_3*~A).valuation(3)
        0

    Therefore `N A^{-1} \in GSp_6(\hat{\ZZ})` as desired.
    """
    if not has_good_precision(M, detM):
        raise ValueError("matrix precision too low!")

    if not M.nrows() % 2 == 0:
        raise ValueError("matrix is not symplectic: it has odd dimension")

    g = M.nrows() // 2
    A_0 = factor_GLQhat(M, detM)
    A_1 = denominator(A_0) * A_0
    E = A_1 * Omega(g) * A_1.transpose()
    _, C = symplectic_basis_over_ZZ(E)
    A = C * A_0
    return A


