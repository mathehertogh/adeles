r"""
We implement a method that returns a random element from
the symplectic group over the rationals.

Let `g \in \ZZ`, `g \geq 1`. Then `GSp_{2g}(\QQ)` is generated
by the matrices

`II = (1, S; 0, 1)` with `S` the zero matrix, except with a 1 at the top-left
entry.

`RR_1 = (U_1, 0; 0, (U_1^t)^{-1})` with `U_1`...
"""


from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.constructor import matrix


def _block_matrix(blocks):
    """
    Construct a matrix out of the blocks ``blocks``

    EXAMPLE::

        sage: A = matrix(ZZ, [[1,2],[3,4]])
        sage: I = matrix(ZZ, [[1,0],[0,1]])
        sage: block_matrix([[A, I, A], [I, I, A]])
        [1 2 1 0 1 2]
        [3 4 0 1 3 4]
        [1 0 1 0 1 2]
        [0 1 0 1 3 4]
        sage: B = matrix(QQ, [[1/2], [1/2]])
        sage: R.<x> = ZZ[]
        sage: C = matrix(R, [[x, -x]])
        sage: D = matrix(ZZ, [[0]])
        sage: E = block_matrix([[A, B], [C, D]]); E
        [  1   2 1/2]
        [  3   4 1/2]
        [  x  -x   0]
        sage: E.parent()
        Full MatrixSpace of 3 by 3 dense matrices over Univariate Polynomial Ring in x over Rational Field
    """
    common_element = ZZ(1)
    for row in blocks:
        for block in row:
            common_element *= block[0,0]
    common_base = common_element.parent()

    result = []
    for row in blocks:
        for i in range(row[0].nrows()):
            long_row = []
            for block in row:
                long_row += list(block[i])
            result.append(long_row)

    return matrix(common_base, result)

def _S_0(g):
    S_0 = MatrixSpace(QQ, g).zero().__copy__()
    S_0[0,0] = QQ(1)
    return S_0

def _II_0(g):
    Z = MatrixSpace(QQ, g).zero()
    I = MatrixSpace(QQ, g).one()
    return _block_matrix([[I, _S_0(g)], [Z, I]])

def _U_1(g):
    U_1 = MatrixSpace(QQ, g).zero().__copy__()
    U_1[0,g-1] = QQ(1)
    for i in range(1, g):
        U_1[i, i-1] = QQ(1)
    return U_1

def _RR_1(g):
    Z = MatrixSpace(QQ, g).zero()
    return _block_matrix([[_U_1(g), Z], [Z, _U_1(g).transpose().inverse()]])

def _U_2(g):
    U_2 = MatrixSpace(QQ, g).one().__copy__()
    U_2[0,1] = QQ(1)
    return U_2

def _RR_2(g):
    Z = MatrixSpace(QQ, g).zero()
    return _block_matrix([[_U_2(g), Z], [Z, _U_2(g).transpose().inverse()]])

def _J_1(g):
    J_1 = MatrixSpace(QQ, g).one().__copy__()
    J_1[0,0] = QQ(0)
    return J_1

def _GG_0(g):
    I = MatrixSpace(QQ, g).one()
    J = _J_1(g)
    return _block_matrix([[J, I-J], [J-I, J]])

def _iota(g, nu):
    Z = MatrixSpace(QQ, g).zero()
    I = MatrixSpace(QQ, g).one()
    return _block_matrix([[I, Z], [Z, nu*I]])

def _random_GSp_generator(g):
    from random import choice
    gen = choice([_II_0, _RR_1, _RR_2, _GG_0])(g)
    invert = choice([True, False])
    if invert:
        gen = gen.inverse()
    return gen

def random_GSp_element(g, n=None):
    if n is None:
        n = -1
        while n < 0:
            n = ZZ.random_element()

    elt = MatrixSpace(QQ, 2*g).one()
    for i in range(n):
        elt *= _random_GSp_generator(g)

    nu = QQ(0)
    while nu == 0:
        nu = QQ.random_element()

    elt *= _iota(g, nu)

    return elt