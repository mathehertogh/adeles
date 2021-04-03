
from sage.rings.infinity import Infinity
from sage.rings.complex_mpfr import ComplexField
from sage.rings.real_mpfi import RIF
from sage.rings.complex_interval_field import ComplexIntervalField
from sage.categories.homset import Hom
from sage.arith.misc import factor
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.padics.factory import Qp

CC = ComplexField()
CIF = ComplexIntervalField()


def completions(K, p, prec=20):
    """
    Return the completions of the number field ``K`` at the primes above ``p``.

    INPUT:

    - ``K`` -- a number field
    - ``p`` -- a (rational) prime number or ``Infinity``
    - ``prec`` -- precision of the ``p``-adic number field we use (default: 20)

    OUTPUT:

    A list of tuples ``(q, L, phi)`` (one for each prime lying above ``p``)
    where: ``q`` is a prime of ``K`` lying above ``p``, ``L`` is the completion
    of ``K`` at ``q``, and ``phi`` is the embedding of ``K`` in ``L``.
    """
    if p is Infinity:
        return infinite_completions(K)
    result = []
    R = PolynomialRing(Qp(p, prec), 'y')
    f = R(K.defining_polynomial())
    for g, e in factor(f):
        if g.degree() == 1:
            L = Qp(p)
            b = g.roots()[0][0]
        else:
            L = R.quotient_by_principal_ideal(g, names=('b',)); (b,) = L._first_ngens(1)
            assert L.is_field(), "L is not a field"
        phi = Hom(K, L)([b])
        hits = []
        for q in K.primes_above(p):
            if all([phi(r).norm().valuation() > 0 for r in q.gens()]):
                hits.append(q)
        assert len(hits) == 1, "No unique q corresponding to g; hits: {}".format(hits)
        q = hits[0]
        result.append((q, L, phi))
    return result

def completion(K, q, prec=20):
    """
    Return the completion of the number field ``K`` at the prime ``q``.

    INPUT:

    - ``K`` -- a number field
    - ``q`` -- a prime ideal of K
    - ``prec`` -- precision of the ``p``-adic number field we use (default: 20)

    OUTPUT:

    A pair ``(L, phi)`` with ``L`` the completion of ``K`` at ``q``, and ``phi``
    the embedding of ``K`` in ``L``.
    """
    if K is QQ:
        L = Qp(q, prec=prec)
        return L, QQ.embeddings(L)[0]
    p = q.gens_two()[0] # the prime below q (i.e. `p = q \cap \ZZ`)
    all_completions = completions(K, p, prec)
    for r, L, phi in all_completions:
        if r == q:
            return L, phi
    raise ValueError("something went wrong in complete...")

def infinite_completions(K):
    """
    Return the infinite completions of ``K``

    INPUT:

    - ``K`` -- a number field

    OUTPUT:

    A list of tuples (K_oo, phi) with K_oo a completion of K and phi: K --> K_oo
    "the" embedding (the same one as K.places() returns, except with codomain
    an interval field).
    """
    completions = []
    for phi in K.places():
        if phi.codomain() is CC:
            K_oo = CIF
        else:
            K_oo = RIF
        im_gen = K_oo(phi(K.gen()))
        psi = Hom(K, K_oo)([im_gen], check=False)
        completions.append((K_oo, psi))
    return completions


