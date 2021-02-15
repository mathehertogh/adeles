

def completions(K, p):
    """
    Return the completions of the number field ``K`` at the primes above ``p``.

    INPUT:

    - ``K`` -- a number field
    - ``p`` -- a (rational) prime number or ``Infinity``

    OUTPUT:

    A list of tuples ``(q, L, phi)`` (one for each prime lying above ``p``)
    where: ``q`` is a prime of ``K`` lying above ``p``, ``L`` is the completion
    of ``K`` at ``q``, and ``phi`` is the embedding of ``K`` in ``L``.
    """
    if p is Infinity:
        return infinite_completions(K)
    result = []
    R.<y> = Qp(p)[]
    f = R(K.defining_polynomial())
    for g, e in factor(f):
        if g.degree() == 1:
            L = Qp(p)
            b = g.roots()[0][0]
        else:
            L.<b> = R.quotient_by_principal_ideal(g)
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

def completion(K, q):
    """
    Return the completion of the number field ``K`` at the prime ``q``.

    INPUT:

    - ``K`` -- a number field
    - ``q`` -- a prime ideal of K

    OUTPUT:

    A pair ``(L, phi)`` with ``L`` the completion of ``K`` at ``q``, and ``phi``
    the embedding of ``K`` in ``L``.
    """
    p = q.gens_two()[0] # the prime below q (i.e. `p = q \cap \ZZ`)
    all_completions = completions(K, p)
    for r, L, phi in all_completions:
        if r == q:
            return L, phi
    raise ValueError("something went wrong in complete...")

def infinite_completions(K):
    """
    Return the infinite completions of ``K`` as a list of RR's/CC's & embeddings

    INPUT:

    - ``K`` -- a number field

    OUTPUT:

    ...
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

