"""
Infinite Completions of Number Fields

This file implements the functionality that `K.places()` claims to provide (for
`K` a number field), but does not provide, in the function
:func:`infinite_completions`.

AUTHORS:

- Mathé Hertogh (2021-07): initial version
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

from sage.categories.homset import Hom
from sage.all import ComplexField
from sage.rings.complex_interval_field import ComplexIntervalField
from sage.rings.real_mpfi import RealIntervalField


def infinite_completions(K, fields_only=False, embeddings_only=False):
    r"""
    Return the infinite completions of the number field ``K``

    INPUT:

    - ``K`` -- a number field
    - ``fields_only`` - boolean (default: ``False``); if ``True``, only return
      the fields, not the embeddings.
    - ``embeddings_only`` - boolean (default: ``False``); if ``True``, only
      return the embeddings, not the fields.

    OUTPUT:

    A list of pairs `(L, \phi)` with `L` equal to ``RIF`` or ``CIF`` and `\phi`
    an embedding `K \to L`. The embeddings returned correspond to
    the infinite primes of `K` and they are returned in the same order as
    ``K.places()``.

    Depending on ``fields_only`` and ``embeddings_only``, only the fields `L` or
    the embeddings `\phi` are returned.
    If they are both set to ``True``, an exception is raised.

    EXAMPLES::

        sage: infinite_completions(QQ)
        [(Real Interval Field with 53 bits of precision,
          Ring morphism:
            From: Rational Field
            To:   Real Interval Field with 53 bits of precision
            Defn: 1 |--> 1)]

    ::

        sage: K.<a> = NumberField(x^3+2)
        sage: infinite_completions(K)
        [(Real Interval Field with 53 bits of precision,
          Ring morphism:
            From: Number Field in a with defining polynomial x^3 + 2
            To:   Real Interval Field with 53 bits of precision
            Defn: a |--> -1.259921049894873?),
         (Complex Interval Field with 53 bits of precision,
          Ring morphism:
            From: Number Field in a with defining polynomial x^3 + 2
            To:   Complex Interval Field with 53 bits of precision
            Defn: a |--> 0.62996052494743671? + 1.0911236359717214?*I)]

    We can obtain only the embeddings as follows::

        sage: K.<sqrt2> = NumberField(x^2-2)
        sage: infinite_completions(K, embeddings_only=True)
        [Ring morphism:
           From: Number Field in sqrt2 with defining polynomial x^2 - 2
           To:   Real Interval Field with 53 bits of precision
           Defn: sqrt2 |--> -1.414213562373095?,
         Ring morphism:
           From: Number Field in sqrt2 with defining polynomial x^2 - 2
           To:   Real Interval Field with 53 bits of precision
           Defn: sqrt2 |--> 1.414213562373095?]

    And we obtain only the fields as follows::

        sage: K.<a> = NumberField(x^8-3*x^5+1)
        sage: infinite_completions(K, fields_only=True)
        [Real Interval Field with 53 bits of precision,
         Real Interval Field with 53 bits of precision,
         Complex Interval Field with 53 bits of precision,
         Complex Interval Field with 53 bits of precision,
         Complex Interval Field with 53 bits of precision]
        sage: K.signature() # the above is consistent with K's signature:
        (2, 3)
    """
    if fields_only and embeddings_only:
        raise ValueError("both fields_only and embeddings_only set to True")

    completions = []
    for psi in K.places():
        if psi.codomain() is ComplexField():
            L = ComplexIntervalField()
        else:
            L = RealIntervalField()
        image_gen = L(psi(K.gen()))
        phi = Hom(K, L)([image_gen], check=False)
        if fields_only:
            completions.append(L)
        elif embeddings_only:
            completions.append(phi)
        else:
            completions.append((L, phi))
    
    return completions

