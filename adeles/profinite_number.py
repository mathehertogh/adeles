r"""
Profinite Numbers of Number Fields

Let `K` be a number field with ring of integers `O`. Let `\hat{O}` be the ring
of profinite `K`-integers (cf. :mod:`profinite_integer`). The ring of profinite
`K`-numbers `\hat{K}` is the ring of fractions of `\hat{O}` with respect to
`\ZZ \setminus \{0\}`. So a profinite `K`-integer is a fraction `a/b` with
`a \in \hat{O}` and `b \in \ZZ \setminus \{0\}`. With
`\{n\hat{\ZZ} \mid n \in \ZZ_{>0}\}` as a basis of open neighborhoods of zero,
`\hat{K}` is a topological ring. It is also a commutative `K`-algebra.

In SageMath profinite numbers are implemented in the class
:class:`ProfiniteNumber` as a pair `a = (n, d)` where `n` is a
:class:`profinite_integer.ProfiniteInteger`, called the *numerator*, and `d` is
a positive integer, called the *denominator*.
We define the *represented subset* of `a` to be the represented subset of `n`
divided by `d`, so an `\alpha \in \hat{K}` is represented by `a` if and only if
`d \cdot \alpha` is represented by `n`.

We define the *value* of `a` to be the value of `n` divided by `d`, which is a
`K`-element. Likewise we define the *modulus* of `a` to be the modulus of `n`
divided by `d`, which is a fractional `O`-ideal in `K`. One can check that the
represented subset of `a` equals `x + m \hat{O}` if `x` and `m` denote the value
and modulus of `a` respectively.

We can approximate any `\alpha \in \hat{K}` arbitrarily closely by a
``ProfiniteNumber`` in the following sense: for any open neighborhood `U` of
`\alpha` in `\hat{K}`, there exists a ``ProfiniteNumber`` representing `\alpha`
whose represented subset is contained in `U`.

In the rest of this file we will refer to instances of ``ProfiniteNumber`` as
profinite numbers. We will refer to exact mathematical profinite `K`-numbers as
elements of `\hat{K}`.

Profinite Numbers over `\QQ`
----------------------------

In SageMath the ring `\hat{\QQ}` is implemented as ``Qhat``::

    sage: Qhat
    Ring of Profinite Numbers of Rational Field
    sage: Qhat is ProfiniteNumbers(QQ)
    True

We can create profinite `\QQ`-integers by giving a numerator, a profinite
`\QQ`-integer, and a denominator, a non-zero integer::

    sage: a = Qhat(Zhat(3, 10), 7)
    sage: a.numerator()
    3 mod 10
    sage: a.denominator()
    7

Such a profinite integer is printed using its value and modulus::

    sage: a.value()
    3/7
    sage: a.modulus()
    10/7
    sage: a
    3/7 mod 10/7

Let's see which elements of `\QQ` of height at most `30` are represented by
``a``::

    sage: X = sorted(set([m/n for m in srange(-30, 31) for n in srange(1, 31)]))
    sage: [x for x in X if a.represents(x)]
    [-21, -11, -27/7, -17/7, -1, 3/7, 13/7, 23/7, 9, 19, 29]

The value and modulus of a profinite number uniquely determine its numerator and
denominator, because we always keep the denominator as small as possible::

    sage: b = Qhat(Zhat(7, 35), 14)
    sage: b.numerator() # this is not 7 mod 35
    1 mod 5
    sage: b.denominator() # this is not 14
    2

We can do the reduction above, since ``(7 mod 35)/14`` and ``(1 mod 5)/2`` have
the same represented subset.

We can also construct profinite numbers by specifying their value and modulus::

    sage: c = Qhat(2/3, 30); c
    2/3 mod 30
    sage: c.numerator(), c.denominator()
    (2 mod 90, 3)

We can perform arithmetic on profinite numbers::

    sage: a + b
    3/14 mod 5/14
    sage: c - b
    1/6 mod 5/2
    sage: a * c
    2/7 mod 10/21
    sage: a / 7
    3/49 mod 10/49

See :meth:`_add_`, :meth:`_sub_`, :meth:`_mul_`, :meth:`_div_` for details.

Recall the natural isomorphism `\hat{\ZZ} \to \prod_p \ZZ_p` (with `p` ranging
over all prime numbers and `\ZZ_p` denoting the ring of `p`-adic integers) induced by the
Chinese Remainder Theorem.
It can be naturally extended to an isomorphism `\hat{\QQ} \to \prod_p' \QQ_p`
with `\QQ_p` denoting the field of `p`-adic numbers and the restricted product
taken with respect to the open subrings `\ZZ_p`.
The ring on the right hand side is also known as the *finite adèle ring* of
`\QQ`. We can take this point of view in SageMath as well::

    sage: d_2 = Qp(2)(3/2, 3); d_2
    2^-1 + 1 + O(2^3)
    sage: d_5 = Qp(5)(10, 2); d_5
    2*5 + O(5^2)
    sage: d_7 = Qp(7)(2/49, -1); d_7
    2*7^-2 + O(7^-1)
    sage: d = Qhat([d_2, d_5, d_7]); d
    2755/98 mod 200/7
    sage: for p in prime_range(10):
    ....:     print(d[p])
    ....:
    2^-1 + 1 + O(2^3)
    O(3^0)
    2*5 + O(5^2)
    2*7^-2 + O(7^-1)

.. NOTE::

    As `p`-adic numbers are currently only fully implemented for rational `p`,
    this functionality is only available for profinite `\QQ`-numbers, not for
    profinite `K`-numbers for non-trivial number fields `K`.

Profinite Integers over `K`
---------------------------

We also have profinite integers over non-trivial number fields `K`. Let's have
a look at an example::

    sage: K.<a> = NumberField(x^3-2)
    sage: Khat = ProfiniteNumbers(K); Khat
    Ring of Profinite Numbers of Number Field in a with defining polynomial x^3 - 2
    sage: b = Khat.an_element(); b
    8/5*a^2 mod (7/5*a^2 - 1/5)

Over general number fields moduli are fractional ideals::

    sage: b.value()
    8/5*a^2
    sage: b.modulus()
    Fractional ideal (7/5*a^2 - 1/5)

We can create profinite numbers as before using numerator/denominator pairs or
value/modulus pairs::

    sage: Ohat = ProfiniteIntegers(K)
    sage: c = Khat(Ohat(a, 10), 2); c # numerator/denominator
    1/2*a mod (5)
    sage: d = Khat(a, 10/3); d # value/modulus
    a mod (10/3)

And we can then perform arithmetic with our profinite numbers::

    sage: c + d
    -1/6*a mod (5/3)
    sage: c - d
    -1/2*a mod (5/3)
    sage: c * d
    1/2*a^2 mod (5/3*a)
    sage: c / a
    1/2 mod (5/2*a^2)

Above, the `K`-element ``a`` is coerced into ``Khat`` as ``a mod 0``, which is a
profinite `K`-number representing only ``a``. The division is then performed
inside ``Khat``.

Profinite `K`-integers coerce to profinite `K`-numbers with denominator `1`::

    sage: Khat(Ohat(0, a))
    0 mod (a)

.. SEEALSO::

    :mod:`profinite_integer`,
    :mod:`adele`

REFERENCES:

[Her2021] Mathé Hertogh, Computing with adèles and idèles, master's thesis,
Leiden University, 2021.

This implementation of profinite numbers is based on [Her2021]. An extensive
exposition of properties, design choices and two applications can be found
there.

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

from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.ring import CommutativeAlgebra
from sage.structure.element import CommutativeAlgebraElement
from sage.arith.misc import factor
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.sets.primes import Primes
from sage.arith.functions import lcm
from sage.modules.free_module_element import vector

from .profinite_integer import ProfiniteIntegers, ProfiniteCompletionFunctor


class ProfiniteNumber(CommutativeAlgebraElement):
    """
    Profinite Number of a Number Field

    .. automethod:: __init__

    .. automethod:: __getitem__

    .. automethod:: _add_

    .. automethod:: _sub_

    .. automethod:: _mul_

    .. automethod:: _div_

    .. automethod:: _richcmp_
    """

    def __init__(self, parent, numerator, denominator=1):
        r"""
        Construct the profinite number ``numerator / denominator``
        
        INPUT:

        - ``parent`` -- A ring of profinite numbers over some number field `K`
        - ``numerator`` -- a profinite `K`-integer
        - ``denominator`` -- a non-zero integral element of `K` (default: 1)

        EXAMPLES::

            sage: Zhat = ProfiniteIntegers(ZZ)
            sage: Qhat = ProfiniteNumbers(QQ)
            sage: ProfiniteNumber(Qhat, Zhat(2, 15), 8)
            1/4 mod 15/8
            sage: ProfiniteNumber(Qhat, -97)
            -97
            sage: ProfiniteNumber(Qhat, Zhat(2, 16), 8)
            1/4 mod 2

        ::

            sage: K.<a> = NumberField(x^4+7)
            sage: Ohat = ProfiniteIntegers(K)
            sage: Khat = ProfiniteNumbers(K)
            sage: ProfiniteNumber(Khat, Ohat(-a^3, a^2+1), 5)
            -1/5*a^3 mod (1/5*a^2 + 1/5)

        TESTS::

            sage: ProfiniteNumber(Qhat, a, 3)
            Traceback (most recent call last):
            ...
            TypeError: numerator should be a profinite integer over the base number field
            sage: ProfiniteNumber(Khat, Ohat(a), 0)
            Traceback (most recent call last):
            ...
            TypeError: denominator should be a non-zero integral element of the base number field
        """
        CommutativeAlgebraElement.__init__(self, parent)
        O = parent.base().maximal_order()
        Ohat = ProfiniteIntegers(O)
        if numerator not in Ohat:
            raise TypeError("numerator should be a profinite integer over the base number field")
        if denominator not in O or denominator == 0:
            raise TypeError("denominator should be a non-zero integral element of the base number field")
        self._numerator = Ohat(numerator)
        self._denominator = O(denominator)
        self._reduce()

    def __getitem__(self, p):
        r"""
        Return the projection of this profinite number to the field of
        ``p``-adic numbers

        INPUT:

        - ``p`` -- a prime number

        Only implemented for profinite numbers over `\QQ`, as no general
        implementation of completions at finite places of a number field
        exists in SageMath at the time of writing.

        This method implements the projection `\hat{\QQ} \to \QQ_p` extending
        the natural projection `\hat{\ZZ} \to \ZZ_p`.

        EXAMPLES::

            sage: a = Qhat(2/35, 2^6 * 3^5 * 5)
            sage: a[2]
            2 + 2^2 + 2^4 + O(2^6)
            sage: a[3]
            1 + 2*3 + O(3^5)
            sage: a[5]
            5^-1 + 2 + O(5)
            sage: a[7]
            6*7^-1 + O(7^0)
            sage: a[11]
            O(11^0)
        """
        if self.parent().base() is not QQ:
            raise NotImplementedError("projection to `p`-adics only implemented over rationals")
        if p not in Primes():
            raise ValueError("p should be a prime number")
        return self.numerator()[p] / self.denominator()

    def _repr_(self):
        """
        Returns a string representation of ``self``

        EXAMPLES::

            sage: Qhat(1/2, 97/5)
            1/2 mod 97/5
            sage: Qhat(1/3)
            1/3

        ::

            sage: K.<a> = NumberField(x^2+5)
            sage: Khat = ProfiniteNumbers(K)
            sage: Khat(a/2, K.ideal(30, 10*a+10))
            1/2*a mod (30, 10*a + 10)
        """
        if self.modulus() == 0:
            return repr(self.value())

        if self.parent().number_field() is QQ:
            modulus = self.modulus()
        else:
            modulus = "(" + ", ".join([str(g) for g in self.modulus().gens()]) + ")"
        
        return "{} mod {}".format(self.value(), modulus)

    def _add_(self, other):
        r"""
        Return the sum of this profinite number and ``other``

        The sum of two profinite numbers `a` and `b` is defined to be the
        profinite number `c` with smallest represented subset (with respect to
        inclusion) containing the sum of the represented subsets of `a` and `b`,
        i.e. containig `\alpha + \beta` for all `\alpha` represented by `a` and
        `\beta` represented by `b`.

        EXAMPLES::
        
            sage: Qhat(Zhat(1, 10), 2) + Qhat(Zhat(2, 15), 3)
            7/6 mod 5
            sage: Qhat(2, 5) + Qhat(1, 10/3)
            4/3 mod 5/3

        ::

            sage: K.<a> = NumberField(x^3-2)
            sage: Ohat = ProfiniteIntegers(K)
            sage: Khat = ProfiniteNumbers(K)
            sage: Khat(Ohat(a, 10), 3) + Khat(Ohat(a^2, 15), 3)
            1/3*a^2 + 1/3*a mod (5/3)
            sage: Khat(Ohat(a+1, 20*a+3), 7) + 5
            -3749/7*a^2 mod (-20/7*a - 3/7)
        """
        numerator = self.numerator() * other.denominator() + other.numerator() * self.denominator()
        denominator = self.denominator() * other.denominator()
        return self.__class__(self.parent(), numerator, denominator)

    def _sub_(self, other):
        r"""
        Return the difference of this profinite number and ``other``

        The difference of two profinite numbers `a` and `b` is defined to be the
        profinite number `c` with smallest represented subset (with respect to
        inclusion) containing the difference of the represented subsets of `a`
        and `b`, i.e. containig `\alpha - \beta` for all `\alpha` represented by
        `a` and `\beta` represented by `b`.

        EXAMPLES::
        
            sage: Qhat(Zhat(4, 12), 3) - 2/3
            2/3 mod 4
            sage: 1 - Qhat(Zhat(3, 6), 2)
            5/2 mod 3
            sage: Qhat(Zhat(4, 12), 3) - Qhat(Zhat(3, 6), 2)
            5/6 mod 1

        ::

            sage: K.<a> = NumberField(x^3-2)
            sage: Khat = ProfiniteNumbers(K)
            sage: Khat(a, 100) - Khat(a, 20)
            0 mod (20)
            sage: Khat(a, 100) - Khat(2, 5*a)
            a - 2 mod (5*a)
        """
        numerator = self.numerator() * other.denominator() - other.numerator() * self.denominator()
        denominator = self.denominator() * other.denominator()
        return self.__class__(self.parent(), numerator, denominator)

    def _mul_(self, other):
        """
        Return the product of this profinite number and ``other``

        The product of two profinite numbers `a` and `b` is defined to be the
        profinite number `c` with smallest represented subset (with respect to
        inclusion) containing the product of the represented subsets of `a`
        and `b`, i.e. containig `\alpha \beta` for all `\alpha` represented by
        `a` and `\beta` represented by `b`.

        EXAMPLES::
        
            sage: 1/2 * Qhat(Zhat(5, 15), 2)
            5/4 mod 15/4
            sage: Qhat(1/3, 30) * Qhat(2, 15)
            2/3 mod 5

        ::

            sage: K.<zeta5> = CyclotomicField(5)
            sage: Khat = ProfiniteNumbers(K)
            sage: zeta5^3 * Khat(2*zeta5^2, 100/3)
            2 mod (100/3*zeta5^3)
            sage: Khat(zeta5, 9) * Khat(1/2, 9*zeta5)
            1/2*zeta5 mod (9/2)
        """
        numerator = self.numerator() * other.numerator()
        denominator = self.denominator() * other.denominator()
        return self.__class__(self.parent(), numerator, denominator)

    def _div_(self, other):
        r"""
        Return the quotient of this profinite number by ``other``

        The profinite number ``other`` must have zero-modulus and non-zero
        value. In other words, it must be the profinite number correspoding to
        a non-zero element of the base maximal order `O`.

        The quotient of a profinite number `a` by an element `b \in O` is the
        unique profinite number whose represented subset equals `R(a)/b`, where
        `R(a)` denotes the represented subset of `a`.

        EXAMPLES::
        
            sage: Qhat(3/2, 5) / 2
            3/4 mod 5/2
            sage: Qhat(4, 8) / Qhat(2/3, 0)
            6 mod 12

        ::

            sage: K.<zeta3> = CyclotomicField(3)
            sage: Khat = ProfiniteNumbers(K)
            sage: Khat(1, 100) / zeta3^2
            zeta3 mod (100*zeta3)

        TESTS::

            sage: Khat(1, 100) / Khat(1, 5)
            Traceback (most recent call last):
            ...
            ValueError: division by profinite number with non-zero modulus
            sage: Khat(1, 100) / Khat(0, 0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: profinite number division by zero
        """
        if other.modulus() != 0:
            raise ValueError("division by profinite number with non-zero modulus")
        if other.value() == 0:
            raise ZeroDivisionError("profinite number division by zero")
        return self * (1/other.value())

    def _richcmp_(self, other, op):
        r"""
        Compare ``self`` and ``other`` based on the relation ``op``

        We only implement equality and non-equality.

        We declare ``a/b`` equal to ``c/d`` if and only if ``a*d == b*c``.

        EXAMPLES::

            sage: Qhat(Zhat(1, 6), 2) == Qhat(Zhat(2, 12), 4)
            True
            sage: Qhat(Zhat(1, 6), 2) == Qhat(Zhat(2, 12), 3)
            False
            sage: Qhat(Zhat(1, 6), 2) != Qhat(Zhat(2, 12), 3)
            True

        ::

            sage: K.<a> = NumberField(x^2+5)
            sage: Khat = ProfiniteNumbers(K)
            sage: Khat(2, 11/3) == Khat(2, 11)
            True
            sage: Khat(1, 11/3) == Khat(2, 11)
            False

        .. WARNING::

            Equality is *not* transitive::

                sage: a, b, c = Khat(2/7, 1), Khat(1/7, 1/7), Khat(1/7, 1/5)
                sage: a == b
                True
                sage: b == c
                True
                sage: a == c
                False
        """
        from sage.structure.richcmp import op_EQ, op_NE
        if op == op_EQ:
            return self.numerator() * other.denominator() == self.denominator() * other.numerator()
        if op == op_NE:
            return not self._richcmp_(other, op_EQ)
        raise NotImplementedError("only equality and inequality are implemented")

    def _reduce(self):
        """
        Adjust the numerator and denominator of this profinite number without
        changing its value or modulus, making the denominator the smallest
        positive integer possible.

        EXAMPLES::

            sage: Qhat = ProfiniteNumbers()
            sage: a = Qhat(Zhat(1, 100), 60)
            sage: a._numerator._value = 10
            sage: a.numerator(), a.denominator()
            (10 mod 100, 60)
            sage: a._reduce()
            sage: a.numerator(), a.denominator()
            (1 mod 10, 6)
            sage: a._numerator._value = 14
            sage: a.numerator(), a.denominator()
            (14 mod 10, 6)
            sage: a._reduce()
            sage: a.numerator(), a.denominator()
            (2 mod 5, 3)

        An example with non-principal ideals.
        This method is called when creating a profinite number::

            sage: K.<a> = NumberField(x^2+5)
            sage: Ohat = ProfiniteIntegers(K)
            sage: Khat = ProfiniteNumbers(K)
            sage: Khat(Ohat(2*a, K.ideal(4, 2*a+2)), 2)
            a mod (2, a + 1)
        """
        K = self.parent().base()
        Ohat = ProfiniteIntegers(K)
        num = self.numerator()

        if K is QQ:
            from sage.arith.misc import gcd
            g = gcd([num.value(), num.modulus(), self.denominator()])
            denominator = abs(self.denominator() // g)
        else:
            g = K.ideal(num.value()) + num.modulus() + K.ideal(self.denominator())
            # We take the smallest positive integer inside self.denominator()/g.
            denominator = (self.denominator() / g).gens_two()[0]
        
        value = denominator * num.value() / self.denominator()
        # Avoid division of the zero-ideal, which is not defined in SageMath.
        if num.modulus() == 0:
            modulus = num.modulus()
        else:
            modulus = denominator * num.modulus() / self.denominator()
        
        self._numerator = Ohat(value, modulus)
        self._denominator = ZZ(denominator)

    def numerator(self):
        """
        Return the numerator of this profinite number

        EXAMPLES::

            sage: Qhat(1/3, 7/5).numerator()
            5 mod 21
            sage: Qhat(Zhat(2, 10), 4).numerator()
            1 mod 5

        ::

            sage: K.<a> = NumberField(x^2+10)
            sage: Khat = ProfiniteNumbers(K)
            sage: Khat(3/2, 10*a).numerator()
            3 mod (20*a)
            sage: Khat(a/2, 10/3).numerator()
            3*a mod (20)
        """
        return self._numerator

    def denominator(self):
        """
        Return the denominator of this profinite number

        EXAMPLES::

            sage: Qhat(1/3, 7/5).denominator()
            15
            sage: Qhat(Zhat(2, 10), 4).denominator()
            2

        ::
        
            sage: K.<a> = NumberField(x^2+10)
            sage: Khat = ProfiniteNumbers(K)
            sage: Khat(a/2, 10/3).denominator()
            6
        """
        return self._denominator

    def value(self):
        """
        Return the value of this profinite number

        EXAMPLES::

            sage: Qhat(7/9, 12).value()
            7/9
            sage: Qhat(Zhat(5, 10), 2).value()
            5/2

        ::

            sage: K.<a> = NumberField(x^2+15)
            sage: Khat = ProfiniteNumbers(K)
            sage: Khat(8*a, 30/7).value()
            -4/7*a
            sage: 8*a - -4/7*a in K.ideal(30/7)
            True
        """
        return self.numerator().value() / self.denominator()

    def modulus(self):
        r"""
        Return the modulus of this profinite integer
        
        Denote our base number field by `K`.
        If `K` is `\QQ`, then the modulus is a non-negative rational number.
        Otherwise the modulus is a fractional `O`-ideal in `K`, with `O` the
        maximal order of `K`.

        EXAMPLES::

            sage: Qhat(1, 79/90).modulus()
            79/90
            sage: Qhat(Zhat(4, 10), 7).modulus()
            10/7

        ::

            sage: K.<a> = NumberField(x^2+15)
            sage: Khat = ProfiniteNumbers(K)
            sage: Khat(3, 4*a).modulus()
            Fractional ideal (4*a)
            sage: Khat(-1, K.prime_above(2)).modulus()
            Fractional ideal (2, 1/2*a - 1/2)
        """
        # Avoid division of the zero ideal, which isn't implemented in SageMath.
        if self.numerator().modulus() == 0:
            return self.numerator().modulus()
        return self.numerator().modulus() / self.denominator()

    def represents(self, element):
        r"""
        Return wether or not this profinite number represents ``element``

        For `x` the value of this profinite number and `m` its modulus, this
        profinite number has represented subset `x + m \hat{O}`, for `\hat{O}`
        the ring of profinite integers over our base number field.
        So this profinite integer represents `\alpha` if and only if
        `\alpha \in x + m \hat{O}`.

        INPUT:

        - ``element`` -- an element of our base number field

        EXAMPLES::

            sage: X = [m/n for m in srange(-30, 30) for n in srange(1, 30)]
            sage: X = sorted(set(X))
            sage: a = Qhat(1/2, 5)
            sage: [x for x in X if a.represents(x)]
            [-29/2, -19/2, -9/2, 1/2, 11/2, 21/2]
            sage: b = Qhat(1/3, 7/2)
            sage: [x for x in X if b.represents(x)]
            [-20/3, -19/6, 1/3, 23/6, 22/3]

        ::

            sage: K.<a> = NumberField(x^2-7)
            sage: Khat = ProfiniteNumbers(K)
            sage: b = Khat(a/2, 10/2)
            sage: Y = [m/n for m in srange(-6, 6) for n in srange(1, 6)]
            sage: Y = sorted(set(Y))
            sage: [y*a+z for y in Y for z in Y if b.represents(y*a+z)] # long time
            [1/2*a - 5, 1/2*a, 1/2*a + 5]
        """
        if self.parent().number_field() is QQ:
            x = (self.value() - element) * self.modulus().denominator()
            if x not in ZZ:
                return False
            return ZZ(x) % self.modulus().numerator() == 0
        return self.value() - element in self.modulus()

    def is_integral(self):
        """
        Return whether or not this profinite number is integral

        A profinite number is integral if its denominator equals one.

        EXAMPLES::

            sage: Qhat(7, 11/3).is_integral()
            False
            sage: Qhat(Zhat(4, 10), 2).is_integral()
            True

        ::

            sage: K.<a> = NumberField(x^2+3)
            sage: Khat = ProfiniteNumbers(K)
            sage: Khat((a+1)/2, 10).is_integral()
            True
            sage: Khat((a+1)/4, 10).is_integral()
            False
        """
        return self.denominator() == 1

    def to_profinite_rational_vector(self):
        r"""
        Convert this profinite number to a vector of profinite `\QQ`-numbers
    
        Let `K` be our base number field, `\alpha` its generator over `\QQ` and
        `n` its degree over `\QQ`. Then every `x \in \hat{K}` can be written
        uniquely as `x = \sum_{i=0}^{n-1} x_i \alpha^i` with
        `x_i \in \hat{\QQ}`. This induces a map `\phi: \hat{K} \to \hat{\QQ}^n`.

        This method implements `\phi`.
        
        OUTPUT:

        A vector `(x_0, ..., x_{n-1})` of profinite `\QQ`-integers such that the
        following holds. Write `\phi_i: \hat{K} \to \hat{\QQ}` for `\phi`
        composed with the projection to the `i`-th `\hat{\QQ}`. Then for any
        `\alpha` that ``self`` represents, `\phi_i(x)` is represented by `x_i`.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2-3)
            sage: Khat = ProfiniteNumbers(K)
            sage: Khat(1/9 + 7*a).to_profinite_rational_vector()
            (1/9, 7)
            sage: b = Khat(a, 6)
            sage: b0, b1 = b.to_profinite_rational_vector(); b0, b1
            (0 mod 6, 1 mod 2)
            sage: b.represents(a)
            True
            sage: b0.represents(a.vector()[0]) and b1.represents(a.vector()[1])
            True
            sage: b.represents(a+6)
            True
            sage: b0.represents((a+6).vector()[0]) and b1.represents((a+6).vector()[1])
            True

        ::

            sage: K.<a> = NumberField(x^3-7)
            sage: Khat = ProfiniteNumbers(K)
            sage: b = Khat(a^2-1/3, 3*a); b
            a^2 - 1/3 mod (3*a)
            sage: c = b.to_profinite_rational_vector(); c
            (8/3 mod 3, 0 mod 3, 1/7 mod 3/7)
            sage: b.represents(a^2-1/3)
            True
            sage: all([c[i].represents((a^2-1/3).vector()[i]) for i in range(3)])
            True
            sage: b.represents(-2*a^2-1/3)
            True
            sage: all([c[i].represents((-2*a^2-1/3).vector()[i]) for i in range(3)])
            True

        TESTS::

            sage: Qhat(1/2, 10).to_profinite_rational_vector()
            (1/2 mod 10,)
        """
        K = self.parent().base()

        if K is QQ:
            return (self,)

        if self.modulus() == 0:
            return vector(Qhat, self.value().vector())
        
        # We compute the set of (rational) prime numbers above which a prime of
        # K exists at which our modulus has non-zero valuation.
        primes = set([ZZ(P.gens_two()[0]) for P in self.modulus().prime_factors()])

        I = self.modulus()
        result = []
        for i in range(K.absolute_degree()):
            modulus = 1
            for p in primes:
                e = min([I.valuation(P) // K(p).valuation(P) for P in K.primes_above(p)])
                modulus *= p**e
            result.append(Qhat(self.value().vector()[i], modulus))
            I /= K.gen()

        return vector(result)


class ProfiniteNumbers(UniqueRepresentation, CommutativeAlgebra):
    """
    Ring of Profinite Number over a Number Field

    .. automethod:: _element_constructor_
    """

    Element = ProfiniteNumber

    def __classcall__(cls, K=QQ):
        """
        Construct the ring of profinite numbers over the number field `K`

        INPUT:

        - ``K`` -- a number field (default: ``QQ``)

        EXAMPLES::

            sage: sage: K.<a> = NumberField(x^2+3)
            sage: ProfiniteNumbers(K)
            Ring of Profinite Numbers of Number Field in a with defining polynomial x^2 + 3

        This method ensures the following::

            sage: ProfiniteNumbers() is ProfiniteNumbers(QQ)
            True
        
        TESTS::

            sage: ProfiniteNumbers("not a number field")
            Traceback (most recent call last):
            ...
            TypeError: K should be a number field
        """
        from sage.rings.number_field.number_field import is_NumberField
        if not is_NumberField(K):
            raise TypeError("K should be a number field")
        return super(ProfiniteNumbers, cls).__classcall__(cls, K)

    def __init__(self, K):
        r"""
        Construct the ring of profinite numbers over the number field ``K``

        INPUT:

        - ``K`` -- a number field

        EXAMPLES::

            sage: ProfiniteNumbers()
            Ring of Profinite Numbers of Rational Field
            sage: K.<a> = NumberField(x^5-3*x+1)
            sage: ProfiniteNumbers(K)
            Ring of Profinite Numbers of Number Field in a with defining polynomial x^5 - 3*x + 1
        
        TESTS:

        Pickling works::

            sage: import __main__
            sage: __main__.ProfiniteNumbers = ProfiniteNumbers
            sage: loads(dumps(ProfiniteNumbers(QQ))) is ProfiniteNumbers(QQ)
            True
            sage: loads(dumps(ProfiniteNumbers(K))) is ProfiniteNumbers(K)
            True
        """
        # Note that the input K is checked by __classcall__().
        CommutativeAlgebra.__init__(self, K)

    def _repr_(self):
        """
        Return a string representation of this ring

        EXAMPLES::

            sage: K = NumberField(x^3+x+1, 'a')
            sage: ProfiniteNumbers(K)
            Ring of Profinite Numbers of Number Field in a with defining polynomial x^3 + x + 1
        """
        return "Ring of Profinite Numbers of {}".format(self.base())

    def _latex_(self):
        r"""
        Return latex-formatted string representation of this ring

        EXAMPLES::

            sage: K.<a> = NumberField(x^2+14)
            sage: Khat = ProfiniteNumbers(K)
            sage: latex(Khat)
             \widehat{ \Bold{Q}[a]/(a^{2} + 14) }
        """
        from sage.misc.latex import latex
        return r" \widehat{" + latex(self.base()) + "} "

    def _element_constructor_(self, x, y=None):
        r"""
        Construct a profinite number

        INPUT:

        Either a value/modulus pair:

        - ``x`` -- an element of our base number field; the value
        - ``y`` -- a fractional ideal of our base number field (default: the
          zero ideal); the modulus

        or a numerator/denominator pair:

        - ``x`` -- a profinite integer over our base number field; the numerator
        - ``y`` -- a non-zero integral element of our base number field
          (default: ``1``); the denominator

        If the base number field is `\QQ`, then we also accept a list of
        `p`-adic numbers, for distinct prime numbers `p`.
        See :meth:`_from_padic_numbers` for details.

        EXAMPLES:

        First we create some profinite numbers by value/modulus pairs::

            sage: Qhat(1/3, 10)
            1/3 mod 10
            sage: Qhat(1000, 7/2)
            5/2 mod 7/2
            sage: Qhat(9/7)
            9/7

        ::

            sage: K.<a> = NumberField(x^3+x+1)
            sage: Khat = ProfiniteNumbers(K)
            sage: Khat(a^2, 7/3*(a+1))
            a^2 mod (7/3*a + 7/3)
            sage: Khat(a)
            a
        
        Now we construct some profinite numbers by numerator/denominator pairs::

            sage: Qhat(Zhat(1, 10), 2)
            1/2 mod 5
            sage: Qhat(Zhat(9, 17))
            9 mod 17

        ::

            sage: Khat(Ohat(a, 6), 2)
            1/2*a mod (3)
            sage: Khat(Ohat(2, 4*a))
            2 mod (4*a)

        For `L/K` a field extension, profinite `K`-numbers can be converted to
        profinite `L`-numbers::

            sage: R.<t> = K[]
            sage: L.<b> = K.extension(t^2-a)
            sage: Lhat = ProfiniteNumbers(L)
            sage: Lhat(Khat(1/2*a, 10))
            1/2*a mod (10)

        Lastly, we can also give a list of `p`-adic numbers to ``Qhat``::

            sage: a_2 = Qp(2)(5/8, 1); a_2
            2^-3 + 2^-1 + O(2)
            sage: a_3 = Qp(3)(2, 2); a_3
            2 + O(3^2)
            sage: a_5 = Qp(5)(7/125, -1); a_5
            2*5^-3 + 5^-2 + O(5^-1)
            sage: a = Qhat([a_2, a_3, a_5]); a
            2081/1000 mod 18/5
            sage: a[2], a[3], a[5], a[7]
            (2^-3 + 2^-1 + O(2), 2 + O(3^2), 2*5^-3 + 5^-2 + O(5^-1), O(7^0))

        TESTS::

            sage: Qhat('bla')
            Traceback (most recent call last):
            ...
            TypeError: Can't construct profinite number from ('bla', None)

        .. automethod:: _from_padic_numbers
        """
        K = self.base()
        O = K.maximal_order()
        Ohat = ProfiniteIntegers(K)

        if y is None:
            P = x.parent() if hasattr(x, "parent") else None

            # Check if x is a profinite K-number for some subfield of K:
            if isinstance(P, ProfiniteNumbers) and K.has_coerce_map_from(P.number_field()):
                return self.element_class(self, x.numerator(), x.denominator())

        if x in K:
            value = K(x)

            if y is None:
                modulus = QQ(0) if K is QQ else K.ideal(0)
            elif y in K.ideal_monoid():
                modulus = QQ(y) if K is QQ else K.ideal(y)
            else:
                raise TypeError("second parameter should be a fractional ideal specifying the modulus")

            modulus_den = modulus.denominator() if K is QQ else modulus.integral_split()[1]
            den = lcm(value.denominator(), modulus_den)
            return self.element_class(self, Ohat(den*value, den*modulus), den)

        if x in Ohat:
            if y is None:
                y = ZZ(1)
            return self.element_class(self, x, y)

        if K is QQ:
            from sage.rings.padics.generic_nodes import is_pAdicRing, is_pAdicField
            def is_pAdic(a):
                """Utility function to check if ``a`` is a p-adic number"""
                if not hasattr(a, "parent"): return False
                return is_pAdicRing(a.parent()) or is_pAdicField(a.parent())
            if isinstance(x, list) and all([is_pAdic(a) for a in x]):
                return self._from_padic_numbers(x)

        raise TypeError("can't construct profinite number from {}".format((x,y)))

    def _from_padic_numbers(self, padics):
        r"""
        Construct a profinite `\QQ`-number from the list of `p`-adic numbers
        ``padics``
        
        Only implemented over `\QQ`, as no general implementation of `p`-adic
        numbers currently exists for general number fields.

        INPUT:

        - ``padics`` -- a list of `p`-adic numbers, for distinct prime numbers
          `p`

        OUTPUT:

        The unique profinite `\QQ` number ``a`` such that

        - for each prime number ``p`` for which a ``p``-adic number ``a_p``
          exists in ``padics``, we have that ``a[p]`` precisely equals ``a_p``;
        - for all other prime numbers ``p``, we have that ``a[p]`` precisely
          equals ``O(p^0)``.

        EXAMPLES::

            sage: a_2 = Qp(2)(3/4, 2); a_2
            2^-2 + 2^-1 + O(2^2)
            sage: a_3 = Qp(3)(1, 1); a_3
            1 + O(3)
            sage: a = Qhat._from_padic_numbers([a_2, a_3]); a
            19/4 mod 12
            sage: a[2]
            2^-2 + 2^-1 + O(2^2)
            sage: a[3]
            1 + O(3)
            sage: a[5]
            O(5^0)

        ::

            sage: b_2 = Qp(2)(1/6, 1)
            sage: b_3 = Qp(3)(3/2, 2)
            sage: b = Qhat._from_padic_numbers([b_3, b_2])
            sage: b[2] == b_2 and b[3] == b_3
            True

        TESTS::

            sage: c_2 = Qp(2)(2^-25)
            sage: c_5 = Qp(5)(0, -3)
            sage: c_7 = Qp(7)(2*7^2, 3)
            sage: c_41 = Qp(41)(11/41, 8)
            sage: c = Qhat._from_padic_numbers([c_41, c_2, c_7, c_5])
            sage: c[2] == c_2 and c[5] == c_5 and c[7] == c_7 and c[41] == c_41
            True

        ::

            sage: Qhat._from_padic_numbers([c_2, c_5, c_7, c_2])
            Traceback (most recent call last):
            ...
            ValueError: multiple 2-adic integers in profinite integer initialization
        """
        from sage.arith.misc import CRT
        from sage.misc.misc_c import prod

        # We compute the denominator of the list of p-adics ``padics``.
        # Each `p`-adic only contributes a power of `p` to the denominator:
        # ``a_p.denominator()`` returns a power of the uniformizer.
        denominator = prod([a_p.denominator().lift() for a_p in padics])

        integral_padics = [denominator * a_p for a_p in padics]
        Zhat = ProfiniteIntegers(QQ)
        numerator = Zhat._from_padic_integers(integral_padics)

        return self.element_class(self, numerator, denominator)

    def _coerce_map_from_(self, S):
        """
        Return whether or not ``S`` coerces into this ring of profinite numbers

        EXAMPLES::

            sage: K.<a> = NumberField(x^3+79)
            sage: Khat = ProfiniteNumbers(K)
            sage: Khat._coerce_map_from_(K)
            True
            sage: Khat._coerce_map_from_(CyclotomicField())
            False
            sage: Ohat = ProfiniteIntegers(K)
            sage: Khat._coerce_map_from_(Ohat)
            True
        """
        if self.base().has_coerce_map_from(S):
            return True
        if ProfiniteIntegers(self.base()).has_coerce_map_from(S):
            return True
        return False

    def _an_element_(self):
        """
        Return a typical element of this ring

        EXAMPLES::

            sage: Qhat.an_element()
            2/5 mod 6/5
            sage: K.<a> = NumberField(x^2+7)
            sage: ProfiniteNumbers(K).an_element()
            8/5*a mod (7/5*a - 1/5)
        """
        Ohat = ProfiniteIntegers(self.base())
        return self.element_class(self, Ohat.an_element(), 5)

    def some_elements(self):
        """
        Return some elements of this ring
        
        EXAMPLES::

            sage: Qhat.some_elements()
            [0, 2/5 mod 6/5, 1, 1/3 mod 100]
        """
        return [self.zero(), self.an_element(), self.one(), self(ZZ(1)/ZZ(3), 100)]

    def random_element(self):
        """
        Return a random profinite number

        EXAMPLES::

            sage: [Qhat.random_element() for i in range(10)] # random
            [0 mod 1/2,
             0 mod 4,
             1/2 mod 2,
             3 mod 7/2,
             1/12 mod 1/6,
             0,
             0 mod 2/9,
             1 mod 2,
             5 mod 6,
             1/2 mod 3/4]

        ::

            sage: K.<a> = NumberField(x^3-3)
            sage: Khat = ProfiniteNumbers(K)
            sage: [Khat.random_element() for i in range(10)] # random
            [-a^2 - 7*a - 3 mod (22*a - 22),
             1/3*a^2 + 1/3*a + 1 mod (7/3*a^2 - 7/3),
             0 mod (1/2),
             0 mod (1),
             0 mod (1/2*a - 1/2),
             -1/2*a^2 + a - 1/2 mod (5*a^2 + 10*a + 5),
             -18993*a^2 + 1 mod (26*a^2 + 80*a + 70),
             -5*a^2 + a mod (2*a^2 + 3*a - 3),
             9/5*a^2 - 1/5 mod (3/5*a^2 - 3/5),
             -167*a^2 + a + 1 mod (2*a^2 + 8*a - 10)]
        """
        Ohat = ProfiniteIntegers(self.base())
        numerator = Ohat.random_element()
        denominator = ZZ.random_element(1, 3+abs(ZZ.random_element()))
        return self.element_class(self, numerator, denominator)

    def gens(self):
        """
        Return a tuple of generators of this ring, which is ``(1,)``

        EXAMPLES::

            sage: Qhat.gens()
            (1,)
        """
        return (self(1), )

    def gen(self, n=0):
        """
        Return the ``n``-th generator of this ring

        As this ring has only one generator, we only accept ``n == 0``.

        INPUT:

        - ``n`` -- the index of the generator to return (default: ``0``); must
          be zero

        EXAMPLES::

            sage: Qhat.gen()
            1
            sage: Qhat.gen(0)
            1

        TESTS::

            sage: Qhat.gen(1)
            Traceback (most recent call last):
            ...
            IndexError: n must be 0
        """
        if n == 0:
            return self(1)
        else:
            raise IndexError("n must be 0")

    def ngens(self):
        """
        Return the number of generators of this ring, which is `1`

        EXAMPLES::

            sage: Qhat.ngens()
            1
        """
        return ZZ(1)

    def is_finite(self):
        """
        Return ``False``, indicating that this ring is not finite

        EXAMPLES::

            sage: Qhat.is_finite()
            False
        """
        return False

    def is_exact(self):
        """
        Return ``False``, indicating that doing arithmetic can lead to precision
        loss

        EXAMPLES::

            sage: Qhat.is_exact()
            False
        """
        return False

    def is_integral_domain(self, proof=None):
        """
        Return ``False``, indicating that this ring is not an integral domain

        EXAMPLES::

            sage: Qhat.is_integral_domain()
            False
        """
        return False

    def is_field(self, proof=True):
        """
        Return ``False``, indicating that this ring is not a field

        EXAMPLES::

            sage: Qhat.is_field()
            False
        """
        return False

    def characteristic(self):
        """
        Return the characteristic of this ring, which is zero

        EXAMPLES::

            sage: Qhat = ProfiniteNumbers()
            sage: Qhat.characteristic()
            0
        """
        return ZZ(0)

    def order(self):
        """
        Return ``Infinity``, indicating that this ring has infinitely many
        elements

        EXAMPLES::

            sage: Qhat.order()
            +Infinity
        """
        from sage.rings.infinity import Infinity
        return Infinity

    def epsilon(self):
        """
        Return the precision error of elements in this ring

        As this depends on the elements, we can give no reasonable answer and
        hence raise a ``NotImplementedError``.

        EXAMPLES::

            sage: Qhat.epsilon()
            Traceback (most recent call last):
            ...
            NotImplementedError: precision error depends on the elements involved
        """
        raise NotImplementedError("precision error depends on the elements involved")

    def construction(self):
        """
        Return a pair ``(functor, parent)`` such that ``functor(parent)``
        returns this profinite integers ring.

        EXAMPLES::

            sage: F, P = Qhat.construction(); F, P
            (ProfiniteCompletionFunctor, Rational Field)
            sage: F(P) is Qhat
            True

        ::

            sage: K.<a> = NumberField(x^7-2)
            sage: Khat = ProfiniteNumbers(K)
            sage: F, P = Khat.construction(); F, P
            (ProfiniteCompletionFunctor,
             Number Field in a with defining polynomial x^7 - 2)
            sage: F(P) is Khat
            True
        """
        return ProfiniteCompletionFunctor(), self.number_field()

    def number_field(self):
        """
        Return the base number field of this ring of profinite numbers

        EXAMPLES::

            sage: K.<a> = NumberField(x^5+7*x+9)
            sage: Khat = ProfiniteNumbers(K)
            sage: Khat.number_field()
            Number Field in a with defining polynomial x^5 + 7*x + 9
        """
        return self.base()


Qhat = ProfiniteNumbers(QQ)

