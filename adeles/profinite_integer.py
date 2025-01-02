r"""
Profinite Integers of Number Fields

This file implements rings of profinite integers over number fields and their
elements.

Let `K` be a number field and `O` its ring of integers. For `O`-ideals `I` and
`J` with `I` dividing `J` we have a canonical map `O/J \to O/I`.
The familie `\{O/I\}_I` with `I` running over the non-zero `O`-ideals together
with these canonical maps form a projective system of topological rings
(where each `O/I` is given the discrete topology).
Taking the projective limit results in the topological ring

.. MATH::

    \hat{O} = \varprojlim_I O/I

of *profinite `K`-integers*. It is naturally a commutative `O`-algebra as well.

Profinite Integers over `\QQ`
-----------------------------

In SageMath the ring of profinite `\QQ`-integers `\hat{\ZZ}` and its elements
look like this::

    sage: Zhat
    Ring of Profinite Integers of Rational Field
    sage: Zhat is ProfiniteIntegers(QQ)
    True
    sage: Zhat.an_element()
    2 mod 6

This element ``2 mod 6`` is an approximation to the integer 2 in `\hat{\ZZ}`.
It represents all profinite integers in the open subset `2 + 6\hat{\ZZ}` of
`\hat{\ZZ}`. Hence ``2 mod 6`` also represents the integers 8 and -4, but not
5::

    sage: a = Zhat(2, 6); a
    2 mod 6
    sage: a.represents(8)
    True
    sage: a.represents(-4)
    True
    sage: a.represents(5)
    False
    sage: [n for n in range(-20, 30) if a.represents(n)]
    [-16, -10, -4, 2, 8, 14, 20, 26]

A sum of a profinite integer in `2+6\hat{\ZZ}` and a profinite integer in
`5+15\hat{\ZZ}` will lie in `1+3\hat{\ZZ}`, while any product will lie in
`10+30\hat{\ZZ}`::

    sage: b = Zhat(5, 15); b
    5 mod 15
    sage: a + b
    1 mod 3
    sage: a * b
    10 mod 30

The Chinese Remainder Theorem induces a natural isomorphism

.. MATH::

    \hat{\ZZ} \to \prod_p \ZZ_p

where `p` ranges over all prime numbers and `\ZZ_p` denotes the ring of `p`-adic
integers. This point of view can be taken in SageMath as well::

    sage: a_2 = Zp(2)(7, 4); a_2
    1 + 2 + 2^2 + O(2^4)
    sage: a_5 = Zp(5)(2, 2); a_5
    2 + O(5^2)
    sage: a = Zhat([a_2, a_5]); a
    327 mod 400
    sage: a[2]
    1 + 2 + 2^2 + O(2^4)
    sage: a[3]
    O(3^0)
    sage: a[5]
    2 + O(5^2)
    sage: print(a.str(style='padic'))
    Profinite integer with values:
      at 2: 1 + 2 + 2^2 + O(2^4)
      at 5: 2 + O(5^2)

As SageMath currently has no implementation of completions of number fields at
finite places, the above `p`-adic functionality only works for profinite
`\QQ`-integers.

The general case
----------------

In general, for `K` a number field with ring of integers `O`, a
``ProfiniteInteger`` consists of an element `x \in O`, called its *value*,
and an ideal `I` of `O`, called its *modulus*. Such a ``ProfiniteInteger``
represents all elements in its *represented subset* `x + I \hat{O}`,
similar to how the ``RealInterval`` ``RIF(1.2, 1.3)`` represents all elements
of the subset `[1.2, 1,3]` of `\RR`. ::

    sage: K.<a> = NumberField(x^2+5)
    sage: Ohat = ProfiniteIntegers(K); Ohat
    Ring of Profinite Integers of Number Field in a with defining polynomial x^2 + 5
    sage: b = Ohat(a+1, K.ideal(4, 2*a+2)); b
    a + 1 mod (4, 2*a + 2)
    sage: b.value()
    a + 1
    sage: b.modulus()
    Fractional ideal (4, 2*a + 2)

We always keep the value HNF-reduced (cf. Algorithm 1.4.12 of [Coh2000])::

    sage: c = Ohat(100000, K.ideal(3, a+1)); c
    -a mod (3, a + 1)
    sage: c.value()
    -a

If ``a`` is a ``ProfiniteInteger`` representing a profinite `K`-integer `\alpha
\in \hat{O}` and ``b`` is a ``ProfiniteInteger`` representing `\beta \in
\hat{O}`, then ``a + b`` represents `\alpha + \beta`. Similar statements holds
for the other arithmetic operations. ::

    sage: 3*b + c
    -a mod (3, a + 1)
    sage: b * c
    0 mod (2, a + 1)

See the arithmetic methods such as :meth:`_add_` for details.

In the rest of this file, we shall write "profinite integer" to mean an instance
of ``ProfiniteInteger``, as opposed to element of `\hat{O}`.

.. SEEALSO::

    To see an application of these profinite integers, see for example
    :class:`~adeles.profinite_graph.ProfiniteGraph`. There a graph of the
    profinite Fibonacci function `\hat{\ZZ} \to \hat{\ZZ}` is created using
    these profinite integers.

.. SEEALSO::

    :mod:`~adeles.profinite_number`

.. NOTE::

    Upon creating a number field in SageMath as follows::

        sage: K.<a> = NumberField(x^2+5)

    the element ``a`` has ``K`` as its parent, not the maximal order of ``K``.
    This can cause arithmetic in a ring of profinite integers to give results
    in a ring of profinite *numbers* when the user might not expect it::

        sage: Ohat = ProfiniteIntegers(K)
        sage: b = Ohat(a, 15) + a; b
        2*a mod (15)
        sage: b.parent() is Ohat
        False
        sage: b.parent() is ProfiniteNumbers(K)
        True

    Although above the element ``a`` lies in the maximal order of ``K``, its
    parent is ``K``. Hence the coercion model will look for a common parent of
    ``Ohat(a, 15)`` and ``a``, which is not ``Ohat``, but the ring of profinite
    *numbers* over ``K``.

REFERENCES:

[Her2021] Mathé Hertogh, Computing with adèles and idèles, master's thesis,
Leiden University, 2021.

This implementation of profinite integers is based on [Her2021]. An extensive
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

from sage.categories.rings import Rings
from sage.categories.integral_domains import IntegralDomains
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.pushout import ConstructionFunctor
from sage.rings.ring import CommutativeAlgebra
from sage.structure.element import CommutativeAlgebraElement
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.arith.misc import gcd
from sage.sets.primes import Primes


class ProfiniteInteger(CommutativeAlgebraElement):
    """
    Profinite Integer over a Number Field

    REFERENCES:

    Section 3.1 of [Her2021].

    .. automethod:: __init__

    .. automethod:: __getitem__

    .. automethod:: _add_

    .. automethod:: _sub_

    .. automethod:: _mul_

    .. automethod:: _div_

    .. automethod:: _floordiv_

    .. automethod:: _richcmp_
    """

    def __init__(self, parent, value, modulus):
        r"""
        Construct the profinite integer "``value`` mod ``modulus``"

        INPUT:

        - ``parent`` -- a ``ProfiniteIntegers`` object for some base number
          field `K`
        - ``value`` -- an integral element of `K`
        - ``modulus`` -- if `K` is `\QQ`: an integer; else: an ideal of `O`,
          the maximal order of `K`

        OUTPUT:

        The profinite integer representing the subset
        ``value`` + ``modulus`` * `\hat{O}` of `\hat{O}`, the profinite
        completion of `O`.

        EXAMPLES::

            sage: ProfiniteInteger(Zhat, 4, 15)
            4 mod 15

        ::

            sage: K.<a> = NumberField(x^7+3)
            sage: Ohat = ProfiniteIntegers(K)
            sage: ProfiniteInteger(Ohat, a^6-a, 30*a^3)
            a^6 - a mod (30*a^3)
            sage: I = K.ideal(a^5+3*a^2, 120)
            sage: ProfiniteInteger(Ohat, 7, I)
            14*a^6 + 1 mod (-a^6 + a^5 + 3*a - 3)

        TESTS::

            sage: Zhat(3, -10)
            3 mod 10
            sage: ProfiniteInteger(Ohat, a^2, 0)
            a^2
            sage: ProfiniteInteger(Ohat, None, a)
            Traceback (most recent call last):
            ...
            TypeError: value must be an element of Maximal Order in Number Field in a with defining polynomial x^7 + 3
            sage: ProfiniteInteger(Ohat, a, None)
            Traceback (most recent call last):
            ...
            TypeError: modulus must be an ideal of Maximal Order in Number Field in a with defining polynomial x^7 + 3
            sage: ProfiniteInteger(Ohat, a, 1/a)
            Traceback (most recent call last):
            ...
            TypeError: modulus must be an ideal of Maximal Order in Number Field in a with defining polynomial x^7 + 3
            sage: ProfiniteInteger(Zhat, 3, 1/25)
            Traceback (most recent call last):
            ...
            TypeError: modulus must be an integer
        """
        CommutativeAlgebraElement.__init__(self, parent)
        O = parent.base()
        K = parent.number_field()
        # if value not in O:
        if value not in K or not K(value).is_integral():
            raise TypeError("value must be an element of {}".format(O))
        self._value = O(value)
        if O is ZZ:
            if modulus not in ZZ:
                raise TypeError("modulus must be an integer")
            if modulus < ZZ(0):
                modulus = -modulus
            self._modulus = ZZ(modulus)
        else:
            if modulus not in K.ideal_monoid() or not modulus.is_integral():
                raise TypeError("modulus must be an ideal of {}".format(O))
            modulus = O.ideal(modulus)
            self._modulus = O.ideal(modulus.gens_reduced())
        self._reduce()

    def __getitem__(self, p):
        r"""
        Return the projection of this profinite integer to the ring of
        ``p``-adic integers

        Only implemented for profinite integers over `\QQ`, as no general
        implementation of completions at finite places of a number field
        exists in SageMath at the time of writing.
        
        The Chinese Remainder Theorem induces a canonical isomorphism
        `\hat{\ZZ} \cong \prod_p \ZZ_p`, where `p` runs over all prime numbers.
        This method computes the image of ``self`` under this isomorphism and
        returns the projection to the ``p``-th coordinate. 

        INPUT:

        - ``p`` -- a prime number

        EXAMPLES::

            sage: a = Zhat(4, 27)
            sage: a[3]
            1 + 3 + O(3^3)
            sage: b = Zhat(5, 2^6 * 5^12)
            sage: b[2]
            1 + 2^2 + O(2^6)
            sage: b[3]
            O(3^0)
            sage: b[5]
            5 + O(5^12)
            sage: c = Zhat(-1, 7^5)
            sage: c[7]
            6 + 6*7 + 6*7^2 + 6*7^3 + 6*7^4 + O(7^5)
            sage: c[7].parent()
            7-adic Ring with capped relative precision 20

        TESTS::

            sage: d[2]
            Traceback (most recent call last):
            ...
            NotImplementedError: projection to `p`-adics only implemented over rationals
            sage: d[K.prime_above(5)]
            Traceback (most recent call last):
            ...
            NotImplementedError: projection to `p`-adics only implemented over rationals
            sage: c[-1]
            Traceback (most recent call last):
            ...
            ValueError: p must be a prime number
            sage: c[6]
            Traceback (most recent call last):
            ...
            ValueError: p must be a prime number
        """
        from sage.rings.padics.factory import Zp
        if self.parent().number_field() is not QQ:
            raise NotImplementedError("projection to `p`-adics only implemented over rationals")
        if p not in Primes():
            raise ValueError("p must be a prime number")
        e = self.modulus().valuation(p)
        return Zp(p)(self.value(), e)

    def _repr_(self):
        """
        Return a string representation of this profinite integer

        EXAMPLES::

            sage: Zhat(6, 20)
            6 mod 20

        ::

            sage: K.<a> = NumberField(x^2+5)
            sage: Ohat = ProfiniteIntegers(K)
            sage: Ohat(a+1, 60)
            a + 1 mod (60)
            sage: I = K.ideal(3, a+2)
            sage: Ohat(a, I)
            a mod (3, a + 2)
        """
        if self.modulus() == 0:
            return repr(self.value())

        if self.parent().base() is ZZ:
            modulus = self.modulus()
        else:
            modulus = "(" + ", ".join([str(g) for g in self.modulus().gens()]) + ")"

        return "{} mod {}".format(self.value(), modulus)

    def _add_(self, other):
        r"""
        Return the sum of this profinite integer and ``other``

        The sum of two profinite integers `a` and `b` is defined to be the
        profinite integer `c` with smallest represented subset (with respect to
        inclusion) containing the sum of the represented subsets of `a` and `b`,
        i.e. containig `\alpha + \beta` for all `\alpha` represented by `a` and
        `\beta` represented by `b`.

        EXAMPLES::
        
            sage: Zhat = ProfiniteIntegers(ZZ)
            sage: a = Zhat(3, 20)
            sage: b = Zhat(-1, 30)
            sage: a + b
            2 mod 10
            sage: a + 3
            6 mod 20

        ::

            sage: K.<a> = NumberField(x^2-7)
            sage: O = K.maximal_order()
            sage: Ohat = ProfiniteIntegers(O)
            sage: b = Ohat(3*a+1, a*13)
            sage: c = Ohat(2*a+1, 7)
            sage: b+c
            2 mod (a)
            sage: O(-a)+c
            a + 1 mod (7)

        REFERENCES:

        Section Arithmetic of Section 3.1 of [Her2021].
        """
        modulus = self._common_modulus(other)
        value = self.value() + other.value()
        return self.__class__(self.parent(), value, modulus)

    def _sub_(self, other):
        r"""
        Return the difference of this profinite integer and ``other``

        The difference of two profinite integers `a` and `b` is defined to be
        the profinite integer with smallest represented subset (with respect to
        inclusion) containing the sum of the represented subsets of `a` and `b`,
        i.e. containig `\alpha + \beta` for all `\alpha` represented by `a` and
        `\beta` represented by `b`.

        EXAMPLES::

            sage: Zhat = ProfiniteIntegers(ZZ)
            sage: a = Zhat(3, 20)
            sage: b = Zhat(-1, 30)
            sage: a-b
            4 mod 10
            sage: b-0
            29 mod 30

        ::

            sage: K.<a> = NumberField(x^2-7)
            sage: Ohat = ProfiniteIntegers(K)
            sage: b = Ohat(3*a+1, a*13)
            sage: c = Ohat(2*a+1, 7)
            sage: b-c
            0 mod (a)
            sage: c-1
            2*a mod (7)

        REFERENCES:

        Section Arithmetic of Section 3.1 of [Her2021].
        """
        modulus = self._common_modulus(other)
        value = self.value() - other.value()
        return self.__class__(self.parent(), value, modulus)

    def _mul_(self, other):
        r"""
        Return the product of this profinite integer and ``other``

        The product of two profinite integers `a` and `b` is defined to be
        the profinite integer with smallest represented subset (with respect to
        inclusion) containing the sum of the represented subsets of `a` and `b`,
        i.e. containig `\alpha + \beta` for all `\alpha` represented by `a` and
        `\beta` represented by `b`.

        EXAMPLES::

            sage: Zhat = ProfiniteIntegers(ZZ)
            sage: a = Zhat(3, 20)
            sage: b = Zhat(-1, 30)
            sage: a*b
            7 mod 10
            sage: 0*a
            0

        ::

            sage: K.<a> = NumberField(x^2-7)
            sage: Ohat = ProfiniteIntegers(K)
            sage: b = Ohat(3*a+1, a*13)
            sage: c = Ohat(2*a+1, 7)
            sage: b*c
            1 mod (a)
            sage: 1*c
            2*a + 1 mod (7)

        REFERENCES:

        Section Arithmetic of Section 3.1 of [Her2021].

        .. TODO::

            Currently this multiplication method implicitly defines powering of
            profinite numbers::

                sage: a = Zhat(1, 2)
                sage: a^2
                1 mod 2

            This is however not optimal: for any `\alpha \in 1 + 2\hat{\ZZ}` we
            have `\alpha^2 \in 1 + 8\hat{\ZZ}`. Hence by implementing powering
            explicitly we could obtain the better behaviour::

                sage: a * a
                1 mod 2
                sage: a^2 # not tested
                1 mod 8
        """
        # We seperate the zero-modulus case because multiplying a non-principal
        # ideal with the zero-ideal is not implemented in Sage.
        if self.modulus().is_zero():
            modulus = self.value() * other.modulus()
        elif other.modulus().is_zero():
            modulus = other.value() * self.modulus()
        else:
            I = self.value() * other.modulus()
            J = self.modulus() * other.value()
            K = self.modulus() * other.modulus()
            if self.parent().base() is ZZ:
                modulus = gcd(I, gcd(J, K))
            else:
                modulus = I + J + K
        value = self.value() * other.value()
        return self.__class__(self.parent(), value, modulus)

    def _div_(self, other):
        """
        Divide this profinite integer by ``other`` in the ring of profinite
        *numbers* and return the result

        Only implemented for ``other`` having zero-modulus and non-zero value.

        EXAMPLES::

            sage: b = Zhat(4, 10) / 7; b
            4/7 mod 10/7
            sage: b.parent()
            Profinite Numbers of Rational Field

        ::

            sage: K.<a> = NumberField(x^2+x-7)
            sage: Ohat = ProfiniteIntegers(K)
            sage: Ohat(a^2, 20) / a
            a mod (20/7*a + 20/7)

        TESTS::

            sage: Zhat(4, 10) / Zhat(2, 10)
            Traceback (most recent call last):
            ...
            NotImplementedError: Cannot divide by profinite integers of non-zero modulus
        """
        if other.modulus() != 0:
            raise ValueError("division by profinite integer with non-zero modulus")
        if other.value() == 0:
            raise ZeroDivisionError("profinite integer division by zero")
        from .profinite_number import ProfiniteNumbers
        Khat = ProfiniteNumbers(self.parent().number_field())
        return Khat(self, other.value())

    def _floordiv_(self, other):
        """
        Return the quotient of this profinite integer by ``other``

        Only implemented when ``other`` has zero modulus and ``self.value()``
        and ``self.modulus()`` are both divisible by ``other.value()``.

        EXAMPLES::

            sage: Zhat(5, 100) // 5
            1 mod 20

        ::

            sage: K.<a> = NumberField(x^2+5)
            sage: Ohat = ProfiniteIntegers(K)
            sage: Ohat(9*a, 75) // O(3*a)
            3 mod (-5*a)

        TESTS::

            sage: b // Zhat(5, 10)
            Traceback (most recent call last):
            ...
            TypeError: can only divide by elements of the base
            sage: b/0
            Traceback (most recent call last):
            ...
            ZeroDivisionError: division by zero
            sage: b/3
            Traceback (most recent call last):
            ...
            NotImplementedError: value is not divisible by 3
            sage: c/9
            Traceback (most recent call last):
            ...
            NotImplementedError: modulus is not divisible by 9
        """
        if other.modulus() != 0:
            raise ValueError("division by profinite integer with non-zero modulus")
        if other.value() == 0:
            raise ZeroDivisionError("profinite integer division by zero")
        if not other.value().divides(self.value()):
            raise ValueError("value is not divisible by {}".format(other))
        if not other.value().divides(self.modulus()):
            raise ValueError("modulus is not divisible by {}".format(other))

        value = self.value() // other.value()
        modulus = 0 if self.modulus() == 0 else self.modulus() / other.value()

        return self.__class__(self.parent(), value, modulus)

    def _richcmp_(self, other, op):
        r"""
        Compare this profinite integer to ``other`` based on the relation ``op``

        We only implement equality and non-equality.

        Two elements ``x mod I`` and ``y mod J`` are considered equal if and
        only if `x \equiv y \mod (I+J)`.

        EXAMPLES::

            sage: K.<a> = NumberField(x^4+17)
            sage: O = K.maximal_order()
            sage: Ohat = ProfiniteIntegers(O)
            sage: b = O(a^2+1)
            sage: c = Ohat(a^2+1, 1000)
            sage: d = Ohat(a^2+11, a*10)
            sage: e = Ohat(a^3+a-5, 170)
            sage: c == c
            True

        Transitivity of equality does *not* hold::

            sage: b == c
            True
            sage: c == d
            True
            sage: b == d
            False

        ::

            sage: c == e
            False
            sage: c != e
            True
            sage: b != d
            True
            sage: c != d
            False

        REFERENCES:

        Section 5.4 of [Her2021].

        TESTS::

            sage: b < d
            Traceback (most recent call last):
            ...
            NotImplementedError: only equality and inequality are implemented
        """
        from sage.structure.richcmp import op_EQ, op_NE
        if op == op_EQ:
            O = self.parent().base()
            common_modulus = O.ideal(self._common_modulus(other))
            return (self.value() - other.value()) in common_modulus
        if op == op_NE:
            return not self._richcmp_(other, op_EQ)
        raise NotImplementedError("only equality and inequality are implemented")

    def _reduce(self):
        """
        HNF-reduce ``self.value()`` modulo ``self.modulus()``

        EXAMPLES::

            sage: a = Zhat(0, 10)
            sage: a._value = 91; a
            91 mod 10
            sage: a._reduce(); a
            1 mod 10

        ::

            sage: K.<a> = NumberField(x^2-6)
            sage: Ohat = ProfiniteIntegers(K)
            sage: b = Ohat(0, a+1)
            sage: b._value = 10*a+11; b
            10*a + 11 mod (a + 1)
            sage: b._reduce(); b
            -a mod (a + 1)

        TESTS::

            sage: b._modulus = K.ideal(0)
            sage: b._reduce(); b
            -a
        """
        if self.parent().base() is ZZ:
            if self.modulus() != 0:
                self._value %= self.modulus()
        else:
            self._value = self.modulus().reduce(self.value())

    def _common_modulus(self, other):
        """
        Return the greatest common divisor of the moduli of this profinite
        integer and ``other``

        INPUT:

        - ``other`` -- a profinite integer with the same parent as ``self``

        EXAMPLES::

            sage: a = Zhat(1, 96)
            sage: b = Zhat(79, 120)
            sage: a._common_modulus(b)
            24

        ::

            sage: K.<a> = NumberField(x^2-3)
            sage: Ohat = ProfiniteIntegers(K)
            sage: b = Ohat(a, 15)
            sage: c = Ohat(0, 25*a)
            sage: b._common_modulus(c)
            Fractional ideal (5*a)

        TESTS::

            sage: b._common_modulus(Ohat(7))
            Fractional ideal (15)
        """
        if self.parent().base() is ZZ:
            return gcd(self.modulus(), other.modulus())
        return self.modulus() + other.modulus()

    def str(self, style="value-modulus"):
        """
        Return a string representation of this profinite integer

        INPUT:

        - ``style`` -- string (default: "value-modulus"); style to use. The
          other options are "factorial" and "padic".

        EXAMPLES::

            sage: a = Zhat(970, 2^7 * 3^3 * 5^2)
            sage: print(a.str(style='value-modulus'))
            970 mod 86400
            sage: print(a.str(style='factorial'))
            2*2! + 1*3! + 2*5! + O(6!)
            sage: print(a.str(style='padic'))
            Profinite integer with values:
              at 2: 2 + 2^3 + 2^6 + O(2^7)
              at 3: 1 + 2*3 + 2*3^2 + O(3^3)
              at 5: 4*5 + O(5^2)

        TESTS::

            sage: print(a.str(style='blah'))
            Traceback (most recent call last):
            ...
            ValueError: unkown style ``blah''
        """
        if style == "value-modulus":
            return self._repr_()
        if style == "factorial":
            d = self.factorial_digits()
            s = "".join(["{}*{}! + ".format(d[i], i+1) for i in range(len(d)) if d[i] != 0])
            s += "O({}!)".format(len(d)+1)
            return s
        if style == "padic":
            padics = ["at {}: {}".format(p, str(self[p])) for p, e in self.modulus().factor()]
            return "Profinite integer with values:\n  " + "\n  ".join(padics)
        raise ValueError("unkown style ``{}''".format(style))

    def value(self):
        """
        Return the value of this profinite integer

        EXAMPLES::

            sage: Zhat(3, 6).value()
            3
            sage: Zhat(7, 6).value() # the value is reduced
            1
        """
        return self._value

    def modulus(self):
        """
        Return the modulus of this profinite integer

        This is an integer if the base number field is ``QQ``.
        Otherwise it is an ideal of the base maximal order.

        EXAMPLES::

            sage: Zhat(50, 100).modulus()
            100

        ::

            sage: K.<a> = NumberField(x^4-17)
            sage: Ohat = ProfiniteIntegers(K)
            sage: Ohat(a^3, a*160).modulus()
            Fractional ideal (2720, 160*a)
        """
        return self._modulus

    def represents(self, element):
        r"""
        Return whether or not this profinite integer represents ``element``

        The *represented subset* of this profinite integer is `x + m \hat{O}`,
        where `x` is our value, `m` is our modulus and `\hat{O}` is the ring
        of profinite numbers over our base field.
        Hence this profinite number represents `\alpha \in \hat{O}` if and only
        if `\alpha \in x + m \hat{O}`.

        INPUT:

        - ``element`` -- an integral element of the base number field

        EXAMPLES::

            sage: b = Zhat(3, 10)
            sage: [n for n in range(-50, 50) if b.represents(n)]
            [-47, -37, -27, -17, -7, 3, 13, 23, 33, 43]

        ::

            sage: K.<a> = NumberField(x^2+2)
            sage: Ohat = ProfiniteIntegers(K)
            sage: c = Ohat(a+1, 3*a)
            sage: [m*a+n for m in range(-7, 7) for n in range(-7, 7) if c.represents(m*a+n)]
            [-5*a - 5, -5*a + 1, -2*a - 5, -2*a + 1, a - 5, a + 1, 4*a - 5, 4*a + 1]
        """
        if self.parent().base() is ZZ:
            return (self.value() - element) % self.modulus() == 0
        return self.value() - element in self.modulus()

    def is_unit(self):
        r"""
        Return whether or not ``self`` could be a unit

        More precisely, we return ``True`` if and only if the subset of 
        profinite integers that ``self`` represents contains a unit.

        If ``self`` is ``x mod m``, this is equivalent to `x \in (O/mO)^*`,
        where `O` denotes our base ring of integers.

        EXAMPLES::

            sage: Zhat(3, 0).is_unit()
            False
            sage: Zhat(-1, 0).is_unit()
            True
            sage: Zhat(3, 25).is_unit()
            True
            sage: Zhat(5, 25).is_unit()
            False

        ::

            sage: K.<a> = NumberField(x^3-2)
            sage: Ohat = ProfiniteIntegers(K)
            sage: Ohat(3, 8).is_unit()
            True
            sage: Ohat(a, 8).is_unit()
            False
        """
        O = self.parent().base()
        x = O.ideal(self.value())
        # We check if x is coprime to our modulus:
        return x + self.modulus() == 1

    def is_integral(self):
        """
        Return ``True``, indicating that this profinite integer is integral

        EXAMPLES::

            sage: Zhat(-5, 12).is_integral()
            True
        """
        return True

    def factorial_digits(self):
        r"""
        Return the factorial digits of this profinite integer

        Only implemented for profinite integers over `\QQ` with non-zero
        modulus.
        
        Let `k` be the largest integer such that `k!` divides our modulus (this
        is called our "factorial precision"). Then our factorial digits are the
        unique integers `d_1, d_2, ..., d_{k-1}` satisfying `0 \leq d_i \leq i`
        such that

        .. MATH::

            x \equiv d_1 \cdot 1! + d_2 \cdot 2! + ... + d_{k-1} \cdot (k-1)! \mod k!
        
        where `x` denotes the value of this profinite integer.

        EXAMPLES::

            sage: digits = Zhat(11, 48).factorial_digits(); digits
            [1, 2, 1]
            sage: sum([digits[i] * factorial(i+1) for i in range(3)])
            11

        REFERENCES:

        Section 7.2 of [Her2021].
        """
        if self.parent().number_field() is not QQ:
            raise NotImplementedError("only implemented for profinite QQ-integers")
        if self.modulus() == ZZ(0):
            raise NotImplementedError("not implemented for zero-moduli")

        # Calculate our factorial precision prec, which is the largest k such
        # that k! divides our modulus.
        prec = ZZ(1)
        prec_factorial = ZZ(1)
        while self.modulus() % (prec_factorial*(prec+1)) == 0:
            prec += 1
            prec_factorial *= prec

        value = self.value()
        digit = 0
        digits = []
        for k in range(1, prec):
            value = (value - digit) // k
            digit = value % (k+1)
            digits.append(digit)

        return digits

    def visual(self):
        r"""
        Return the smallest closed interval within the unitinterval `[0,1]`
        in which the image of this profinite integer under the visualization map
        is contained

        The visualization function is defined by
        
        .. MATH::

            \phi(a) = \sum_{i=1}^\infty \frac{d_i}{(1+i)!}

        for `a \in \hat{\ZZ}` with factorial digit sequence
        `(d_i)_{i=1}^\infty`.

        EXAMPLES:

        The subset `1 + 2\hat{\ZZ}` of `\hat{\ZZ}` maps onto `[1/2, 1]` by
        `\phi`. Hence we have ::

            sage: Zhat(1, 2).visual()
            (1/2, 1)

        The subset `2 + 3\hat{\ZZ}` is mapped onto `[1/6, 1/3] \cup [5/6, 1]` by
        `\phi`. So we get ::

            sage: Zhat(2, 3).visual()
            (1/6, 1)

        REFERENCES:

        Section 7.3 of [Her2021].
        """
        from sage.functions.other import factorial

        # We compute the smallest positive integer prec such that self.modulus()
        # divides prec!:
        prec = ZZ(1)
        while not self.modulus().divides(factorial(prec)):
            prec += 1

        lefts, rights = [], []
        for k in range(0, factorial(prec)/self.modulus()):
            a = Zhat(self.value() + k * self.modulus(), factorial(prec))
            digits = a.factorial_digits()
            left = sum([digits[i-1] / factorial(i+1) for i in range(1, prec)])
            right = left + 1/factorial(prec)
            lefts.append(left)
            rights.append(right)

        return min(lefts), max(rights)


class ProfiniteIntegers(UniqueRepresentation, CommutativeAlgebra):
    """
    Ring of profinite integers over a number field

    REFERENCES:

    Section 3.1 of [Her2021].

    .. automethod:: _element_constructor_
    """

    Element = ProfiniteInteger

    @staticmethod
    def __classcall__(cls, R=ZZ):
        """
        Construct the ring of profinite integers over ``R``

        INPUT:

        - ``R`` -- a number field or a maximal order in a number field (default:
          ``ZZ``)

        This method makes sure input to ``UniqueRepresentation`` (in particular,
        ``CachedRepresentation``) is normalized such that a number field and its
        maximal order return the same object::

            sage: ProfiniteIntegers(ZZ) is ProfiniteIntegers(QQ)
            True
        """
        from sage.rings.number_field.number_field import is_NumberField
        if is_NumberField(R):
            O = R.maximal_order()
        elif R is ZZ:
            O = ZZ
        else:
            try:
                K = R.ambient()
            except AttributeError:
                raise TypeError("R must be (the maximal order of) a number field")
            if not is_NumberField(K) or R != K.maximal_order():
                raise TypeError("R must be (the maximal order of) a number field")
            O = R
        return super(ProfiniteIntegers, cls).__classcall__(cls, O)

    def __init__(self, O):
        """
        Construct the ring of profinite ``O``-integers

        INPUT:

        - ``O`` -- the maximal order of a number field

        EXAMPLES::

            sage: ProfiniteIntegers()
            Profinite Integers of Rational Field
            sage: K.<a> = NumberField(x^5-3*x+1)
            sage: ProfiniteIntegers(K)
            Profinite Integers of Number Field in a with defining polynomial x^5 - 3*x + 1
            sage: ProfiniteIntegers("not a number field")
            Traceback (most recent call last):
            ...
            TypeError: R should be (the maximal order of) a number field

        TESTS:

        Pickling works::

            sage: import __main__
            sage: __main__.ProfiniteIntegers = ProfiniteIntegers
            sage: loads(dumps(ProfiniteIntegers(QQ))) is ProfiniteIntegers(QQ)
            True
            sage: loads(dumps(ProfiniteIntegers(K))) is ProfiniteIntegers(K)
            True
        """
        # Note that the input O is checked by __classcall__().
        CommutativeAlgebra.__init__(self, O)

    def _repr_(self):
        """
        Return a string representation of this ring of profinite integers

        EXAMPLES::

            sage: K = NumberField(x^3+x+1, 'a')
            sage: ProfiniteIntegers(K)
            Ring of Profinite Integers of Number Field in a with defining polynomial x^3 + x + 1
        """
        return "Ring of Profinite Integers of {}".format(self.number_field())

    def _latex_(self):
        r"""
        Return a latex-formatted string representation of this ring of profinite
        integers

        EXAMPLES::

            sage: K.<a> = NumberField(x^2+14)
            sage: Ohat = ProfiniteIntegers(K)
            sage: latex(Ohat)
             \widehat{ \mathcal{O}_{ \Bold{Q}[a]/(a^{2} + 14) } }
        """
        from sage.misc.latex import latex
        return r" \widehat{ \mathcal{O}_{" + latex(self.number_field()) + "} } "

    def _element_constructor_(self, x, y=None):
        r"""
        Construct a profinite integer

        INPUT:

        We accept many input formats. The most common one is the following:

        - ``x`` -- element of the base maximal order; the value
        - ``y`` -- ideal of the base maximal order; the modulus

        All other formats only accept one argument, which can be one of the
        following:

        - a profinite `K`-integer for `K` a subfield of our base number field
        - an integral profinite `K`-number for `K` a subfield of our base number
          field
        - an element of a quotient of our base maximal order
        - an element of our base maximal order

        For constructing profinite `\QQ`-integers, even more formats are
        accepted, namely:

        - a factorial digit list, i.e. a list of non-negative integers such that
          the `i`-th entry is at most `i+1` (see
          :meth:`~adeles.profinite_integer.ProfiniteInteger.factorial_digits`).
        - a `p`-adic integer, for some prime number `p`
        - a list of `p`-adic integers for distinct prime numbers `p`

        EXAMPLES:

        We start with standard ``(value, modulus)`` input::

            sage: Zhat(-4, 17)
            13 mod 17
            sage: K.<a> = NumberField(x^3-5*x^2+1)
            sage: Ohat = ProfiniteIntegers(K)
            sage: Ohat(a^2+1, 4*a)
            a^2 + 1 mod (4*a)

        Upon giving a ``ProfiniteInteger`` as input we get::

            sage: Zhat(Zhat(3, 5))
            3 mod 5
            sage: Ohat(Zhat(4, 6))
            -2 mod (6)
            sage: Ohat(Ohat(a, 8))
            a mod (8)

        Integral profinite numbers can be converted to profinite integers::

            sage: Zhat(Qhat(-1, 5))
            4 mod 5
            sage: Khat = ProfiniteNumbers(K)
            sage: Ohat(Khat(a+1, 6))
            a + 1 mod (6)

        Quotient ring elements can be given as input as well::

            sage: Zhat(Zmod(40)(7))
            7 mod 40
            sage: R = K.maximal_order().quotient(30*a, 'b')
            sage: Ohat(R(2*a+3))
            2*a + 3 mod (30*a)

        Elements of the base maximal order a coerced to profinite integers with
        zero modulus::

            sage: Zhat(97)
            97
            sage: Ohat(79*a)
            79*a

        Profinite integers over `\QQ` can be constructed from factorial digits::

            sage: Zhat([1, 0, 0, 2, 0])
            49 mod 720

        And profinite integers over `\QQ` can also be constructed from `p`-adic
        integers::

            sage: Zhat(Zp(5)(-1))
            95367431640624 mod 95367431640625
            sage: Zhat([Zp(2)(20, 5), Zp(3)(7, 2)])
            52 mod 288
        """
        from .profinite_number import ProfiniteNumbers
        from sage.rings.quotient_ring import is_QuotientRing
        import sage.rings.abc

        is_pAdic = lambda P: isinstance(P, (sage.rings.abc.pAdicRing, sage.rings.abc.pAdicField))

        if y is None:
            P = x.parent() if hasattr(x, "parent") else None

            # Check if x is a profinite K-integer for some subfield K of our base field:
            if isinstance(P, ProfiniteIntegers) and self.number_field().has_coerce_map_from(P.number_field()):
                return self.element_class(self, x.value(), x.modulus())

            # Check if x is a profinite K-number for some subfield K of our base field:
            if isinstance(P, ProfiniteNumbers) and self.number_field().has_coerce_map_from(P.number_field()):
                return self._from_profinite_number(x)
            
            # Check if x is an element of a quotient of our base maximal order:
            if is_QuotientRing(P) and P.ambient() in [self.base(), self.number_field()]:
                    return self._from_modulo_element(x)

            if self.number_field() is QQ:
                # Check if x is a list of factorial digits:
                if isinstance(x, list) and all([x[i] in ZZ and 0 <= x[i] <= i+1 for i in range(len(x))]):
                    return self._from_factorial_digits(x)
                
                # Check if x is a `p`-adic number:
                if is_pAdic(P):
                    return self._from_padic_integers([x])
                
                # Check if x is a list of `p`-adic numbers:
                if isinstance(x, list) and all([hasattr(a, "parent") and is_pAdic(a.parent()) for a in x]):
                        return self._from_padic_integers(x)

            # No conversions seem to fit. Hence we just expect x to be an
            # element of our base maximal order and we coerce it to a profinite
            # integer with modulus zero.
            y = ZZ(0) if self.base() is ZZ else self.base().ideal(0)

        return self.element_class(self, x, y)

    def _from_modulo_element(self, element):
        """
        Construct a profinite integer from the element ``element`` of a quotient
        of our base maximal order

        EXAMPLES::

            sage: b = Zmod(25)(10)
            sage: Zhat._from_modulo_element(b)
            10 mod 25
            sage: Zhat(b)
            10 mod 25

        ::

            sage: K.<a> = NumberField(x^3-2)
            sage: O = K.maximal_order()
            sage: I = O.ideal(5*a)
            sage: OmodI = O.quotient_ring(I, 'z')
            sage: c = OmodI(3-a)
            sage: Ohat = ProfiniteIntegers(O)
            sage: Ohat._from_modulo_element(c)
            -a + 3 mod (5*a)
            sage: Ohat(c)
            -a + 3 mod (5*a)

        REFERENCES:

        Section 4.2 of [Her2021].
        """
        OmodI = element.parent()
        modulus = OmodI.defining_ideal()
        if self.base() is ZZ:
            modulus = modulus.gens_reduced()[0]
        value = element.lift()
        return self.element_class(self, value, modulus)

    def _from_profinite_number(self, number):
        """
        Construct a profinite integer from the profinite number ``number``

        INPUT:

        - ``number`` -- an integral profinite number over a subfield of our base
          number field

        EXAMPLES::

            sage: K.<a> = NumberField(x^2-5)
            sage: Khat = ProfiniteNumbers(K)
            sage: Ohat = ProfiniteIntegers(O)
            sage: b = Khat((a+1)/2, 10)
            sage: b.denominator()
            1
            sage: Ohat._from_profinite_number(b)
            1/2*a + 1/2 mod (10)

        If ``number`` is not integral then we throw an exception::

            sage: Ohat._from_profinite_number(Khat(a, 1/2))
            Traceback (most recent call last):
            ...
            ValueError: can't convert non-integral profinite number to profinite integer

        REFERENCES:

        Section 4.1 of [Her2021].
        """
        if not number.is_integral():
            raise ValueError("can't convert non-integral profinite number to profinite integer")
        return self.element_class(self, number.value(), number.modulus())

    def _from_padic_integers(self, padics):
        """
        Construct the profinite integer corresponding to ``padics``

        INPUT:

        - ``padics`` -- a list of `p`-adic integers for distinct prime number
          `p`

        EXAMPLES::

            sage: a_2 = Zp(2)(17, 5); a_2
            1 + 2^4 + O(2^5)
            sage: a_5 = Zp(5)(17, 3); a_5
            2 + 3*5 + O(5^3)
            sage: a = Zhat([a_2, a_5]); a
            17 mod 4000
            sage: factor(a.modulus())
            2^5 * 5^3
        
        ::

            sage: b_2 = Qp(2)(-1, 4); b_2
            1 + 2 + 2^2 + 2^3 + O(2^4)
            sage: b_3 = Zp(3)(3, 2); b_3
            3 + O(3^2)
            sage: b_7 = Qp(7)(8, 3); b_7
            1 + 7 + O(7^3)
            sage: b = Zhat([b_7, b_2, b_3]); b
            16815 mod 49392
            sage: b[2]
            1 + 2 + 2^2 + 2^3 + O(2^4)
            sage: b[3]
            3 + O(3^2)
            sage: b[7]
            1 + 7 + O(7^3)

        TESTS::

            sage: b_5 = Qp(5)(1/5)
            sage: Zhat([b_2, b_3, b_5, b_7])
            Traceback (most recent call last):
            ...
            ValueError: non-integral p-adic in profinite integer initialization
            sage: Zhat([b_2, b_3, b_2])
            Traceback (most recent call last):
            ...
            ValueError: multiple 2-adic integers in profinite integer initialization
        """
        if any(not a_p.is_integral() for a_p in padics):
            raise ValueError("non-integral p-adic in profinite integer initialization")
        primes = set()
        for a_p in padics:
            p = a_p.parent().prime()
            if p in primes:
                raise ValueError("multiple {}-adic integers in profinite integer initialization".format(p))
            primes.add(p)

        from sage.arith.misc import CRT
        from sage.misc.misc_c import prod

        values, moduli = [], []
        for a_p in padics:
            p = a_p.parent().prime()
            values.append(a_p.lift())
            moduli.append(p**a_p.precision_absolute())
        value = CRT(values, moduli)
        modulus = prod(moduli)
        return self.element_class(self, value, modulus)

    def _from_factorial_digits(self, digits):
        """
        Construct the profinite integer corresponding to the factorial digits
        ``digits``

        This is the ``ProfiniteInteger`` with smallest represented subset such
        that the start of the factorial digit sequence of each profinite integer
        it represents is given by ``digits`` (cf.
        :meth:`~adeles.profinite_integer.ProfiniteInteger.factorial_digits`).

        EXAMPLES::

            sage: Zhat._from_factorial_digits([1, 0, 2])
            13 mod 24
            sage: Zhat._from_factorial_digits([1, 0, 2, 0])
            13 mod 120
            sage: Zhat._from_factorial_digits([])
            0 mod 1

        REFERENCES:

        Section 7.2 of [Her2021].
        """
        from sage.functions.other import factorial
        value = sum([digits[i] * factorial(i+1) for i in range(len(digits))])
        modulus = factorial(len(digits)+1)
        return self.element_class(self, value, modulus)

    def _coerce_map_from_(self, S):
        """
        Return whether or not ``S`` coerces into this ring of profinite integers

        EXAMPLES::

            sage: K.<a> = NumberField(x^3+79)
            sage: Ohat = ProfiniteIntegers(K)
            sage: Ohat._coerce_map_from_(K.maximal_order())
            True
            sage: Ohat._coerce_map_from_(ZZ)
            True
            sage: Ohat._coerce_map_from_(Zhat)
            True
            sage: R = K.maximal_order().quotient(7*a^2, 'b')
            sage: Ohat._coerce_map_from_(R)
            True
            sage: Ohat._coerce_map_from_(K)
            False
            sage: Ohat._coerce_map_from_("blah")
            False
        """
        if self.base().has_coerce_map_from(S):
            return True
        if isinstance(S, ProfiniteIntegers) and self.number_field().has_coerce_map_from(S.number_field()):
            return True
        from sage.rings.quotient_ring import is_QuotientRing
        if is_QuotientRing(S) and S.ambient() in [self.base(), self.number_field()]:
            return True
        return False

    def _an_element_(self):
        """
        Return a typical element of this ring

        EXAMPLES::

            sage: Zhat.an_element()
            2 mod 6
            sage: K.<a> = NumberField(x^3-2)
            sage: ProfiniteIntegers(K).an_element()
            8*a^2 mod (7*a^2 - 1)
        """
        g = self.base().gens()[-1]
        return self.element_class(self, g+1, 7*g-1)

    def some_elements(self):
        """
        Return some elements of this ring
        
        EXAMPLES::

            sage: Zhat.some_elements()
            [0, 2 mod 6, 1, 66 mod 79, 20, 10 mod 100]
        """
        return [self.zero(), self.an_element(), self.one(), self(66, 79), self(20), self(10, 100)]

    def random_element(self):
        """
        Return a random element of this ring

        EXAMPLES::

            sage: Ohat.random_element() # random
            0 mod (1)
            sage: Ohat.random_element() # random
            2/3*a^2 mod (a - 1)
            sage: Ohat.random_element() # random
            1/3*a^2 + a mod (-1/3*a^2 + a - 1)
            sage: Ohat.random_element() # random
            4/3*a^2 + a mod (-10773*a^2 - 75978*a - 229392)
            sage: Ohat.random_element() # random
            0 mod (1)
            sage: Ohat.random_element() # random
            a^2 + 3*a - 1 mod (104*a^2 + 52*a + 364)
        """
        O = self.base()

        value = O.random_element()
        
        if O is ZZ:
            modulus = ZZ.random_element()
            while modulus == -1:
                modulus = ZZ.random_element()
        else:
            from sage.misc.prandom import choice
            from sage.arith.misc import factor
            n_primes = ZZ.random_element()
            while n_primes < 0:
                n_primes = ZZ.random_element()

            modulus = O.ideal(1)
            for i in range(n_primes):
                I = O.ideal(O.random_element())
                while I.is_one() or I.is_zero():
                    I = O.ideal(O.random_element())

                # We pick individual primes, so that we can end up with
                # non-principal ideals as modulus.
                p, e = choice(factor(I))  
                modulus *= p**e

        multiplier = ZZ.random_element()
        while multiplier <= 0:
            multiplier = ZZ.random_element()

        modulus *= multiplier
            
        return self.element_class(self, value, modulus)

    def gens(self):
        """
        Return a tuple of generators of this ring, which is ``(1,)``

        EXAMPLES::

            sage: Zhat.gens()
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

            sage: Zhat.gen()
            1
            sage: Zhat.gen(0)
            1

        TESTS::

            sage: Zhat.gen(1)
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

            sage: Zhat.ngens()
            1
        """
        return ZZ(1)

    def is_finite(self):
        """
        Return ``False``, indicating that this ring is not finite

        EXAMPLES::

            sage: Zhat.is_finite()
            False
        """
        return False

    def is_exact(self):
        """
        Return ``False``, indicating that doing arithmetic can lead to precision
        loss

        EXAMPLES::

            sage: Zhat.is_exact()
            False
        """
        return False

    def is_integral_domain(self, proof=None):
        r"""
        Return ``False``, indicating that this ring is not an integral domain

        Note that this ring of profinite integers is (canonically isomorphic to)
        the product of all completions of `K` at the finite primes of `K`, where
        `K` denotes our base number field.
        As a product of non-trivial rings, this ring is clearly not an integral
        domain.

        EXAMPLES::

            sage: Zhat = ProfiniteIntegers()
            sage: Zhat.is_integral_domain()
            False
        """
        return False

    def is_field(self, proof=True):
        """
        Return ``False``, indicating that this ring is not a field

        Note that this ring of profinite integers is (canonically isomorphic to)
        the product of all completions of `K` at the finite primes of `K`, where
        `K` denotes our base number field.
        As a product of non-trivial rings, this ring is clearly not a field.

        EXAMPLES::

            sage: Zhat.is_field()
            False
        """
        return False

    def is_noetherian(self):
        """
        Return ``False``, indicating that this ring is not Noetherian

        Note that this ring of profinite integers is (canonically isomorphic to)
        the product of all completions of `K` at the finite primes of `K`, where
        `K` denotes our base number field.
        As a product of infinitely many non-trivial rings, this ring is clearly
        not Noetherian.

        EXAMPLES::

            sage: Zhat.is_noetherian()
            False
        """
        return False

    def krull_dimension(self):
        """
        Return ``1``, indiciting that the Krull dimension of this ring is one

        Note that this ring of profinite integers is (canonically isomorphic to)
        the product of all completions of `K` at the finite primes of `K`, where
        `K` denotes our base number field. Each such completion has Krull
        dimension one and therefore this ring has Krull dimension one as well.

        EXAMPLES::

            sage: Zhat.krull_dimension()
            1
        """
        return ZZ(1)

    def characteristic(self):
        """
        Return the characteristic of this ring, which is zero

        EXAMPLES::

            sage: Zhat = ProfiniteIntegers()
            sage: Zhat.characteristic()
            0
        """
        return ZZ(0)

    def order(self):
        """
        Return ``Infinity``, indicating that this ring has infinitely many
        elements

        EXAMPLES::

            sage: Zhat.order()
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

            sage: Zhat.epsilon()
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

            sage: F, P = Zhat.construction(); F, P
            (ProfiniteCompletionFunctor, Integer Ring)
            sage: F(P) is Zhat
            True

        ::

            sage: K.<a> = NumberField(x^3-2)
            sage: Ohat = ProfiniteIntegers(K)
            sage: F, P = Ohat.construction(); F, P
            (ProfiniteCompletionFunctor,
             Maximal Order in Number Field in a with defining polynomial x^3 - 2)
            sage: F(P) is Ohat
            True
        """
        return ProfiniteCompletionFunctor(), self.base()

    def number_field(self):
        """
        Return the base number field of this ring of profinite integers

        The base number field equals the fraction field of ``self.base()``.

        EXAMPLES::

            sage: Zhat = ProfiniteIntegers()
            sage: Zhat.number_field()
            Rational Field
            sage: K.<a> = NumberField(x^2+x-7)
            sage: Ohat = ProfiniteIntegers(K)
            sage: Ohat.number_field()
            Number Field in a with defining polynomial x^2 + x - 7
        """
        return self.base().fraction_field()


class ProfiniteCompletionFunctor(ConstructionFunctor):
    """
    The functor sending (the maximal order of) a number field to its profinite
    completion

    Due to this functor, we have the following functionality. ::

        sage: Zhat(3, 6) + 1/2
        7/2 mod 6

    ::

        sage: R.<X> = ZZ['X']
        sage: f = X^2 + 1
        sage: f + Zhat(5, 20)
        X^2 + 6 mod 20
    """
    # The rank below should be at least 6, since the FractionField functor has
    # rank 5 and we want the pushout of Zhat and QQ to be Qhat (the other way
    # around is undefined: we cannot create the fraction field of Zhat, which is
    # not an integral domain).
    # The rank below should also be at most 8, since the PolynomialFunctor
    # Poly[X] has rank 9 and we want the pushout of ZZ[X] and Zhat to be
    # Zhat[X] (as the profinite completion of ZZ[X] is not implemented).
    # Based on the possible values 6, 7 and 8 for the rank, we choose 7 without
    # any good reason over 6 and 8.
    rank = 7

    def __init__(self, args=None, kwds=None):
        """
        Create a ProfiniteCompletionFunctor

        EXAMPLES::

            sage: ProfiniteCompletionFunctor()
            ProfiniteCompletionFunctor
        """
        self.args = args or ()
        self.kwds = kwds or {}
        ConstructionFunctor.__init__(self, IntegralDomains(), Rings())

    def _apply_functor(self, R):
        """
        Apply this functor to ``R``

        INPUT:

        - ``R`` -- a number field or the maximal order of a number field

        OUPUT:

        If ``R`` is a number field, then we return the ring of profinite
        `R`-numbers.
        If ``R`` is the maximal order of a number field, then we return the ring
        of profinite `R`-integers.

        EXAMPLES::

            sage: F = ProfiniteCompletionFunctor()
            sage: F(ZZ)
            Ring of Profinite Integers of Rational Field
            sage: F(QQ)
            Ring of Profinite Numbers of Rational Field

        ::

            sage: K.<a> = NumberField(x^4+2)
            sage: F(K.maximal_order())
            Ring of Profinite Integers of Number Field in a with defining polynomial x^4 + 2
            sage: F(K)
            Ring of Profinite Numbers of Number Field in a with defining polynomial x^4 + 2
        """
        if R.is_field():
            from .profinite_number import ProfiniteNumbers
            return ProfiniteNumbers(R, *self.args, **self.kwds)
        return ProfiniteIntegers(R, *self.args, **self.kwds)


Zhat = ProfiniteIntegers(ZZ)

