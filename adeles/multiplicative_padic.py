r"""
Multiplicative `p`-adics over Number Fields

Let `K` be a number field and let `p` be a finite prime of `K`, i.e. a prime
ideal of its ring of integers. Write `K_p` for the field of `p`-adic numbers,
i.e. the completion of `K` at `p`. This file implementes the unit group `K_p^*`
of `K_p`, i.e. the group of multiplicative `p`-adics, in the class
:class:`MultiplicativePAdics`. ::

    sage: M = MultiplicativePAdics(QQ, 2); M
    Group of multiplicative 2-adics of Rational Field
    sage: M.prime()
    2
    sage: M.number_field()
    Rational Field

The corresponding elements are implemented as instances of
:class:`MultiplicativePAdic`. Such an instance `a` consists of a *center*
`x \in K^*` and a *precision* `n \in \ZZ_{\geq 0} \cup \{\infty\}`. The
represented subset of `a` is the subgroup `x \cdot U_p(n)` of `K_p^*`, where we
write `U_p(n)` for the `n`-th *multiplicative subgroup* at `p`:

.. MATH::

  U_p(n) = \begin{cases}
               O_p^*          & \text{ if } n = 0; \\
               1 + p^n O_p    & \text{ if } n \in \ZZ_{>0}; \\
               \{1\}          & \text{ if } n = \infty,
           \end{cases}

where `O_p` denotes the ring of `p`-adic integers, i.e. the valuation ring of
`K_p`. ::

    sage: a = M(3, 4); a
    3 * U(4)
    sage: a.center()
    3
    sage: a.precision()
    4
    sage: [n for n in range(-50, 50) if a.represents(n)]
    [-45, -29, -13, 3, 19, 35]

The following multiplicative `p`-adic only representes `1/7`::

    sage: b = M(1/7, oo); b
    1/7
    sage: b.precision()
    +Infinity

We keep the center relatively small using HNF-reduction::

    sage: c = M(-100, 1); c
    4 * U(1)
    sage: c.center()
    4

Multiplication and division are implemented::

    sage: a * b
    3/7 * U(4)
    sage: b / c
    1/28 * U(1)
    sage: c / a
    4/3 * U(1)

These operations satisfy that the represented subset of the product/quotient
equals the product/quotient of the represented subsets of the inputs.

Each multiplicative `p`-adic only represents elements of `K_p^*` of equal
`p`-adic valuation. ::

    sage: a.valuation()
    0
    sage: b.valuation()
    0
    sage: c.valuation()
    2

Let us demonstrate the same functionality as above over a non-trivial number
field.

We create a number field, take a prime `p` of it above `11` and create the
corresponding group of multiplicative `p`-adics::

    sage: K.<a> = NumberField(x^3+x+1)
    sage: p = K.prime_above(11)
    sage: M = MultiplicativePAdics(K, p); M
    Group of multiplicative (11, a^2 - 4)-adics of Number Field in a with defining polynomial x^3 + x + 1
    sage: M.prime()
    Fractional ideal (a - 2)
    sage: M.number_field()
    Number Field in a with defining polynomial x^3 + x + 1

We create some multiplicative `p`-adics and perform arithmetic. Note that the
centers are being "reduced" automatically (although for very small input
values the reduced center may actually be bigger). ::

    sage: b = M(a^2-4, 0); b
    55*a^2 * U(0)
    sage: b.valuation()
    1
    sage: c = M(97, 2); c
    -39*a^2 * U(2)
    sage: b * c
    11*a^2 * U(0)
    sage: b / c
    (-47*a^2 - 16/39) * U(0)

As the the open multiplicative subgroups `\{U_p(n) \mid n \in \ZZ_{\geq 0}\}`
form a basis of open neighborhoods of `1` in `K_p^*` and `K^*` lies dense in
`K_p^*`, we can approximate any `\alpha \in K_p^*` arbitrarily closely by a
multiplicative `p`-adic in the following sense: for every neighborhood `U` of
`\alpha` in `K_p^*`, there exists a multiplicative `p`-adic representing
`\alpha` and with represented subset contained in `U`.

.. SEEALSO::

    :mod:`~adeles.idele`

REFERENCES:

[Her2021] Mathé Hertogh, Computing with adèles and idèles, master's thesis,
Leiden University, 2021.

This implementation of multiplicative `p`-adics is based on [Her2021]. An
extensive exposition of properties, design choices and two applications can be
found there.

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

from sage.structure.element import MultiplicativeGroupElement
from sage.structure.unique_representation import UniqueRepresentation
from sage.groups.group import Group
from sage.categories.groups import Groups
from sage.rings.infinity import Infinity
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.arith.functions import lcm
from sage.sets.primes import Primes


def is_finite_prime(p, K):
    r"""
    Return whether or not ``p`` is a finite prime of the number field ``K``
    
    By a finite prime of `\QQ` we mean a prime number. For ``K`` a non-trivial
    number field, a finite prime of ``K`` is a prime ideal of the ring of
    integers of ``K``

    INPUT:

    - ``p`` -- any object; ``p`` will be checked to see if it is a finite prime
      of ``K``
    - ``K`` -- a number field

    EXAMPLES::

        sage: is_finite_prime(2, QQ)
        True
        sage: is_finite_prime(10, QQ)
        False
        sage: is_finite_prime(ZZ.ideal(5), QQ)
        False

    ::

        sage: K.<a> = NumberField(x^2+5)
        sage: is_finite_prime(2, K) # 2 ramifies
        False
        sage: is_finite_prime(K.ideal(2, a+1), K)
        True
        sage: is_finite_prime(11, K) # 11 is inert
        True

    TESTS::

        sage: is_finite_prime("blah", QQ)
        False
        sage: is_finite_prime([], K)
        False
    """
    if K is QQ:
        return p in Primes()
    else:
        return p in K.ideal_monoid() and K.ideal(p).is_prime()


class MultiplicativePAdic(MultiplicativeGroupElement):
    r"""
    Multiplicative `p`-adic, for `p` a finite prime of a number field

    REFERENCES:

    Section 3.5 of [Her2021].

    .. automethod:: _mul_
    
    .. automethod:: _div_
    
    .. automethod:: _richcmp_
    """

    def __init__(self, parent, center, prec):
        r"""
        Construct a multiplicative `p`-adic

        INPUT:

        - ``parent`` -- a ring of multiplicative `p`-adics
        - ``center`` -- a non-zero element of the base number field; the center
        - ``prec`` -- a non-negative integer or ``Infinity``; the precision

        OUTPUT:
        A multiplicative `p`-adic representing the subset `x \cdot U_p^n` of
        `K_p^*`, where

        - `K` denotes our base number field;
        - `K_p^*` denotes the unit group of the completion of `K` at `p`, i.e.
          the multiplicative group of the field of `p`-adic numbers;
        - `x \in K^*` denotes the ``center``;
        - `n \in \ZZ_{\geq 0} \cup \{\infty\}` denotes ``prec``;
        - `U_p(n)` denotes the `n`-th multiplicative subgroup at `p`, i.e. the
          subgroup

          .. MATH::

              U_p(n) = \begin{cases}
                           O_p^*          & \text{ if } n = 0; \\
                           1 + p^n O_p    & \text{ if } n \in \ZZ_{>0}; \\
                           \{1\}          & \text{ if } n = \infty,
                       \end{cases}

          with `O_p` denoting the ring of `p`-adic integers.

        EXAMPLES::

            sage: M = MultiplicativePAdics(QQ, 3)
            sage: MultiplicativePAdic(M, 1/2, 10)
            1/2 * U(10)
            sage: MultiplicativePAdic(M, 6, oo)
            6

        ::

            sage: K.<a> = NumberField(x^3-2)
            sage: M = MultiplicativePAdics(K, K.ideal(a))
            sage: MultiplicativePAdic(M, a^2+1, 40)
            (a^2 + 1) * U(40)
            sage: MultiplicativePAdic(M, a-1, oo)
            a - 1
        """
        MultiplicativeGroupElement.__init__(self, parent)
        
        K = parent.number_field()
        if center not in K or center == 0:
            raise ValueError("center must lie in K^*")
        self._center = K(center)

        if prec is Infinity:
            self._prec = prec
        elif prec in ZZ and prec >= 0:
            self._prec = ZZ(prec)
            self._reduce()
        else:
            raise ValueError("prec must be a non-negative integer or Infinity")

    def _repr_(self):
        """
        Return a string representation of this multiplicative `p`-adic

        EXAMPLE::

            sage: K.<a> = NumberField(x^7-10)
            sage: M = MultiplicativePAdics(K, K.prime_above(7))
            sage: M(a^3, 8)
            a^3 * U(8)
            sage: M(a^2+a+1, 10)
            (a^2 + a + 1) * U(10)
            sage: M(-a^6 + 1/2, oo)
            -a^6 + 1/2
        """
        if self.prec() is Infinity:
            return str(self.center())

        need_parentheses = 1 < len([c for c in self.center().list() if not c.is_zero()])
        if need_parentheses:
            return "({}) * U({})".format(self.center(), self.prec())
        else:
            return "{} * U({})".format(self.center(), self.prec())

    def _mul_(self, other):
        """
        Return the product of this multiplicative `p`-adic and ``other``

        The product of two multiplicative `p`-adics `a` and `b` is a
        multiplicative `p`-adic with represented subset equal to the product of
        the represented subsets of `a` and `b`.

        EXAMPLES::

            sage: M = MultiplicativePAdics(QQ, 5)
            sage: M(1/2, 3) * M(1/3, 10)
            1/6 * U(3)
            sage: M(5/9, 11) * M(3)
            5/3 * U(11)

        ::

            sage: K.<a> = NumberField(x^2+1)
            sage: M = MultiplicativePAdics(K, 3)
            sage: M(1/3, 4) * M(a, 8)
            1/3*a * U(4)
            sage: M(1, 0) * M(-1, 3)
            -1 * U(0)

        REFERENCES:

        Section Arithmetic in Section 3.5 of [Her2021].
        """
        center = self.center() * other.center()
        prec = min(self.prec(), other.prec())
        return self.__class__(self.parent(), center, prec)

    def _div_(self, other):
        """
        Return the quotient of this multiplicative `p`-adic by ``other``

        This quotient is defined to be ``self * other.inverse()``.

        EXAMPLES::

            sage: M = MultiplicativePAdics(QQ, 2)
            sage: M(2, 5) / M(4, 7)
            1/2 * U(5)
            sage: M(1)/M(4,7)
            1/4 * U(7)
        
        ::

            sage: K.<a> = NumberField(x^7-10)
            sage: M = MultiplicativePAdics(K, K.prime_above(23))
            sage: M(2*a) / M(a)
            2
            sage: M(a^5, 2) / M(a^4, 3)
            49/976 * U(2)

        REFERENCES:

        Section Arithmetic in Section 3.5 of [Her2021].
        """
        return self * other.inverse()

    def _richcmp_(self, other, op):
        """
        Return the result of operator ``op`` applied to ``self`` and ``other``

        Only equality and inequality are implented.
    

        Two representations of multipicative `p`-adics are considered equal if
        their represented subsets have non-empty intersection. This is
        equivalent to one of them being included in the other.

        Inequality is defined as not being equal, i.e. having disjoint
        represented subsets.

        EXAMPLES::

            sage: M = MultiplicativePAdics(QQ, 2)
            sage: M(17, 10) == M(1, 4)
            True
            sage: M(17, 10) == M(1, 5)
            False
            sage: M(17, 10) != M(1, 5)
            True

        ::

            sage: K.<a> = NumberField(x^2-x-1)
            sage: M = MultiplicativePAdics(K, 7)
            sage: M(1/2, 2) == M(99/2, 3)
            True
            sage: M(1/2, 3) == M(99/2, 3)
            False
            sage: M(a, 0) == 3*a+1
            True

        REFERENCES:

        Section 5.4 of [Her2021].
        """
        from sage.structure.richcmp import op_EQ, op_NE
        if op == op_EQ:
            if self.prec() > other.prec():
                self, other = other, self
            return self.represents(other.center())
        if op == op_NE:
            return not (self == other)
        raise NotImplementedError()

    def _reduce(self):
        r"""
        HNF-reduce our center modulo `p^{\max(1, n)+\ord_p(x)}` where `n`
        denotes our precision and `x` denotes our current center.

        EXAMPLES::

            sage: M = MultiplicativePAdics(QQ, 5)
            sage: u = M(-10, 2); u
            115 * U(2)
            sage: u.represents(-10)
            True
            sage: u.center()
            115

        ::

            sage: K.<a> = NumberField(x^3-2)
            sage: M = MultiplicativePAdics(K, K.ideal(a))
            sage: u = M(a^10, 3); u
            8*a * U(3)
            sage: u.represents(a^10)
            True
            sage: u.center()
            8*a
        """
        P = self.parent().prime()
        e = max(1, self.prec()) + self.center().valuation(P)
        I = P**e
        if self.parent().number_field() is QQ:
            den = lcm(self.center().denominator(), I.denominator())
            num_new_center = (den*self.center()) % (den*I)
        else:
            den = lcm(self.center().denominator(), I.integral_split()[1])
            num_new_center = (den*I).reduce(den*self.center())
        new_center = num_new_center / den
        self._center = new_center

    def inverse(self):
        """
        Return the product of this multiplicative `p`-adic and ``other``

        The inverse of this multiplicative `p`-adic is the multiplicative
        `p`-adic with the same precision and whose center is inversed.

        EXAMPLES::

            sage: M = MultiplicativePAdics(QQ, 7)
            sage: M(1/2, 5).inverse()
            2 * U(5)
            sage: M(2/11, oo).inverse()
            11/2

        ::

            sage: K.<a> = NumberField(x^3+10)
            sage: M = MultiplicativePAdics(K, K.ideal(11, a-1))
            sage: M(a^2, 0).inverse()
            -1/10*a * U(0)
            sage: M(-a/3, oo).inverse()
            3/10*a^2

        REFERENCES:

        Section Arithmetic in Section 3.5 of [Her2021].
        """
        return self.__class__(self.parent(), 1/self.center(), self.prec())

    def center(self):
        """
        Return the center of this multiplicative `p`-adic

        EXAMPLES::

            sage: M = MultiplicativePAdics(QQ, 97)
            sage: M(-1, 1).center()
            96
        """
        return self._center

    def prec(self):
        """
        Return the precision of this multiplicative `p`-adic

        EXAMPLES::

            sage: M = MultiplicativePAdics(QQ, 97)
            sage: M(79, 0).precision()
            0
            sage: M(79, oo).precision()
            +Infinity
        """
        return self._prec

    def precision(self):
        """
        Return the precision of this multiplicative `p`-adic

        EXAMPLES::

            sage: M = MultiplicativePAdics(QQ, 97)
            sage: M(79, 20).precision()
            20
        """
        return self._prec

    def valuation(self):
        """
        Return the valuation at `p` of this multipicative `p`-adic

        EXAMPLES::

            sage: M = MultiplicativePAdics(QQ, 5)
            sage: M(1/25, 10).valuation()
            -2

        ::

            sage: K.<a> = NumberField(x^2+1)
            sage: M = MultiplicativePAdics(K, a+1)
            sage: M(4, 100).valuation()
            4
        """
        P = self.parent().prime()
        return self.center().valuation(P)

    def represents(self, element):
        r"""
        Return whether or not this multipicative `p`-adic represents ``element``
        
        This multiplicative `p`-adic represents all elements in its *represented
        subset*, which is the subset subset `x \cdot U_p^n` of `K_p^*`, where

        - `K` denotes our base number field;
        - `K_p^*` denotes the unit group of the completion of `K` at `p`, i.e.
          the multiplicative group of the field of `p`-adic numbers;
        - `x \in K^*` denotes ``self.center()``;
        - `n \in \ZZ_{\geq 0} \cup \{\infty\}` denotes ``self.precision()``;
        - `U_p(n)` denotes the `n`-th multiplicative subgroup at `p`, i.e. the
          subgroup

          .. MATH::

              U_p(n) = \begin{cases}
                           O_p^*          & \text{ if } n = 0; \\
                           1 + p^n O_p    & \text{ if } n \in \ZZ_{>0}; \\
                           \{1\}          & \text{ if } n = \infty,
                       \end{cases}

          with `O_p` denoting the ring of `p`-adic integers.

        EXAMPLES::

            sage: M = MultiplicativePAdics(QQ, 3)
            sage: a = M(1, 3)
            sage: X = sorted(set([m/n for m in srange(-18, 18) for n in srange(1, 18)]))
            sage: [x for x in X if a.represents(x)]
            [-17/10, -16/11, -14/13, -13/14, -11/16, -10/17, 1]
            sage: -16/11 - 1 # does this really have valuation at least 3 at 3? (yes)
            -27/11
            sage: b = M(2/3, 3)
            sage: [x for x in X if b.represents(x)]
            [-17/15, 2/3]
            sage: (-17/15)/(2/3)-1 # this should have valuation at least 3 at 3
            -27/10
            sage: M(9, 0).represents(9/100)
            True

        ::

            sage: K.<a> = NumberField(x^2+1)
            sage: M = MultiplicativePAdics(K, a+1)
            sage: b = M(a, 6)
            sage: Y = sorted(set([m/n for m in srange(-7, 7) for n in srange(1, 7)]))
            sage: R = [y*a+z for y in Y for z in Y if b.represents(y*a+z)]; R
            [-7*a, -5/3*a, -3/5*a, a]
            sage: all([(r/a-1).valuation(a+1) >= 6 for r in R])
            True
        """
        if self.prec() is Infinity:
            return element == self.center()

        P = self.parent().prime()
        if self.prec() == 0:
            return element.valuation(P) == self.valuation()

        # self.prec() is a positive integer, so our represented subset is
        # self.center() * (1 + P^self.prec() O_P)
        # for O_P the ring of P-adic integers.
        return (element / self.center() - 1).valuation(P) >= self.prec()


class MultiplicativePAdics(UniqueRepresentation, Group):
    """
    Group of Multiplicative `p`-adics, for `p` a finite prime of a number field

    REFERENCES:

    Section 3.5 of [Her2021].

    .. automethod:: _element_constructor_
    """

    Element = MultiplicativePAdic

    def __classcall__(cls, K, p):
        """
        Construct the group of multiplicative `p`-adics

        INPUT:

        - ``K`` -- a number field
        - ``p`` -- a finite prime of ``K``, as described by
          :func:`is_finite_prime`

        OUPUT:

        The multiplicative group of the field of ``p``-adic numbers, i.e. the
        unit group of the completion of ``K`` at ``p``.

        EXAMPLES::

            sage: MultiplicativePAdics(QQ, 5)
            Group of multiplicative 5-adics of Rational Field

        This method ensures the following::

            sage: K.<a> = NumberField(x^2-x-1)
            sage: MultiplicativePAdics(K, K.ideal(2)) is MultiplicativePAdics(K, 2)
            True
        
        TESTS::

            sage: MultiplicativePAdics("not a number field", "not a prime")
            Traceback (most recent call last):
            ...
            TypeError: K should be a number field
            sage: MultiplicativePAdics(QQ, "not a prime")
            Traceback (most recent call last):
            ...
            ValueError: p should be a finite prime of Rational Field
        """
        from sage.rings.number_field.number_field import is_NumberField
        if not is_NumberField(K):
            raise TypeError("K should be a number field")
        if not is_finite_prime(p, K):
            raise ValueError("p should be a finite prime of {}".format(K))
        p = ZZ(p) if K is QQ else K.ideal(p)
        return super(MultiplicativePAdics, cls).__classcall__(cls, K, p)

    def __init__(self, K, p):
        """
        Construct the group of multiplicative `p`-adics

        INPUT:

        - ``K`` -- a number field
        - ``p`` -- a finite prime of ``K``, as described by
          :func:`is_finite_prime`

        OUPUT:

        The multiplicative group of the field of ``p``-adic numbers, i.e. the
        unit group of the completion of ``K`` at ``p``.

        EXAMPLES::

            sage: MultiplicativePAdics(QQ, 3)
            Group of multiplicative 3-adics of Rational Field

        ::

            sage: K.<a> = NumberField(x^9+2)
            sage: MultiplicativePAdics(K, K.ideal(3, a-1))
            Group of multiplicative (3, a - 1)-adics of Number Field in a with defining polynomial x^9 + 2
        """
        # Note that input vericifaction is done in __classcall__().
        self._number_field = K
        self._prime = p
        Group.__init__(self, category=Groups().Commutative())

    def _repr_(self):
        """
        Return a string representation of this group

        EXAMPLES::

            sage: K.<a> = NumberField(x^2-3)
            sage: MultiplicativePAdics(K, 5)
            Group of multiplicative (5)-adics of Number Field in a with defining polynomial x^2 - 3
        """
        K, P = self.number_field(), self.prime_name()
        return "Group of multiplicative {}-adics of {}".format(P, K)

    def _latex_(self):
        r"""
        Return a latex-formatted string representation of this group

        EXAMPLES::

            sage: latex(MultiplicativePAdics(QQ, 5))
             { \Bold{Q} }_{ 5 }^*

        ::

            sage: K.<a> = NumberField(x^2+5)
            sage: latex(MultiplicativePAdics(K, K.prime_above(3)))
             { \Bold{Q}[a]/(a^{2} + 5) }_{ \left(3, a + 1\right) }^*
        """
        from sage.misc.latex import latex
        return " { " + latex(self.number_field()) + " }_{ " + latex(self.prime()) + " }^* "

    def _element_constructor_(self, center, prec=Infinity):
        r"""
        Construct a multiplicative `p`-adic

        INPUT:

        - ``center`` -- non-zero element of our base number field; the center
        - ``prec`` -- non-negative integer of ``Infinity``; the precision

        EXAMPLES::

            sage: M = MultiplicativePAdics(QQ, 2)
            sage: M(3/4, 10)
            3/4 * U(10)
            sage: M(79)
            79

        ::

            sage: K.<a> = NumberField(x^2+21)
            sage: M = MultiplicativePAdics(K, K.ideal(5, a+3))
            sage: M(-a, 0)
            -a * U(0)
            sage: M(a+1, 5)
            -661*a * U(5)
            sage: M(a/3)
            1/3*a
        """
        return self.element_class(self, center, prec)

    def _coerce_map_from_(self, domain):
        """
        Return whether or not ``domain`` coerces into this group

        At the moment nothing coerces into this group, except the group itself.

        EXAMPLES::

            sage: M = MultiplicativePAdics(QQ, 5)
            sage: M._coerce_map_from_(M)
            True
            sage: M._coerce_map_from_(QQ)
            False
        """
        # We may not define the conversion from our base field `K` to this group
        # to be a coercion, since it is not defined on the zero element of `K`.
        # if self.number_field().has_coerce_map_from(domain):
        #     return True
        return domain is self

    def _an_element_(self):
        """
        Return a typical element of this group

        EXAMPLES::

            sage: MultiplicativePAdics(QQ, 5).an_element()
            2/3 * U(6)

        ::

            sage: K.<a> = NumberField(x^4+2)
            sage: MultiplicativePAdics(K, K.prime_above(7)).an_element()
            (-39734*a^3 + 21488*a^2 - 1/3*a) * U(6)
        """
        center = 2 * self.number_field().gen(0) / 3
        return self.element_class(self, center, 6)

    def some_elements(self):
        """
        Return some elements of this group
        
        EXAMPLES::

            sage: MultiplicativePAdics(QQ, 5).some_elements()
            [2/3 * U(6), 1, 1 * U(0), 2, 121/5 * U(3)]

        ::

            sage: K.<a> = CyclotomicField(17)
            sage: MultiplicativePAdics(K, K.prime_above(17)).some_elements()
            [(5*a^15 - 6*a^14 - 4*a^13 - 8*a^12 + a^11 - 4*a^10 - 1/3*a) * U(6),
             1,
             a^15 * U(0),
             a + 1,
             (-6*a^15 + 8*a^14 - 3*a^13 + 1/5*a) * U(3)]
        """
        g = self.number_field().gen(0)
        return [self.an_element(), self.one(), self(g, 0), self(g+1), self(g/5-1, 3)]

    def random_element(self):
        """
        Return a random element of this group

        EXAMPLES::

            sage: [M.random_element() for n in range(8)] # random
            [59048 * U(10),
             5/2 * U(0),
             8 * U(15),
             2 * U(1),
             1/15 * U(1),
             7,
             1/2 * U(5),
             8 * U(2)]

        ::

            sage: K.<a> = NumberField(x^3-7)
            sage: M = MultiplicativePAdics(K, K.prime_above(2))
            sage: [M.random_element() for n in range(8)] # random
            [a^2 * U(1),
             a^2 + 1/13*a - 3,
             -15/44*a^2 * U(2),
             -1/40*a^2,
             1/2*a * U(0),
             a^2 * U(1),
             -1/8*a * U(3),
             1/2*a - 1]
        """
        center = self.number_field().random_element()
        while center == 0:
            center = self.number_field().random_element()

        prec = ZZ.random_element()
        if prec < 0:
            prec = ZZ.random_element()
            if prec < 0:
                # This happens approximately 16% of the time.
                prec = Infinity
            
        return self.element_class(self, center, prec)

    def gens(self):
        """
        Return a tuple of generators of this group, which is ``(1,)``

        EXAMPLES::

            sage: MultiplicativePAdics(QQ, 3).gens()
            (1,)
        """
        return (self(1), )

    def gen(self, n=0):
        """
        Return the ``n``-th generator of this group

        As this group has only one generator, we only accept ``n == 0``.

        INPUT:

        - ``n`` -- the index of the generator to return (default: ``0``); must
          be zero

        EXAMPLES::

            sage: MultiplicativePAdics(QQ, 7).gen()
            1
            sage: MultiplicativePAdics(QQ, 7).gen(0)
            1

        TESTS::

            sage: MultiplicativePAdics(QQ, 7).gen(1)
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
        Return the number of generators of this group, which is `1`

        EXAMPLES::

            sage: MultiplicativePAdics(QQ, 97).ngens()
            1
        """
        return ZZ(1)

    def is_abelian(self):
        """
        Return ``True``, indicating that this group is abelian

        EXAMPLES::

            sage: MultiplicativePAdics(QQ, 2).is_abelian()
            True
        """
        return True

    def is_commutative(self):
        """
        Return ``True``, indicating that this group is commutative

        EXAMPLES::

            sage: MultiplicativePAdics(QQ, 5).is_commutative()
            True
        """
        return True

    def is_exact(self):
        """
        Return ``False``, indicating that doing arithmetic can lead to precision
        loss

        EXAMPLES::

            sage: MultiplicativePAdics(QQ, 2).is_exact()
            False
        """
        return False

    def is_finite(self):
        """
        Return ``False``, indicating that this group is not finite

        EXAMPLES::

            sage: MultiplicativePAdics(QQ, 17).is_finite()
            False
        """
        return False

    def order(self):
        """
        Return the order of this group, which is ``Infinity``

        The order of a group is its number of elements.

        EXAMPLES::

            sage: MultiplicativePAdics(QQ, 5).order()
            +Infinity
        """
        from sage.rings.infinity import Infinity
        return Infinity

    def number_field(self):
        """
        Return the base number field of this group of multiplicative `p`-adics

        EXAMPLES::

            sage: K.<a> = NumberField(x^4+5)
            sage: M = MultiplicativePAdics(K, K.ideal(7, a+2))
            sage: M.number_field()
            Number Field in a with defining polynomial x^4 + 5
        """
        return self._number_field

    def prime(self):
        r"""
        Return the prime `p` of this group of multiplicative `p`-adics

        This group is the unit group of the completion of a number field `K` at
        a finite prime `p` of `K`. We return this `p`.
        
        If our base number field is `\QQ`, then `p` is a prime number in `\ZZ`.
        Otherwise `p` is a prime ideal of the ring of integers of `K`.

        EXAMPLES::

            sage: MultiplicativePAdics(QQ, 17).prime()
            17
            sage: K.<a> = NumberField(x^2-x-1)
            sage: MultiplicativePAdics(K, 3).prime()
            Fractional ideal (3)
        """
        return self._prime

    def prime_name(self):
        """
        Return a string representing our prime, as returned by :meth:`prime`

        EXAMPLES::

            sage: MultiplicativePAdics(QQ, 5).prime_name()
            '5'
            sage: K.<a> = NumberField(x^2-71)
            sage: MultiplicativePAdics(K, 3).prime_name()
            '(3)'
            sage: MultiplicativePAdics(K, -a+8).prime_name()
            '(7, a + 6)'
        """
        if self.number_field() is QQ:
            return str(self.prime())
        else:
            gens = self.prime().gens_two()
            if gens[1] == 0:
                return "({})".format(gens[0])
            else:
                return str(gens)


def MulPAdic(p, data):
    """
    Create a multiplicative ``p``-adic based on ``data``

    INPUT:

    - ``p`` -- a finite prime of some number field, as described by
      :func:`is_finite_prime`; determines the prime `p` of the multiplicative
      `p`-adic we create
    - ``data`` -- data from which to construct the multiplicative ``p``-adic.
      This can be:

      - a single value from which a multiplicative ``p``-adic can be created,
        such as a non-zero element of the base number field or a multiplicative
        ``p``-adic itself;
      - an iterable of length 2, where the first value denotes a center and the
        second value denotes a precision, i.e. ``data[0]`` is a non-zero element
        of the base number field and ``data[1]`` is a non-negative integer or
        ``Infinity``.
    
    EXAMPLES::

        sage: a = MulPAdic(2, (2/3, 5))
        sage: a, a.parent()
        (2/3 * U(5), Group of multiplicative 2-adics of Rational Field)
        sage: b = MulPAdic(7, 5/7)
        sage: b, b.parent()
        (5/7, Group of multiplicative 7-adics of Rational Field)
        sage: c = MulPAdic(7, b)
        sage: c, c.parent()
        (5/7, Group of multiplicative 7-adics of Rational Field)

    ::

        sage: K.<a> = NumberField(x^2-5)
        sage: b = MulPAdic(a, (a+1, 3))
        sage: b, b.parent()
        ((a + 1) * U(3),
         Group of multiplicative (5, 1/2*a + 5/2)-adics of Number Field in a with defining polynomial x^2 - 5)
        sage: c = MulPAdic(K.ideal(3), [a/10, 0])
        sage: c, c.parent()
        (1/10*a * U(0),
         Group of multiplicative (3)-adics of Number Field in a with defining polynomial x^2 - 5)

    TESTS:

        sage: MulPAdic(4, 1)
        Traceback (most recent call last):
        ...
        ValueError: p must be a finite prime of a number field
        sage: MulPAdic(5, oo)
        Traceback (most recent call last):
        ...
        TypeError: can't construct a multiplicative p-adic from +Infinity
        sage: MulPAdic(3, ["blah"])
        Traceback (most recent call last):
        ...
        TypeError: can't construct a multiplicative p-adic from ['blah']
    """
    from sage.rings.number_field.number_field_ideal import is_NumberFieldIdeal
    from sage.rings.number_field.number_field_element import is_NumberFieldElement
    if p in Primes():
        K = QQ
    elif is_NumberFieldIdeal(p):
        K = p.number_field()
    elif is_NumberFieldElement(p):
        K = p.parent()
    else:
        raise ValueError("p must be a finite prime of a number field")
    parent = MultiplicativePAdics(K, p)

    if data in parent or data in K:
        return parent(data)

    try:
        center, prec = data
    except:
        raise TypeError("can't construct a multiplicative p-adic from {}".format(data))
    return parent(center, prec)

