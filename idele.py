r"""
Idele Groups of Number Fields

This implementation of ideles is similar to how the Real Interval Field ``RIF``
implements the field of real numbers. Storing a single real number takes in
general an infinite amount of memory, which we don't have. Hence ``RIF`` chooses
to only represent intervals of the real line, instead of actual real numbers.
Working with very small intervals corresponds to working with high precision.

Similarly, we do not represent individual ideles, but rather open subsets of
the idele group. A user can work with arbitrary high precision "ideles" by using
very small open subsets in its calculations.

Let `K` be a number field. For a place `p` of `K`, write `K_p` for the
completion of `K` at `p`. Write `Z_p` for the valuation ring of `K_p`. Write
`Z_p^*` for the unit group of `Z_p` and `M_p` for its maximal ideal.

The idele group `J_K` of `K` is is the restricted product
`\prod_p'(K_p^*, Z_p^*)` of `K_p^*`'s with respect to `Z_p`'s, where `p` ranges
over all places (finite and infinite) of `K`.
Concretely, that means

.. MATH::

    J_K = \left\{ (x_p)_p \in \prod_p K_p^* : x_p \in Z_p^*
                  \text{ for all but finitely many } p \right\}

This is a multiplicative group with component-wise multiplication.

It is even a topological group, with basic open sets being translates of open
subgroups of the form

.. MATH::

    \prod_{p \in S} U_p \times \prod_{p \not\in S} Z_p^*

where `S` is a finite set of places of `K` and each `U_p` is of the form
`1+M_p^i` with `i \in \NN`. These basic open sets our what our
``Idele``'s  will represent.

.. TODO::

    Check all idele-code for edge cases:
        - finite empty
        - exact vs non-exact
    Also check idele-code in other files (i.e. Adeles._from_idele())
"""
import sys # TODO erase these two lines
sys.path.append('/home/mathe/adeles/src')

from sage.structure.element import MultiplicativeGroupElement
from sage.structure.unique_representation import UniqueRepresentation
from sage.groups.group import Group
from sage.arith.misc import factor
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.real_mpfi import RIF
from sage.rings.complex_interval_field import ComplexIntervalField
from sage.rings.infinity import Infinity
from sage.sets.primes import Primes

from completion import completions, infinite_completions
from ray_class_group import Modulus, ray_class_group

CIF = ComplexIntervalField()
oo = Infinity


class Idele(MultiplicativeGroupElement):
    r"""
    An idele: an element of an :class:`Ideles` of a number field

    .. automethod:: __init__
    """

    def __init__(self, parent, exact, infinite, finite):
        r"""
        Construct an ``Idele``

        We denote the number field to which this idele belongs by ``K``.

        INPUT:

        - ``parent`` -- instance of :class:`Ideles`; parent idele group
        - ``exact`` -- an element of ``K`` or ``None``; determines the exact
          value of this idele at all unstored finite primes
        - ``infinite`` -- ``None`` or a list of elements of ``RIF`` and ``CIF``;
          the list should have as many entries as ``K`` has infinite places. The
          ``i``th entry specifies the  value of this idele at the infinite place
          corresponding to ``K.places()[i]``. The entry should lie in ``RIF``
          or ``CIF``, depending on the place being real or complex.
        - ``finite`` -- a ``dict`` or ``None``; ``None`` is the same as passing
          an empty ``dict``. The ``dict`` should contain (``key``, ``value``)
          pairs such that:
            * ``key`` is a prime ``p`` of ``K``; if ``K`` is `\QQ`, this is a
              prime number, else it's a prime fractional ideal of ``K``.
            * ``value`` should be a pair (``x``, ``i``) with ``x`` a non-zero
              element of ``K`` and ``i`` a non-negative integer.
          Such an entry means: at place ``p`` of ``K``, this idele is equal to
          (the image of) ``x`` (in `K_p`), up to elements in the open subgroup
          `1+M_p^i` of `K_p`.

        EXAMPLES:

        We begin with the easiest case: `K = \QQ`::

            sage: J = Ideles(QQ)
            sage: Idele(J, None, [3.14], {2: (5/3, 2), 5: (1, 9)})
            Idele with values:
              infinity_0:   3.1400000000000002?
              2:            5/3 * U(2)
              5:            1 * U(9)
            sage: Idele(J, None, [-1], {7: (1/10, 100)})
            Idele with values:
              infinity_0:   -1
              7:            1/10 * U(100)
            sage: Idele(J, 1/3, [1], {7: (1/10, 100)})
            Idele with values:
              infinity_0:   1
              7:            1/10 * U(100)
              elsewhere:    1/3

        Now let's take a non-trivial number field::

            sage: K.<a> = NumberField(x^4-17)
            sage: Jk = Ideles(K)
            sage: Idele(Jk, a/2, [3.14, -1, 1+I], {K.prime_above(3): (a^3+a+1, 10)})
            Idele with values:
              infinity_0:   3.1400000000000002?
              infinity_1:   -1
              infinity_2:   1 + 1*I
              (3, 1/2*a^2 - a - 1/2):       (a^3 + a + 1) * U(10)
              elsewhere:    1/2*a
        """
        self._exact = exact
        self._infinite = infinite
        self._finite = finite
        MultiplicativeGroupElement.__init__(self, parent)
        self._validate_exact()
        self._validate_infinite()
        self._validate_finite()

    def _validate_exact(self):
        """
        Check if the ``exact`` attribute is valid and adjust it if feasable

        Valid means: ``exact`` is None, or is a non-zero element of `K`.

        We throw an exception if exact is invalid.

        If ``exact`` can be coerces into `K`, but is not an actual element of
        `K`, we cast it into an actual element of `K`.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3-2)
            sage: J = Ideles(K)
            sage: u = J(1)
            sage: u._exact = None
            sage: u._validate_exact()
            sage: u._exact = QQ(3/2)
            sage: u._validate_exact()
            sage: u._exact.parent()
            Number Field in a with defining polynomial x^3 - 2
            sage: u._exact = []
            sage: u._validate_exact()
            Traceback (most recent call last):
            ...
            TypeError: exact should be an element of the number field
        """
        if self.exact() is None:
            return
        K = self.parent().number_field()
        if self.exact() not in K:
            raise TypeError("exact should be an element of the number field")
        self._exact = K(self.exact())
        if self.exact().is_zero():
            raise ValueError("exact must be a unit (i.e. non-zero)")

    def _validate_infinite(self):
        """
        Check that the ``infinite`` attribute is valid and adjust if faesable

        Valid means: ``infinite`` is a list of length ``len(K.places())`` with
        the `i`-th entry an element of ``RIF`` or ``CIF``, depending on the
        place being real or complex, or ``infinite`` is ``None``.

        If the `i`-th entry's parent is not equal to ``RIF``/``CIF``, but does
        lie in it, than we cast the element to an actual element of such an
        interval field.

        Upon invalid (and uncastable) ``infinite``, we throw an exception.

        EXAMPLES:

        We look at a number field with one real infinite place and one complex
        one::

            sage: K.<a> = NumberField(x^3 + x + 1)
            sage: Jk = Ideles(K)
            sage: u = Jk(1)
            sage: u._infinite = [1.2, I]
            sage: u._validate_infinite()

        Below the ``Integer 3`` will be cast into ``RIF``::

            sage: u._infinite = [ZZ(3), I]
            sage: u._validate_infinite()
            sage: u._infinite[0].parent()
            Real Interval Field with 53 bits of precision

        And here we cast ``None`` to a list of unknown intervals::

            sage: u._infinite = None
            sage: u._validate_infinite()
            sage: u._infinite
            [[-infinity .. +infinity],
             [-infinity .. +infinity] + [-infinity .. +infinity]*I]

        The exceptions that we may throw:

            sage: u._infinite = [3.14, I, I]
            sage: u._validate_infinite()
            Traceback (most recent call last):
            ...
            ValueError: infinite should have length 2
            sage: u._infinite = {CIF: I+1}
            sage: u._validate_infinite()
            Traceback (most recent call last):
            ...
            TypeError: infinite should be a list
            sage: u._infinite = [I, I]
            sage: u._validate_infinite()
            Traceback (most recent call last):
            ...
            TypeError: 0th infinite value (I) should lie in Real Interval Field with 53 bits of precision
        """
        K = self.parent().number_field()
        K_oo = infinite_completions(K)
        t = len(K_oo)
        if self.infinite() is None:
            r1, r2 = K.signature()
            RR = RIF(-oo, oo)
            CC = CIF(RR, RR)
            self._infinite = [RR for i in range(r1)] + [CC for i in range(r2)]
            return
        if not isinstance(self.infinite(), list):
            raise TypeError("infinite should be a list")
        if len(self.infinite()) != t:
            raise ValueError("infinite should have length {}".format(t))
        for i in range(t):
            val = self.infinite(i)
            if val not in K_oo[i][0]:
                raise TypeError("{}th infinite value ({}) should lie in {}".format(i, val, K_oo[i][0]))
            if val == 0:
                raise ValueError("{}th infinite value ({}) must be non-zero".format(i, val))
            self._infinite[i] = K_oo[i][0](val)

    def _validate_finite(self):
        r"""
        Check if the ``finite`` attribute is valid and adjust if feasable

        Valid means: ``None`` or a ``dict`` with keys primes of `K` and values
        pairs in `(K, \NN)`.

        We throw an exception upon invalid ``finite``.

        EXAMPLES:

        A few examples of the small adjustments we can do::

            sage: K.<a> = NumberField(x^4+x+1)
            sage: J = Ideles(K)
            sage: u = J(1)
            sage: u._finite = None
            sage: u._validate_finite()
            sage: u._finite
            {}
            sage: u._finite = {K.prime_above(2): (1, 1)}
            sage: u._validate_finite()
            sage: u._finite[K.prime_above(2)][0].parent()
            Number Field in a with defining polynomial x^4 + x + 1

        And some exceptions we may throw::

            sage: u._finite = "blah"
            sage: u._validate_finite()
            Traceback (most recent call last):
            ...
            TypeError: finite should be a dict
            sage: u._finite = {40: (1, 1)}
            sage: u._validate_finite()
            Traceback (most recent call last):
            ...
            TypeError: keys of finite should be prime ideals of the number field
            sage: u._finite = {K.prime_above(2): (1, -1)}
            sage: u._validate_finite()
            Traceback (most recent call last):
            ...
            TypeError: values at finite primes should be pairs in K x NN
        """
        if self.finite() is None:
            self._finite = {}
            return
        if not isinstance(self.finite(), dict):
            raise TypeError("finite should be a dict")
        K = self.parent().number_field()
        already_seen = set()
        new_finite = {}
        for q, val in self.finite().items():
            if q not in K.ideal_monoid() or not q.is_prime():
                if K is QQ:
                    raise TypeError("keys of finite should be prime numbers")
                else:
                    raise TypeError("keys of finite should be prime ideals of the number field")
            try:
                v, i = val
            except TypeError:
                raise TypeError("values at finite primes should be pairs in K x NN")
            if v not in K or i not in ZZ or i < 0:
                raise TypeError("values at finite primes should be pairs in K x NN")
            if v.is_zero():
                raise ValueError("values at finite primes must be units (i.e. non-zero)")
            if q in already_seen:
                raise TypeError("multiple values specified for finite prime {}".format(q))
            already_seen.add(q)

            if K is QQ:
                new_q = ZZ(q)
            else:
                new_q = K.ideal(q)
            if not new_q.is_prime():
                raise ValueError("{} is not prime as an ideal of {}".format(q, K))
            new_finite[new_q] = (K(v), ZZ(i))
        self._finite = new_finite

    def _repr_(self):
        """
        Return a representation of this idele.

        EXAMPLES::

            sage: J = Ideles(QQ)
            sage: u = J(1/2, None, {10000079: (7/9, 20)}); u
            Idele with values:
              infinity_0:  RR^*
              10000079:     7/9 * U(20)
              elsewhere:    1/2
            sage: J(None, [7.9], {10000079: (7/9, 20)})
            Idele with values:
              infinity_0:   7.9000000000000004?
              10000079:     7/9 * U(20)
        """
        K = self.parent().number_field()
        r1, r2 = K.signature()
        RR = RIF(-oo, oo)
        CC = CIF(RR, RR)
        rep = "Idele with values:"
        for i in range(len(K.places())):
            x = self.infinite(i)
            if x.endpoints() == RR.endpoints():
                rep += "\n  infinity_{}:\tRR^*".format(i)
            elif x.endpoints() == CC.endpoints():
                rep += "\n  infinity_{}:\tCC^*".format(i)
            else:
                rep += "\n  infinity_{}:\t{}".format(i, self.infinite(i))

        finite = list(self.finite().items())
        if K is QQ:
            # sort on the (rational) prime itself
            finite.sort(key=lambda x: x[0])
        else:
            # sort on rational prime lying below the stored prime:
            finite.sort(key=lambda x: x[0].gens_two()[0])

        for q, val in finite:
            x_q, i_q = val
            q_name = str(q) if K is QQ else str(q.gens_two())
            tab = "\t\t" if len(q_name) < 5 else "\t"
            if len([c for c in x_q.list() if not c.is_zero()]) > 1:
                x_q = "(" + str(x_q) + ")"
            rep += "\n  {}:{}{} * U({})".format(q_name, tab, x_q, i_q)
        if self._has_exact():
            rep += "\n  elsewhere:\t{}".format(self.exact())
        return rep

    def _equals(self, other):
        """
        Return whether or not ``self`` equals ``other``

        We define equality of ideles very loosely: if two ideles *could* be
        equal, that is if the open subsets they represent have non-empty
        intersection, than we declare them equal.

        EXAMPLES::

            sage: J = Ideles(QQ)
            sage: u = J(None, None, {2: (3, 7), 3: (1/2, 6)})
            sage: v = J(4, [-1.23], {3: (5, 10)})
            sage: w = J(None, None, {2: (3, 7), 3: (1/2, 6), 5: (1, 0)})
            sage: u._equals(u)
            True
            sage: (u/u)._equals(J(1))
            True
            sage: u._equals(v)
            False
            sage: u._equals(w)
            True
            sage: v._equals(w)
            False

        If an idele ``x`` represents a strict subset of another idele ``v``,
        they are equal::

            sage: x = J(4, [-1.23], {3: (5, 11)})
            sage: v._equals(x)
            True

        Also in this case the result is not equal, due to the non-exactness of
        ``u``::
        """
        K = self.parent().number_field()

        if self._has_exact() and other._has_exact():
            if self.exact() != other.exact():
                return False
        elif self._has_exact() or other._has_exact():
            ex, nex = (self, other) if self._has_exact() else (other, self)
            if not nex._contains(ex.exact()):
                return False

        t = len(self.infinite())
        for i in range(t):
            try:
                intersection = self.infinite(i).intersection(other.infinite(i))
                if intersection.is_zero():
                    return False
            except ValueError:
                # This indicates the intersection is empty
                return False

        for q, v in self.finite().items():
            if q in other.finite():
                w = other.finite(q)
                small, big = (self, other) if v[1] >= w[1] else (other, self)
                if not big._contains_at(small.finite(q)[0], q):
                    # the small open subset does not lie in the big open subset
                    # so they are disjoint
                    return False
            else: # other[q] not stored
                if other._has_exact():
                    if not self._contains_at(other.exact(), q):
                        return False
                else:
                    x = v[0]
                    # other[q] is Z_q*
                    # So x must lie in Z_q*, i.e. have valuation zero
                    if K is not QQ:
                        x = K.ideal(x)
                    if  x.valuation(q) != 0:
                        return False
        for q, v in other.finite().items():
            if q not in self.finite():
                if self._has_exact():
                    # self[q] is exact, while other[q] is an approximation.
                    return False
                else:
                    x = v[0]
                    # self[q] is Z_q*
                    # So x must lie in Z_q*, i.e. have valuation zero
                    if K is not QQ:
                        x = K.ideal(x)
                    if x.valuation(q) != 0:
                        return False

        # Passed all checks. We are equal!
        return True

    def _richcmp_(self, other, op):
        """
        Return the result of operator ``op`` applied to ``self`` and ``other``

        Only equality and inequality are implented.

        EXAMPLES::

            sage: K.<a> = NumberField(x^5-x+2)
            sage: J = Ideles(K)
            sage: u = J(a^4+1, None, {K.prime_above(11): (a^3-a, 20)})
            sage: v = J(None, None, {K.prime_above(97): (a, 20)})
            sage: v == u
            False
            sage: v != u
            True
            sage: u == u
            True
            sage: u != u
            False
        """
        from sage.structure.richcmp import op_EQ, op_NE
        if op == op_EQ:
            return self._equals(other)
        if op == op_NE:
            return not self._equals(other)
        raise NotImplementedError()

    def _has_exact(self):
        """
        Return whether or not our ``exact`` attribute is set

        EXAMPLES::

            sage: J = Ideles(QQ)
            sage: u = J(None, [3.1415])
            sage: u._has_exact()
            False
            sage: u._exact = 1/2
            sage: u._has_exact()
            True
        """
        return self.exact() is not None

    def exact(self):
        """
        Return the exact part of ``self`` as an element of our number field (or
        None)
        """
        return self._exact

    def infinite(self, index=None):
        """
        Return the infinite part of ``self`` as a list, or an entry of it

        INPUT:

        - ``index`` -- integer (optional) such that ``0 <= index < t`` where
          ``t`` denotes the number of infinite places of our base number field

        OUTPUT:

        If ``index`` is not specified, return the infinite part of ``self`` as
        a list of RIF/CIF/None elements.
        Else, return ``self.infinite()[i]``.

        .. WARNING::

            The input ``index`` is not checked for validity.
        """
        if index is None:
            return self._infinite
        return self._infinite[index]

    def finite(self, prime=None):
        """
        Return the finite part of ``self`` as a dictionary, or an entry of it

        INPUT:

        - ``prime`` -- a (finite) prime of our base number field (optional)

        OUTPUT:

        If ``index`` is not specified, return the finite part of ``self`` as
        a dictionary with (finite) primes as keys and pairs `(x, i)` as values,
        indiciating that we are equals at that prime to the `K`-element `x`, up
        to multiplication by `U(i)`.
        Else, return ``self.finite()[prime]``.

        .. WARNING::

            The input ``prime`` is not checked for validity.

        .. TODO::

            Improve this docstring
        """
        if prime is None:
            return self._finite
        if prime not in self._finite:
            K = self.parent().number_field()
            return (K(1), ZZ(0))
        return self._finite[prime]

    
    def _mul_(self, other):
        """
        Multiply ``self`` and ``other`` together

        .. NOTE::

            Precision may be lost. See the examples below.

        EXAMPLES:

        First a few multiplications of rational ideles::

            sage: J = Ideles(QQ)
            sage: u = J(None, None, {2: (1/2, 7), 3: (2/5, 8)})
            sage: v = J(None, [1.2345], None)
            sage: w = J(1/7, [-1.0], {3: (-1, 7)})
            sage: u*v
            Idele with values:
              infinity_0:  RR^*
              2:            1/2 * U(0)
              3:            2/5 * U(0)
            sage: u*w
            Idele with values:
              infinity_0:  RR^*
              2:            1/14 * U(7)
              3:            -2/5 * U(7)
              7:            1/7 * U(0)
            sage: w*u
            Idele with values:
              infinity_0:       RR^*
              2:                1/14 * U(7)
              3:                -2/5 * U(7)
              7:                1/7 * U(0)
            sage: v*w
            Idele with values:
              infinity_0:   -1.2345000000000000?
              3:            -1 * U(0)
              7:            1/7 * U(0)

        And now a few over a non-trivial number field::

            sage: K.<i> = NumberField(x^2+1)
            sage: J = Ideles(K)
            sage: p2, p3 = K.prime_above(2), K.prime_above(3)
            sage: u = J(i+1)
            sage: v = J(None, [I+1], {p3: (i/2, 7)})
            sage: w = J(i-1, None, {p2: (2, 5), p3: (3*i, 20)})
            sage: u*v
            Idele with values:
              infinity_0:       2*I
              (2, i + 1):       (i + 1) * U(0)
              (3, 0):   (1/2*i - 1/2) * U(7)
            sage: u*w
            Idele with values:
              infinity_0:  CC^*
              (2, i + 1):   (2*i + 2) * U(5)
              (3, 0):       (3*i - 3) * U(20)
              elsewhere:    -2
            sage: v*w
            Idele with values:
              infinity_0:  CC^*
              (2, i + 1):   2 * U(0)
              (3, 0):       -3/2 * U(7)
        """
        K = self.parent().number_field()

        exact = None
        if self._has_exact() and other._has_exact():
            exact = self.exact() * other.exact()

        infinite = self.infinite().copy()
        for i in range(len(infinite)):
            infinite[i] *= other.infinite(i)

        finite = self.finite().copy()
        for q, val in finite.items():
            if q in other.finite():
                x, i = val
                y, j = other.finite(q)
                finite[q] = (x*y, min(i, j))
            else:
                if other._has_exact():
                    finite[q] = (val[0] * other.exact(), val[1])
                else: # other is unknown (i.e. Z_q^*) at q
                    # U_q^i is a subgroup of Z_q^*, so x*U_q^i * Z_q^* = x*Z_q^*
                    finite[q] = (val[0], ZZ(0)) 
        for q, val in other.finite().items():
            if q not in finite:
                if self._has_exact():
                    finite[q] = (val[0] * self.exact(), val[1])
                else: # self is unknown (i.e. Z_q^*) at q
                    # U_q^i is a subgroup of Z_q^*, so x*U_q^i * Z_q^* = x*Z_q^*
                    finite[q] = (val[0], ZZ(0))
        if ((self._has_exact() and not other._has_exact())
                or (not self._has_exact() and other._has_exact())):
            # One has an exact value and one is Z_q^* almost everywhere.
            value = self.exact() if self._has_exact() else other.exact()
            # Where value has valuation 0, we have value * Z_q^* = Z_q^*, so we
            # don't have to store this.
            # But at the finitely many places where value has non-zero
            # valuation, we must store "value * Z_q^*".
            if K is QQ:
                I = value
            else:
                I = K.ideal(value)
            for q, e in factor(I):
                if q not in finite:
                    finite[q] = (value, ZZ(0))

        return self.__class__(self.parent(), exact, infinite, finite)

    def inverse(self):
        """
        Return the inverse of ``self`` in its idele group (i.e. "1/``self``")

        EXAMPLES::

            sage: K.<a> = NumberField(x^2+x+1)
            sage: p2, p3 = K.prime_above(2), K.prime_above(3)
            sage: J = Ideles(K)
            sage: u = J(None, [I], {p2: (a, 10), p3: (1/2, 7)})
            sage: u.inverse()
            Idele with values:
              infinity_0:   -1*I
              (2, 0):       (-a - 1) * U(10)
              (3, a + 2):   2 * U(7)
            sage: v = J(5, [1/2-I], {p3: (a+1, 97)})
            sage: v.inverse()
            Idele with values:
              infinity_0:   0.4000000000000000? + 0.80000000000000005?*I
              (3, a + 2):   -a * U(97)
              elsewhere:    1/5
        """
        exact = None
        if self._has_exact():
            exact = ZZ(1) / self.exact()

        infinite = self.infinite().copy()
        for i in range(len(infinite)):
            infinite[i] = ZZ(1) / infinite[i]

        finite = self.finite().copy()
        for q, val in finite.items():
            finite[q] = (ZZ(1) / val[0], val[1])

        return self.__class__(self.parent(), exact, infinite, finite)

    def _div_(self, other):
        """
        Divide ``self`` by ``other``

        EXAMPLES::

            sage: K.<a> = NumberField(x^2+8)
            sage: p2, q3 = K.prime_above(2), K.primes_above(3)[1]
            sage: J = Ideles(K)
            sage: u = J(None, [I+1], {p2: (a-1, 7), q3: (a/2, 9)})
            sage: v = J(a, [3], {q3: (2, 10)})
            sage: u/v
            Idele with values:
              infinity_0:   0.3333333333333334? + 0.3333333333333334?*I
              (2, 1/2*a):   (1/8*a + 1) * U(7)
              (3, 1/2*a + 2):       1/4*a * U(9)
            sage: v/u
            Idele with values:
              infinity_0:   1.5000000000000000? - 1.5000000000000000?*I
              (2, 1/2*a):   (-1/9*a + 8/9) * U(7)
              (3, 1/2*a + 2):       -1/2*a * U(9)
            sage: u/u
            Idele with values:
              infinity_0:   1
              (2, 1/2*a):   1 * U(7)
              (3, 1/2*a + 2):       1 * U(9)
        """
        return self * other.inverse()

    def integral_split(self):
        """
        Return a tuple ``(u, d)`` with ``u` an "integral" idele and ``d`` an
        integer such that ``self = u/d``

        The "integrality" of ``u`` is strong, in the following sense: at every
        finite prime `q` at which we stored a value `(x, i)` (representing
        `x*U_q^i`), the `K`-element ``x`` is integral in `K`; and if ``u`` is
        exact, its exact value is also integral in `K`.
        This is stronger than ``u`` being locally everywhere integral.

        EXAMPLES::

            sage: J = Ideles(QQ)
            sage: u = J(None, None, {2: (1/4, 1), 3: (3/5, 2)})
            sage: u.integral_split()
            (Idele with values:
               infinity_0: RR^*
               2:           5 * U(1)
               3:           12 * U(2)
               5:           20 * U(0),
             20)
            sage: v = J(3/14, [1], {2: (3/4, 5), 3: (1, 3)})
            sage: v.integral_split()
            (Idele with values:
               infinity_0:  28
               2:           21 * U(5)
               3:           28 * U(3)
               elsewhere:   6,
             28)
            sage: K.<a> = NumberField(x^2+10)
            sage: p3, p5 = K.prime_above(3), K.prime_above(5)
            sage: Jk = Ideles(K)
            sage: u = Jk(None, [I], {p3: (1/a, 3), p5: (7/6, 8)})
            sage: u.integral_split()
            (Idele with values:
               infinity_0:  30*I
               (2, a):      30 * U(0)
               (3, 0):      -3*a * U(3)
               (5, a):      35 * U(8),
             30)
        """
        from sage.arith.functions import lcm
        K = self.parent().number_field()

        denominators = []
        factors = factor(ZZ(1))
        for q, val in self.finite().items():
            d_q = val[0].denominator()
            denominators.append(d_q)
            if K is not QQ:
                d_q = K.ideal(d_q)
            factors *= factor(d_q)
        if self.exact() is not None:
            d_exact = self.exact().denominator()
            denominators.append(d_exact)
            if K is not QQ:
                d_exact = K.ideal(d_exact)
            factors *= factor(d_exact)
        d = lcm(denominators)

        # Now returning (d*self, d) would be correct. The multiplication will
        # usually need to factor d, which may be expensive. Less expensive is
        # to factor all the individual d_q's (as we did above) and build up our
        # ``u`` directly. This is what we will do below.

        finite = {}
        for q, val in self.finite().items():
            finite[q] = (d*val[0], val[1])
        exact = None
        if self._has_exact():
            exact = self.exact() * d
        else:
            for q, e in factors:
                if q not in finite:
                    finite[q] = (d, ZZ(0))
        infinite = self.infinite().copy()
        for i in range(len(infinite)):
            infinite[i] = d*infinite[i]
        u = self.__class__(self.parent(), exact, infinite, finite)

        return (u, d)

    def to_modulo_element(self):
        r"""
        Convert this idele to its image in O/I, with O the maximal order of K
        and I the (biggest) ideal where ``self`` is defined

        This implements the canonical homomorphism `\hat{O} \to O/I`, where
        we view `\hat{O}` as the subring of the idele group of ideles which are
        everywhere integral.

        Throws an exception if ``self`` does not lie in this subring.

        INPUT:

        - ``I`` -- an ideal of the maximal order ``O`` of our number field ``K``

        EXAMPLES:
            
        sage: J = Ideles(QQ)
        sage: u = J(None, None, {2: (1/3, 3), 5: (7, 1)})
        sage: u_bar = u.to_modulo_element()
        sage: u_bar, u_bar.parent()
        (27, Ring of integers modulo 40)
        """
        from adele import Adeles
        A = Adeles(self.parent().number_field())
        return A(self).to_modulo_element()

    def to_ray_class(self, modulus):
        r"""
        Return the image of ``self`` in the ray class group modulo ``modulus``

        Let `K` denote our number field and let `O` denote its maximal order.
        Let `I(modulus)` be the group of fractional `O`-ideals coprime to
        ``modulus`` and let `R(modulus)` denote the ray modulo ``modulus``.
        Hence the ray class group of `K` modulo ``modulus`` is
        `I(modulus)/R(modulus)`.
        Write `J_K` for the group of ideles of `K` (to which ``self`` belongs).

        This method implementst the homomorphism

        .. MATH::

            \phi: J_K \to I(modulus)/R(modulus)

        that sends a prime element at a finite prime `q` (not dividing `m`) to
        `q \mod R(modulus)`.

        INPUT:

        - ``modulus`` -- a :class:`Modulus` of `K`

        OUTPUT:

        The image of ``self`` under `\phi`.

        EXAMPLES::

            sage: Q = NumberField(x-1, "one")
            sage: J = Ideles(Q)
            sage: u = J(None, [-9], {2: (10, 1), 3: (3, 2), 7: (1/2, 4)})
            sage: m = Modulus(Q.ideal(18), [0])
            sage: u.to_ray_class(m).ideal()
            Fractional ideal (7)

        ::

            sage: K.<a> = NumberField(x^3-x-1)
            sage: Jk = Ideles(K)
            sage: p7, q7 = K.primes_above(7)
            sage: u = Jk(None, [2, I], {p7: (1/2, 1), q7: (7*a, 2)})
            sage: m = Modulus(K.ideal(7), [0])
            sage: u.to_ray_class(m).ideal()
            Fractional ideal (-2*a^2 - 3*a + 3)

        TESTS:

            sage: K = NumberField(x^2 + 5*x - 3, 'a')
            sage: J = Ideles(K)
            sage: G = ray_class_group(K, Modulus(K.ideal(14), [0])); G
            Ray class group of order 18 with structure C6 x C3 of Number Field in a with defining polynomial x^2 + 5*x - 3 of modulus (Fractional ideal (14)) * infinity_0
            sage: c0, c1 = G.gens()
            sage: for e in range(6):
            ....:     for f in range(3):
            ....:         r = c0^e * c1^f
            ....:         assert G(J(r)) == r, "bug!"  # long

        .. TODO::

            Implement the case that ``self`` has exact.
        
        ALGORITHM:

        We construct an idele `v` satisfying:
            - `v \equiv self \mod K^*`
            - `v \equiv 1 \mod^* modulus`
        We do this by starting of with `v = self` and then improve `v` in three
        steps:
            1. Make `v` integral at the primes dividing ``modulus``.
            2. Make `v \equiv 1 \mod^*` ``modulus.finite_part()``
            3. Make `v \equiv 1 \mod^*` ``modulus`` (fix the infinite part)
        In every step we only change `v` by multiplying it with principal ideles
        (i.e. elements of K^*) as to never violate the first desired condition
        on `v`.
        For Step 2 we use the number field version of the Chinese Remainder
        Theorem (cf. :func:`solve_CRT`). For Step 3 we use
        :meth:`get_one_mod_star_finite_with_fixed_signs`.

        Once we have such a `v`, we return the image of the ideal
        `\prod_{p} p^{\ord_p(v)}` in ``G``, where `p` ranges over the finite
        primes of `K`.
        """
        if self._has_exact():
            raise NotImplementedError("to_ray_class() not implemented yet for exact ideles")

        J = self.parent()
        K = J.number_field()
        G = ray_class_group(K, modulus)

        # First we check if the precision of this idele is high enough to have
        # a well-defined image in the ray class group ``G``.
        for i in modulus.infinite_part():
            if (not (self.infinite(i) <= 0 or self.infinite(i) >= 0)
                    or self.infinite(i).is_zero()):
                raise ValueError("idele has no well-defined sign at infinite prime {}".format(i))
        for q, e in modulus.finite_factors():
            if q not in self.finite() or self.finite(q)[1] < e:
                raise ValueError("idele must be known up to at least U_q^{} at q={}".format(e, q))

        v = self
        # Step 1. We make `v` integral at the primes dividing modulus.
        for q, e in modulus.finite_factors():
            x, i = self.finite(q)
            if x not in K.maximal_order():
                v *= x.denominator()

        # Step 2. We find an element y in the maximal order O_K of K such that
        # y \equiv v \mod modulus.finite_part() using the Chinese Remainder
        # Theorem.
        values = []
        moduli = []
        for q, e in modulus.finite_factors():
            x, i = v.finite(q)
            f = i + x.valuation(q)
            # v represents x*U_q^i at q, which equals x+q^f\Z_q because
            # i >= modulus.valuation(q) >= 1 and x.valuation(q) >= 0 (so we
            # do not need to worry about the case i==0 representing x*Z_q^*)
            values.append(x)
            moduli.append(q**f)
        y = K.solve_CRT(values, moduli)
        if not y.is_zero():
            # y is zero if and only if values and moduli are empty, so when
            # modulus is infinite. In that case we can just leave v as it is.
            v /= J(y)
        # Now we have v \equiv 1 \mod^* modulus.finite_part().
        

        # Step 3. We address the infinite part using the Modulus-method
        # get_one_mod_star_finite_with_fixed_signs().
        positive = []
        negative = []
        for i in modulus.infinite_part():
            if v.infinite(i) >= 0:
                positive.append(i)
            else:
                negative.append(i)
        t = modulus.get_one_mod_star_finite_with_fixed_signs(positive, negative)
        v *= t

        # TODO do not do the redundant check below to save its costs
        assert v.is_one_mod_star(modulus), r"Assertion error in Idele.to_ray_class(): v \not\equiv 1 \mod^* m"

        # Our `v` is finished. We can now build up an ideal representing the
        # image of `v` (which is equal to the image of ``self``) in ``G``.
        I = K.unit_ideal()
        for q, val in v.finite().items():
            I *= q**(val[0].valuation(q))

        return G(I)


    def is_one_mod_star(self, modulus):
        r"""
        Return whether or not ` ``self`` \equiv 1 \mod^\ast ``modulus`` ` holds

        EXAMPLES::

            sage: Q = NumberField(x-1, "one")
            sage: J = Ideles(Q)
            sage: u = J(None, None, {5: (6, 2)})
            sage: u.is_one_mod_star(Modulus(Q.ideal(3)))
            False
            sage: u.is_one_mod_star(Modulus(Q.ideal(5)))
            True
            sage: u.is_one_mod_star(Modulus(Q.ideal(25)))
            False
            sage: u.is_one_mod_star(Modulus(Q.ideal(5), [0]))
            False
            sage: u._infinite[0] = RIF(1.2345)
            sage: u.is_one_mod_star(Modulus(Q.ideal(5), [0]))
            True
            sage: u._infinite[0] = RIF(-1.2345)
            sage: u.is_one_mod_star(Modulus(Q.ideal(5), [0]))
            False
            sage: u._exact = Q.one()
            sage: u.is_one_mod_star(Modulus(Q.ideal(2*3*5*7*11)))
            True
        """
        for q, e in modulus.finite_factors():
            if not q in self.finite():
                if not self._has_exact():
                    return False
                x = self.exact()
            else:
                x, i = self.finite(q)
                if i < e:
                    return False
            if (x-1).valuation(q) < e:
                return False
        for i in modulus.infinite_part():
            if not (self.infinite(i) >= 0):
                return False
        return True

    def _contains_at(self, x, q):
        r"""
        Check if this idele contais the `K`-element ``x`` at the prime ``q``

        Suppose that the open subset this idele represents is `U_q^i` at ``q``.
        Then we check whether or not `x \in U_q^i` holds.

        INPUT:

        - ``x`` -- element of `K`
        - ``q`` - prime (finite or infinite) of `K`; in the infinite case, it
          should be one of the embbedings returned by ``K.places()``.

        EXAMPLES:
        
        First some examples with `\QQ`::

            sage: J = Ideles(QQ)
            sage: u = J(None, None, {2: (5, 3), 3: (1/3, 0)})
            sage: u._contains_at(5, 2)
            True
            sage: u._contains_at(2, 2)
            False
            sage: u._contains_at(5/3, 3)
            True
            sage: u._contains_at(5/9, 3)
            False
            sage: u._contains_at(13, 5)
            True
            sage: u._contains_at(130, 5)
            False
            sage: u._contains_at(3.1415, QQ.places()[0])
            True

        And some examples with a non-trivial number field::

            sage: K.<a> = NumberField(x^3+x+1)
            sage: Jk = Ideles(K)
            sage: v = Jk(a, [RIF(-0.68, -0.69), I], {K.prime_above(7): (a+1, 8)})
            sage: v._contains_at(a, K.places()[0])
            True
            sage: v._contains_at(a+1, K.places()[0])
            False
            sage: v._contains_at(a, K.places()[1])
            False
            sage: v._contains_at(a, K.prime_above(7))
            False
            sage: v._contains_at((a+1)*(7^8+1), K.prime_above(7))
            True

        And lastly some errors that may occur:

            sage: u._contains_at(1, 8)
            Traceback (most recent call last):
            ...
            TypeError: ``q`` should be a prime of `K`, but got 8
            sage: v._contains_at(a^2, None)
            Traceback (most recent call last):
            ...
            TypeError: ``q`` should be a prime of `K`, but got None
        """
        K = self.parent().number_field()
        if q in K.places():
            i = K.places().index(q)
            phi = completions(K, oo)[i][1] # phi: K --> RIF/CIF
            return phi(x) in self.infinite(i)

        if (K is QQ and q in Primes()) or (K is not QQ and q in K.ideal_monoid() and q.is_prime()):
            if q in self.finite():
                y, i = self.finite(q)
                if i == ZZ(0):
                    return (x/y).valuation(q) == ZZ(0)
                else:
                    return (x/y - 1).valuation(q) >= i
            elif self._has_exact():
                return x == self.exact()
            else:
                # Not stored represents Z_q^*: the elements of valuation 0.
                return x.valuation(q) == 0

        raise TypeError("``q`` should be a prime of `K`, but got {}".format(q))

    def _contains(self, x):
        """
        Check if the subset this idele represents contains the `K`-element ``x``

        INPUT:

        - ``x`` -- element of `K`, the number field to which this idele belongs

        EXAMPLES::

            sage: K.<a> = NumberField(x^7-4/7*x^3+1)
            sage: J = Ideles(K)
            sage: u = J(a^5-a)
            sage: v_oo = completions(K, oo)[0][1](a^3)
            sage: C = CIF(RIF(-oo, oo), RIF(-oo, oo))
            sage: v = J(None, [v_oo, C, C, C], {K.prime_above(5): (a^3, 4)})
            sage: u._contains(a^5-a)
            True
            sage: u._contains(1)
            False
            sage: v._contains(a^3)
            True
            sage: v._contains(a)
            False
        """
        if self._has_exact() and x != self.exact():
            return False
        K = self.parent().number_field()
        for phi in K.places():
            if not self._contains_at(x, phi):
                return False
        for q in self.finite():
            if not self._contains_at(x, q):
                return False
        return True

    def is_principal(self):
        """
        Return whether or not this idele is principal

        We check if there exists an "original" rational ``r`` whose image lies in
        the open subset op the idele group that this idele represents.

        ALGORITHM:

        If this idele has an exact value, that is our only candidate for ``r``.
        Else we create a generator by ... TODO

        EXAMPLES::

            sage: K.<a> = NumberField(x^2+3)
            sage: J = Ideles(K)
            sage: u = J(a, None, {K.prime_above(2): (-a+1, 10)})
            sage: v = J(None, None, {K.prime_above(2): (-a+1, 4), K.prime_above(5): (-a-1, 0)})
            sage: w = J(None, None, {K.prime_above(2): (-a+1, 4), K.prime_above(5): (1/5, 0)})
            sage: u.is_principal()
            False
            sage: v.is_principal()
            True
            sage: w.is_principal()
            False
        """
        K = self.parent().number_field()
        if self._has_exact():
            # The only possible ``r`` is ``self.exact()``.
            return self._contains(self.exact())

        # we are not exact.
        r = K.ideal(1)
        for q, val in self.finite().items():
            x, i = val
            r *= q**(x.valuation(q))
        if not r.is_principal():
            # There does not exist an element of `K` with the right valuations.
            return False
        # The ideal ``r`` is principal, so we can take a generator.
        r = r.gens_reduced()[0]
        # This will be our candidate element of `K` that could lie in ``self``.
        # It is only uniquely determined up to a unit of `K`.
        U = K.unit_group()
        if not U.is_finite():
            raise NotImplementedError("K has infinite unit group, we can't deal with that yet")
        if len(U.gens()) > 1:
            raise NotImplementedError("K has non-cyclic unit group, we can't deal with that yet")
        u = U.gen().value()
        for e in range(U.order()):
            if self._contains(u**e * r):
                return True
        return False
        r"""
        old:
        if all([u_oo is None for u_oo in self.infinite()]):
            # Only finite values stored.
            if K is QQ and len(self.finite()) == 1:
                p = list(self.finite().keys())[0]
                x, i = self.finite(p)
                # The only candidate element is p^e, where:
                e = x.valuation(p)
                # Need to find some `z \in \ZZ_p` such that
                #     `p^e = x*(1+p^i*z)`.
                # Defining:
                y = x / p^e
                # Gives as only option:
                #     `z = (p^e-x)/(x*p^i) = (1-y)/(y*p^i)`
                # which lies in `\QQ`, but we should check that it lies in
                # `\ZZ_p \cap \QQ = \ZZ_{(p)}` (localization at `p`). 
                # That means that the valuation of `z` at `p` should be at least
                # zero. Since the valuation of `y` at `p` is zero, this boils
                # down to:
                return (1-y).valuation(p) >= i
        """

    def equivalent_modulo_principal_ideles(self, v):
        r"""
        Return whether or not `self \equiv v \mod K*`

        INPUT:

        - v -- an idele
        """
        return (self/v).is_principal()


    def increase_precision(self, primes, prec_increment=1):
        """
        Increase the precision of ``self`` at the primes given in ``primes``

        INPUT:

        - ``primes`` -- an iterable containing prime ideals of our number field
          and/or rational prime numbers, or a single prime
        - ``prec_increment`` -- integer (default = 1); the amount by which we
          increase the precision at each prime in ``primes``

        Let `p` be a prime ideal in ``primes``. Suppose ``self`` represents the
        open subset `x * U_p^i` at `p`. Then after calling this method, ``self``
        will represent `x * U_p^(i + prec_increment)` at `p`.

        If `p` in ``primes`` is a rational prime number, then the above is done
        for each prime `q` lying above `p`.

        .. NOTE::

            If ``self`` has an exact known value at a prime `p`, then nothing
            changes. If one sees exactness as having infinite precision, this
            just corresponds to ``oo + prec_increment == oo``.

        .. NOTE::

            Setting ``prec_increment`` to a negative value will decrease the
            precision of ``self``. If the precision drops below zero anywhere,
            we throw an exception.

        EXAMPLE::

            sage: K.<a> = NumberField(x^3-2)
            sage: J = Ideles(K)
            sage: p2 = K.prime_above(2)
            sage: p5 = K.prime_above(5)
            sage: u = J(None, None, {p2: (a, 5), p5: (1/3, 1)})

        Let's increase the precision of ``p2`` and *both* prime ideals above 5
        by 3::

            sage: u.increase_precision([p2, 5], 3)
            sage: u
            Idele with values:
              infinity_0:  RR^*
              infinity_1:  CC^*
              (2, a):       a * U(8)
              (5, a + 2):   1/3 * U(4)
              (5, a^2 - 2*a - 1):   1 * U(3)

        We can also decrease precision::

            sage: u.increase_precision(p2, -1)
            sage: u
            Idele with values:
              infinity_0:  RR^*
              infinity_1:  CC^*
              (2, a):       a * U(7)
              (5, a + 2):   1/3 * U(4)
              (5, a^2 - 2*a - 1):   1 * U(3)

        The precision of exact values is unchanged::
        
            sage: u._exact = a+1
            sage: p3 = K.prime_above(3)
            sage: u.increase_precision([p2, p3, p5])
            sage: u
            Idele with values:
              infinity_0:  RR^*
              infinity_1:  CC^*
              (2, a):       a * U(8)
              (5, a + 2):   1/3 * U(5)
              (5, a^2 - 2*a - 1):   1 * U(3)
              elsewhere:    a + 1
        """
        K = self.parent().number_field()

        if primes in K.ideal_monoid() and K.ideal(primes).is_prime():
            # primes is just a single prime ideal, let's call it p.
            p = ZZ(primes) if K is QQ else K.ideal(primes)
            if p in self.finite():
                x, i = self.finite(p)
            elif self._has_exact():
                # We know the value at p exactly, so we don't change
                return
            else:
                x, i = K(1), ZZ(0)
            new_prec = i + prec_increment
            if new_prec < 0:
                raise ValueError("Trying to give idele negative precision")
            self._finite[p] = (x, new_prec)
            return

        if primes in Primes():
            p = primes
            for q in K.primes_above(p):
                self.increase_precision(q, prec_increment)
            return

        for p in primes:
            if p in K.ideal_monoid() and K.ideal(p).is_prime():
                self.increase_precision(p, prec_increment)
            elif p in Primes():
                for q in K.primes_above(p):
                    self.increase_precision(q, prec_increment)
            else:
                raise TypeError("primes should be a list of primes")



class Ideles(UniqueRepresentation, Group):
    Element = Idele

    def __init__(self, K):
        # TODO: implement default K=QQ, via __classcall__()
        from sage.misc.functional import is_field
        if not is_field(K) or not K.absolute_degree() in ZZ:
            raise TypeError("K should be a number field")
        self._number_field = K
        Group.__init__(self)

    def _repr_(self):
        """
        Return a string representation of ``self``
        """
        return "Idele Group of {}".format(self.number_field())

    def _latex_(self):
        from sage.misc.latex import latex
        return r"\Bold{A}_{" + latex(self.number_field()) + "}^*"

    def _element_constructor_(self, exact, infinite=None, finite=None):
        if infinite is None and finite is None:
            x = exact
            if x is None:
                raise TypeError("No arguments supplied to Idele-constructor")
            K = self.number_field()
            if x in K:
                infinite = [phi(x) for L, phi in infinite_completions(K)]
                return self.element_class(self, x, infinite, {})
            from adele import Adeles
            Ak = Adeles(K)
            if Ak.has_coerce_map_from(x.parent()): # conversion A_K --> J_K
                return self._from_adele(Ak(x))
            if hasattr(x.parent(), "_bnr"): # conversion Cl_m --> J_K
                # TODO make the check above less hacky and more robust
                # for some reason checking isinstance(x, RayClassGroupElement) fails
                return self._from_ray_class(x)
        return self.element_class(self, exact, infinite, finite)

    def _from_adele(self, adele):
        r"""
        Convert the adele ``adele`` to an idele

        Let U be the subset of the adele ring that ``adele`` represents. Write J
        for the idele group. Then the idele returned represents a subset that
        *contains* `U \cap J`. Note that in general this looses precision.
        TODO: example of precision loss

        EXAMPLES::
            
            sage: from profinite_number import ProfiniteNumbers
            sage: from adele import Adeles
            sage: J = Ideles(QQ)
            sage: A = Adeles(QQ)
            sage: Qhat = ProfiniteNumbers(QQ)
            sage: a = A([-1.5], Qhat(7, 24, 5))
            sage: J._from_adele(a)
            Idele with values:
              infinity_0:   -1.5000000000000000?
              2:            7 * U(3)
              3:            7 * U(1)
            sage: b = A([RIF(-1, 1)], Qhat(5, 25, 2))
            sage: J._from_adele(b)
            Idele with values:
              infinity_0:   0.?
              5:            5 * U(1)
            sage: c = A([pi.n()], 20/3)
            sage: J._from_adele(c)
            Idele with values:
              infinity_0:   3.1415926535897932?
              elsewhere:    20/3

        If the given adele has value zero modulo at one of its stored primes
        (i.e. ``4 mod 12`` has value zero modulo 2^2), then there is no idele
        that represents the adele. For example: the given adele then represents
        multiple element which have different valuation at that prime. But an
        idele always has a unique valuation at every prime. In this case, we
        throw an exception::

            sage: d = A([1], Qhat(4, 12, 3))
            sage: J._from_adele(d)
            Traceback (most recent call last):
            ...
            ValueError: adele is zero at the prime 2
        """
        K = self.number_field()
        value = adele.finite().numerator().value()
        modulus = adele.finite().numerator().modulus()
        denominator = adele.finite().denominator()

        exact = None
        if modulus.is_zero():
            exact = value / denominator
            finite = None
        else:
            x = value
            if K is not QQ:
                x = K.ideal(x)

            finite = {}
            for q, e in factor(modulus):
                v = x.valuation(q)
                if v >= e:
                    # value is zero modulo q^e, which we cannot represent as an
                    # idele
                    raise ValueError("adele is zero at the prime {}".format(q))
                finite[q] = (value, e-v)

        return self.element_class(self, exact, adele.infinite().copy(), finite)

    def _from_ray_class(self, r):
        r"""
        Convert the ray class group element ``r`` to an idele

        INPUT:

        - ``r`` -- a ray class group element

        OUPUT:
        
        Denote the ray class group to which ``r`` belongs by `G` and denote the
        modulus of `G` by `m`. So `G = I(m)/R(m)`.
        Consider the homomorphism `\phi: ` ``self`` `\to` `G` that sends a prime
        element at a finite prime `q` (not dividing `m`) to `q \mod R(m)`.
        The kernel of `\phi` is `K^* W_m` where `W_m = \prod_p U_p^{\ord_p(m)}`.

        Given ``r``, let `H` be the inverse image of ``r`` under `\phi`. We can
        try to find an idele that represents the subset `H` of ``self``. This is
        however not precisely possible. We can exactly represent `W_m`, but 
        we can not represent `K^*`. Hence what we do is the following: we find
        some `x \in H` en return an idele that represents `x \cdot W_m`.

        Although we do always return the same idele for equal inputs, the user
        should be aware that from a mathematical perspective, the output is only
        defined up to multiplication by a principal idele.
        

        EXAMPLES::

            sage: Q = NumberField(x-1, "one")
            sage: J = Ideles(Q)
            sage: G = ray_class_group(Q, Modulus(Q.ideal(10), [0]))
            sage: r = G(Q.ideal(9))
            sage: factor(r.ideal())
            (Fractional ideal (3)) * (Fractional ideal (163))
            sage: J._from_ray_class(r)
            Idele with values:
              infinity_0:   [0.0000000000000000 .. +infinity]
              (2, 0):       1 * U(1)
              (5, 0):       1 * U(1)
              (3, 0):       3 * U(0)
              (163, 0):     163 * U(0)
            sage: s = G(Q.ideal(7))
            sage: s.ideal()
            Fractional ideal (67)
            sage: J._from_ray_class(s)
            Idele with values:
              infinity_0:   [0.0000000000000000 .. +infinity]
              (2, 0):       1 * U(1)
              (5, 0):       1 * U(1)
              (67, 0):      67 * U(0)
           
        :: 

            sage: K.<a> = NumberField(x^2-6)
            sage: G = ray_class_group(K, Modulus(K.ideal(10*a), [1]))
            sage: r = G([3, 0, 1])
            sage: factor(r.ideal())
            (Fractional ideal (25*a + 19)) * (Fractional ideal (28*a - 25)) * (Fractional ideal (-67*a + 109)) * (Fractional ideal (-1507*a - 5011))
            sage: Jk = Ideles(K)
            sage: Jk(r)
            Idele with values:
              infinity_0:  RR^*
              infinity_1:   [0.0000000000000000 .. +infinity]
              (2, a):       1 * U(3)
              (3, a):       1 * U(1)
              (5, a + 1):   1 * U(1)
              (5, a + 4):   1 * U(1)
              (3389, a + 543):      (a + 543) * U(0)
              (4079, a + 2767):     (a - 1312) * U(0)
              (15053, a + 13254):   (a - 1799) * U(0)
              (11483827, a + 1242116):      (a + 1242116) * U(0)
        """
        K = self.number_field()
        G = r.parent()  # ray class group of r
        exact = None
        r1, r2 = K.signature()
        RR = RIF(-oo, oo)
        CC = CIF(RR, RR)
        infinite = [RR for i in range(r1)] + [CC for i in range(r2)]
        for i in G.modulus().infinite_part():
            infinite[i] = RIF(0, oo)

        finite = {}
        for q, e in G.modulus().finite_factors():
            finite[q] = (K(1), e)

        if not r.is_one():
            for q, e in factor(r.ideal()):
                if self.number_field() is QQ:
                    finite[q] = (q**e, ZZ(0))
                else:
                    pi = K.uniformizer(q)
                    finite[q] = (pi**e, ZZ(0))
        
        return self.element_class(self, exact, infinite, finite)

    def cardinality(self):
        return Infinity

    def _coerce_map_from_(self, S):
        if self.number_field().has_coerce_map_from(S):
            return True
        return False

    def number_field(self):
        """
        Return the base number field of ``self``
        """
        return self._number_field

    def _an_element_(self):
        """
        Return a typical element of this group

        EXAMPLE::

            sage: Ideles(QQ).an_element()
            Idele with values:
              infinity_0:   3.1415926535897901?
              2:            3/2 * U(7)
              3:            2 * U(9)
              elsewhere:    -1
            sage: K.<a> = NumberField(x^3-2)
            sage: Ideles(K).an_element()
            Idele with values:
              infinity_0:   1
              infinity_1:   1
              (2, a):       (a + 1) * U(7)
              (3, a + 1):   1/2*a^2 * U(9)
              elsewhere:    -a
        """
        K = self.number_field()
        if K is QQ:
            infinite = [3.14159265358979]
            p2, p3 = ZZ(2), ZZ(3)
        else:
            r1, r2 = K.signature()
            infinite = [i for i in range(1,r1+1)] + [1-i*CIF.gen() for i in range(r2)]
            p2, p3 = K.prime_above(2), K.prime_above(3)
        finite = {
            p2: (1+K.an_element(), 7),
            p3: (1/K.an_element(), 9)
        }
        return self.element_class(self, -K.gens()[-1], infinite, finite)

