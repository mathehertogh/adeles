import sys # TODO erase these two lines
sys.path.append('/home/mathe/adeles/src')

from sage.structure.element import MultiplicativeGroupElement
from sage.structure.unique_representation import UniqueRepresentation
from sage.groups.group import Group
from sage.categories.groups import Groups
from sage.arith.misc import factor
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.real_mpfi import RIF
from sage.rings.complex_interval_field import ComplexIntervalField
from sage.rings.infinity import Infinity
from sage.sets.primes import Primes
from sage.arith.functions import lcm

from completion import completions, infinite_completions
from multiplicative_padic import is_finite_prime, multPAdic
from ray_class_group import Modulus, ray_class_group

CIF = ComplexIntervalField()
oo = Infinity







class Idele(MultiplicativeGroupElement):
    r"""
    An idele: an element of an :class:`Ideles` of a number field

    .. automethod:: __init__
    """

    def __init__(self, parent, infinite=None, finite=None, check=True):
        r"""
        Construct an idèle

        We denote the base number field by `K` and its maximal order by `O`.

        INPUT:

        - ``parent`` -- instance of :class:`Ideles`; parent idèle group
        - ``infinite`` -- ``None`` or a list or ``RIF`` or ``CIF`` elements with
          as many entries as `K` has infinite places. The ``i``-th entry
          specifies the value of this idèle at the infinite place corresponding
          to ``K.places()[i]``. The entry should be non-zero and lie in ``RIF``
          or ``CIF``, depending on the place being real or complex.
          ``None`` means we do not know anything about the values of this idele
          at infinite primes.
        - ``finite`` -- a dictionary or an element of `K^*` (default: empty
          dictionary);

          If ``finite`` is a dictionary, it should contain (``key``, ``value``)
          pairs such that:

              - ``key`` is a finite prime `p` of ``K``; if ``K`` is `\QQ`, this
                is a prime number, else it's a prime ideal of ``O``.
              - ``value`` should be a pair (``x``, ``i``) with ``x`` an element
                of `K^*` and ``i`` a non-negative integer.

          Such an entry means: at prime `p` of `K`, this idele is equal to (the
          image of) ``x`` (in `K_p`), up to multiplication with elements of the
          ``n``-th standard subgroup of `K_p^*`.

          If ``finite`` lies in `K^*`, then this specifies the exact
          finite value of this idele. Then the value of this idele is exactly
          equal to ``finite`` at each finite prime of `K`.

        - ``check`` -- boolean (default: ``True``): whether or not to validate
          the input data

        EXAMPLES:

        We begin with the easiest case: `K = \QQ`::

            sage: J = Ideles(QQ)
            sage: Idele(J, [3.14], {2: (5/3, 2), 5: (1, 9)})
            Idele with values:
              infinity_0:   3.1400000000000002?
              2:            5/3 * U(2)
              5:            1 * U(9)
              other primes: 1 * U(0)
            sage: Idele(J, [RIF(-1.122, -1.119)], 1/3)
            Idele with values:
              infinity_0:   -1.12?
              other primes: 1/3
            sage: Idele(J, None, None)
            Idele with values:
              infinity_0:   RR^*
              other primes: 1 * U(0)

        Now let's take a non-trivial number field::

            sage: K.<a> = NumberField(x^4-17)
            sage: Jk = Ideles(K)
            sage: p3, p7 = K.prime_above(3), K.prime_above(7)
            sage: Idele(Jk, [3.14, -1, 1+I], {p3: (a^3, 7), p7: (a^2, 4)})
            Idele with values:
              infinity_0:   3.1400000000000002?
              infinity_1:   -1
              infinity_2:   1 + 1*I
              (3, 1/2*a^2 - a - 1/2):       a^3 * U(7)
              (7, 1/2*a^2 + a - 5/2):       a^2 * U(4)
              other primes: 1 * U(0)
            sage: Idele(Jk, None, a^3-a^2+2*a-5)
            Idele with values:
              infinity_0:   RR^*
              infinity_1:   RR^*
              infinity_2:   CC^*
              other primes: a^3 - a^2 + 2*a - 5 

        TESTS::

            sage: Idele(J, [1], 0)
            Traceback (most recent call last):
            ...
            ValueError: exact finite value must be non-zero
            sage: Idele(J, [1], "blah")
            Traceback (most recent call last):
            ...
            TypeError: finite must be a dictionary or a non-zero number field element

        .. TODO::

            Implement parameter ``check=True``
        """
        MultiplicativeGroupElement.__init__(self, parent)

        self._infinite = infinite
        if check:
            self._validate_infinite()

        self._finite = finite
        if check:
            self._validate_finite()

    def _validate_infinite(self):
        """
        Raise an exception if ``self._infinite`` is not valid

        ``self._infinite`` is valid if it's a list of length ``len(K.places())``
        with the `i`-th entry a non-zero element of ``RIF`` if ``K.places()[i]``
        is real and a non-zero element of ``CIF`` otherwise.
        Note that only the intervals `[0, 0]` in ``RIF`` and `[0, 0] + [0,0]*I`
        in ``RIC`` are zero. For example, `[-1, 5]` is non-zero although it does
        represent zero.

        If the `i`-th entry's parent is not equal to ``RIF``/``CIF``, but does
        lie in it, than we cast the element to an actual element of such an
        interval field.

        EXAMPLES:

        We look at a number field with one real and one complex prime::

            sage: K.<a> = NumberField(x^3 + x + 1)
            sage: Jk = Ideles(K)
            sage: u = Jk(1)
            sage: u._infinite = [1.2, I]
            sage: u._validate_infinite()

        Below the ``Integer 3`` will be cast into ``RIF`` and ``I`` will be cast
        into ``CIF``::

            sage: u._infinite = [ZZ(3), I]
            sage: u._validate_infinite()
            sage: u._infinite[0].parent()
            Real Interval Field with 53 bits of precision
            sage: u._infinite[1].parent()
            Complex Interval Field with 53 bits of precision

        If ``_infinite`` is ``None``, we make a list of ``RR``s and ``CC``s of
        if::

            sage: u._infinite = None
            sage: u._validate_infinite()
            sage: u
            Idele with values:
              infinity_0:       RR^*
              infinity_1:       CC^*
              other primes:     1

        The exceptions that we may throw:

            sage: u._infinite = [3.14, I, I]
            sage: u._validate_infinite()
            Traceback (most recent call last):
            ...
            ValueError: infinite must have length 2
            sage: u._infinite = {CIF: I+1}
            sage: u._validate_infinite()
            Traceback (most recent call last):
            ...
            TypeError: infinite must be a list
            sage: u._infinite = [I, I]
            sage: u._validate_infinite()
            Traceback (most recent call last):
            ...
            TypeError: 0th infinite value (I) must lie in Real Interval Field with 53 bits of precision
            sage: u._infinite = [1, CIF(0)]
            sage: u._validate_infinite()
            Traceback (most recent call last):
            ...
            ValueError: 1th infinite value (0) must be non-zero
        """
        K = self.parent().number_field()

        if self._infinite is None:
            r1, r2 = K.signature()
            RR = RIF(-oo, oo)
            CC = CIF(RR, RR)
            self._infinite = [RR for i in range(r1)] + [CC for i in range(r2)]
            return

        K_oo = infinite_completions(K, fields_only=True)
        t = len(K_oo)

        try:
            if t == 1 and self._infinite in K_oo[0]:
                self._infinite = [self._infinite]
            else:
                self._infinite = list(self._infinite)
        except TypeError:
            raise TypeError("infinite must be iterable")

        if len(self._infinite) != t:
            raise ValueError("infinite must have length {}".format(t))

        for i in range(t):
            val = self._infinite[i]
            if val not in K_oo[i]:
                raise TypeError("{}th infinite value ({}) must lie in {}".format(i, val, K_oo[i]))
            if val == 0:
                raise ValueError("{}th infinite value ({}) must be non-zero".format(i, val))
            self._infinite[i] = K_oo[i](val)

    def _validate_finite(self):
        r"""
        Raise an exception if ``self._finite`` is not valid
    
        ``self._finite`` is valid if it is a dictionary with keys finite primes of
        `K` (our base number field) and values pairs in `K^* \times \NN`.
        A finite prime of `K` is a prime number if `K` is `\QQ` and a prime
        ideal of the ring of integers of `K` otherwise.

        This method also casts the entries of each value pair ``(x, n)`` into
        `K` and ``ZZ``.

        EXAMPLES:

        A few examples of the small adjustments we can do::

            sage: K.<a> = NumberField(x^4+x+1)
            sage: J = Ideles(K)
            sage: u = J(1)
            sage: u._finite = {a+3: (int(3), int(1))}
            sage: u._validate_finite()
            sage: u._finite # a+3 is turned into an ideal:
            {Fractional ideal (a + 3): (35*a^3, 1)}
            sage: u[a+3][0].parent() # 3 into a number field element:
            Number Field in a with defining polynomial x^4 + x + 1
            sage: u[a+3][1].parent() # 1 into a Sage Integer
            Integer Ring

        Note that above ``3`` was turned into the number field element
        ``35*a^3``. That is because ``35*a^3`` is the HNF-reduction of ``3``
        modulo ``p7``. Hence ``3 * U(1)`` is the same set as ``35*a^3 * U(1)``
        at ``p7``.

        A ``None`` finite value is turned into an empty dictionary::

            sage: u._finite = None
            sage: u._validate_finite()
            sage: u._finite
            {}

        And some exceptions we may throw::

            sage: u._finite = 0
            sage: u._validate_finite()
            Traceback (most recent call last):
            ...
            ValueError: exact finite value must be non-zero
            sage: u._finite = ["blah"]
            sage: u._validate_finite()
            Traceback (most recent call last):
            ...
            TypeError: finite must be a dictionary or a non-zero number field element

        Keys must be primes of K. Although 5 is a prime number, it is not a
        prime of `K` as 5 splits into two primes in `K`::

            sage: u._finite = {5: (1, 1)}
            sage: u._validate_finite()
            Traceback (most recent call last):
            ...
            TypeError: keys of finite must be prime ideals of the number field

        Values are checked for correctness as well::

            sage: u._finite = {K.prime_above(2): (1, -1)}
            sage: u._validate_finite()
            Traceback (most recent call last):
            ...
            TypeError: values at finite primes must be pairs in K^* x NN
            sage: u._finite = {K.prime_above(2): ([], 3)}
            sage: u._validate_finite()
            Traceback (most recent call last):
            ...
            TypeError: values at finite primes must be pairs in K^* x NN
        """
        K = self.parent().number_field()

        if self._finite is None:
            self._finite = {}
        elif self._finite in K:
            if self._finite == 0:
                raise ValueError("exact finite value must be non-zero")
            self._finite = K(self._finite)
        elif isinstance(self._finite, dict):
            # We will replace self._finite will our own dictionary ``finite`` in
            # which we ensure all keys and values have the correct parents.
            finite = {}
            for P, data in self._finite.items():
                if not is_finite_prime(P, K):
                    raise TypeError("keys of finite must be finite primes")
                P = ZZ(P) if K is QQ else K.ideal(P)

                finite[P] = multPAdic(P, data)

            self._finite = finite
        else:
            raise TypeError("finite must be a dictionary or a non-zero number field element")

    def _repr_(self):
        """
        Return a representation of this idele.

        EXAMPLES::

            sage: J = Ideles(QQ)
            sage: u = J(None, {97: (7/9, 20)}); u
            Idele with values:
              infinity_0:       RR^*
              97:               7/9 * U(20)
              other primes:     1 * U(0)
            sage: v = J([-7.123456], 5/8); v
            Idele with values:
              infinity_0:   -7.1234560000000001?
              other primes: 5/8
        """
        K = self.parent().number_field()
        rep = "Idele with values:"

        RR = RIF(-oo, oo)
        CC = CIF(RR, RR)
        for i in range(len(self.infinite_part())):
            x_i = self[oo, i]
            if x_i.endpoints() == RR.endpoints():
                rep += "\n  infinity_{}:\tRR^*".format(i)
            elif x_i.endpoints() == CC.endpoints():
                rep += "\n  infinity_{}:\tCC^*".format(i)
            else:
                rep += "\n  infinity_{}:\t{}".format(i, x_i)

        for P in self.stored_primes():
            P_name = self[P].parent().prime_name()
            tab = "\t\t" if len(P_name) < 5 else "\t"
            rep += "\n  {}:{}{}".format(P_name, tab, self[P])

        if self.has_exact_finite_part():
            rep += "\n  other primes:\t{}".format(self.finite_part())
        else:
            rep += "\n  other primes:\t1 * U(0)"

        return rep

    def __getitem__(self, prime):
        """
        Get the value of this idèle at the prime ``prime``

        TODO INPUT

        EXAMPLES::

            sage: J = Ideles(QQ)
            sage: u = J([3.1], {2: (1,2)})
            sage: u[oo]
            3.1000000000000001?
            sage: u[2]
            (1, 2)
            sage: u[3]
            (1, 0)
            sage: u[4]
            Traceback (most recent call last):
            ...
            KeyError: 'You can only index by a prime ideal or (Infinity, index)'

        ::

            sage: K.<a> = NumberField(x^2+7)
            sage: Jk = Ideles(K)
            sage: v = Jk([I], {K.prime_above(5): (a, 3)})
            sage: v[oo]
            1*I
            sage: v[oo, 0]
            1*I
            sage: v[K.prime_above(2)]
            (1, 0)
            sage: v[K.prime_above(5)]
            (a, 3)
            sage: v[7]
            Traceback (most recent call last):
            ...
            KeyError: 'You can only index by a prime ideal or (Infinity, index)'
        """
        K = self.parent().number_field()
        if is_finite_prime(prime, K):
            if self.has_exact_finite_part():
                return self.finite_part()
            try:
                return self.finite_part()[prime]
            except KeyError:
                return multPAdic(prime, (K(1), ZZ(0)))
        if prime is Infinity and len(self.infinite_part()) == 1:
            return self.infinite_part()[0]
        try:
            inf, index = prime
            if inf is Infinity:
                return self.infinite_part()[index]
        except TypeError:
            pass
        raise KeyError("You can only index by a prime ideal or (Infinity, index)")

    def infinite_part(self):
        """
        Return the infinite part of this idele as a list

        EXAMPLES::

            sage: K.<a> = NumberField(x^5-3*x+1)
            sage: J = Ideles(K)
            sage: K.signature()
            (3, 1)
            sage: u = J([-1, 7.9, 1, 1+I], a^2)
            sage: u.infinite_part()
            [-1, 7.9000000000000004?, 1, 1 + 1*I]
        """
        return self._infinite

    def finite_part(self):
        """
        Return the finite part of this idèle

        This could be either a number field element or a dictionary.

        EXAMPLES::

            sage: K.<a> = NumberField(x^5-3*x+1)
            sage: J = Ideles(K)
            sage: u = J([-1, 7.9, 1, 1+I], a^2)
            sage: u.finite_part()
            a^2
            sage: p3, q3 = K.primes_above(3)
            sage: v = J(None, {p3: (a,10), q3: (-1-a^2,7)})
            sage: v.finite_part()
            {Fractional ideal (-2*a^4 - a^3 - 2*a^2 - a + 4): (-920*a^4 - 31*a^3 - 90*a^2 - 961*a, 7),
             Fractional ideal (-a^2 + 1): (-11107*a^4, 10)}
        """
        return self._finite

    def has_exact_finite_part(self):
        """
        Return whether or not this idele has exact finite part

        EXAMPLES::

            sage: J = Ideles(QQ)
            sage: u = J([-1], 7)
            sage: u.has_exact_finite_part()
            True
            sage: v = J([-1], {2: (1, 1)})
            sage: v.has_exact_finite_part()
            False
        """
        K = self.parent().number_field()
        return self._finite in K

    def stored_primes(self):
        r"""
        Return the set of stored primes of this idele as an ordered list.

        A stored prime is a finite prime of the base number field `K` for which
        this idele has stored a value in its ``_finite`` attribute.

        If this idele has exact finite part, then it has no stored primes.
        
        The order is determined by the implemented order on number field ideals,
        or by the natural order on `\NN` if `K = \QQ` (in which case the primes
        are prime numbers).

        EXAMPLES::

            sage: K.<a> = NumberField(x^2+17)
            sage: J = Ideles(K)
            sage: p2, p5 = K.prime_above(2), K.prime_above(5)
            sage: p7, q7 = K.primes_above(7)
            sage: u = J([I], {q7: (a, 1), p2: (8, 0), p7: (1,0), p5: (a, 3)})
            sage: u.stored_primes()
            [Fractional ideal (2, a + 1),
             Fractional ideal (5),
             Fractional ideal (7, a + 2),
             Fractional ideal (7, a + 5)]

        ::

            sage: v = J([-1], a+1)
            sage: v.stored_primes()
            []
        """
        if self.has_exact_finite_part():
            return []
        return sorted(list(self.finite_part().keys()))

    def valuation(self, prime):
        r"""
        Return the valuation of this idèle at the finite prime ``prime``

        INPUT:

        - ``prime`` -- a finite prime

        EXAMPLES::

            sage: J = Ideles(QQ)
            sage: u = J([1], 1/25)
            sage: u.valuation(2)
            0
            sage: u.valuation(5)
            -2
            sage: v = J([1], {2: (8, 10)})
            sage: v.valuation(2)
            3
            sage: v.valuation(3)
            0

        ::

            sage: K.<a> = NumberField(x^2+5)
            sage: J = Ideles(K)
            sage: p3, q3 = K.primes_above(3)
            sage: u = J([I], {p3: (3, 0)})
            sage: u.valuation(p3)
            1
            sage: u.valuation(q3)
            0

        TESTS::

            sage: Q.<one> = NumberField(x-1)
            sage: J = Ideles(Q)
            sage: u = J([-1], {7: (7^20, 3)})
            sage: u.valuation(Q.prime_above(7))
            20
        """
        if self.has_exact_finite_part():
            return self.finite_part().valuation(prime)
        return self[prime].valuation()

    def increase_precision(self, primes, prec_increment=1):
        """
        Increase the precision of this idele at the primes given in ``primes``

        INPUT:

        - ``primes`` -- an iterable containing finite primes and/or (rational)
          prime numbers, or a single prime
        - ``prec_increment`` -- integer (default = 1); the amount by which we
          increase the precision at each prime in ``primes``

        Let `p` be a finite prime in ``primes``. Suppose ``self`` represents the
        open subset `x * U_p^n` at `p`. Then after calling this method, ``self``
        will represent `x * U_p^(n + prec_increment)` at `p`.

        If `p` in ``primes`` is a rational prime number, then the above is done
        for each prime `q` lying above `p` with ``prec_increment`` multiplied
        by the ramification index of `q` over `p`.

        .. NOTE::

            If this idele has exact finite part, then this method does not do
            anything. If one sees exactness as having infinite precision, this
            just corresponds to ``oo + prec_increment == oo``.

        .. NOTE::

            Setting ``prec_increment`` to a negative value will decrease the
            precision, but not below zero.

        EXAMPLE::

            sage: K.<a> = NumberField(x^3-2)
            sage: J = Ideles(K)
            sage: p2, p5 = K.prime_above(2), K.prime_above(5)
            sage: u = J(None, {p2: (a, 5), p5: (1/3, 1)})

        Let's increase the precision of ``p2`` and *both* prime ideals above 5
        by 3::

            sage: u.increase_precision([p2, 5], 3); u
            Idele with values:
              infinity_0:   RR^*
              infinity_1:   CC^*
              (2, a):       a * U(8)
              (5, a^2 - 2*a - 1):   (36*a^2 + 33*a) * U(3)
              (5, a + 2):   1/3 * U(4)
              other primes: 1 * U(0)

        We can also decrease precision::

            sage: u.increase_precision(p2, -1); u
            Idele with values:
              infinity_0:   RR^*
              infinity_1:   CC^*
              (2, a):       a * U(7)
              (5, a^2 - 2*a - 1):   (36*a^2 + 33*a) * U(3)
              (5, a + 2):   1/3 * U(4)
              other primes: 1 * U(0)

        As ``p2`` has ramification index 3, the following will increase the
        precision of ``u`` at ``p2`` by 3::

            sage: u.increase_precision(2); u
            Idele with values:
              infinity_0:   RR^*
              infinity_1:   CC^*
              (2, a):       a * U(10)
              (5, a^2 - 2*a - 1):   (36*a^2 + 33*a) * U(3)
              (5, a + 2):   1/3 * U(4)
              other primes: 1 * U(0)

        Nothing changes for ideles with exact finite part::
        
            sage: v = J(None, a^2)
            sage: v.increase_precision([p2, p5]); v
            Idele with values:
              infinity_0:   RR^*
              infinity_1:   CC^*
              other primes: a^2

        TESTS::

            sage: u.increase_precision([QQ])
            Traceback (most recent call last):
            ...
            ValueError: primes should be a (list of) prime(s)
        """
        if self.has_exact_finite_part():
            return

        K = self.parent().number_field()

        if is_finite_prime(primes, K):
            P = primes # primes is a single finite prime
            new_prec = max(ZZ(0), self[P].prec() + prec_increment)
            self._finite[P] = multPAdic(P, (self[P].center(), new_prec))
            return

        if primes in Primes():
            p = primes # primes is a single rational prime number
            for P in K.primes_above(p):
                self.increase_precision(P, prec_increment*P.ramification_index())
            return

        for prime in primes:
            if not (is_finite_prime(prime, K) or prime in Primes()):
                raise ValueError("primes should be a (list of) prime(s)")
            self.increase_precision(prime, prec_increment)

    def integral_split(self):
        """
        TODO
        """
        from sage.arith.functions import lcm
        K = self.parent().number_field()

        if self.has_exact_finite_part():
            d = self.finite_part().denominator()
        else:
            d = lcm([self[P].center().denominator() for P in self.stored_primes()])

        return (self.parent()(d)*self, d)

    def _mul_(self, other):
        """
        Multiply this idele with ``other``

        EXAMPLES::

            sage: J = Ideles(QQ)
            sage: u = J(-1, {2: (1/2, 7), 3: (2/5, 8)}); u
            Idele with values:
              infinity_0:   -1
              2:            1/2 * U(7)
              3:            2/5 * U(8)
              other primes: 1 * U(0)
            sage: v = J(2.5, {3: (-1, 4)}); v
            Idele with values:
              infinity_0:   2.5000000000000000?
              3:            80 * U(4)
              other primes: 1 * U(0)
            sage: u*v
            Idele with values:
              infinity_0:   -2.5000000000000000?
              2:            1/2 * U(0)
              3:            32 * U(4)
              other primes: 1 * U(0)
        """
        from sage.modules.free_module_element import vector

        # At the infinite primes we perform component-wise multiplication.
        infinite = vector(self.infinite_part()) * vector(other.infinite_part())

        if self.has_exact_finite_part() and other.has_exact_finite_part():
            finite = self.finite_part() * other.finite_part()
        elif not self.has_exact_finite_part() and not other.has_exact_finite_part():
            stored_primes = set(self.stored_primes() + other.stored_primes())
            finite = dict([(P, self[P]*other[P]) for P in stored_primes])
        else:
            if not self.has_exact_finite_part():
                self, other = other, self
            # Now self has exact finite part while other does not.
            stored_primes = set(self.finite_part().support() + other.stored_primes())
            finite = dict([(P, self[P]*other[P]) for P in stored_primes])

        return self.__class__(self.parent(), infinite, finite)

    def inverse(self):
        """
        Return the inverse of this idèle
        """
        infinite = self.infinite_part().copy()
        for i in range(len(infinite)):
            infinite[i] = 1 / infinite[i]

        if self.has_exact_finite_part():
            finite = 1 / self.finite_part()
        else:
            finite = dict([(P, self[P].inverse()) for P in self.stored_primes()])

        return self.__class__(self.parent(), infinite, finite)

    def _div_(self, other):
        """"
        TODO
        """
        return self * other.inverse()

    def _equals(self, other):
        """
        TODO
        """
        raise NotImplementedError("_equals() not implemtend yet TODO")

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


class Ideles(UniqueRepresentation, Group):
    Element = Idele

    def __init__(self, K):
        # TODO: docstring
        # TODO: check is_NumberField mooier maken
        # TODO: implement default K=QQ, via __classcall__()
        from sage.misc.functional import is_field
        if not is_field(K) or not K.absolute_degree() in ZZ:
            raise TypeError("K should be a number field")
        self._number_field = K
        Group.__init__(self, category=Groups().Commutative())

    def _repr_(self):
        """
        Return a string representation of this idèle group

        EXAMPLES::

            sage: K.<a> = NumberField(x^2-3)
            sage: Ideles(K)
            Idele Group of Number Field in a with defining polynomial x^2 - 3
        """
        return "Idele Group of {}".format(self.number_field())

    def _latex_(self):
        r"""
        Return a latex string representation of this idèle group

        EXAMPLES::

            sage: K.<a> = NumberField(x^2-5)
            sage: latex(Ideles(K))
            J_{ \Bold{Q}[a]/(a^{2} - 5) }
        """
        from sage.misc.latex import latex
        return "J_{" + latex(self.number_field()) + "}"

    def number_field(self):
        """
        Return the base number field of ``self``

        EXAMPLES::

            sage: K.<a> = NumberField(x^7-1000*x+1)
            sage: J = Ideles(K)
            sage: J.number_field()
            Number Field in a with defining polynomial x^7 - 1000*x + 1
        """
        return self._number_field

    def _element_constructor_(self, x, y=None):
        # TODO: docstring
        if y is None:
            K = self.number_field()
            if x in K:
                infinite = [phi(x) for L, phi in infinite_completions(K)]
                return self.element_class(self, infinite, x)
            if hasattr(x, "parent") and hasattr(x.parent(), "_bnr"):
                # TODO make the check above less hacky and more robust
                # for some reason checking isinstance(x, RayClassGroupElement) fails
                return self._from_ray_class(x)
        return self.element_class(self, x, y)

    def _from_modulo_element(self, x):
        """
        Create an idele from ``x`` in `(O/I)^*`
        """
        I = x.parent().defining_ideal()
        finite_part = {}
        for P, e in I.factor():
            finite_part[P] = (x.lift(), e)
        return self.element_class(self, None, finite_part)

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
              (3, 0):       3 * U(0)
              (5, 0):       1 * U(1)
              (163, 0):     163 * U(0)
              other primes: 1 * U(0)
            sage: s = G(Q.ideal(7))
            sage: s.ideal()
            Fractional ideal (67)
            sage: J._from_ray_class(s)
            Idele with values:
              infinity_0:   [0.0000000000000000 .. +infinity]
              (2, 0):       1 * U(1)
              (5, 0):       1 * U(1)
              (67, 0):      67 * U(0)
              other primes: 1 * U(0)

        :: 

            sage: K.<a> = NumberField(x^2-6)
            sage: G = ray_class_group(K, Modulus(K.ideal(10*a), [1]))
            sage: r = G([3, 0, 1])
            sage: factor(r.ideal())
            (Fractional ideal (25*a + 19)) * (Fractional ideal (28*a - 25)) * (Fractional ideal (-67*a + 109)) * (Fractional ideal (-1507*a - 5011))
            sage: Jk = Ideles(K)
            sage: Jk(r)
            Idele with values:
              infinity_0:   RR^*
              infinity_1:   [0.0000000000000000 .. +infinity]
              (2, a):       1 * U(3)
              (3, a):       1 * U(1)
              (5, a + 1):   -a * U(1)
              (5, a + 4):   a * U(1)
              (3389, a + 543):      2846760*a * U(0)
              (4079, a + 2767):     -2916485*a * U(0)
              (15053, a + 13254):   -94683370*a * U(0)
              (11483827, a + 1242116):      -22108284774109*a * U(0)
              other primes: 1 * U(0)
        """
        K = self.number_field()
        G = r.parent()  # ray class group of r
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
        
        return self.element_class(self, infinite, finite)