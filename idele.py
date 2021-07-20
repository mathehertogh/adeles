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
from ray_class_group import Modulus, ray_class_group

CIF = ComplexIntervalField()
oo = Infinity


def _associated_multPAdic(center, precision, prime):
    r"""
    Return the associated multiplicative ``prime``-adic of the pair
    (``center``, ``precision``) in `K \times \ZZ_{\geq 0}`

    EXAMPLES::

        sage: _associated_multPAdic(-3, 3, 2)
        (5, 3)
        sage: _associated_multPAdic(-3/5, 4, 5)
        (622/5, 4)
        sage: -3/5 - 622/5
        -125

    ::

        sage: K.<a> = NumberField(x^2+5)
        sage: p3 = K.prime_above(3)
        sage: _associated_multPAdic(-1+2*a, 4, p3)
        (-75*a, 4)
        sage: (-1+2*a - (-75*a)).valuation(p3)
        5
    """
    e = max(1, precision) + center.valuation(prime)
    I = prime**e
    if I in ZZ:
        den = lcm(center.denominator(), I.denominator())
        num_new_center = (den*center) % (den*I)
    else:
        den = lcm(center.denominator(), I.integral_split()[1])
        num_new_center = (den*I).reduce(den*center)
    new_center = num_new_center / den
    return (new_center, precision)



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
        K = parent.number_field()

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

        if not isinstance(self._infinite, list):
            raise TypeError("infinite must be a list")

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
            self._exact = K(self._finite)
        elif isinstance(self._finite, dict):
            # We will replace self._finite will our own dictionary ``finite`` in
            # which we ensure all keys and values have the correct parents.
            finite = {}
            for P, val in self._finite.items():
                if K is QQ:
                    if P not in Primes():
                        raise TypeError("keys of finite must be prime numbers")
                    P = ZZ(P)
                else:
                    if P not in K.ideal_monoid() or not K.ideal(P).is_prime():
                        raise TypeError("keys of finite must be prime ideals of the number field")
                    P = K.ideal(P)
                
                try:
                    center, prec = val
                except TypeError:
                    raise TypeError("values at finite primes must be pairs in K^* x NN")
                if center not in K or center == 0 or prec not in ZZ or prec < 0:
                    raise TypeError("values at finite primes must be pairs in K^* x NN")

                finite[P] = _associated_multPAdic(K(center), ZZ(prec), P)

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
            center, prec = self[P]
            p_name = str(P) if K is QQ else str(P.gens_two())
            tab = "\t\t" if len(p_name) < 5 else "\t"
            if len([c for c in center.list() if not c.is_zero()]) > 1:
                center = "(" + str(center) + ")"
            rep += "\n  {}:{}{} * U({})".format(p_name, tab, center, prec)

        if self.has_exact_finite_part():
            rep += "\n  other primes:\t{}".format(self.finite_part())
        else:
            rep += "\n  other primes:\t1 * U(0)"

        return rep

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

    def infinite_part(self, index=None):
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

    def finite_part(self, prime=None):
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

    def __getitem__(self, prime):
        """
        Get the value of this idèle at the prime ``prime``

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
        if ((K is QQ and prime in Primes()) or
                (prime in K.ideal_monoid() and K.ideal(prime).is_prime())):
            if self.has_exact_finite_part():
                return self.finite_part()
            try:
                return self.finite_part()[prime]
            except KeyError:
                return (K(1), ZZ(0))
        if prime is Infinity and len(self.infinite_part()) == 1:
            return self.infinite_part()[0]
        try:
            inf, index = prime
            if inf is Infinity:
                return self.infinite_part()[index]
        except TypeError:
            pass
        raise KeyError("You can only index by a prime ideal or (Infinity, index)")


class Ideles(UniqueRepresentation, Group):
    Element = Idele

    def __init__(self, K):
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

    def _element_constructor_(self, x, y=None):
        if y is None:
            K = self.number_field()
            if x in K:
                infinite = [phi(x) for L, phi in infinite_completions(K)]
                return self.element_class(self, infinite, x)
            from adele import Adeles
            Ak = Adeles(K)
            if x in Ak:
                return self._from_adele(Ak(x))
            if hasattr(x, "parent") and hasattr(x.parent(), "_bnr"):
                # TODO make the check above less hacky and more robust
                # for some reason checking isinstance(x, RayClassGroupElement) fails
                return self._from_ray_class(x)
        return self.element_class(self, x, y)

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
