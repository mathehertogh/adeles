"""
Profinite Integers of Number Fields

.. TODO::

    Write documentation for this module
"""
import sys # TODO erase these two lines
sys.path.append('/home/mathe/adeles/src')

from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.ring import CommutativeAlgebra
from sage.structure.element import CommutativeAlgebraElement
from sage.arith.misc import factor
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.sets.primes import Primes

from profinite_integer import ProfiniteIntegers


class ProfiniteNumber(CommutativeAlgebraElement):
    """
    Profinite Integer of a Number Field

    We try to keep gcd(value, modulus, denominator) small (in, say, norm). But
    this is difficult for number fields with non-trivial class group. See
    :meth:`_reduce`.
    
    .. TODO::
    
        Write documentation for this class definition
    
    :automethod:`_repr_`
    """

    def __init__(self, parent, numerator, denominator=None):
        r"""
        Construct the profinite number ``numerator / denominator``
        
        INPUT:

        - ``parent`` -- A ring of profinite numbers over some number field `K`
        - ``numerator`` -- an element of ``ProfiniteIntegers(K)``
        - ``denominator`` -- an element of ``K.maximal_order()``

        OUTPUT:

        The profinite number ``numerator``/``denominator``.

        EXAMPLES::

            sage: Zhat = ProfiniteIntegers(ZZ)
            sage: Qhat = ProfiniteNumbers(QQ)
            sage: ProfiniteNumber(Qhat, Zhat(2, 15), 8)
            (2 mod 15)/8
            sage: ProfiniteNumber(Qhat, -97)
            -97
            sage: ProfiniteNumber(Qhat, Zhat(2, 16), 8)
            (1 mod 8)/4

        ::

            sage: K.<a> = NumberField(x^4+7)
            sage: Ohat = ProfiniteIntegers(K)
            sage: Khat = ProfiniteNumbers(K)
            sage: ProfiniteNumber(Khat, Ohat(a^3+a, a^3), 2*a-1)
            (a^3 + a mod (a^3))/(2*a - 1)
            sage: ProfiniteNumber(Khat, Ohat(a+2, a^3), 2*a-1)
            (-3/4*a^3 - 7/4*a^2 - 3/4*a + 1/4 mod (a^3))/(2*a - 1)

        TODO:

        The following SHOULD work but does not::

            sage: ProfiniteNumber(Khat, a)
            a

        This is because the parent of a is K, not O. And there is not coercion
        from K to Ohat. Hence a in Ohat returns False.
        We should fix this by the following from sage manual:
            "Note x∈P holds if and only if the test x==P(x) does not raise an
             error and evaluates as true."
        Currently a==Ohat(a) returns False, but I think if we make sure this
        is evaluated inside Khat, we get True. Hence we shouuld make pushouts
        work. See also:
        https://doc.sagemath.org/html/en/thematic_tutorials/coercion_and_categories.html

        TESTS::

            sage: ProfiniteNumber(Qhat, a, 3)
            Traceback (most recent call last):
            ...
            TypeError: numerator should be an element of Profinite Integers of Rational Field
            sage: ProfiniteNumber(Qhat, 3, a)
            Traceback (most recent call last):
            ...
            TypeError: denominator should be a non-zero element of Integer Ring
            sage: ProfiniteNumber(Khat, Ohat(a), 0)
            Traceback (most recent call last):
            ...
            TypeError: denominator should be a non-zero element of Maximal Order in Number Field in a with defining polynomial x^4 + 7
        """
        CommutativeAlgebraElement.__init__(self, parent)
        O = parent.base().maximal_order()
        Ohat = ProfiniteIntegers(O)
        if denominator is None:
            denominator = O.one()
        if numerator not in Ohat:
            raise TypeError("numerator should be an element of {}".format(Ohat))
        if denominator not in O or denominator == 0:
            raise TypeError("denominator should be a non-zero element of {}".format(O))
        self._numerator = Ohat(numerator)
        self._denominator = O(denominator)
        self._reduce()

    def _reduce(self):
        """
        Try to divide out the greatest common divisor of ``value``, ``modulus``
        and ``denominator``, where ``self == (value mod modulus)/denominator``.

        Note that this is not possible if the greatest common divisor is a
        non-principal ideal. In this case we leave ``self`` unchanged.

        EXAMPLES::

            sage: Qhat = ProfiniteNumbers()
            sage: a = Qhat(1, 100, 60); a
            (1 mod 100)/60
            sage: a._numerator._value = 10; a
            (10 mod 100)/60
            sage: a._reduce(); a
            (1 mod 10)/6
            sage: a._numerator._value = 14; a
            (14 mod 10)/6
            sage: a._reduce(); a
            (2 mod 5)/3

        ::

            sage: K.<a> = NumberField(x^2+5)
            sage: Khat = ProfiniteNumbers(K)
            sage: Khat(a, 5, a)
            (-1 mod (a))/-1
            sage: Khat(3*a, 25, 7*a)
            (-3 mod (5*a))/-7
            sage: Khat(3*a, 75, 21*a)
            (-1 mod (5*a))/-7

        And lastly some examples with non-principal ideals, in which case we
        cannot always reduce::

            sage: O = K.maximal_order()
            sage: I = O.ideal(2, a+1)
            sage: I.divides(2)
            True
            sage: I.divides(a+1)
            True
            sage: I.is_principal()
            False
            sage: Khat(a+1, 2*I, 2)
            (a + 1 mod (4, 2*a + 2))/2

        Above the greatest common divisor of the three entries is the ideal
        ``I``. But since ``I`` is non-principal, we cannot divide it out of the
        elements ``a+1`` and ``2`` of ``O``.

        .. TODO::

            There is still improvement possible in this case of non-principal
            GCD. For example::

                sage: Khat(2*(a+1), 5*I^3, 4)
                (2*a + 2 mod (20, 10*a + 10))/4

            Here the GCD is ``I^3``, which is not principal. But ``I^2 = (2)``
            is. So we could divide out a factor 2.

            But:
                - how do we find these principal divisors?
                - is there a unique maximal one or should we chooose one?
                - what to base this choice on?
        """
        O = self.parent().base().maximal_order()
        if O is ZZ:
            from sage.arith.misc import gcd
            common = gcd(self.numerator().value(), self.numerator().modulus())
            common = gcd(common, self.denominator())
            self._numerator /= common
            self._denominator = ZZ(self.denominator() // common)
        else:
            common = O.ideal(self.numerator().value()) + O.ideal(self.numerator().modulus())
            common += O.ideal(self.denominator())
            gens = common.gens_reduced()
            if len(gens) == 1:
                self._numerator /= O(gens[0])
                self._denominator = O(self.denominator() // gens[0])

    def _repr_(self):
        """
        Returns a string representation of ``self``

        EXAMPLES::

            sage: K.<a> = NumberField(x^2+5)
            sage: Khat = ProfiniteNumbers(K)
            sage: Khat(a-1, K.ideal(4, 2*a+2), 2*a-5)
            (-a + 1 mod (4, 2*a + 2))/(2*a - 5)
        """
        if self.denominator() == 1:
            return repr(self.numerator())
        numerator = repr(self.numerator())
        if (not self.numerator().modulus().is_zero() or
                sum([not c.is_zero() for c in self.numerator().value().list()]) > 1):
            numerator = "(" + numerator + ")"
        denominator = repr(self.denominator())
        if not (self.denominator() in QQ or
                (sum([not c.is_zero() for c in self.denominator().list()]) == 1
                    and sum([c for c in self.denominator().list()]) == 1)):
            denominator = "(" + denominator + ")"
        return numerator + "/" + denominator

    def _richcmp_(self, other, op):
        r"""
        Compare ``self`` and ``other`` based on the relation ``op``

        We only implement equality and non-equality.

        We declare ``a/b`` equal to ``c/d`` if and only if ``a*d == b*c``.

        EXAMPLES::

            sage: Qhat = ProfiniteNumbers()
            sage: Qhat(1, 6, 2) == Qhat(2, 12, 4)
            True
            sage: Qhat(1, 6, 2) == Qhat(2, 12, 3)
            False
            sage: Qhat(1, 6, 2) != Qhat(2, 12, 3)
            True

        ::

            sage: K.<a> = NumberField(x^2+5)
            sage: Khat = ProfiniteNumbers(K)
            sage: I = K.ideal(2, a+1)
            sage: b = Khat(2*(a+1), 5*I^3, 4); b
            (2*a + 2 mod (20, 10*a + 10))/4
            sage: c = Khat(a+1, 5*I, 2); c
            (a + 1 mod (10, 5*a + 5))/2
            sage: b == c
            True
            sage: b != c
            False
        """
        from sage.structure.richcmp import op_EQ, op_NE
        if op == op_EQ:
            return self.numerator() * other.denominator() == self.denominator() * other.numerator()
        if op == op_NE:
            return not self._richcmp_(other, op_EQ)
        raise NotImplementedError("only equality and inequality are implemented")

    def _add_(self, other):
        """
        Add ``other`` to ``self`` and return the result

        EXAMPLES::
        
            sage: Qhat = ProfiniteNumbers()
            sage: Qhat(1, 10, 2) + Qhat(2, 15, 3)
            (7 mod 30)/6
            sage: Qhat(2, 5) + Qhat(3, 11, 3)
            (0 mod 1)/3

        ::

            sage: K.<a> = NumberField(x^3-2)
            sage: Khat = ProfiniteNumbers(K)
            sage: Khat(a, 10, 3) + Khat(a^2, 15, 3)
            (a^2 + a mod (5))/3
            sage: Khat(a+1, 20*a+3, 7) + 5
            (-3749*a^2 mod (-20*a - 3))/7
        """
        numerator = self.numerator() * other.denominator() + other.numerator() * self.denominator()
        denominator = self.denominator() * other.denominator()
        return self.__class__(self.parent(), numerator, denominator)
        
    def _sub_(self, other):
        """
        Subtract ``other`` from ``self`` and return the result

        EXAMPLES::
        
            sage: Qhat = ProfiniteNumbers(QQ)
            sage: Qhat(4, 12, 3) - 2/3
            (2 mod 12)/3
            sage: 1 - Qhat(3, 6, 2)
            (5 mod 6)/2
            sage: Qhat(4, 12, 3) - Qhat(3, 6, 2)
            (5 mod 6)/6

        ::

            sage: K.<a> = NumberField(x^3-2)
            sage: Khat = ProfiniteNumbers(K)
            sage: 4 - Khat(a, 7, 2) == Khat(8-a, 7, 2)
            True
            sage: Khat(a^2, 21, 2*(a+1)) - a^2/(2*(a+1))
            (0 mod (7*a^2 + 14*a + 7))/(2*a^2 + 2*a + 2)
            sage: Khat(a^2, 21, 2*(a+1)) - Khat(a, 7, 2)
            (-a mod (7*a + 7))/(2*a + 2)
        """
        numerator = self.numerator() * other.denominator() - other.numerator() * self.denominator()
        denominator = self.denominator() * other.denominator()
        return self.__class__(self.parent(), numerator, denominator)
        
    def _mul_(self, other):
        """
        Multiply ``self`` and ``other`` and return the result

        EXAMPLES::
        
            sage: Qhat = ProfiniteNumbers(QQ)
            sage: Qhat(5, 25, 3) * Qhat(0, 75, 2)
            (0 mod 125)/2
            sage: Qhat(5, 25, 3) * Qhat(3, 75, 2)
            (5 mod 25)/2
            sage: 3/4 * Qhat(3, 10, 5)
            (9 mod 30)/20

        ::

            sage: K.<a> = CyclotomicField(5)
            sage: Khat = ProfiniteNumbers(K)
            sage: a^3 * Khat(a^2, 21, 5)
            (1 mod (21*a^3))/5
            sage: Khat(a^2, 21, 5) * Khat(a+1, 21, 10)
            (a^3 + a^2 mod (21))/50
        """
        numerator = self.numerator() * other.numerator()
        denominator = self.denominator() * other.denominator()
        return self.__class__(self.parent(), numerator, denominator)
        
    def _div_(self, other):
        r"""
        Divide ``self`` and ``other`` and return the result

        Only implemented for exact ``other`` (i.e. ``other`` a number field
        element).

        EXAMPLES::
        
            sage: Qhat = ProfiniteNumbers(QQ)
            sage: Qhat(3, 10, 2) / 2
            (3 mod 10)/4
            sage: Qhat(3, 10, 2) / (2/3)
            (9 mod 30)/4
            sage: Qhat(3, 10, 2) / Qhat(3, 0, 2)
            (3 mod 10)/3

        ::

            sage: K.<a> = CyclotomicField(5)
            sage: Khat = ProfiniteNumbers(K)
            sage: Khat(a+1, 7, 2) / a
            (a + 1 mod (7))/(2*a)
            sage: Khat(a+1, 7, 2) / (3*a/2)
            (a + 1 mod (7))/(3*a)
        """
        if not other.numerator().modulus().is_zero():
            raise NotImplementedError("Division of profinite numbers only implemented for exact denominators")
        numerator = self.numerator() * other.denominator()
        denominator = self.denominator() * other.numerator().value()
        return self.__class__(self.parent(), numerator, denominator)

    def __getitem__(self, p):
        r"""
        Return `x_p` where ``self`` `= \prod_p x_p \in \prod_p' \QQ_p`

        INPUT:

        - ``p`` -- a prime number

        .. NOTE:

            Only implemented over `\QQ`, since we currently only have `p`-adics
            for rational `p` in Sage.
        """
        if self.parent().base() is not QQ:
            raise NotImplementedError("Projection to `p`-adics only implemented over rationals")
        if p not in Primes():
            raise ValueError("p should be a prime number")

        return self.numerator()[p] / self.denominator()

    def to_profinite_rational_vector(self, enclosure=True):
        r"""
        Return ``self`` as a profinite rational vector

        Denote our base field by `K`.
        Then we have a canonical isomorphism of topological rings
        `\phi: \hat{K} \to \hat{\QQ} \otimes K`.
        The codomain is a free `\hat{\QQ}`-module with basis
        `B = \{1, a, a^2, ... a^{n-1}\}`, where `a` is the algebraic number
        adjoined to `\QQ` to obtain `K` and `n = \deg(K/\QQ)`.

        This method returns the image of ``self`` under `\phi` in the form of an
        `\hat{\QQ}`-vector relative to the basis `B`.

        INPUT:

        - ``enclosure`` -- boolean (default: ``True``); whether or not to return
          an enclosure of the exact result `u`. We cannot in general represent
          `u`. Hence we must choose to either return an enclosure `v` of `u`,
          meaning `v` represents at least all profinite numbers that `u`
          represents; or return a vector with "too much" precision: some `w`
          that represents at most all profinite numbers that `u` represents.

        EXAMPLES::

            sage: K.<a> = NumberField(x^4-2)
            sage: Khat = ProfiniteNumbers(K)
            sage: b = Khat(a, 12*(a-1), 5)
            sage: b.to_profinite_rational_vector()
            ((0 mod 12)/5, (1 mod 6)/5, (0 mod 6)/5, (0 mod 6)/5)
            sage: b.to_profinite_rational_vector(enclosure=False)
            ((0 mod 12)/5, (1 mod 12)/5, (0 mod 12)/5, (0 mod 12)/5)
        """
        from sage.modules.free_module_element import vector
        from sage.arith.functions import lcm
        K = self.parent().base()
        n = K.absolute_degree()
        Qhat = ProfiniteNumbers(QQ)

        if self.numerator().modulus().is_zero():
            # self is (exactly) the following element x of K:
            x = self.numerator().value() / self.denominator()
            return vector(Qhat, x.vector())
        
        x = self.numerator().value() / self.denominator()
        I = self.numerator().modulus() / self.denominator()
        # self is `x mod I` with x in K, I a fractional ideal of K

        values = x.vector()
        moduli = []
        if enclosure:
            e_p = {}
            for i in range(n):
                for q, e in factor(I):
                    p = q.gens_two()[0]  # p = q \cap ZZ
                    e_q = p.valuation(q)
                    if p not in e_p:
                        e_p[p] = e // e_q
                    else:
                        e_p[p] = min(e_p[p], e // e_q)
                modulus = QQ(1)
                for p in e_p:
                    modulus *= p**e_p[p]
                moduli.append(modulus)
                I /= K.gen()
        else:
            for i in range(n):
                m = I.gens_two()[0]  # m = I \cap QQ
                moduli.append(m)
                I /= K.gen()

        profinite_rationals = []
        for i in range(n):
            denominator = lcm(values[i].denominator(), moduli[i].denominator())
            value = values[i] * denominator
            modulus = moduli[i] * denominator
            profinite_rationals.append(Qhat(value, modulus, denominator))

        return vector(profinite_rationals)

    def numerator(self):
        """
        Return the numerator of ``self`` as a profinite integer
        """
        return self._numerator

    def denominator(self):
        """
        Return the denominator of ``self`` as an integral element of our base
        field
        """
        return self._denominator

    def modulus(self):
        """
        Return the modulus of self as a fractional ideal of our number field
        """
        return self.numerator().modulus() / self.denominator()

    def value(self):
        """
        Return the value of self as an element of our number field
        """
        return self.numerator().value() / self.denominator()

    def is_integral(self):
        """
        Return whether or not ``self`` is integral
        """
        return self.value().is_integral() and self.modulus().is_integral()


class ProfiniteNumbers(UniqueRepresentation, CommutativeAlgebra):
    Element = ProfiniteNumber

    def __classcall__(cls, K=QQ):
        """
        Construct the ring of profinite numbers over the number field `K`

        INPUT:

        - ``K`` -- a number field (default: ``QQ``)

        EXAMPLES::

            sage: ProfiniteNumbers("not a number field")
            Traceback (most recent call last):
            ...
            TypeError: K should be a number field

        We make sure this returns True::

            sage: ProfiniteNumbers() is ProfiniteNumbers(QQ)
            True
        """
        try:
            from sage.misc.functional import is_field
            if not is_field(K) or not K.absolute_degree() in ZZ:
                raise TypeError("K should be a number field")
        except AttributeError:
            raise TypeError("K should be a number field")
        return super(ProfiniteNumbers, cls).__classcall__(cls, K)

    def __init__(self, K):
        r"""
        Construct the ring of profinite numbers over the number field ``K``

        INPUT:

        - ``K`` -- a number field

        EXAMPLES::

            sage: ProfiniteNumbers()
            Profinite Numbers of Rational Field
            sage: K.<a> = NumberField(x^5-3*x+1)
            sage: ProfiniteNumbers(K)
            Profinite Numbers of Number Field in a with defining polynomial x^5 - 3*x + 1
        """
        CommutativeAlgebra.__init__(self, K)

    def _repr_(self):
        """
        Return a string representation of ``self``

        EXAMPLES::

            sage: K = NumberField(x^3+x+1, 'a')
            sage: ProfiniteNumbers(K)
            Profinite Numbers of Number Field in a with defining polynomial x^3 + x + 1
        """
        return "Profinite Numbers of {}".format(self.base())

    def _latex_(self):
        r"""
        Return latex-formatted string representation of ``self``

        EXAMPLES::

            sage: K.<a> = NumberField(x^2+14)
            sage: Khat = ProfiniteNumbers(K)
            sage: latex(Khat)
             \widehat{ \Bold{Q}[a]/(a^{2} + 14) }
        """
        from sage.misc.latex import latex
        return r" \widehat{" + latex(self.base()) + "} "

    def characteristic(self):
        """
        Return the characteristic of this ring, which is zero

        EXAMPLES::

            sage: Qhat = ProfiniteNumbers()
            sage: Qhat.characteristic()
            0
        """
        return ZZ(0)

    def _element_constructor_(self, value, modulus=None, denominator=None):
        """
        Construct the element "(``value`` mod ``modulus``)/``denominator``" of
        ``self``

        EXAMPLES::

            sage: Qhat = ProfiniteNumbers(QQ)
            sage: Qhat(7, 25, 8)
            (7 mod 25)/8
            sage: Qhat(3, 10, 1)
            3 mod 10
            sage: Qhat(3, 10)
            3 mod 10
            sage: Qhat(3/2, 7/2)
            (3 mod 7)/2

        ::

            sage: K.<a> = NumberField(x^3+x+1)
            sage: Khat = ProfiniteNumbers(K)
            sage: Khat(a, 18, a^2+2)
            (a mod (18))/(a^2 + 2)
            sage: Khat(a^2+a, 9*a^2, a^3-1)
            (a^2 + a mod (9*a^2))/(-a - 2)
            sage: Khat(a^2+a, 9*a^2, a^3+1)
            (a^2 + a mod (9*a^2))/(-a)
            sage: Khat(3*a-2, 4*a^2+16, a)
            (-24*a^2 - a + 2 mod (4*a^2 + 16))/a

        The number field which is our base coerces in::

            sage: Qhat(7/9)
            7/9
            sage: Khat(1/(-4*a^2 - 2*a - 6))
            1/(-4*a^2 - 2*a - 6)
            sage: Khat((a-1)/(2*a^2+10))
            1/(-4*a^2 - 2*a - 6)

        As well as the corresponding profinite integer rings::

            sage: Zhat = ProfiniteIntegers(QQ)
            sage: Ohat = ProfiniteIntegers(K)
            sage: Qhat(Zhat(6, 24))
            6 mod 24
            sage: Khat(Ohat(a^2, 3*a+5))
            a^2 mod (3*a + 5)
            sage: Khat(Ohat(a+1, a^2+18*a+20))
            1436*a^2 mod (a^2 + 18*a + 20)
        """
        K = self.base()
        O = K.maximal_order()
        Ohat = ProfiniteIntegers(O)
        if modulus is None and denominator is None:
            if value in K:  # K --> Khat
                d = K(value).denominator()
                return self.element_class(self, Ohat(value*d), d)
            if value in Ohat:  # Ohat --> Khat
                return self.element_class(self, Ohat(value))
        if denominator is None:
            if value in K and modulus in K.ideal_monoid():
                from sage.arith.functions import lcm
                if K is QQ:
                    den = QQ(modulus).denominator()
                else:
                    _, den = K.ideal(modulus).integral_split()
                d = lcm(K(value).denominator(), den)
                return self.element_class(self, Ohat(d*value, d*modulus), d)
        return self.element_class(self, Ohat(value, modulus), denominator)

    def _coerce_map_from_(self, S):
        """
        Return a coerce map ``self`` --> ``S`` (or ``True`` to use
        ``_element_constructor_``) if it exists, ``False`` otherwise.

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

    def is_integral_domain(self, proof=None):
        """
        Return ``False``, indicating that ``self`` is not an integral domain

        EXAMPLES::

            sage: Qhat = ProfiniteNumbers()
            sage: Qhat.is_integral_domain()
            False
        """
        return False

    def number_field(self):
        """
        Return the base number field of ``self``

        EXAMPLES::

            sage: K.<a> = NumberField(x^5+7*x+9)
            sage: Khat = ProfiniteNumbers(K)
            sage: Khat.number_field()
            Number Field in a with defining polynomial x^5 + 7*x + 9
        """
        return self.base()

Qhat = ProfiniteNumbers(QQ)