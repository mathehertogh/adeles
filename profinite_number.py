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

    def __init__(self, parent, numerator, denominator):
        r"""
        Construct the profinite number ``numerator / denominator``
        
        INPUT:

        - ``parent`` -- A ring of profinite numbers over some number field `K`
        - ``numerator`` -- an element of ``ProfiniteIntegers(K)``
        - ``denominator`` -- an element of ``K.maximal_order()``

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
            TypeError: numerator should be a profinite integer
            sage: ProfiniteNumber(Qhat, 3, a)
            Traceback (most recent call last):
            ...
            TypeError: denominator should be a non-zero integer
            sage: ProfiniteNumber(Khat, Ohat(a), 0)
            Traceback (most recent call last):
            ...
            TypeError: denominator should be a non-zero integer
        """
        CommutativeAlgebraElement.__init__(self, parent)
        Ohat = ProfiniteIntegers(parent.base())
        if numerator not in Ohat:
            raise TypeError("numerator should be a profinite integer")
        if denominator not in ZZ or denominator == 0:
            raise TypeError("denominator should be a non-zero integer")
        self._numerator = Ohat(numerator)
        self._denominator = ZZ(denominator)
        self._reduce()

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
            b = gcd([num.value(), num.modulus(), self.denominator()])
        else:
            O = K.maximal_order()
            # The greatest common divisor of ideals is their sum.
            I = O.ideal(num.value()) + num.modulus() + O.ideal(self.denominator())
            # The element ``(1/I).gens_two()[0].abs()`` is the positive
            # generator of `1/I \cap \QQ`. Write this element as `a/b` with `a`
            # and `b` coprime and both positive.
            # We have `(a/b)I \subset O`, i.e. `aI \subset bO`, i.e. `bO`
            # divides `aI`. As `a` and `b` are coprime this implies `bO` divides
            # `I`. In turn this is equivalent to `(1/b)I \subset O` and so
            # `1/b \in 1/I \cap \QQ`.
            # As `a/b` was the generator of `1/I \cap \QQ` we must have
            # `a/b = 1/b`, i.e. `a = 1` (note that `a/b` and `b` are positive).
            # This means `1/b` is the largest integer satisfying `I/b \subset O`
            # i.e. `I \subset bO`.
            # We conclude: `b` is the largest integer dividing `I`.
            b = ZZ(1/(1/I).gens_two()[0])
        # Avoid division of the zero-ideal, which is not defined in SageMath.
        modulus = num.modulus()/b if num.modulus() != 0 else num.modulus()
        self._numerator = Ohat(num.value()//b, modulus)
        self._denominator = self.denominator() // b

    def _repr_(self):
        """
        Returns a string representation of ``self``

        EXAMPLES::

            sage: Qhat(1/2, 97/5)
            1/2 mod 97/5

        ::

            sage: K.<a> = NumberField(x^2+5)
            sage: Ohat = ProfiniteIntegers(K)
            sage: Khat = ProfiniteNumbers(K)
            sage: n = Ohat(2*a, K.ideal(4, 2*a+2))
            sage: Khat(n, 6)
            1/3*a mod (2/3, 1/3*a + 1/3)
        """
        K = self.parent().base()
        if K is QQ:
            modulus = self.modulus()
        else:
            modulus = self.modulus().gens()
            if len(modulus) == 1:
                modulus = "(" + str(modulus[0]) + ")"
        return str(self.value()) + " mod " + str(modulus)

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
        """
        Subtract ``other`` from ``self`` and return the result

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
        # Division of the zero ideal is not allowed, so we handle that case
        # separately.
        if self.numerator().modulus() == 0:
            return self.numerator().modulus()
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

    def _element_constructor_(self, x, y=None):
        """
        Construct the element "(``value`` mod ``modulus``)/``denominator``" of
        ``self``

        Write `K` for our base number field.

        INPUT:

        One of the following options:

        - ``x`` an element of `K`, ``y`` a fractional ideal of `K` (default:
          zero); gives ``x mod y``
        - ``x`` a profinite `K`-integer, ``y`` a non-zero integer (default:
          ``1``); gives ``x/y``

        EXAMPLES::

            Zhat = ProfiniteIntegers(QQ)
            Qhat = ProfiniteNumbers(QQ)
            n = Zhat(3, 10)
            Qhat(n)
            3 mod 10
            sage: Qhat(n, 7)
            3/7 mod 10/7
            sage: Qhat(3/2)
            3/2 mod 0
            sage: Qhat(3/2, 11/10)
            2/5 mod 11/10

        ::

            sage: K.<a> = NumberField(x^3+x+1)
            sage: Ohat = ProfiniteIntegers(K)
            sage: Khat = ProfiniteNumbers(K)
            sage: n = Ohat (a^2, a^2+a+1)
            sage: Khat(n)
            a^2 mod (a^2 + a + 1)
            sage: Khat(n, 17)
            1/17*a^2 mod (1/17*a^2 + 1/17*a + 1/17)
            sage: Khat(a/3)
            1/3*a mod (0)
            sage: Khat(a/3, (a^2-2)/9)
            8/9*a^2 mod (1/9*a^2 - 2/9)

        TESTS::

            sage: Qhat('bla')
            Traceback (most recent call last):
            ...
            TypeError: Can't construct profinite number from ('bla', None)
        """
        K = self.base()
        Ohat = ProfiniteIntegers(K)

        if x in K:
            value = K(x)

            if y is None:
                modulus = QQ(0) if K is QQ else K.ideal(0)
            elif y in K.ideal_monoid():
                modulus = QQ(y) if K is QQ else K.ideal(y)
            else:
                raise TypeError("second parameter should be a fractional ideal specifying the modulus")

            from sage.arith.functions import lcm
            modulus_den = modulus.denominator() if K is QQ else modulus.integral_split()[1]
            den = lcm(value.denominator(), modulus_den)
            return self.element_class(self, Ohat(den*value, den*modulus), den)

        if x in Ohat:
            num = Ohat(x)

            if y is None:
                den = ZZ(1)
            elif y in ZZ and y != 0:
                den = ZZ(y)
            else:
                raise TypeError("second parameter should be a non-zero integer specifying the denominator")

            if den < 0:
                num, den = -num, -den

            return self.element_class(self, num, den)

        raise TypeError("Can't construct profinite number from {}".format((x,y)))
        # if modulus is None and denominator is None:
        #     if value in K:  # K --> Khat
        #         d = K(value).denominator()
        #         return self.element_class(self, Ohat(value*d), d)
        #     if value in Ohat:  # Ohat --> Khat
        #         return self.element_class(self, Ohat(value))
        # if denominator is None:
        #     if value in K and modulus in K.ideal_monoid():
        #         from sage.arith.functions import lcm
        #         if K is QQ:
        #             den = QQ(modulus).denominator()
        #         else:
        #             _, den = K.ideal(modulus).integral_split()
        #         d = lcm(K(value).denominator(), den)
        #         return self.element_class(self, Ohat(d*value, d*modulus), d)
        # return self.element_class(self, Ohat(value, modulus), denominator)

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

    def random_element(self):
        """
        Return a random profinite number
        """
        K = self.base()

        value = K.random_element()
        
        if K is QQ:
            modulus = QQ.random_element()
            while modulus == -1:
                modulus = QQ.random_element()
        else:
            from sage.misc.prandom import choice
            from sage.arith.misc import factor
            n_primes = ZZ.random_element()
            while n_primes < 0:
                n_primes = ZZ.random_element()

            modulus = K.ideal(1)
            for i in range(n_primes):
                I = K.ideal(K.random_element())
                while I.is_one() or I.is_zero():
                    I = K.ideal(K.random_element())

                # We pick individual primes, so that we can end up with
                # non-principal ideals as modulus.
                p, e = choice(factor(I))  
                modulus *= p**e

        multiplier = ZZ.random_element()
        while multiplier <= 0:
            multiplier = ZZ.random_element()

        modulus *= multiplier
            
        return self(value, modulus)

Qhat = ProfiniteNumbers(QQ)