"""
Profinite Integers of Number Fields

.. TODO::

    Write documentation for this module
"""
load("completions.sage")

from sage.structure.element import CommutativeRingElement

class ProfiniteInteger(CommutativeRingElement):
    """
    Profinite Integer of a Number Field
    
    .. TODO::
    
        Write documentation for this class definition
    
    :automethod:`_repr_`
    """

    def __init__(self, parent, value, modulus):
        r"""
        Construct the profinite integer "``value`` mod ``modulus``"

        Write `O` for the order such that the profinite completion `\hat{O}` of
        `O` is the ring this profinite integer belongs to. Write `K` for the
        ambient number field of `O`.
        
        INPUT:

        - ``parent`` -- `\hat{O}` as an ``ProfiniteIntegers`` object
        - ``value`` -- an element of `O`
        - ``modulus`` -- if `K` is `\QQ`: an integer; else: an ideal of `O`.

        OUTPUT:

        The profinite integer representing the open subset
        ``value`` + ``modulus``*`\hat{O}` of `\hat{O}` (or the point ``value``
        if ``modulus`` is zero).

        EXAMPLES::

            sage: Zhat = ProfiniteIntegers()
            sage: ProfiniteInteger(Zhat, 4, 15)
            4 mod 15
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
        CommutativeRingElement.__init__(self, parent)
        O = parent.base()
        K = parent.number_field()
        if value not in O:
            raise TypeError("value must be an element of {}".format(O))
        self.value = O(value)
        if O is ZZ:
            if modulus not in ZZ:
                raise TypeError("modulus must be an integer")
            if modulus < ZZ(0):
                modulus = -modulus
            self.modulus = ZZ(modulus)
        else:
            if modulus not in K.ideal_monoid() or not modulus.is_integral():
                raise TypeError("modulus must be an ideal of {}".format(O))
            modulus = O.ideal(modulus)
            self.modulus = O.ideal(modulus.gens_reduced())
        self._reduce()

    def _repr_(self):
        """
        Returns a string representation of ``self``

        EXAMPLES::

            sage: Zhat = ProfiniteIntegers()
            sage: Zhat(6, 20)
            6 mod 20
            sage: K.<a> = NumberField(x^2+5)
            sage: Ohat = ProfiniteIntegers(K)
            sage: Ohat(a+1, 60)
            a + 1 mod (60)
            sage: I = K.ideal(3, a+2)
            sage: Ohat(a, I)
            a mod (3, a + 2)
        """
        if self.modulus == 0:
            return repr(self.value)
        modulus = self.modulus
        if self.parent().base() is not ZZ:
            t = len(self.modulus.gens())
            modulus = "("
            for i in range(t):
                modulus += str(self.modulus.gens()[i])
                if i < t-1:
                    modulus += ", "
            modulus += ")"
        return "{} mod {}".format(self.value, modulus)

    def _richcmp_(self, other, op):
        r"""
        Compare ``self`` and ``other`` based on the relation ``op``

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
            return (self.value - other.value) in common_modulus
        if op == op_NE:
            return not self._richcmp_(other, op_EQ)
        raise NotImplementedError("only equality and inequality are implemented")

    def _reduce(self):
        """
        Reduce ``self.value`` modulo ``self.modulus``

        .. TODO::

            Use modulus.small_residue() instead of modulus.reduce()?
            Which one is better for our purposes?

        EXAMPLES::

            sage: Zhat = ProfiniteIntegers()
            sage: a = Zhat(0, 10)
            sage: a.value = 91; a
            91 mod 10
            sage: a._reduce(); a
            1 mod 10

        ::

            sage: K.<a> = NumberField(x^2-6)
            sage: Ohat = ProfiniteIntegers(K)
            sage: b = Ohat(0, a+1)
            sage: b.value = 10*a+11; b
            10*a + 11 mod (a + 1)
            sage: b._reduce(); b
            -a mod (a + 1)

        TESTS::

            sage: b.modulus = K.ideal(0)
            sage: b._reduce(); b
            -a
        """
        if self.parent().base() is ZZ:
            if self.modulus != ZZ(0):
                self.value %= self.modulus
        else:
            self.value = self.modulus.reduce(self.value)

    def _common_modulus(self, other):
        """
        Return the biggest modulus modulo which both ``self`` and ``other`` are
        defined

        For base ring ZZ, this is the GCD of the moduli.
        Otherwise, it is the sum of the moduli.

        INPUT:

        - ``other`` -- a profinite integer with the same parent as ``self``

        EXAMPLES::

            sage: Zhat = ProfiniteIntegers(ZZ)
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
            return gcd(self.modulus, other.modulus)
        return self.modulus + other.modulus

    def _add_(self, other):
        """
        Add ``other`` to ``self`` and return the result

        EXAMPLES::
        
            sage: Zhat = ProfiniteIntegers(ZZ)
            sage: a = Zhat(3, 20)
            sage: b = Zhat(-1, 30)
            sage: a+b
            2 mod 10
            sage: a+3
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
        """
        modulus = self._common_modulus(other)
        value = self.value + other.value
        return self.__class__(self.parent(), value, modulus)
        
    def _sub_(self, other):
        """
        Subtract ``other`` from ``self`` and return the result

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
        """
        modulus = self._common_modulus(other)
        value = self.value - other.value
        return self.__class__(self.parent(), value, modulus)
        
    def _mul_(self, other):
        """
        Multiply ``self`` and ``other`` and return the result

        EXAMPLES::
        
            sage: load("profinite_integers.sage")
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
        """
        I = self.value * other.modulus
        J = self.modulus * other.value
        K = self.modulus * other.modulus
        if self.parent().base() is ZZ:
            modulus = gcd(I, gcd(J, K))
        else:
            modulus = I + J + K
        value = self.value * other.value
        return self.__class__(self.parent(), value, modulus)

    def prime_repr(self):
        """
        TODO fix projection() first, then write this
        """
        from sage.arith.misc import factor
        factorization = factor(self.modulus)
        rep = "("
        for p, e in factorization:
            rep += str(Qp(p)(self.value, e)) + ", "
        rep += "...)"
        return rep

    def projection(self, p):
        r"""
        Return the projection of ``self`` to the ``p``-adic integers

        Write `O` for ``self.parent().base()``, `K` for the ambient number field
        of `O` and `\ZZ_q` for the ring of `q`-adic integers.

        This implements the composite map

        .. MATH::

            \hat{O} \to \prod_q \ZZ_q \to \ZZ_p

        where the product is taken over all (finite) primes `q` of `O`, the
        first map is the natural isomorphism induced by the Chinese Remainder
        Theorem and the second map is the natural projection to the `p`th
        coordinate.

        INPUT:

        - ``p`` -- a prime of `O`; if `O` is `\ZZ` this means a prime number,
          else this means a prime ideal of `O`

        OUTPUT:

        The image of ``self`` under the map described above, as an element of
        the field returned by :func:`completion <completions.sage>`.

        .. TODO::
            
            Get output in the right field ((a finite exteions of) Qp).
            Some initial attempt is written down already (commented out).
            But there for non-trivial number fields, we calculate the precision
            of this projection wrong...

        EXAMPLES::

            TODO
        """
        K = self.parent().number_field()
        if K is QQ:
            if p not in ZZ or not p.is_prime():
                raise ValueError("p should be a prime number")
        else:
            if p not in K.ideal_monoid() or not p.is_prime():
                raise ValueError("p should be a prime ideal of {}".format(K))
        #prec = self.modulus.valuation(p)
        #Kp, phi = completion(K, p, prec) # phi is the natural embedding K --> Kp
        #return phi(self.value)
        O = self.parent().base()
        I = p^(self.modulus.valuation(p))
        R = O.quotient(I, 'a')
        return R(self.value)


from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.ring import CommutativeAlgebra
class ProfiniteIntegers(UniqueRepresentation, CommutativeAlgebra):
    Element = ProfiniteInteger

    def __init__(self, R=ZZ):
        """
        Construct a profinite completion of a ring of integers

        INPUT:

        - ``R`` -- a number field or a maximal order in a number field
          (default: ``ZZ``)

        OUTPUT:

        The profinite completion of the maximal order in the ambient field of
        `R`.

        EXAMPLES::

            sage: ProfiniteIntegers()
            Profinite Completion of Integer Ring
            sage: K.<a> = NumberField(x^5-3*x+1)
            sage: ProfiniteIntegers(K)
            Profinite Completion of Maximal Order in Number Field in a with defining polynomial x^5 - 3*x + 1
            sage: ProfiniteIntegers("not a number field")
            Traceback (most recent call last):
            ...
            TypeError: R should be (the maximal order of) a number field
        """
        try:
            if R is ZZ or R is QQ:
                O = ZZ
            else:
                K = R
                if not is_field(R):
                    K = R.ambient()
                if not K.absolute_degree() in ZZ:
                    raise TypeError("R should be (the maximal order of) a number field")
                O = K.maximal_order()
        except AttributeError:
            raise TypeError("R should be (the maximal order of) a number field")
        CommutativeAlgebra.__init__(self, O)

    def _repr_(self):
        """
        Return a string representation of ``self``

        EXAMPLES::

            sage: K = NumberField(x^3+x+1, 'a')
            sage: ProfiniteIntegers(K)
            Profinite Completion of Maximal Order in Number Field in a with defining polynomial x^3 + x + 1
        """
        return "Profinite Completion of {}".format(self.base())

    def _latex_(self):
        r"""
        Return latex-formatted string representation of ``self``

        EXAMPLES::

            sage: K.<a> = NumberField(x^2+14)
            sage: Ohat = ProfiniteIntegers(K)
            sage: latex(Ohat)
             \widehat{ \mathcal{O}_{ \Bold{Q}[a]/(a^{2} + 14) } }
        """
        K = self.number_field()
        return r" \widehat{ \mathcal{O}_{" + latex(K) + "} } "

    def number_field(self):
        """
        Return the ambient number field of our base
        
        EXAMPLES:

            sage: Zhat = ProfiniteIntegers()
            sage: Zhat.number_field()
            Rational Field
            sage: K.<a> = NumberField(x^2+x-7)
            sage: Ohat = ProfiniteIntegers(K)
            sage: Ohat.number_field()
            Number Field in a with defining polynomial x^2 + x - 7
        """
        if self.base() is ZZ:
            return QQ
        return self.base().ambient()

    def characteristic(self):
        """
        Return the characteristic of this ring, which is zero

        EXAMPLES::

            sage: Zhat = ProfiniteIntegers()
            sage: Zhat.characteristic()
            0
        """
        return ZZ(0)

    def _element_constructor_(self, value, modulus=None):
        """
        Construct an element "``value`` mod ``modulus``" of ``self``

        EXAMPLES::

            sage: Zhat = ProfiniteIntegers()
            sage: Zhat(-4, 17)
            13 mod 17
            sage: K.<a> = NumberField(x^3-5*x^2+1)
            sage: Ohat = ProfiniteIntegers(K)
            sage: Ohat(a^2+1)
            a^2 + 1
            sage: Ohat(a^2+1, 4*a)
            a^2 + 1 mod (4*a)
        """
        if modulus is None:
            modulus = ZZ(0) if self.base() is ZZ else self.base().ideal(0)
        return self.element_class(self, value, modulus)

    def _coerce_map_from_(self, S):
        """
        Return a coerce map ``self`` --> ``S`` (or ``True`` to use
        ``_element_constructor_``) if it exists, ``False`` otherwise.

        .. TODO::

            This does not work very nice right now...
            We would like to make whole K coerce into self, so that things like
            ``a*Ohat(a^2+1, 120)`` or ``a^3+1 == Ohat(a^3+1, 100)`` work.
            They do not right now, since the parent of a^3+1 is not O, but K...
            ``O(a)*Ohat(a^2+1, 120)`` does work, but yeah...

            But could also go for an output in Khat instead of Ohat of course

        EXAMPLES::

            sage: K.<a> = NumberField(x^3+79)
            sage: Ohat = ProfiniteIntegers(K)
            sage: Ohat._coerce_map_from_(K)
            False
            sage: Ohat._coerce_map_from_(K.maximal_order())
            True
        """
        if S is self.base() or S is ZZ:
            return True
        return False

    def is_integral_domain(self, proof=None):
        """
        Return ``False``, indicating that ``self`` is not an integral domain

        EXAMPLES::

            sage: Zhat = ProfiniteIntegers()
            sage: Zhat.is_integral_domain()
            False
        """
        return False

