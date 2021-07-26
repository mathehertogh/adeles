
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
    """
    Return ``true`` if and only if ``p`` is a finite prime of the number field
    ``K``
    """
    if K is QQ:
        return p in Primes()
    else:
        return p in K.ideal_monoid() and K.ideal(p).is_prime()


class MultiplicativePAdic(MultiplicativeGroupElement):
    r"""
    A multiplicative `P`-adic number, for `P` a finite prime of a number field
    """

    def __init__(self, parent, center, prec):
        r"""
        TODO
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
            raise ValueError("prec must lie in NN or be Infinity")

    def _reduce(self):
        """
        TODO

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

    def _repr_(self):
        """
        Return a representation of this multiplicative `P`-adic.
        """
        if self.prec() is Infinity:
            return str(self.center())

        need_parentheses = 1 < len([c for c in self.center().list() if not c.is_zero()])
        if need_parentheses:
            return "({}) * U({})".format(self.center(), self.prec())
        else:
            return "{} * U({})".format(self.center(), self.prec())

    def center(self):
        return self._center

    def prec(self):
        return self._prec

    def precision(self):
        return self._prec

    def valuation(self):
        """
        Return the valuation of this multipicative p-adic
        """
        P = self.parent().prime()
        return self.center().valuation(P)

    def _mul_(self, other):
        """
        TODO
        """
        center = self.center() * other.center()
        prec = min(self.prec(), other.prec())
        return self.__class__(self.parent(), center, prec)

    def inverse(self):
        return self.__class__(self.parent(), 1/self.center(), self.prec())

    def _div_(self, other):
        """
        TODO:

        sage: x = M(7, 3)
        sage: 1/x
        ERROR
        """
        return self * other.inverse()

    def represents(self, x):
        """
        Return whether or not this multipicative p-adic represents ``x``
        """
        if self.prec() is Infinity:
            return x == self.center()

        P = self.parent().prime()
        if self.prec() == 0:
            return x.valuation(P) == self.valuation()

        # self.prec() is a positive integer, so our represented subset is
        # self.center() * (1 + P^self.prec() O_P)
        # for O_P the ring of P-adic integers.
        return (x / self.center() - 1).valuation(P) >= self.prec()

    def _equals(self, other):
        """
        Return whether or not this multipicative p-adic equals ``other``

        Equality here means: have non-empty intersection of represented subsets.
        """
        if self.prec() > other.prec():
            self, other = other, self

        return other.represents(self.center())

    def _richcmp_(self, other, op):
        """
        Return the result of operator ``op`` applied to ``self`` and ``other``

        Only equality and inequality are implented.
        """
        from sage.structure.richcmp import op_EQ, op_NE
        if op == op_EQ:
            return self._equals(other)
        if op == op_NE:
            return not self._equals(other)
        raise NotImplementedError()


class MultiplicativePAdics(UniqueRepresentation, Group):
    Element = MultiplicativePAdic

    def __init__(self, K, prime):
        from sage.misc.functional import is_field
        if not is_field(K) or not K.absolute_degree() in ZZ:
            raise TypeError("K should be a number field")
        self._number_field = K

        if not is_finite_prime(prime, K):
            raise TypeError("prime should be a finite prime of the number field")
        self._prime = ZZ(prime) if K is QQ else K.ideal(prime)

        Group.__init__(self, category=Groups().Commutative())

    def _repr_(self):
        """
        Return a string representation of this group

        EXAMPLES::

            sage: K.<a> = NumberField(x^2-3)
            sage: Ideles(K)
            Idele Group of Number Field in a with defining polynomial x^2 - 3
        """
        K, P = self.number_field(), self.prime_name()
        return "Group of multiplicative {}-adics of {}".format(P, K)

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

    def prime(self):
        """
        Return the finite prime `P` for which this group is the unit group of
        the completion of our base number field at `P`
        """
        return self._prime

    def prime_name(self):
        if self.number_field() is QQ:
            return str(self.prime())
        else:
            return str(self.prime().gens_two())

    def _element_constructor_(self, center, prec=Infinity):
        """
        TODO
        """
        return self.element_class(self, center, prec)

    def _coerce_map_from_(self, domain):
        if self.number_field().has_coerce_map_from(domain):
            return True
        return False


def multPAdic(prime, data):
    """
    Shortcut to create multiplicative p-adics
    """
    if prime in Primes():
        K = QQ
    else:
        K = prime.number_field()
    parent = MultiplicativePAdics(K, prime)

    if data in parent:
        return parent(data)

    try:
        center, prec = data
        return parent(center, prec)
    except TypeError:
        raise TypeError("Can't construct a multiplicative p-adic from {}".format(data))

