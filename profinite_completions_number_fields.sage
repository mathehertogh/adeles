


class ProfiniteCompletionNumberFieldElement(CommutativeRingElement):
    # TODO reduce value modulo the modulus

    def __init__(self, parent, value, modulus):
        """
        Construct an element of the profinite completion

        EXAMPLES::

            sage: Khat = ProfiniteCompletionNumberField(K)
            sage: K.<a> = NumberField(x^2+7)
            sage: Khat = ProfiniteCompletionNumberField(K)
            sage: ProfiniteCompletionNumberFieldElement(Khat, a+1, K.ideal(5*a))
            a + 1 mod Fractional ideal (5*a)

        TESTS::

            sage: Khat(3, 0)
            3 mod Ideal (0) of Number Field in a with defining polynomial x^2 + 7
            sage: ProfiniteCompletionNumberFieldElement(Khat, None, K.ideal(5*a))
            Traceback (most recent call last):
            ...
            TypeError: value must be an element of Number Field in a with defining polynomial x^2 + 7
            sage: ProfiniteCompletionNumberFieldElement(Khat, a+1, None)
            Traceback (most recent call last):
            ...
            TypeError: modulus must be an ideal of Number Field in a with defining polynomial x^2 + 7
        """
        K = parent.number_field
        if value not in K:
            raise TypeError("value must be an element of {}".format(K))
        if modulus not in K.ideal_monoid():
            raise TypeError("modulus must be an ideal of {}".format(K))
        self.value = K(value)
        self.modulus = K.ideal(modulus)
        CommutativeRingElement.__init__(self, parent)

    def _repr_(self):
        """
        Returns a string representation of ``self``

        EXAMPLES::

            sage: Khat(3*a, 121)
            3*a mod Fractional ideal (121)
        """
        return "{} mod {}".format(self.value, self.modulus)

    def _richcmp_(self, other, op):
        r"""
        Compare ``self`` and ``other`` based on the relation ``op``

        We only implement equality and non-equality.

        Two elements `(x, I)` and `(y, J)` are considered equal if and only if
        `x \equiv y \mod (I+J)`.

        EXAMPLES::

            sage: K.<a> = NumberField(x^4+17)
            sage: Khat = ProfiniteCompletion(K)
            sage: b = Khat(a^2+1, 1000)
            sage: c = Khat(a^2+11, a*10)
            sage: d = Khat(a^3+a-5, 170)
            sage: b == b
            True
            sage: b == c
            True
            sage: b == d
            False
            sage: c == d
            False
            sage: b != c
            False
            sage: c != d
            True

        TESTS::

            sage: b < d
            Traceback (most recent call last):
            ...
            NotImplementedError: only equality and inequality are implemented
        """
        from sage.structure.richcmp import op_EQ, op_NE
        if op == op_EQ:
            common_modulus = self.modulus + other.modulus
            return (self.value - other.value) in common_modulus
        if op == op_NE:
            return not self._richcmp_(other, op_EQ)
        raise NotImplementedError("only equality and inequality are implemented")

    def _add_(self, other):
        val = self.value + other.value
        a, b = self.modulus.numerator(), self.modulus.denominator()
        c, d = other.modulus.numerator(), other.modulus.denominator()
        mod = gcd(a*d, b*c) / (b*d)
        return self.__class__(self.parent(), val, mod)
        
    def _sub_(self, other):
        val = self.value - other.value
        a, b = self.modulus.numerator(), self.modulus.denominator()
        c, d = other.modulus.numerator(), other.modulus.denominator()
        mod = gcd(a*d, b*c) / (b*d)
        return self.__class__(self.parent(), val, mod)
        
    def _mul_(self, other):
        val = self.value * other.value
        x_n, x_d = self.value.numerator(), self.value.denominator()
        m_n, m_d = self.modulus.numerator(), self.modulus.denominator()
        y_n, y_d = other.value.numerator(), other.value.denominator()
        n_n, n_d = other.modulus.numerator(), other.modulus.denominator()
        den = x_d * m_d * y_d * n_d
        num = gcd(x_n*n_n*y_d*m_d, gcd(y_n*m_n*x_d*n_d, m_n*n_n*x_d*y_d)) # TODO: pull factor m_n*x_d out of second gcd
        mod = num / den
        return self.__class__(self.parent(), val, mod)

    def prime_repr(self):
        from sage.arith.misc import factor
        from sage.rings.padics.factory import Qp
        factorization = factor(self.modulus)
        rep = "("
        for p, e in factorization:
            rep += str(Qp(p)(self.value, e)) + ", "
        rep += "...)"
        return rep

    def to_Qp(self, p):
        if p not in ZZ or not p.is_prime():
            raise ValueError("p should be a rational prime number")
        return Qp(p)(self.value, self.modulus.valuation(p))


from sage.structure.unique_representation import UniqueRepresentation
class ProfiniteCompletionNumberField(UniqueRepresentation, CommutativeRing):
    Element = ProfiniteCompletionNumberFieldElement

    def __init__(self, K):
        """
        Construct the profinite completion of the number field ``K``

        EXAMPLES::

            sage: ProfiniteCompletionNumberField(QQ)
            Profinite Completion of Rational Field
            sage: K.<a> = NumberField(x^5-3*x+1)
            sage: ProfiniteCompletionNumberField(K)
            Profinite Completion of Number Field in a with defining polynomial x^5 - 3*x + 1
            sage: ProfiniteCompletionNumberField("not a number field")
            Traceback (most recent call last):
            ...
            TypeError: K should be a number field
        """
        try:
            if not is_field(K) or not K.absolute_degree() in ZZ:
                raise TypeError("K should be a number field")
        except AttributeError:
            raise TypeError("K should be a number field")
        self.number_field = K
        # We choose K as our base ring.
        CommutativeRing.__init__(self, K)

    def _repr_(self):
        """
        Return a string representation of ``self``

        EXAMPLES::

            sage: ProfiniteCompletionNumberField(K)
            Profinite Completion of Number Field in a with defining polynomial x^4 + x^2 + 3
        """
        return "Profinite Completion of {}".format(self.number_field)

    def _latex_(self):
        return r"\hat{K}" # TODO

    def characteristic(self):
        """
        Return the characteristic of this ring, which is zero

        EXAMPLES::

            sage: Qhat = ProfiniteCompletionNumberField(QQ)
            sage: Qhat.characteristic()
            0
        """
        return 0

    def _element_constructor_(self, value, modulus=None):
        """
        Construct an element "``value`` mod ``modulus``" of ``self``
        """
        if modulus is None:
            modulus = 0
        return self.element_class(self, value, modulus)

    def _coerce_map_from_(self, S):
        if S is self.number_field:
            return True
        return False

    def is_integral_domain(self, proof=None):
        return False

    
def ProfiniteCompletion(K):
    """
    Return the profinite completion of the number field ``K``
    """
    return ProfiniteCompletionNumberField(K)