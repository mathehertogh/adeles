


class PfRat(CommutativeRingElement):
    def __init__(self, parent, val, mod=None):
        if mod is None:
            mod = 0
        if val not in QQ or mod not in QQ:
            raise ValueError("value and modulus must be rationals")
        if mod < 0:
            mod = -mod
        self.value = QQ(val)
        self.modulus = QQ(mod)
        self._reduce_value()
        CommutativeRingElement.__init__(self, parent)

    def _repr_(self):
        mode = PfRats.print_mode
        if mode == "mod":
            return "{} mod {}".format(self.value, self.modulus)
        if mode == "factorial":
            return self.factorial_repr()
        if mode == "prime":
            return self.prime_repr()
        raise ValueError("Error: unknown print_mode '{}'".format(mode))

    def _reduce_value(self):
        """Make sure 0 <= value < modulus when modulus!=0"""
        if self.modulus != ZZ(0):
            quo = (self.value/self.modulus).floor()
            rem = self.value - self.modulus*quo
            self.value = rem

    def _richcmp_(self, other, op):
        from sage.structure.richcmp import op_EQ, op_NE
        if op == op_EQ:
            common_modulus = gcd(self.modulus, other.modulus)
            if common_modulus == 0:
                return self.value == other.value
            return self.value % common_modulus == other.value % common_modulus
        if op == op_NE:
            return not self._richcmp_(other, op_EQ)
        raise NotImplementedError()

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
class PfRats(UniqueRepresentation, CommutativeRing):
    Element = PfRat
    print_mode = "mod"

    def __init__(self):
        # We choose QQ as our base ring.
        CommutativeRing.__init__(self, QQ)

    def _repr_(self):
        return "Profinite Rational Ring"

    def _latex_(self):
        return r"\hat{\QQ}"

    def set_print_mode(self, mode):
        modes = ["mod", "prime", "factorial"]
        if mode in modes:
            PfRats.print_mode = mode
        else:
            print("mode '{}' not recognized; choose one of the following:".format(mode))
            print(modes)

    def characteristic(self):
        return 0

    def _element_constructor_(self, val, mod=None):
        return self.element_class(self, val, mod)

    def _coerce_map_from_(self, S):
        if S is QQ:
            return True
        return False

    def is_integral_domain(self, proof=None):
        return False

Qhat = PfRats()
    
