


class PfInt(CommutativeRingElement):
    def __init__(self, parent, val, mod=None):
        if mod is None:
            mod = 0
        if val not in ZZ or mod not in ZZ:
            raise ValueError("value and modulus must be integers")
        if mod < 0:
            mod = -mod
        if mod != 0:
            val = val % mod
        self.value = ZZ(val)
        self.modulus = ZZ(mod)
        CommutativeRingElement.__init__(self, parent)

    def _repr_(self):
        mode = ProfiniteIntegers.print_mode
        if mode == "mod":
            return "{} mod {}".format(self.value, self.modulus)
        if mode == "factorial":
            return self.factorial_repr()
        if mode == "prime":
            return self.prime_repr()
        raise ValueError("Error: unknown print_mode '{}'".format(mode))

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
        mod = gcd(self.modulus, other.modulus)
        val = self.value + other.value
        return self.__class__(self.parent(), val, mod)
        
    def _sub_(self, other):
        mod = gcd(self.modulus, other.modulus)
        val = self.value - other.value
        return self.__class__(self.parent(), val, mod)
        
    def _mul_(self, other):
        mod = gcd(self.value * other.modulus, self.modulus * other.value)
        mod = gcd(mod, self.modulus * other.modulus)
        val = self.value * other.value
        return self.__class__(self.parent(), val, mod)

    def factorial_digits(self):
        val = self.value
        # Handle exact integer case.
        if self.modulus == ZZ(0):
            raise ValueError("not implemented yet!")
        # Calculate precision, which is the largest k such that k! divides
        # self.modulus.
        prec = ZZ(1)
        prec_fac = ZZ(1)
        while self.modulus % (prec_fac*(prec+1)) == 0:
            prec += 1
            prec_fac *= prec
        # Calculate the projection of val to ZZ/(prec!)ZZ.
        val = val % prec_fac
        # Calculate the factorial digits of val.
        digits = {}
        k_fac = prec_fac // prec
        for k in range(prec-1, ZZ(0), ZZ(-1)):
            digit, val = divmod(val, k_fac)
            digits[k] = digit
            k_fac //= k
        return digits, prec

    def factorial_repr(self):
        digits, prec = self.factorial_digits()
        rep = "("
        for k in range(ZZ(1), prec):
            rep += "{}|".format(digits[k])
        rep += "...)"
        return rep

    def prime_repr(self):
        from sage.arith.misc import factor
        from sage.rings.padics.factory import Qp
        factorization = factor(self.modulus)
        rep = "("
        for p, e in factorization:
            rep += str(Qp(p)(self.value, e)) + ", "
        rep += "...)"
        return rep

    def visual(self):
        """
        Return an interval within the unitinterval `[0,1]` representing ``self``
        """
        digits, prec = self.factorial_digits()
        v = ZZ(0)
        k_fac = ZZ(1)
        for k in range(ZZ(1), prec):
            k_fac *= k + 1
            v += digits[k] / k_fac
        w = v + ZZ(1)/k_fac
        return v, w

from sage.structure.unique_representation import UniqueRepresentation
class PfInts(UniqueRepresentation, CommutativeRing):
    Element = PfInt
    print_mode = "mod"

    def __init__(self):
        # We choose ZZ as our base ring.
        CommutativeRing.__init__(self, ZZ)

    def _repr_(self):
        return "Profinite Integer Ring"

    def _latex_(self):
        return r"\hat{\ZZ}"

    def set_print_mode(self, mode):
        modes = ["mod", "prime", "factorial"]
        if mode in modes:
            ProfiniteRationals.print_mode = mode
        else:
            print("mode '{}' not recognized; choose one of the following:".format(mode))
            print(modes)

    def characteristic(self):
        return 0

    def _element_constructor_(self, val, mod=None):
        return self.element_class(self, val, mod)

    def _coerce_map_from_(self, S):
        if S is ZZ:
            return True
        return False

    def is_integral_domain(self, proof=None):
        return False

Zhat = PfInts()
    
