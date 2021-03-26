""""
Profinite Integers
"""
from sage.rings.integer_ring import ZZ
from sage.categories.rings import Rings
from sage.rings.ring import Ring
from sage.structure.element import RingElement
from sage.arith.functions import lcm
from sage.arith.misc import gcd

class ProfiniteInteger(RingElement):
    r"""
    v mod m represents the open subset v + m*Zhat

    .. TODO::

        implement __contains__() (hence the ``in`` keyword)
    """
    
    def __init__(self, val, mod, print_mode="mod"):
        if val not in ZZ:
            raise ValueError("value must be an integer (provided: {})".format(val))
        if mod not in ZZ:
            raise ValueError("modulus must be an integer (provided: {})".format(val) )
        self.modulus = ZZ(mod)
        self.value = ZZ(val)
        if self.modulus != ZZ(0):
            self.value = self.value % self.modulus

    def parent(self):
        return ProfiniteIntegers() # TODO is this really the way this should work?
    
    def _repr_(self):
        mode = ProfiniteIntegers.print_mode
        if mode == "mod":
            return "{} mod {}".format(self.value, self.modulus)
        if mode == "factorial":
            return self.factorial_repr()
        if mode == "prime":
            return self.prime_repr()
        raise ValueError("Error: don't know how print_mode '{}' works".format(mode))
    
    def _richcmp_(self, other, op):
        from sage.structure.richcmp import richcmp
        
    def _add_(self, other):
        mod = gcd(self.modulus, other.modulus)
        val = self.value + other.value
        if mod != ZZ(0):
            val = val % mod
        return self.__class__(val, mod)
        
    def _sub_(self, other):
        mod = gcd(self.modulus, other.modulus)
        val = self.value - other.value
        if mod != ZZ(0):
            val = val % mod
        return self.__class__(val, mod)
        
    def _mul_(self, other):
        mod = gcd(self.value * other.modulus, self.modulus * other.value)
        mod = gcd(mod, self.modulus * other.modulus)
        val = self.value * other.value
        if mod != ZZ(0):
            val = val % mod
        return self.__class__(val, mod)

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




class ProfiniteIntegers(Ring):
    Element = ProfiniteInteger
    print_mode = "mod"
    def __init__(self, base=ZZ, category=Rings(), print_mode="mod"):
        # We must provide a base... What should base be? What is base?
        Ring.__init__(self, base, category=category) 
        ProfiniteIntegers.print_mode = print_mode
    
    def _repr_(self):
        return "Ring of profinite integers"
    
    def _element_constructor_(self, *args, **kwds):
        if len(args) < 2:
            return self.element_class(args[0], 0)
        return self.element_class(*args, **kwds)
    
    def _coerce_map_from_(self, ZZ):
        # use the above conversion ("_element_constructor_") as coercion from ZZ
        return True

    def set_print_mode(self, mode):
        modes = ["mod", "prime", "factorial"]
        if mode in modes:
            ProfiniteRationals.print_mode = mode
        else:
            print("mode '{}' not recognized; choose one of the following:".format(mode))
            print(modes)

    def is_finite(self):
        return False

Zhat = ProfiniteIntegers()
