""""
Profinite Rationals
"""
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.categories.rings import Rings
from sage.rings.ring import Ring
from sage.structure.element import RingElement
from sage.arith.functions import lcm
from sage.arith.misc import gcd

class ProfiniteRational(RingElement):
    r"""
    v mod m represents the open subset v + m*Zhat

    .. TODO::

        implement __contains__() (hence the ``in`` keyword)

    .. TODO::
    
        garanderen dat altijd geldt: 0 <= value < modulus of modulus=0? Of niet?

    .. TODO::

        latex macros \Zhat and \Qhat maken?

    """
    
    def __init__(self, val, mod, print_mode="mod"):
        if val not in QQ:
            raise ValueError("value must be a rational number (provided: {})".format(val))
        if mod not in QQ:
            raise ValueError("modulus must be a rational number (provided: {})".format(val) )
        self.value = QQ(val)
        self.modulus = QQ(mod)

    def parent(self):
        return ProfiniteRationals() # TODO is this really the way this should work?
    
    def _repr_(self):
        mode = ProfiniteRationals.print_mode
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
        val = self.value + other.value
        a, b = self.modulus.numerator(), self.modulus.denominator()
        c, d = other.modulus.numerator(), other.modulus.denominator()
        mod = gcd(a*d, b*c) / (b*d)
        return self.__class__(val, mod)
        
    def _sub_(self, other):
        val = self.value - other.value
        a, b = self.modulus.numerator(), self.modulus.denominator()
        c, d = other.modulus.numerator(), other.modulus.denominator()
        mod = gcd(a*d, b*c) / (b*d)
        return self.__class__(val, mod)
        
    def _mul_(self, other):
        val = self.value * other.value
        x_n, x_d = self.value.numerator(), self.value.denominator()
        m_n, m_d = self.modulus.numerator(), self.modulus.denominator()
        y_n, y_d = other.value.numerator(), other.value.denominator()
        n_n, n_d = other.modulus.numerator(), other.modulus.denominator()
        den = x_d * m_d * y_d * n_d
        num = gcd(x_n*n_n*y_d*m_d, gcd(y_n*m_n*x_d*n_d, m_n*n_n*x_d*y_d)) # TODO: pull factor m_n*x_d out of second gcd
        mod = num / den
        return self.__class__(val, mod)

    def denominator(self):
        return lcm(self.value.denominator(), self.modulus.denominator())

    def integral_part(self):
        d = self.denominator()
        val = self.value.numerator() * (d / self.value.denominator())
        mod = self.modulus.numerator() * (d / self.modulus.denominator())
        return self.__class__(val, mod)

    def decompose(self):
        return 1/self.denominator(), self.integral_part()

    def factorial_digits(self):
        n = self.integral_part()
        mod = n.modulus.numerator()
        val = n.value.numerator()
        # handle exact integer case:
        if mod == 0:
            raise ValueError("not implemented yet!")
        # calculate precision, which is the largest k such that k! divides mod:
        prec = ZZ(1)
        prec_fac = ZZ(1)
        while mod % (prec_fac*(prec+1)) == 0:
            prec += 1
            prec_fac *= prec
        # calculate the projection of val to ZZ/(prec!)ZZ:
        val = val % prec_fac
        # calculate the factorial digits of val:
        digits = {}
        k_fac = prec_fac // prec
        for k in range(prec-1, ZZ(0), ZZ(-1)):
            digit, val = divmod(val, k_fac)
            digits[k] = digit
            k_fac //= k
        return digits, prec

    def factorial_repr(self):
        digits, prec = self.factorial_digits()
        # build a string to print:
        rep = "1/{} * (".format(self.denominator())
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
        #returns an interval within the unitinterval [0,1]
        digits, prec = self.factorial_digits()
        v = ZZ(0)
        k_fac = ZZ(1)
        for k in range(ZZ(1), prec):
            k_fac *= k + 1
            v += digits[k] / k_fac
        w = v + ZZ(1)/k_fac
        return v, w




class ProfiniteRationals(Ring):
    Element = ProfiniteRational
    print_mode = "mod"
    def __init__(self, base=ZZ, category=Rings(), print_mode="mod"):
        # We must provide a base... What should base be? What is base?
        Ring.__init__(self, base, category=category) 
        ProfiniteRationals.print_mode = print_mode
    
    def _repr_(self):
        return "Ring of profinite rationals"
    
    def _element_constructor_(self, *args, **kwds):
        if len(args) < 2:
            return self.element_class(args[0], 0)
        return self.element_class(*args, **kwds)
        '''x = args[0]
        if not x in QQ:
            return self.element_class(x, **kwds)
        return self.element_class(x, 0)'''
    
    def _coerce_map_from_(self, QQ):
        return True # use the above conversion ("_element_constructor_") as coercion

    def set_print_mode(self, mode):
        modes = ["mod", "prime", "factorial"]
        if mode in modes:
            ProfiniteRationals.print_mode = mode
        else:
            print("mode {} not recognized; choose one of the following:".format(mode))
            print(modes)

    def is_finite(self):
        return False

Qhat = ProfiniteRationals()
