INFINITY = 10000000000000

def factorial(n):
    if n <= 1:
        return 1
    return n * factorial(n-1)


class ProfiniteInteger(RingElement):
    # zero is represented by the empty digits sequence
    default_prec = 10
    
    def __init__(self, parent, digits):
        for i in range(0, len(digits)):
            if digits[i] < 0:
                raise ValueError("{}th factorial digit is negative".format(i+1))
            if digits[i] > i+1:
                raise ValueError("{}th factorial digit is bigger than {}".format(i+1, i+1))
        self.digits = digits
        RingElement.__init__(self, parent)
    
    def _repr_(self):
        if not self.digits:
            return "0"
        s = ""
        all_zero = True
        for i in range(0, len(self.digits)):
            if self.digits[i] != 0:
                s += "{}*{}! + ".format(self.digits[i], i+1)
                all_zero = False
        if all_zero:
            s += "0 + "
        s += "O({}!)".format(len(self.digits)+1)
        return s
    
    def _richcmp_(self, other, op):
        from sage.structure.richcmp import richcmp
        
    def _add_(self, other):
        prec = min(len(self.digits), len(other.digits))
        digits = []
        carry = 0
        for i in range(0, prec):
            d = self.digits[i] + other.digits[i] + carry
            if d > i + 1:
                d -= i + 2
                carry = 1
            else:
                carry = 0
            digits.append(d)
        return self.__class__(self.parent(), digits)
        
    def _sub_(self, other):
        prec = min(len(self.digits), len(other.digits))
        digits = []
        carry = 0
        for i in range(0, prec):
            d = self.digits[i] - other.digits[i] - carry
            if d < 0:
                d += i + 2
                carry = 1
            else:
                carry = 0
            digits.append(d)
        return self.__class__(self.parent(), digits)
        
    def _mul_(self, other):
        prec = min(len(self.digits), len(other.digits))
        a = self.factorial_sum(prec)
        b = other.factorial_sum(prec)
        print(prec)
        return self.parent()(a*b, prec=prec) # note: only the first prec digits are (always) correct
    
    def factorial_sum(self, end=None):
        if end is None:
            end = len(self.digits)
        fac_sum = 0
        k = 1
        k_fac = 1 # k!
        for i in range(0, end):
            fac_sum += k_fac * self.digits[i]
            k += 1
            k_fac *= k
        return fac_sum



from sage.categories.rings import Rings
class ProfiniteIntegers(Ring):
    Element = ProfiniteInteger
    def __init__(self, base=ZZ, category=Rings()):
        # We must provide a base... What should base be? What is base?
        Ring.__init__(self, base, category=category) 
    
    def _repr_(self):
        return "Ring of profinite integers"
    
    def _element_constructor_(self, *args, **kwds):
        if not args:
            return self.element_class(self, *args, **kwds)
        n = args[0]
        if not n in ZZ:
            return self.element_class(self, n, **kwds)
        if n < 0:
            n = factorial(self.__class__.Element.default_prec + 1) + n
        if "prec" in kwds:
            k_max = kwds["prec"]
        else:
            k_max = INFINITY
        k = 1
        k_fac = 1 # k!
        kp1_fac = 2 # (k+1)!
        digits = []
        while n > 0 and k <= k_max:
            dkfac = n % kp1_fac # d_k * k!
            digits.append(dkfac//k_fac)
            n -= dkfac
            k += 1
            k_fac = kp1_fac
            kp1_fac *= k + 1
        return self.element_class(self, digits)
    
    def _coerce_map_from_(self, ZZ):
        return True # use the above conversion ("_element_constructor_") as coercion


Zhat = ProfiniteIntegers()
