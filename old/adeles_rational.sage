
load("pfrats.sage")

class RationalAdele(CommutativeRingElement):
    def __init__(self, parent, x, y=None):
        if y is None: 
            if x not in QQ:
                raise TypeError("no conversion from {} to RationalAdele".format(type(x)))
            # We embed QQ diagonally into AQ.
            y = Qhat(x)
            x = RR(x)
        if x not in RR or y not in Qhat:
            raise TypeError("real should be in RR and profinite_rational in Qhat")
        self.real = x
        self.profinite_rational = y
        CommutativeRingElement.__init__(self, parent)

    def _repr_(self):
        return "({}, {})".format(self.real, self.profinite_rational)

    def _richcmp_(self, other, op):
        #print("DEBUG - RationalAdele.__rechcmp__({}, {}, {})".format(self, other, op))
        from sage.structure.richcmp import op_EQ, op_NE
        if op == op_EQ:
            return self.real == other.real and self.profinite_rational == other.profinite_rational
        if op == op_NE:
            return not self._richcmp_(other, op_EQ)
        raise NotImplementedError()

    def _add_(self, other):
        real = self.real + other.real
        profinite_rational = self.profinite_rational + other.profinite_rational
        return self.__class__(self.parent(), real, profinite_rational)
        
    def _sub_(self, other):
        real = self.real - other.real
        profinite_rational = self.profinite_rational - other.profinite_rational
        return self.__class__(self.parent(), real, profinite_rational)
        
    def _mul_(self, other):
        real = self.real * other.real
        profinite_rational = self.profinite_rational * other.profinite_rational
        return self.__class__(self.parent(), real, profinite_rational)

    def at_place(self, p):
        if p is Infinity:
            return self.real
        if p not in ZZ or not p.is_prime():
            raise ValueError("p should be a rational prime number or Infinity")
        return self.profinite_rational.to_Qp(p)

from sage.structure.unique_representation import UniqueRepresentation
class RationalAdeles(UniqueRepresentation, CommutativeRing):
    Element = RationalAdele

    def __init__(self):
        # We choose QQ as our base ring.
        CommutativeRing.__init__(self, QQ)

    def _repr_(self):
        return "Rational Adele Ring"

    def _latex_(self):
        return r"\Bold{A}_\QQ"

    def characteristic(self):
        return 0

    def _element_constructor_(self, x, y=None):
        return self.element_class(self, x, y)

    def _coerce_map_from_(self, S):
        if S is QQ:
            return True
        return False

    def is_integral_domain(self, proof=None):
        return False

AQ = RationalAdeles()
    
