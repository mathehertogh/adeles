attach("profinite_integers.sage")

class Adele_factorial:
    def __init__(self, real, denominator, pfint):
        if real not in RR:
            raise ValueError("real coordinate not in RR")
        if denominator not in NN:
            raise ValueError("denominator is not a natural number")
        if pfint not in ProfiniteIntegers():
            pass #raise ValueError("pfint is not a profinite number") TODO
        self.real = real
        self.denominator = denominator
        self.pfint = pfint
    
    def __repr__(self):
        return "({}, 1/{} * ({}))".format(self.real, self.denominator, self.pfint)
    
    def __add__(self, other):
        real = self.real + other.real
        denominator = self.denominator * other.denominator
        pfint = other.denominator*self.pfint + self.denominator*other.pfint
        return Adele_factorial(real, denominator, pfint)

    def __mul__(self, other):
        real = self.real * other.real
        denominator = self.denominator * other.denominator
        pfint = self.pfint * other.pfint
        return Adele_factorial(real, denominator, pfint)

x = Adele_factorial(0.324, 6, Zhat(25))
y = Adele_factorial(0.615, 10, Zhat(-7))