class Adele_per_place:
    def __init__(self, real, coords):
        if not real in RR:
            raise ValueError("real coordinate not in RR")
        for p, x_p in coords.items():
            if not p in ZZ or not p.is_prime():
                raise ValueError("{} is not a prime number".format(p))
            if not x_p in Qp(p):
                raise ValueError("{} is not a {}-adic number".format(x_p, p))
        self.real = real
        self.coords = coords
    
    def __repr__(self):
        s = "([{}] at infty, ".format(self.real)
        for p, x_p in self.coords.items():
            s += "[" + str(x_p) + "] at " + str(p) + ", "
        s = s[0:-2] + ")"
        return s
    
    def __add__(self, other):
        real = self.real + other.real
        coords = {}
        for p, x_p in self.coords.items():
            if p in other.coords:
                coords[p] = x_p + other.coords[p]
        return self.__class__(real, coords)
    
    def __mul__(self, other):
        real = self.real * other.real
        coords = {}
        for p, x_p in self.coords.items():
            if p in other.coords:
                coords[p] = x_p * other.coords[p]
        return self.__class__(real, coords)


Q2 = Qp(2, prec=5)
Q3 = Qp(3, prec=5)
Q5 = Qp(5, prec=5)
Q7 = Qp(7, prec=5)

a = Adele_per_place(0.9797, {2: Q2(8), 5: Q5(1/10)})
b = Adele_per_place(1.5, {3: Q3(12), 5: Q5(250)})
a+b