class FiniteAdele_marco:
    def __init__(self, x, m):
        self.x = x
        self.m = m
    
    def __str__(self):
        return str(self.x) + " mod " + str(self.m)
    
    def __repr__(self):
        return str(self.x) + " mod " + str(self.m)
    
    def __add__(self, other):
        x = self.x + other.x
        a, b = self.m.as_integer_ratio()
        c, d = other.m.as_integer_ratio()
        m = gcd(a*d, b*c) / (b*d)
        return self.__class__(x, m)
    
    def __mul__(self, other):
        x = self.x * other.x
        x_n, x_d = self.x.as_integer_ratio()
        m_n, m_d = self.m.as_integer_ratio()
        y_n, y_d = other.x.as_integer_ratio()
        n_n, n_d = other.m.as_integer_ratio()
        den = x_d * m_d * y_d * n_d
        num = gcd(x_n*n_n*y_d*m_d, gcd(y_n*m_n*x_d*n_d, m_n*n_n*x_d*y_d))
        m = num / den
        return self.__class__(x, m)

    
a = FiniteAdele_marco(5/17, 20) # 5/17 + 20Z^
b = FiniteAdele_marco(6, 35) # 6 + 35Z^
a+b