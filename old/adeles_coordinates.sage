
load("adeles_rational.sage")
# TODO check dat 0-moduli juist verwerkt worden overal

class Adele(CommutativeRingElement):
    """
    An adele
    """

    def __init__(self, parent, x):
        if x in parent.number_field:
            coords = [AQ(xi) for xi in x.list()]
            self.coordinates = vector(coords) 
        elif x in AQ:
            n = parent.number_field.absolute_degree()
            # Length n all-zero vector in AQ:
            self.coordinates = vector(AQ, n)
            self.coordinates[0] = x
        else:
            r"""
            Let ``K`` be the associated number field.
            Let ``a`` be the generator of ``K`` over `\QQ` (i.e.
            K.absolute_generator()).
            Let ``n`` be the degree of ``K`` over `\QQ`.
            Let ``coordinates`` be the vector/list/...  [c_0, c_1, c_2, ..., c_{n-1}],
            where each c_i is a rational adèle.

            Create the adèle
                ``c_0 + c_1*a + c_2*a^2 + ... + c_{n-1}*a^(n-1)``
            of ``K``.
            """
            coordinates = x
            #print("DEBUG: start Adele.__init__")
            if len(coordinates) != parent.number_field.absolute_degree():
                raise ValueError("there should be as many coordinates as the degree of the number field")
            #print("DEBUG: len checked")
            for c in coordinates:
                if c not in AQ:
                    raise TypeError("coordinates should be rational adeles (but c={} of type {}".format(c, type(c)))
            #print("DEBUG: types checked")
            self.coordinates = vector(coordinates)
        CommutativeRingElement.__init__(self, parent)

    def _repr_(self):
        n = self.parent().number_field.degree()
        gen = self.parent().number_field.gen()
        rep = str(self.coordinates[0])
        if n > 1:
            rep += " + {}*{}".format(self.coordinates[1], gen)
        for i in range(2, n):
            rep += " + {}*{}^{}".format(self.coordinates[i], gen, i)
        return rep

    def _richcmp_(self, other, op):
        from sage.structure.richcmp import op_EQ, op_NE
        if op == op_EQ:
            n = self.parent().number_field.degree()
            c = self.coordinates
            d = other.coordinates
            return all([c[i] == d[i] for i in range(n)])
        if op == op_NE:
            return not self._richcmp_(other, op_EQ)
        raise NotImplementedError()

    def _add_(self, other):
        coordinates = self.coordinates + other.coordinates
        return self.__class__(self.parent(), coordinates)

    def _sub_(self, other):
        coordinates = self.coordinates - other.coordinates
        return self.__class__(self.parent(), coordinates)
        
    def _mul_(self, other):
        n = self.parent().number_field.degree()
        power = self.parent().power_table
        c = self.coordinates
        d = other.coordinates
        result = vector(AQ, n)  # (AQ(0), AQ(0), ..., AQ(0)) of length ``n``
        for i in range(n):
            for j in range(n):
                if i+j < n:
                    result[i+j] += c[i] * d[j]
                else:
                    for k in range(n):
                        result[k] += power[i+j][k] * (c[i] * d[j])
        return self.__class__(self.parent(), result)

    def prime_repr(self):
        K = self.parent().number_field
        n = K.degree()
        result = []

        # We start of with the infinite primes of `K`. These are given by the
        # embeddings `K \to \CC` up to complex conjugation. # TODO kies de juiste geconjugeerde!
        oo_primes = set()
        images = set()
        for phi in K.embeddings(CC):
            if phi(K.gen()).conjugate() not in images:
                images.add(phi(K.gen()))
                oo_primes.add(phi)
        # Now we compute the images of ``self`` in the completions of K at the
        # infinite primes (which are RR or CC).
        for phi in oo_primes:
            im_gen = phi(K.gen())
            # im_gen is either in RR or in CC. The python int 0 below coerces
            # into both.
            im = 0  
            for i in range(n):
                im += self.coordinates[i].real * im_gen^i
            result.append(im)

        # Next up are the finite primes of `K`. We only have to consider those
        # at which not all moduli of ``self`` are units. # TODO moet "at which none of the moduli" zijn toch? ook code veranderen
        primes = set()
        for i in range(n):
            F = factor(self.coordinates[i].profinite_rational.modulus)
            for p, e in F:
                primes.add(p)
        primes = sorted(primes)
        # Now, for each prime `p`, we factor the defining polynomial of `K`
        # over the `p`-adics. For each irreducible factor `g` of `f` in Qp[x],
        # we construct the completion `L` of `K` at the prime of `K` above `p`
        # corresponding to `g`.
        for p in primes:
            S.<x> = PolynomialRing(Qp(p))
            im = S(0)
            for i in range(n):
                im += self.coordinates[i].at_place(p) * x^i
            f = S(K.defining_polynomial())
            for g, e in factor(f):
                L = S.quotient_by_principal_ideal(g, str(K.gen()))
                result.append(L(im))
        return str(result)




from sage.structure.unique_representation import UniqueRepresentation
class AdeleRing(UniqueRepresentation, CommutativeRing):
    Element = Adele

    def __init__(self, K):
        """
        Create the adèle ring of the number field ``K``
        """
        self.number_field = K
        self.power_table = []
        n = K.degree()
        for i in range(0, 2*n-1):
            self.power_table.append(vector(K.gen()^i))
        # We choose K as our base ring.
        CommutativeRing.__init__(self, K)

    def _repr_(self):
        return "Adele Ring of {}".format(self.number_field)

    def _latex_(self):
        return r"\Bold{A}_{ " + latex(self.number_field) + " }"

    def characteristic(self):
        return 0

    def _element_constructor_(self, coordinates):
        return self.element_class(self, coordinates)

    def _coerce_map_from_(self, S):
        if S is self.number_field:
            return True
        if S is AQ:
            return True
        return False

    def is_integral_domain(self, proof=None):
        return False

def AA(K):
    return AdeleRing(K)