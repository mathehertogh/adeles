
load("completions.sage")
load("profinite_numbers.sage")
#load("ideles.sage")

from sage.structure.element import CommutativeAlgebraElement
class Adele(CommutativeAlgebraElement):
    """
    TODO docstring
    """

    def __init__(self, parent, infinite, finite):
        r"""
        Construct the adele determined by ``infinite`` and ``finite``

        INPUT:

        - ``parent`` -- an adele ring of some number field `K`
        - ``infinite`` -- a list of elements in real/complex interval fields,
          corresponding to ``completions(K, oo)``
        - ``finite`` -- a profinite number over `K`

        EXAMPLES::

            sage: A = Adeles(QQ)
            sage: Adele(A, [3.14], 7)
            (3.1400000000000002?, 7)
            sage: Qhat = ProfiniteNumbers(QQ)
            sage: Adele(A, [7.9], Qhat(7, 9, 2))
            (7.9000000000000004?, (7 mod 9)/2)

        ::

            sage: K.<a> = NumberField(x^2+15)
            sage: Ak = Adeles(K)
            sage: Khat = ProfiniteNumbers(K)
            sage: Adele(Ak, [I], Khat(a, 21, a+1))
            (1*I, (a mod (21))/(a + 1))
            sage: Adele(Ak, [0.0], Khat(2*a-1, 100, a))
            (0, (2*a - 1 mod (100))/a)

        TESTS::

            sage: Adele(A, [3.0, -1.0], 7)
            Traceback (most recent call last):
            ...
            TypeError: infinite should have length 1
            sage: Adele(A, (1.0,), 7)
            Traceback (most recent call last):
            ...
            TypeError: infinite should be a list
            sage: Adele(A, [I], 7)
            Traceback (most recent call last):
            ...
            TypeError: 0th infinite value (I) should lie in Real Interval Field with 53 bits of precision
            sage: Adele(A, [1.0], a+1)
            Traceback (most recent call last):
            ...
            TypeError: finite should lie in Profinite Numbers of Rational Field
        """
        debug("Adele.__init__({}, {})".format(infinite, finite))
        CommutativeAlgebraElement.__init__(self, parent)
        K = parent.base()
        K_oo = infinite_completions(K)
        t = len(K_oo)
        if not isinstance(infinite, list):
            raise TypeError("infinite should be a list")
        if len(infinite) != t:
            raise TypeError("infinite should have length {}".format(t))
        self.infinite = infinite
        for i in range(t):
            val = infinite[i]
            if val not in K_oo[i][0]:
                raise TypeError("{}th infinite value ({}) should lie in {}".format(i, val, K_oo[i][0]))
            self.infinite[i] = K_oo[i][0](val)
        Khat = ProfiniteNumbers(K)
        if finite not in Khat:
            raise TypeError("finite should lie in {}".format(Khat))
        self.finite = Khat(finite)

    def _repr_(self):
        """
        Return a string representation of ``self``

        EXAMPLES::

            sage: K.<zeta> = CyclotomicField(5)
            sage: A = Adeles(K)
            sage: Khat = ProfiniteNumbers(K)
            sage: A([-1, I], Khat(zeta-1, zeta^3+zeta, 10))
            (-1, 1*I, (0 mod (zeta^3 + zeta))/10)
            sage: inf = CIF(RIF(-oo, oo), RIF(-oo, oo))
            sage: A([-97*I, inf], Khat(zeta, 100, 1+zeta^2))
            (-97*I, CC, (zeta mod (100))/(zeta^2 + 1))
        """
        debug("Adele._repr_()")
        rif_oo = RIF(-oo, oo)
        cif_oo = CIF(rif_oo, rif_oo)
        rep = "("
        for i in range(len(self.infinite)):
            x_oo = self.infinite[i]
            if x_oo.endpoints() == rif_oo.endpoints():
                rep += "RR, "
            elif x_oo.endpoints() == cif_oo.endpoints():
                rep += "CC, "
            else:
                rep += repr(x_oo) + ", "
        rep += repr(self.finite) + ")"
        return rep

    def _richcmp_(self, other, op):
        r"""
        Compare ``self`` and ``other`` based on the relation ``op``

        We only implement equality and non-equality.

        We declare two adeles equal if the *could* be equal: the open subsets
        they represents should have non-empty intersection.

        EXAMPLES::

            sage: A = Adeles()
            sage: Qhat = ProfiniteNumbers(QQ)
            sage: b = A([RIF(2.0, 3.0)], Qhat(3, 10, 7))
            sage: c = A([RIF(2.5, 3.5)], Qhat(13, 20, 7))
            sage: d = A([2.1], Qhat(13, 20, 7))
            sage: e = A([2.1], Qhat(3, 20, 7))
            sage: b == c
            True
            sage: b == d
            True
            sage: c == d
            False
            sage: b != e
            False
            sage: d == e
            False
            sage: c != d
            True

        ::

            sage: K.<i> = NumberField(x^2+1)
            sage: Ak = Adeles(K)
            sage: Khat = ProfiniteNumbers(K)
            sage: b = Ak([I], Khat(i, 10+i, 3))
            sage: c = Ak([RIF(-1, 4)*I], Khat(-10, 20+2*i, 3))
            sage: b == c
            True
            sage: b != c
            False
        """
        debug("Adele._richcmp_({}, {})".format(self, other))
        from sage.structure.richcmp import op_EQ, op_NE
        if op == op_EQ:
            for i in range(len(self.infinite)):
                try:
                    self.infinite[i].intersection(other.infinite[i])
                except ValueError:
                    # This indicates the intersection is empty
                    return False
            return self.finite == other.finite
        if op == op_NE:
            return not self._richcmp_(other, op_EQ)
        raise NotImplementedError("only equality and inequality are implemented")

    def _add_(self, other):
        """
        Add ``other`` to ``self`` and return the result

        EXAMPLES::

            sage: A = Adeles(QQ)
            sage: Qhat = ProfiniteNumbers(QQ)
            sage: b = A([1.5], Qhat(-1, 12, 6))
            sage: c = A([2.5], Qhat(3, 24, 2))
            sage: b + c
            (4, (4 mod 6)/3)
        
        ::

            sage: K.<a> = CyclotomicField(5)
            sage: Ak = Adeles(K)
            sage: Khat = ProfiniteNumbers(K)
            sage: b = Ak([I, -I], Khat(a+5, 30*a^2, a))
            sage: c = Ak([3.5, I], Khat(1, 50, a^3))
            sage: b + c
            (3.5000000000000000? + 1*I, 0, (4*a^3 - a^2 - 1 mod (10))/(-a^3 - a^2 - a - 1))
            sage: b + a
            (-0.80901699437494746? + 1.5877852522924732?*I, 0.30901699437494746? - 0.048943483704846358?*I, (a^2 + a + 5 mod (30))/a)
        """
        debug("Adele._add_({}, {})".format(self, other))
        infinite = []
        for i in range(len(self.infinite)):
            infinite.append(self.infinite[i] + other.infinite[i])
        finite = self.finite + other.finite
        return self.__class__(self.parent(), infinite, finite)

    def _sub_(self, other):
        """
        Subtract ``other`` from ``self`` and return the result

        EXAMPLES::

            sage: A = Adeles(QQ)
            sage: Qhat = ProfiniteNumbers(QQ)
            sage: b = A([-10.375], Qhat(4, 15, 3))
            sage: c = A([5], Qhat(2, 25, 5))
            sage: b - c
            (-15.375000000000000?, (14 mod 75)/15)
            sage: c - b
            (15.375000000000000?, (61 mod 75)/15)

        ::

            sage: K.<a> = NumberField(x^2+3)
            sage: Ak = Adeles(K)
            sage: Khat = ProfiniteNumbers(K)
            sage: b = Ak([I], Khat(a-1, 5*a, 2))
            sage: c = Ak([1], Khat(2*a, 10*a+20, a+1))
            sage: b - c
            (-1 + 1*I, (1/2*a + 1/2 mod (5))/(a + 1))
            sage: c - b
            (1 - 1*I, (-1/2*a - 1/2 mod (5))/(a + 1))
            sage: b - 1/a
            (1.577350269189626?*I, (13/2*a + 5/2 mod (15))/2*a)
        """
        debug("Adele._sub_({}, {})".format(self, other))
        infinite = []
        for i in range(len(self.infinite)):
            infinite.append(self.infinite[i] - other.infinite[i])
        finite = self.finite - other.finite
        return self.__class__(self.parent(), infinite, finite)

    def _mul_(self, other):
        """
        Multiply ``self`` and ``other`` and return the result

        EXAMPLES::

            sage: A = Adeles(QQ)
            sage: Qhat = ProfiniteNumbers(QQ)
            sage: b = A([-1], Qhat(7, 9, 4))
            sage: c = A([79], Qhat(70, 90, 3))
            sage: b * c
            (-79, (20 mod 45)/6)

        ::

            sage: K.<a> = NumberField(x^2+x+1)
            sage: Ak = Adeles(K)
            sage: Khat = ProfiniteNumbers(K)
            sage: b = Ak([I], Khat(a-1, 8*a, 3))
            sage: c = Ak([CIF(2.3+I)], Khat(7*a+2, 30*a-3, a))
            sage: b * c
            (-1 + 2.2999999999999999?*I, (0 mod (1))/(-a - 2))
            sage: a * b
            (-0.86602540378443860? - 0.50000000000000000?*I, (-2*a - 1 mod (-8*a - 8))/3)
            sage: b * (3/(a+1))
            (2.598076211353316? + 1.5000000000000000?*I, 2*a + 1 mod (8*a + 8))
        """
        debug("Adele._mul_({}, {})".format(self, other))
        infinite = []
        for i in range(len(self.infinite)):
            infinite.append(self.infinite[i] * other.infinite[i])
        finite = self.finite * other.finite
        return self.__class__(self.parent(), infinite, finite)

    def to_modulo_element(self):
        """
        Convert this *integral* adele to an element of `O/I`, with `O` the
        maximal order in our number feild and `I` the ideal up to which our
        finite part is defined

        If ``self`` is not integral, throw an exception.

        EXAMPLES::

            sage: A = Adeles(QQ)
            sage: Qhat = ProfiniteNumbers(QQ)
            sage: b = A([1], Qhat(3, 10))
            sage: b_bar = b.to_modulo_element()
            sage: b_bar, b_bar.parent()
            (3, Ring of integers modulo 10)
            sage: c = A([pi.n()], Qhat(97, 790))
            sage: c_bar = c.to_modulo_element()
            sage: c_bar, c_bar.parent()
            (97, Ring of integers modulo 790)

        ::

            sage: K.<a> = NumberField(x^2+5)
            sage: Ak = Adeles(K)
            sage: Khat = ProfiniteNumbers(K)
            sage: b = Ak([I], Khat(a, K.ideal(2, a+1)))
            sage: b_bar = b.to_modulo_element()
            sage: b_bar, b_bar.parent()
            (a,
             Quotient of Maximal Order in Number Field in a with defining polynomial x^2 + 5 by the ideal (2, a + 1))
            sage: c = Ak([I], Khat(-3, K.ideal(9, 3*a+3), 3))
            sage: c_bar = c.to_modulo_element()
            sage: c_bar, c_bar.parent()
            (a,
             Quotient of Maximal Order in Number Field in a with defining polynomial x^2 + 5 by the ideal (3, a + 1))

        TESTS::

            sage: d = Ak([1], Khat(a, 7, 3))
            sage: d.to_modulo_element()
            Traceback (most recent call last):
            ...
            ValueError: non-integral adele can't be converted to O/I
        """
        if self.finite.denominator != 1:
            raise ValueError("non-integral adele can't be converted to O/I")
        K = self.parent().base()
        O = K.maximal_order()
        I = self.finite.numerator.modulus
        x = self.finite.numerator.value
        if I.is_zero():
            return O(x)
        name = K.variable_name() + "_bar"
        OmodI = O.quotient(I, name)
        return OmodI(x)

    def to_rational_adele_vector(self):
        r"""
        Return ``self`` as a rational adele vector

        Write `A_K` for ``self.parent()``, with base the number field `K` and
        denote the rational adeles by `A_\QQ`.
        Then we have a canonical isomorphism of topological rings
        `\phi: A_K \to A_\QQ \otimes K`.
        The codomain is a free `A_\QQ`-module with basis
        `B = \{1, a, a^2, ... a^{n-1}\}`, where `a` is the algebraic number
        adjoined to `\QQ` to obtain `K` and `n = \deg(K/\QQ)`.

        This method returns the image of ``self`` under `\phi` in the form of an
        `A_\QQ`-vector relative to the basis `B`.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3-2)
            sage: Khat = ProfiniteNumbers(K)
            sage: AK = Adeles(K)
            sage: b = AK([-1, CIF(7.9*I)], Khat(a, 10*(a+1), 3))
            sage: c = b.to_rational_adele_vector()
            sage: d = AK._from_rational_adele_vector(c)
            sage: c
            ((-0.33333333333333?, (0 mod 10)/3), (3.35555453543495?, (1 mod 5)/3), (-3.08327908304135?, (0 mod 5)/3))
            sage: b == d
            True
            sage: e = AK([2, -1], Khat(a/3))
            sage: e.to_rational_adele_vector()
            ((0.?e-14, 0), (0.79370052598410?, 1/3), (0.629960524947437?, 0))
        """
        K = self.parent().base()
        n = K.absolute_degree()
        
        finites = self.finite.to_profinite_rational_vector()

        r, s = K.signature()
        K_oo = infinite_completions(K)
        A = []
        for i in range(r+s):
            phi = K_oo[i][1] # phi: K --> \RR
            gen_im = phi(K.gen())
            row = []
            if i < r: # Real place
                for j in range(n):
                    row.append(gen_im^j)
            else: # Complex place
                for j in range(n):
                    row.append((gen_im^j).real())
                A.append(row)
                row = []
                for j in range(n):
                    row.append((gen_im^j).imag())
            A.append(row)
        A = matrix(RIF, A)

        Y = []
        for i in range(r+s):
            if i < r: # Real place
                Y.append(self.infinite[i])
            else: # Complex place
                Y.append(self.infinite[i].real())
                Y.append(self.infinite[i].imag())
        Y = vector(Y)

        infinites = A.solve_right(Y)

        A_Q = Adeles(QQ)
        result = []
        for i in range(n):
            a = A_Q([infinites[i]], finites[i])
            result.append(a)

        return vector(result)


from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.ring import CommutativeAlgebra
class Adeles(UniqueRepresentation, CommutativeAlgebra):
    Element = Adele

    def __classcall__(cls, K=QQ):
        """
        Construct the adele ring of the number field `K`

        INPUT:

        - ``K`` -- a number field (default: ``QQ``)

        EXAMPLES::

            sage: Adeles("not a number field")
            Traceback (most recent call last):
            ...
            TypeError: K should be a number field

        We make sure this returns True::

            sage: Adeles() is Adeles(QQ)
            True
        """
        debug("Adeles.__classcall__({})".format(K))
        try:
            if not is_field(K) or not K.absolute_degree() in ZZ:
                raise TypeError("K should be a number field")
        except AttributeError:
            raise TypeError("K should be a number field")
        return super(Adeles, cls).__classcall__(cls, K)

    def __init__(self, K=QQ):
        r"""
        Construct the adele ring of the number field ``K``

        INPUT:

        - ``K`` -- a number field (default: `\QQ`)

        EXAMPLES::

            sage: Adeles()
            Adele Ring of Rational Field
            sage: K.<a> = NumberField(x^5-3*x+1)
            sage: Adeles(K)
            Adele Ring of Number Field in a with defining polynomial x^5 - 3*x + 1
        """
        debug("Adeles.__init__({})".format(K))
        CommutativeAlgebra.__init__(self, K)

    def _repr_(self):
        """
        Return a string representation of ``self``

        EXAMPLES::

            sage: K = NumberField(x^3+x+1, 'a')
            sage: Adeles(K)
            Adele Ring of Number Field in a with defining polynomial x^3 + x + 1
        """
        debug("Adeles._repr_()")
        return "Adele Ring of {}".format(self.base())

    def _latex_(self):
        r"""
        Return latex-formatted string representation of ``self``

        EXAMPLES::

            sage: K.<a> = NumberField(x^2+14)
            sage: latex(Adeles(K))
             \Bold{A}_{ \Bold{Q}[a]/(a^{2} + 14) }
        """
        debug("Adeles._latex_()")
        return r" \Bold{A}_{" + latex(self.base()) + "} "

    def characteristic(self):
        """
        Return the characteristic of this ring, which is zero

        EXAMPLES::

            sage: Adeles().characteristic()
            0
        """
        debug("Adeles.characteristic()")
        return ZZ(0)

    def _element_constructor_(self, x, y=None):
        """
        Construct the adele determined by ``x`` and ``y``

        EXAMPLES::

            sage: A = Adeles(QQ)
            sage: Qhat = ProfiniteNumbers(QQ)
            sage: A([-1.2345], Qhat(4, -10, 9))
            (-1.2345000000000000?, (4 mod 10)/9)
            sage: A([97], 97)
            (97, 97)
            sage: A(1/3)
            (0.3333333333333334?, 1/3)

        ::

            sage: K.<a> = NumberField(x^3-7)
            sage: Ak = Adeles(K)
            sage: Khat = ProfiniteNumbers(K)
            sage: Ak([3.14, -I], Khat(a, 9*a^2, a+1))
            (3.1400000000000002?, -1*I, (a mod (9*a^2))/(a + 1))
            sage: Ak(1/a)
            (0.522757958574711?, -0.261378979287356? - 0.452721672156193?*I, 1/a)
        """
        debug("Adeles._element_constructor_({}, {})".format(x, y))
        K = self.base()
        J_K = IdeleGroup(K)
        if y is None:
            if x in K:  # coercion K --> A_K
                infinite = [phi(x) for L, phi in completions(K, oo)]
                return self.element_class(self, infinite, x)
            if x.parent() is self: # make copy
                return self.element_class(self, x.infinite, x.finite)
            if x in J_K:  # coercion J_K --> A_K
                return self._from_idele(x)
        return self.element_class(self, x, y)

    def _from_idele(self, idele):
        """
        Construct the adele corresponding to the idele ``idele``

        This implements the natural embedding of the ideles into the adeles. Due
        to the ways we store ideles/adeles, this may lose precision: the returend
        ``Adele`` represents a subset of the adeles that contains all ideles
        ``self`` represents. TODO happens when p^i*Z_p* (i.e. x*U_q^0 w/
        v_q(x)>0) explain+example. Also if exact set?

        EXAMPLES::

            sage: A = Adeles(QQ)
            sage: J = IdeleGroup(QQ)
            sage: u = J(None, None, {2: (1/4, 1), 3: (3/5, 2)})
            sage: v = J(None, [-1], {2: (4, 2), 3: (1/4, 1), 5: (3, 0)})
            sage: w = J(None, [RIF(7, 9)], {2: (4, 2), 3: (1/4, 1), 5: (5, 0)})
            sage: A(u)
            (RR, (51 mod 54)/4)
            sage: A(v)
            (-1, 4 mod 48)
            sage: A(w)
            (8.?, 100 mod 240)
            sage: A(u*v)
            (RR, 15 mod 18)
            sage: A(u)*A(v)
            (RR, 15 mod 18)

        ::

            sage: K.<a> = NumberField(x^2+5)
            sage: Ak = Adeles(K)
            sage: Jk = IdeleGroup(K)
            sage: u = Jk(None, [I], {})
            sage: Ak(u)
            (1*I, 0 mod (1))
            sage: p2, p3 = K.prime_above(2), K.prime_above(3)
            sage: v = Jk(None, None, {p2: (a, 3), p3: (1/3-a, 3)})
            sage: Ak(v)
            (CC, (a + 2 mod (108, 6*a + 42))/3)

        This defines a coercion::

            sage: J = IdeleGroup(QQ)
            sage: A = Adeles(QQ)
            sage: Qhat = ProfiniteNumbers(QQ)
            sage: a = A([2.5], Qhat(3, 10, 2))
            sage: u = J(None, [-0.5], {2: (2, 2), 5: (1/2, 1)})
            sage: a+u
            (2, (9 mod 10)/2)
            sage: u - a
            (-3, (3 mod 10)/2)
            sage: a * u
            (-1.2500000000000000?, 7 mod 10)

        .. TODO::
            
            - empty finite
            - non-matching exact and finite part
        """
        K = self.base()
        Khat = ProfiniteNumbers(K)

        if idele._has_exact() and not idele._contains(idele.exact):
            raise NotImplementedError("non-matching exact and finite part not implemented yet...")

        K_oo = completions(K, oo)
        infinite = idele.infinite.copy()
        for i in range(len(infinite)):
            if infinite[i] is None:
                if idele._has_exact():
                    phi = K_oo[i][1]  # phi: K --> RIF/CIF
                    infinite[i] = phi(idele.exact)
                elif K_oo[i][0] is RIF:
                    infinite[i] = RIF(-oo, oo)
                else:
                    infinite[i] = CIF(RIF(-oo, oo), RIF(-oo, oo))

        integral, denominator = idele.integral_split()
        values = []
        moduli = []
        for q, val in integral.finite.items():
            x, i = val
            if i == ZZ(0):
                values.append(0)
            else:
                values.append(x)
            if K is not QQ:
                x = K.ideal(x)
            e = x.valuation(q)
            moduli.append(q^(i+e))
        if K is QQ:
            value = CRT(values, moduli)
        else:
            value = K.solve_CRT(values, moduli)
        modulus = prod(moduli)
        finite = Khat(value, modulus, denominator)

        if len(moduli) == 0 and idele._has_exact():
            finite = idele.exact

        return self.element_class(self, infinite, finite)

    def _from_rational_adele_vector(self, v):
        r"""
        Build a `K`-adele out of the rational adele vector ``v``.

        Let `K = \QQ(a)` be our number field of degree `n` over `\QQ`.

        INPUT:

        - ``v`` -- a vector of rational adeles of length `n`

        OUPUT:

        The `K`-adele `v[0] + v[1]*a + v[2]*a^2 + ... + v[n-1]*a^{n-1}`.

        .. TODO::

            Move the computation on the finite part to a separate method of
            :class:`ProfiniteNumber`.
        """
        K = self.base()
        n = K.absolute_degree()

        infinite = []
        places = [phi for L, phi in infinite_completions(K)]
        for k in range(len(places)):
            phi = places[k]
            x_oo = sum([v[i].infinite[0] * phi(K.gen())^i for i in range(n)])
            infinite.append(x_oo)
        denominator = lcm([v[i].finite.denominator for i in range(n)])
        v = denominator * v
        value = sum([v[i].finite.numerator.value * K.gen()^i for i in range(n)])
        modulus = sum([K.ideal(v[i].finite.numerator.modulus * K.gen()^i) for i in range(n)])
        finite = ProfiniteNumbers(K)(value, modulus, denominator)

        return self.element_class(self, infinite, finite)




    def _coerce_map_from_(self, S):
        r"""
        Return a coerce map ``self`` --> ``S`` (or ``True`` to use
        ``_element_constructor_``) if it exists, ``False`` otherwise.

        EXAMPLES::

            sage: K.<a> = NumberField(x^9+x+1)
            sage: A = Adeles(K)
            sage: A._coerce_map_from_(K)
            True
            sage: A._coerce_map_from_(CyclotomicField())
            False
            sage: A._coerce_map_from_(ZZ)
            True

        .. TODO::

            For a field extension K --> L, define coercion A_K --> A_L
            Also Ohat_K --> Ohat_L and Khat --> Lhat

            Define Coercion `\QQ_p` --> A_Q; even K_p --> A_K
            Maybe via Khat? So two coercions: K_p --> Khat --> A_K
        """
        debug("Adeles._coerce_map_from_({})".format(S))
        if self.base().has_coerce_map_from(S):
            return True
        J = IdeleGroup(self.base())
        if J.has_coerce_map_from(S):
            return True
        return False

    def is_integral_domain(self, proof=None):
        """
        Return ``False``, indicating that ``self`` is not an integral domain

        EXAMPLES::

            sage: Adeles().is_integral_domain()
            False
        """
        debug("Adeles.is_integral_domain()")
        return False





