
load("completions.sage")
load("profinite_numbers.sage")



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
        """
        debug("Adele._repr_()")
        rep = "("
        for x_oo in self.infinite:
            rep += repr(x_oo) + ", "
        rep += repr(self.finite) + ")"
        return rep

    def _richcmp_(self, other, op):
        r"""
        Compare ``self`` and ``other`` based on the relation ``op``

        We only implement equality and non-equality.

        We declare adeles equal if...?. TODO

        EXAMPLES::

        """
        debug("Adele._richcmp_({}, {})".format(self, other))
        from sage.structure.richcmp import op_EQ, op_NE
        if op == op_EQ:
            return self.numerator * other.denominator == self.denominator * other.numerator
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
            sage: ProfiniteNumbers(K)
            Adele Ring of Number Field in a with defining polynomial x^3 + x + 1
        """
        debug("Adeles._repr_()")
        return "Adele Ring of {}".format(self.base())

    def _latex_(self):
        r"""
        Return latex-formatted string representation of ``self``

        EXAMPLES::

            sage: K.<a> = NumberField(x^2+14)
            sage: Khat = ProfiniteNumbers(K)
            sage: latex(Khat)
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
        if y is None:
            if x in K:  # coercion from base()
                infinite = [phi(x) for L, phi in completions(K, oo)]
                return self.element_class(self, infinite, x)
        return self.element_class(self, x, y)

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

            Define Coercion `\QQ_p` --> A_Q
            Maybe even K_p --> A_K?
        """
        debug("Adeles._coerce_map_from_({})".format(S))
        if self.base().has_coerce_map_from(S):
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





