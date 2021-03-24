r"""
Idele Groups of Number Fields

This implementation of ideles is similar to how the Real Interval Field ``RIF``
implements the field of real numbers. Storing a single real number takes in
general an infinite amount of memory, which we don't have. Hence ``RIF`` chooses
to only represent intervals of the real line, instead of actual real numbers.
Working with very small intervals corresponds to working with high precision.

Similarly, we do not represent individual ideles, but rather open subsets of
the idele group. A user can work with arbitrary high precision "ideles" by using
very small open subsets in its calculations.

Let `K` be a number field. For a place `p` of `K`, write `K_p` for the
completion of `K` at `p`. Write `Z_p` for the valuation ring of `K_p`. Write
`Z_p^*` for the unit group of `Z_p` and `M_p` for its maximal ideal.

The idele group `J_K` of `K` is is the restricted product
`\prod_p'(K_p^*, Z_p^*)` of `K_p^*`'s with respect to `Z_p`'s, where `p` ranges
over all places (finite and infinite) of `K`.
Concretely, that means

.. MATH::

    J_K = \left\{ (x_p)_p \in \prod_p K_p^* : x_p \in Z_p^*
                  \text{ for all but finitely many } p \right\}

This is a multiplicative group with component-wise multiplication.

It is even a topological group, with basic open sets being translates of open
subgroups of the form

.. MATH::

    \prod_{p \in S} U_p \times \prod_{p \not\in S} Z_p^*

where `S` is a finite set of places of `K` and each `U_p` is of the form
`1+M_p^i` with `i \in \NN`. These basic open sets our what our
``Idele``'s  will represent.

.. TODO::

    Check all idele-code for edge cases:
        - finite empty
        - exact vs non-exact
        - infinite: None vs value
    Also check idele-code in other files (i.e. Adeles._from_idele())
"""

load("completions.sage")
load("ray_class_groups.sage")
load("adeles.sage")


from sage.structure.element import MultiplicativeGroupElement
class Idele(MultiplicativeGroupElement):
    r"""
    An idele: an element of an :class:`IdeleGroup` of a number field

    .. automethod:: __init__
    """

    def __init__(self, parent, exact, infinite, finite):
        r"""
        Construct an ``Idele``

        We denote the number field to which this idele belongs by ``K``.

        INPUT:

        - ``parent`` -- the :class:`IdeleGroup` to which this idele belongs
        - ``exact`` -- an element of ``K`` or ``None``; determines the exact
          value of this idele at all placed unspecified in ``infinite`` and
          ``finite``.
        - ``infinite`` -- a ``list`` or ``None``; the list should have as many
          entries as ``K`` has infinite places. The ``i``th entry specifies the
          value of this idele at the infinite place corresponding to
          ``K.places()[i]``. The entry should lie in ``RIF`` or ``CIF``
          (depending on the place being real/complex), or be ``None`` for
          "unspecified".
          Setting ``list`` to None is the same as passing an all-``None``'s
          list.
        - ``finite`` -- a ``dict`` or ``None``; ``None`` is the same as passing
          an empty ``dict``. The ``dict`` should contain (``key``, ``value``)
          pairs such that:
            * ``key`` is a prime ``p`` of ``K``; if ``K`` is `\QQ`, this is a
              prime number, else it's a prime fractional ideal of ``K``.
            * ``value`` should be a pair (``x``, ``i``) with ``x`` a non-zero
              element of ``K`` and ``i`` a non-negative integer.
          Such an entry means: at place ``p`` of ``K``, this idele is equal to
          (the image of) ``x`` (in `K_p`), up to elements in the open subgroup
          `1+M_p^i` of `K_p`.

        EXAMPLES:

        We begin with the easiest case: `K = \QQ`::

            sage: J = IdeleGroup(QQ)
            sage: Idele(J, None, [3.14], {2: (5/3, 2), 5: (1, 9)})
            (3.1400000000000002?, 5/3*(1+M_2^2), Z_3*, 1*(1+M_5^9), ...)
            sage: Idele(J, None, None, {7: (1/10, 100)})
            (RR*, Z_2*, Z_3*, Z_5*, 1/10*(1+M_7^100), ...)
            sage: Idele(J, 1/3, None, {7: (1/10, 100)})
            (0.3333333333333334?, 1/3, 1/3, 1/3, 1/10*(1+M_7^100), 1/3, 1/3, 1/3, ...)

        Now let's take a non-trivial number field::

        sage: K.<a> = NumberField(x^4-17)
        sage: Jk = IdeleGroup(K)
        sage: Idele(Jk, None, None, {K.prime_above(3): (a^3+a+1, 10)})
        (RR*, RR*, CC*, Z_p2*, Z_q2*, Z_r2*, (a^3 + a + 1)*(1+M_p3^10), Z_q3*, ...)
        where:
            p2 = Fractional ideal (2, 1/4*a^3 + 1/4*a^2 + 1/4*a - 3/4)
            q2 = Fractional ideal (2, 1/4*a^3 - 1/4*a^2 + 1/4*a + 3/4)
            r2 = Fractional ideal (2, -1/2*a^2 + a + 1/2)
            p3 = Fractional ideal (3, 1/2*a^2 - a - 1/2)
            q3 = Fractional ideal (3, 1/2*a^2 + a - 1/2)
        """
        #print("DEBUG Idele.__init__({}, {}, {}, {})".format(parent, exact, infinite, finite))
        self.exact = exact
        self.infinite = infinite
        self.finite = finite
        MultiplicativeGroupElement.__init__(self, parent)
        self._validate_exact()
        self._validate_infinite()
        self._validate_finite()

    def _validate_exact(self):
        """
        Check if the ``exact`` attribute is valid and adjust it if feasable

        Valid means: ``exact`` is None, or is a non-zero element of `K`.

        We throw an exception if exact is invalid.

        If ``exact`` can be coerces into `K`, but is not an actual element of
        `K`, we cast it into an actual element of `K`.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3-2)
            sage: J = IdeleGroup(K)
            sage: u = J(1)
            sage: u.exact = None
            sage: u._validate_exact()
            sage: u.exact = QQ(3/2)
            sage: u._validate_exact()
            sage: u.exact.parent()
            Number Field in a with defining polynomial x^3 - 2
            sage: u.exact = []
            sage: u._validate_exact()
            Traceback (most recent call last):
            ...
            TypeError: exact should be an element of the number field
        """
        if self.exact is None:
            return
        K = self.parent().number_field
        if self.exact not in K:
            raise TypeError("exact should be an element of the number field")
        self.exact = K(self.exact)
        if self.exact.is_zero():
            raise ValueError("exact must be a unit (i.e. non-zero)")

    def _validate_infinite(self):
        """
        Check that the ``infinite`` attribute is valid and adjust if faesable

        Valid means: ``infinite`` is ``None``, or a list of length
        ``len(K.places())`` with the `i`-th entry an element of
        ``RIF`` or ``CIF`` (depending on the place being real/complex) or
        ``None``.

        If the `i`-th entry's parent is not equal to ``RIF``/``CIF``, but does
        lie in it, than we cast the element to an actual element of such an
        interval field.

        Upon invalid (and uncastable) ``infinite``, we throw an exception.

        EXAMPLES:

        We look at a number field with one real infinite place and one complex
        one::

            sage: K.<a> = NumberField(x^3 + x + 1)
            sage: Jk = IdeleGroup(K)
            sage: u = Jk(1)
            sage: u.infinite = [1.2, I]
            sage: u._validate_infinite()

        Below the ``Integer 3`` will be cast into ``RIF``::

            sage: u.infinite = [ZZ(3), I]
            sage: u._validate_infinite()
            sage: u.infinite[0].parent()
            Real Interval Field with 53 bits of precision

        And here we cast ``None`` to a list of ``None``s::

            sage: u.infinite = None
            sage: u._validate_infinite()
            sage: u.infinite
            [None, None]

        The exceptions that we may throw:

            sage: u.infinite = [3.14, I, I]
            sage: u._validate_infinite()
            Traceback (most recent call last):
            ...
            TypeError: infinite should have length 2
            sage: u.infinite = {CIF: I+1}
            sage: u._validate_infinite()
            Traceback (most recent call last):
            ...
            TypeError: infinite should be a list
            sage: u.infinite = [I, I]
            sage: u._validate_infinite()
            Traceback (most recent call last):
            ...
            TypeError: 0th infinite value (I) should lie in Real Interval Field with 53 bits of precision
        """
        K = self.parent().number_field
        K_oo = infinite_completions(K)
        t = len(K_oo)
        if self.infinite is None:
            self.infinite = [None for i in range(t)]
            return
        if not isinstance(self.infinite, list):
            raise TypeError("infinite should be a list")
        if len(self.infinite) != t:
            raise TypeError("infinite should have length {}".format(t))
        for i in range(t):
            val = self.infinite[i]
            if val is not None:
                if val not in K_oo[i][0]:
                    raise TypeError("{}th infinite value ({}) should lie in {}".format(i, val, K_oo[i][0]))
                if val.is_zero():
                    raise ValueError("{}th infinite value ({}) should be a unit (i.e. non-zero)".format(i, val))
                self.infinite[i] = K_oo[i][0](val)

    def _validate_finite(self):
        r"""
        Check if the ``finite`` attribute is valid and adjust if feasable

        Valid means: ``None`` or a ``dict`` with keys primes of `K` and values
        pairs in `(K, \NN)`.

        We throw an exception upon invalid ``finite``.

        EXAMPLES:

        A few examples of the small adjustments we can do::

            sage: K.<a> = NumberField(x^4+x+1)
            sage: J = IdeleGroup(K)
            sage: u = J(1)
            sage: u.finite = None
            sage: u._validate_finite()
            sage: u.finite
            {}
            sage: u.finite = {K.prime_above(2): (1, 1)}
            sage: u._validate_finite()
            sage: u.finite[K.prime_above(2)][0].parent()
            Number Field in a with defining polynomial x^4 + x + 1

        And some exceptions we may throw::

            sage: u.finite = "blah"
            sage: u._validate_finite()
            Traceback (most recent call last):
            ...
            TypeError: finite should be a dict
            sage: u.finite = {40: (1, 1)}
            sage: u._validate_finite()
            Traceback (most recent call last):
            ...
            TypeError: keys of finite should be prime ideals of the number field
            sage: u.finite = {K.prime_above(2): (1, -1)}
            sage: u._validate_finite()
            Traceback (most recent call last):
            ...
            TypeError: values at finite primes should be pairs in K x NN
        """
        if self.finite is None:
            self.finite = {}
            return
        if not isinstance(self.finite, dict):
            raise TypeError("finite should be a dict")
        K = self.parent().number_field
        already_seen = set()
        new_finite = {}
        for q, val in self.finite.items():
            if q not in K.ideal_monoid() or not q.is_prime():
                if K is QQ:
                    raise TypeError("keys of finite should be prime numbers")
                else:
                    raise TypeError("keys of finite should be prime ideals of the number field")
            try:
                v, i = val
            except TypeError:
                raise TypeError("values at finite primes should be pairs in K x NN")
            if v not in K or i not in NN:
                raise TypeError("values at finite primes should be pairs in K x NN")
            if v.is_zero():
                raise ValueError("values at finite primes must be units (i.e. non-zero)")
            if q in already_seen:
                raise TypeError("multiple values specified for finite prime {}".format(q))
            already_seen.add(q)

            if K is QQ:
                new_q = ZZ(q)
            else:
                new_q = K.ideal(q)
            if not new_q.is_prime():
                raise ValueError("{} is not prime as an ideal of {}".format(q, K))
            new_finite[new_q] = (K(v), ZZ(i))
        self.finite = new_finite

    def _repr_(self):
        """
        Return a representation of this idele.

        If ``self`` only stored values at primes that lie above rational primes
        smaller than 50, we return a dense representation (see
        :meth:`_repr_dense`).
        Else, we return a sparse representation (see :meth:`_repr_sparse`).
        """
        K = self.parent().number_field
        for q in self.finite:
            p = q if K is QQ else q.gens_two()[0] # the rational prime below q
            if p >= 50:
                return self._repr_sparse()
        return self._repr_dense()

    def _repr_sparse(self):
        """
        Return a sparse representation of this idele.

        EXAMPLES::

            sage: J = IdeleGroup(QQ)
            sage: u = J(1/2, None, {10000079: (7/9, 20)})
            sage: 
            Idele over Rational Field with values
                    10000079: 7/9 * U_q^20
            and which equals exactly 1/2 at all other primes
            sage: J(None, [7.9], {10000079: (7/9, 20)})
            Idele over Rational Field with values
                    infinity_0: 7.9000000000000004?
                    10000079: 7/9 * U_q^20
        """
        K = self.parent().number_field
        rep = "Idele over {} with values\n".format(K)
        for i in range(len(K.places())):
            x_oo_i = self.infinite[i]
            if x_oo_i is not None:
                rep += "\tinfinity_{}: {}\n".format(i, x_oo_i)
        for q in self.finite:
            x_q, i_q = self.finite[q]
            q_name = q if K is QQ else q.gens_reduced()
            rep += "\t{}: {} * U_q^{}\n".format(q_name, x_q, i_q)
        if self._has_exact():
            rep += "and which equals exactly {} at all other primes.".format(self.exact)
        return rep[0:-1]


    def _repr_dense(self):
        """
        Return a dense representation of this idele

        We represent the idele "densely": first we print (the values at) all
        infinite primes, then at the primes above 2, then at the primes above
        3, then at the primes above 5, etc. Each time we use the order as
        returned by :meth:`NumberField.primes_above`. We keep printing until
        we printed all stored values. All not-stored values in between are
        printed as "RR*", "CC*" or "Z_p*".
        
        EXAMPLES::

            sage: K.<a> = NumberField(x^2+5)
            sage: J = IdeleGroup(K)
            sage: J(a, [2*I], {K.prime_above(3): (a/2+1, 8)})
            (2*I, a, (1/2*a + 1)*(1+M_p3^8), a, a, a, a, ...)
            where:
                    p3 = Fractional ideal (3, a + 1)
            sage: J(None, None, {K.prime_above(5): (-a-7, 0)})
            (CC*, Z_p2*, Z_p3*, Z_q3*, (-a - 7)*Z_p5*, ...)
            where:
                    p2 = Fractional ideal (2, a + 1)
                    p3 = Fractional ideal (3, a + 1)
                    q3 = Fractional ideal (3, a + 2)
                    p5 = Fractional ideal (-a)
        """
        if not self.finite and all([x is None for x in self.infinite]):
            return str(self.exact)
        rep = "("
        K = self.parent().number_field
        K_oo = infinite_completions(K)
        for i in range(len(K_oo)):
            val = self.infinite[i]
            if val is None:
                if self._has_exact():
                    phi = K_oo[i][1] # phi: K --> RIF/CIF
                    rep += "{}, ".format(phi(self.exact))
                elif K_oo[i][0] is CIF:
                    rep += "CC*, "
                else:
                    rep += "RR*, "
            else:
                rep += "{}, ".format(val)
        start = 0
        step = 100
        to_do = set(self.finite)
        if K is QQ:
            while to_do:
                for p in prime_range(start, start+step):
                    if p in to_do:
                        v, i = self.finite[p]
                        if i > 0:
                            rep += "{}*(1+M_{}^{}), ".format(v, p, i)
                        else:
                            rep += "{}*Z_{}*, ".format(v, p)
                        to_do.remove(p)
                        if not to_do:
                            break
                    elif self.exact is None:
                        rep += "Z_{}*, ".format(p)
                    else:
                        rep += "{}, ".format(self.exact)
                start += step
        else:
            legend = "\nwhere:"
            while to_do:
                for p in prime_range(start, start+step):
                    qs = K.primes_above(p)
                    for i in range(len(qs)):
                        q_str = self._prime_name(p, i)
                        if qs[i] in to_do:
                            v, w = self.finite[qs[i]]
                            if sum([c != 0 for c in v.list()]) > 1:
                                v = "(" + str(v) + ")"
                            if w > 0:
                                rep += "{}*(1+M_{}^{}), ".format(v, q_str, w)
                            else:
                                rep += "{}*Z_{}*, ".format(v, q_str)
                            legend += "\n\t{} = {}".format(q_str, qs[i])
                            to_do.remove(qs[i])
                        elif self.exact is None:
                            rep += "Z_{}*, ".format(q_str)
                            legend += "\n\t{} = {}".format(q_str, qs[i])
                        else:
                            rep += "{}, ".format(self.exact)
                    if not to_do:
                        break
                start += step
        if self._has_exact():
            rep += "{}, {}, {}, ".format(self.exact, self.exact, self.exact)
        rep += "...)"
        if K is not QQ:
            rep += legend
        return rep

    def _prime_name(self, p, i):
        r"""
        Return our name for the ``i``-th prime above the rational prime ``p``

        INPUT:

        - `p` -- a prime number in `\ZZ`
        - `i` -- a non-negative integer

        OUTPUT:

        A short name for the prime ``K.primes_above(p)[i]``.

        For `X` a prime number, we call the 0th, 1th, 2th, etc. primes above `X`
        as follows: "pX", "qX", "rX", "sX", "tX", "ppX", "pqX", "prX", etc.
        
        We do not validate the input. For example, we don't check if ``p`` is a
        prime number and we do not check if ``self.parent().number_field``
        actually has at least ``i+1`` primes above ``p``.

        EXAMPLES::

            sage: K.<a> = NumberField(x^6+x+1)
            sage: J = IdeleGroup(K)
            sage: u = J(1)
            sage: u._prime_name(2, 0)
            'p2'
            sage: for i in range(12):
            ....:     u._prime_name(5, i)
            'p5'
            'q5'
            'r5'
            's5'
            't5'
            'pp5'
            'pq5'
            'pr5'
            'ps5'
            'pt5'
            'qp5'
            'qq5'
            sage: u._prime_name(100, 0)
            'p100'
        """
        if i == 0:
            return "p" + str(p)
        i = ZZ(i)
        result = ""
        digits = i.digits(base=5)
        if len(digits) > 1:
            digits[-1] -= 1
        for d in digits:
            result = chr(d + ord("p")) + result
        return result + str(p)

    def _equals(self, other):
        """
        Return whether or not ``self`` equals ``other``

        We define equality of ideles very loosely: if two ideles *could* be
        equal, that is if the open subsets they represent have non-empty
        intersection, than we declare them equal.

        EXAMPLES::

            sage: J = IdeleGroup(QQ)
            sage: u = J(None, None, {2: (3, 7), 3: (1/2, 6)})
            sage: v = J(4, [-1.23], {3: (5, 10)})
            sage: w = J(None, None, {2: (3, 7), 3: (1/2, 6), 5: (1, 0)})
            sage: u._equals(u)
            True
            sage: (u/u)._equals(J(1))
            True
            sage: u._equals(v)
            False
            sage: u._equals(w)
            True
            sage: v._equals(w)
            False

        If an idele ``x`` represents a strict subset of another idele ``v``,
        they are equal::

            sage: x = J(4, [-1.23], {3: (5, 11)})
            sage: v._equals(x)
            True

        Also in this case the result is not equal, due to the non-exactness of
        ``u``::
        """
        K = self.parent().number_field

        if self._has_exact() and other._has_exact():
            if self.exact != other.exact:
                return False
        elif self._has_exact() or other._has_exact():
            ex, nex = (self, other) if self._has_exact() else (other, self)
            if not nex._contains(ex.exact):
                return False

        t = len(self.infinite)
        for i in range(t):
            if self.infinite[i] is None:
                if other.infinite[i] is not None:
                    return False
            else:
                if other.infinite[i] is None:
                    return False
                try:
                    self.infinite[i].intersection(other.infinite[i])
                except ValueError:
                    # This indicates the intersection is empty
                    return False

        for q, v in self.finite.items():
            if q in other.finite:
                w = other.finite[q]
                small, big = (self, other) if v[1] >= w[1] else (other, self)
                if not big._contains_at(small.finite[q][0], q):
                    # the small open subset does not lie in the big open subset
                    # so they are disjoint
                    return False
            else: # other[q] not stored
                if other._has_exact():
                    if not self._contains_at(other.exact, q):
                        return False
                else:
                    x = v[0]
                    # other[q] is Z_q*
                    # So x must lie in Z_q*, i.e. have valuation zero
                    if K is not QQ:
                        x = K.ideal(x)
                    if  x.valuation(q) != 0:
                        return False
        for q, v in other.finite.items():
            if q not in self.finite:
                if self._has_exact():
                    # self[q] is exact, while other[q] is an approximation.
                    return False
                else:
                    x = v[0]
                    # self[q] is Z_q*
                    # So x must lie in Z_q*, i.e. have valuation zero
                    if K is not QQ:
                        x = K.ideal(x)
                    if x.valuation(q) != 0:
                        return False

        # Passed all checks. We are equal!
        return True

    def _richcmp_(self, other, op):
        """
        Return the result of operator ``op`` applied to ``self`` and ``other``

        Only equality and inequality are implented.

        EXAMPLES::

            sage: K.<a> = NumberField(x^5-x+2)
            sage: J = IdeleGroup(K)
            sage: u = J(a^4+1, None, {K.prime_above(11): (a^3-a, 20)})
            sage: v = J(None, None, {K.prime_above(97): (a, 20)})
            sage: v == u
            False
            sage: v != u
            True
            sage: u == u
            True
            sage: u != u
            False
        """
        from sage.structure.richcmp import op_EQ, op_NE
        if op == op_EQ:
            return self._equals(other)
        if op == op_NE:
            return not self._equals(other)
        raise NotImplementedError()

    def _has_exact(self):
        """
        Return whether or not our ``exact`` attribute is set

        EXAMPLES::

            sage: J = IdeleGroup(QQ)
            sage: u = J(None, [3.1415])
            sage: u._has_exact()
            False
            sage: u.exact = 1/2
            sage: u._has_exact()
            True
        """
        return self.exact is not None
    
    def _mul_(self, other):
        """
        Multiply ``self`` and ``other`` together

        .. NOTE::

            Precision may be lost. See the examples below.

        EXAMPLES:

        First a few multiplications of rational ideles::

            sage: J = IdeleGroup(QQ)
            sage: u = J(None, None, {2: (1/2, 7), 3: (2/5, 8)})
            sage: v = J(None, [1.2345], None)
            sage: w = J(1/7, [-1.0], {3: (-1, 7)})
            sage: u*v
            (RR*, 1/2*Z_2*, 2/5*Z_3*, ...)
            sage: u*w
            (RR*, 1/14*(1+M_2^7), -2/5*(1+M_3^7), Z_5*, 1/7*Z_7*, ...)
            sage: w*u
            (RR*, 1/14*(1+M_2^7), -2/5*(1+M_3^7), Z_5*, 1/7*Z_7*, ...)
            sage: v*w
            (-1.2345000000000000?, Z_2*, -1*Z_3*, Z_5*, 1/7*Z_7*, ...)

        And now a few over a non-trivial number field::

            sage: K.<i> = NumberField(x^2+1)
            sage: J = IdeleGroup(K)
            sage: p2, p3 = K.prime_above(2), K.prime_above(3)
            sage: u = J(i+1)
            sage: v = J(None, [I+1], {p3: (i/2, 7)})
            sage: w = J(i-1, None, {p2: (2, 5), p3: (3*i, 20)})
            sage: u*v
            (2*I, (i + 1)*Z_p2*, (1/2*i - 1/2)*(1+M_p3^7), ...)
            where:
                    p2 = Fractional ideal (i + 1)
                    p3 = Fractional ideal (3)
            sage: u*w
            (-2, (2*i + 2)*(1+M_p2^5), (3*i - 3)*(1+M_p3^20), -2, -2, -2, ...)
            where:
                    p2 = Fractional ideal (i + 1)
                    p3 = Fractional ideal (3)
            sage: v*w
            (-2, 2*Z_p2*, -3/2*(1+M_p3^7), ...)
            where:
                    p2 = Fractional ideal (i + 1)
                    p3 = Fractional ideal (3)
        """
        exact = None
        if self._has_exact() and other._has_exact():
            exact = self.exact * other.exact

        infinite = self.infinite.copy()
        K = self.parent().number_field
        K_oo = infinite_completions(K)
        for i in range(len(K_oo)):
            if infinite[i] is None:
                if other.infinite[i] is None:
                    pass # keep None, depends on exact which is already set correctly
                else: # other has infinite value
                    if self._has_exact():
                        phi = K_oo[i][1] # phi: K --> RIF/CIF
                        infinite[i] = phi(self.exact) * other.infinite[i]
                    else:
                        pass
            else: # infinite[i] set
                if other.infinite[i] is None:
                    if other._has_exact():
                        phi = K_oo[i][1] # phi: K --> RIF/CIF
                        infinite[i] = infinite[i] * phi(other.exact)
                    else: # other.infinite[i] unknown
                        infinite[i] = None
                else: # both values set
                    infinite[i] *= other.infinite[i]

        finite = self.finite.copy()
        for q, val in finite.items():
            if q in other.finite:
                x, i = val
                y, j = other.finite[q]
                finite[q] = (x*y, min(i, j))
            else:
                if other._has_exact():
                    finite[q] = (val[0] * other.exact, val[1])
                else: # other is unknown (i.e. Z_q^*) at q
                    # U_q^i is a subgroup of Z_q^*, so x*U_q^i * Z_q^* = x*Z_q^*
                    finite[q] = (val[0], ZZ(0)) 
        for q, val in other.finite.items():
            if q not in finite:
                if self._has_exact():
                    finite[q] = (val[0] * self.exact, val[1])
                else: # self is unknown (i.e. Z_q^*) at q
                    # U_q^i is a subgroup of Z_q^*, so x*U_q^i * Z_q^* = x*Z_q^*
                    finite[q] = (val[0], ZZ(0))
        if ((self._has_exact() and not other._has_exact())
                or (not self._has_exact() and other._has_exact())):
            # One has an exact value and one is Z_q^* almost everywhere.
            value = self.exact if self._has_exact() else other.exact
            # Where value has valuation 0, we have value * Z_q^* = Z_q^*, so we
            # don't have to store this.
            # But at the finitely many places where value has non-zero
            # valuation, we must store "value * Z_q^*".
            if K is QQ:
                I = value
            else:
                I = K.ideal(value)
            for q, e in factor(I):
                if q not in finite:
                    finite[q] = (value, ZZ(0))

        return self.__class__(self.parent(), exact, infinite, finite)

    def inverse(self):
        """
        Return the inverse of ``self`` in its idele group (i.e. "1/``self``")

        EXAMPLES::

            sage: K.<a> = NumberField(x^2+x+1)
            sage: p2, p3 = K.prime_above(2), K.prime_above(3)
            sage: J = IdeleGroup(K)
            sage: u = J(None, [I], {p2: (a, 10), p3: (1/2, 7)})
            sage: u.inverse()
            (-1*I, (-a - 1)*(1+M_p2^10), 2*(1+M_p3^7), ...)
            where:
                    p2 = Fractional ideal (2)
                    p3 = Fractional ideal (-2*a - 1)
            sage: v = J(5, None, {p3: (a+1, 97)})
            sage: v.inverse()
            (0.2000000000000000?, 1/5, -a*(1+M_p3^97), 1/5, 1/5, 1/5, ...)
            where:
                    p3 = Fractional ideal (-2*a - 1)
        """
        exact = None
        if self._has_exact():
            exact = 1/self.exact

        infinite = self.infinite.copy()
        for i in range(len(infinite)):
            if infinite[i] is not None:
                infinite[i] = 1/infinite[i]

        finite = self.finite.copy()
        for q, val in finite.items():
            finite[q] = (1/val[0], val[1])

        return self.__class__(self.parent(), exact, infinite, finite)

    def _div_(self, other):
        """
        Divide ``self`` by ``other``

        EXAMPLES::

            sage: K.<a> = NumberField(x^2+8)
            sage: p2, q3 = K.prime_above(2), K.primes_above(3)[1]
            sage: J = IdeleGroup(K)
            sage: u = J(None, [I+1], {p2: (a-1, 7), q3: (a/2, 9)})
            sage: v = J(a, None, {q3: (2, 10)})
            sage: u/v
            (0.35355339059327374? - 0.35355339059327374?*I, (1/8*a + 1)*(1+M_p2^7), Z_p3*, 1/4*a*(1+M_q3^9), ...)
            where:
                    p2 = Fractional ideal (-1/2*a)
                    p3 = Fractional ideal (1/2*a + 1)
                    q3 = Fractional ideal (-1/2*a + 1)
            sage: v/u
            (1.4142135623730950? + 1.4142135623730950?*I, (-1/9*a + 8/9)*(1+M_p2^7), Z_p3*, -1/2*a*(1+M_q3^9), ...)
            where:
                    p2 = Fractional ideal (-1/2*a)
                    p3 = Fractional ideal (1/2*a + 1)
                    q3 = Fractional ideal (-1/2*a + 1)
            sage: u/u
            (1, 1*(1+M_p2^7), Z_p3*, 1*(1+M_q3^9), ...)
            where:
                    p2 = Fractional ideal (-1/2*a)
                    p3 = Fractional ideal (1/2*a + 1)
                    q3 = Fractional ideal (-1/2*a + 1)

        .. TODO::
            
            Fix the following bug::

                sage: J = IdeleGroup(QQ)
                sage: u = J(None, None, {2: (2, 3)})
                sage: u/J(2)
                (RR*, 1*(1+M_2^3), ...)
                sage: u/2
                Traceback (most recent call last):
                ...
                KeyError: (Idele Group of Rational Field, Integer Ring, <built-in function truediv>)

                During handling of the above exception, another exception occurred:

                Traceback (most recent call last):
                ...
                KeyError: 'gens'

                During handling of the above exception, another exception occurred:

                Traceback (most recent call last):
                ...
                AttributeError: 'IdeleGroup_with_category' object has no attribute 'gens'

        """
        return self * other.inverse()

    def integral_split(self):
        """
        Return a tuple ``(u, d)`` with ``u` an "integral" idele and ``d`` an
        integer such that ``self = u/d``

        The "integrality" of ``u`` is strong, in the following sense: at every
        finite prime `q` at which we stored a value `(x, i)` (representing
        `x*U_q^i`), the `K`-element ``x`` is integral in `K`; and if ``u`` is
        exact, its exact value is also integral in `K`.
        This is stronger than ``u`` being locally everywhere integral.

        EXAMPLES::

            sage: J = IdeleGroup(QQ)
            sage: u = J(None, None, {2: (1/4, 1), 3: (3/5, 2)})
            sage: u.integral_split()
            ((RR*, 5*(1+M_2^1), 12*(1+M_3^2), 5*Z_5*, ...), 20)
            sage: v = J(3/14, [1], {2: (3/4, 5), 3: (1, 3)})
            sage: v.integral_split()
            ((28, 21*(1+M_2^5), 28*(1+M_3^3), 6, 6, 6, ...), 28)

            sage: K.<a> = NumberField(x^2+10)
            sage: p3, p5 = K.prime_above(3), K.prime_above(5)
            sage: Jk = IdeleGroup(K)
            sage: u = Jk(None, [I], {p3: (1/a, 3), p5: (7/6, 8)})
            sage: u.integral_split()
            ((30*I, 30*Z_p2*, -3*a*(1+M_p3^3), 35*(1+M_p5^8), ...)
             where:
                    p2 = Fractional ideal (2, a)
                    p3 = Fractional ideal (3)
                    p5 = Fractional ideal (5, a),
             30)
        """
        K = self.parent().number_field

        denominators = []
        factors = factor(1)
        for q, val in self.finite.items():
            d_q = val[0].denominator()
            denominators.append(d_q)
            if K is not QQ:
                d_q = K.ideal(d_q)
            factors *= factor(d_q)
        if self.exact is not None:
            d_exact = self.exact.denominator()
            denominators.append(d_exact)
            if K is not QQ:
                d_exact = K.ideal(d_exact)
            factors *= factor(d_exact)
        d = lcm(denominators)

        # Now returning (d*self, d) would be correct. The multiplication will
        # usually need to factor d, which may be expensive. Less expensive is
        # to factor all the individual d_q's (as we did above) and build up our
        # ``u`` directly. This is what we will do below.

        finite = {}
        for q, val in self.finite.items():
            finite[q] = (d*val[0], val[1])
        exact = None
        if self._has_exact():
            exact = self.exact * d
        else:
            for q, e in factors:
                if q not in finite:
                    finite[q] = (d, 0)
        infinite = self.infinite.copy()
        for i in range(len(infinite)):
            if infinite[i] is not None:
                infinite[i] = d*infinite[i]
        u = self.__class__(self.parent(), exact, infinite, finite)

        return (u, d)



    def to_modulo_element(self):
        r"""
        Convert this idele to its image in O/I, with O the maximal order of K
        and I the (biggest) ideal where ``self`` is defined

        This implements the canonical homomorphism `\hat{O} \to O/I`, where
        we view `\hat{O}` as the subring of the idele group of ideles which are
        everywhere integral.

        Throws an exception if ``self`` does not lie in this subring.

        INPUT:

        - ``I`` -- an ideal of the maximal order ``O`` of our number field ``K``

        EXAMPLES:
            
        sage: J = IdeleGroup(QQ)
        sage: u = J(None, None, {2: (1/3, 3), 5: (7, 1)})
        sage: u_bar = u.to_modulo_element()
        sage: u_bar, u_bar.parent()
        (27, Ring of integers modulo 40)
        """
        A = Adeles(self.parent().number_field)
        return A(self).to_modulo_element()

    def to_ray_class(self, modulus):
        r"""
        Return the image of ``self`` in the ray class group modulo ``modulus``

        Let `K` denote our number field and let `O` denote its maximal order.
        Let `I(modulus)` be the group of fractional `O`-ideals coprime to
        ``modulus`` and let `R(modulus)` denote the ray modulo ``modulus``.
        Hence the ray class group of `K` modulo ``modulus`` is
        `I(modulus)/R(modulus)`.
        Write `J_K` for the group of ideles of `K` (to which ``self`` belongs).

        This method implementst the homomorphism

        .. MATH::

            \phi: J_K \to I(modulus)/R(modulus)

        that sends a prime element at a finite prime `q` (not dividing `m`) to
        `q \mod R(modulus)`.

        INPUT:

        - ``modulus`` -- a :class:`Modulus` of `K`

        OUTPUT:

        The image of ``self`` under `\phi`.

        EXAMPLES::

            sage: Q = NumberField(x-1, "one")
            sage: J = IdeleGroup(Q)
            sage: u = J(None, [-9], {2: (10, 1), 3: (3, 2), 7: (1/2, 4)})
            sage: m = Modulus(Q.ideal(18), [0])
            sage: u.to_ray_class(m).ideal()
            Fractional ideal (7)

        ::

            sage: K.<a> = NumberField(x^3-x-1)
            sage: Jk = IdeleGroup(K)
            sage: p7, q7 = K.primes_above(7)
            sage: u = Jk(None, [2, I], {p7: (1/2, 1), q7: (7*a, 2)})
            sage: m = Modulus(K.ideal(7), [0])
            sage: u.to_ray_class(m).ideal()
            Fractional ideal (-a^2 - 3*a + 2)

        TESTS:

            sage: K = NumberField(x^2 + 5*x - 3, 'a')
            sage: J = IdeleGroup(K)
            sage: G = ray_class_group(K, Modulus(K.ideal(14), [0])); G
            Ray class group of order 18 with structure C6 x C3 of Number Field in a with defining polynomial x^2 + 5*x - 3 of modulus (Fractional ideal (14)) * infinity_0
            sage: c0, c1 = G.gens()
            sage: for e in range(6):
            ....:     for f in range(3):
            ....:         r = c0^e * c1^f
            ....:         assert G(J(r)) == r, "bug!"  # long

        .. TODO::

            Implement the case that ``self`` has exact.
        
        ALGORITHM:

        We construct an idele `v` satisfying:
            - `v \equiv self \mod K^*`
            - `v \equiv 1 \mod^* modulus`
        We do this by starting of with `v = self` and then improve `v` in three
        steps:
            1. Make `v` integral at the primes dividing ``modulus``.
            2. Make `v \equiv 1 \mod^*` ``modulus.finite_part()``
            3. Make `v \equiv 1 \mod^*` ``modulus`` (fix the infinite part)
        In every step we only change `v` by multiplying it with principal ideles
        (i.e. elements of K^*) as to never violate the first desired condition
        on `v`.
        For Step 2 we use the number field version of the Chinese Remainder
        Theorem (cf. :func:`solve_CRT`). For Step 3 we use
        :meth:`get_one_mod_star_finite_with_fixed_signs`.

        Once we have such a `v`, we return the image of the ideal
        `\prod_{p} p^{\ord_p(v)}` in ``G``, where `p` ranges over the finite
        primes of `K`.
        """
        if self._has_exact():
            raise NotImplementedError("to_ray_class() not implemented yet for exact ideles")

        J = self.parent()
        K = J.number_field
        G = ray_class_group(K, modulus)

        # First we check if the precision of this idele is high enough to have
        # a well-defined image in the ray class group ``G``.
        for i in modulus.infinite_part():
            if (self.infinite[i] is None
                    or not (self.infinite[i] <= 0 or self.infinite[i] >= 0)
                    or self.infinite[i].is_zero()):
                raise ValueError("idele has no well-defined sign at infinite prime {}".format(i))
        for q, e in modulus.finite_factors():
            if q not in self.finite or self.finite[q][1] < e:
                raise ValueError("idele must be known up to at least U_q^{} at q={}".format(e, q))

        v = self
        # Step 1. We make `v` integral at the primes dividing modulus.
        for q, e in modulus.finite_factors():
            x, i = self.finite[q]
            if x not in K.maximal_order():
                v *= x.denominator()

        # Step 2. We find an element y in the maximal order O_K of K such that
        # y \equiv v \mod modulus.finite_part() using the Chinese Remainder
        # Theorem.
        values = []
        moduli = []
        for q, e in modulus.finite_factors():
            x, i = v.finite[q]
            f = i + x.valuation(q)
            # v represents x*U_q^i at q, which equals x+q^f\Z_q because
            # i >= modulus.valuation(q) >= 1 and x.valuation(q) >= 0 (so we
            # do not need to worry about the case i==0 representing x*Z_q^*)
            values.append(x)
            moduli.append(q^f)
        y = K.solve_CRT(values, moduli)
        if not y.is_zero():
            # y is zero if and only if values and moduli are empty, so when
            # modulus is infinite. In that case we can just leave v as it is.
            v /= J(y)
        # Now we have v \equiv 1 \mod^* modulus.finite_part().
        

        # Step 3. We address the infinite part using the Modulus-method
        # get_one_mod_star_finite_with_fixed_signs().
        positive = []
        negative = []
        for i in modulus.infinite_part():
            if v.infinite[i] >= 0:
                positive.append(i)
            else:
                negative.append(i)
        t = modulus.get_one_mod_star_finite_with_fixed_signs(positive, negative)
        v *= t

        # TODO do not do the redundant check below to save its costs
        assert v.is_one_mod_star(modulus), r"Assertion error in Idele.to_ray_class(): v \not\equiv 1 \mod^* m"

        # Our `v` is finished. We can now build up an ideal representing the
        # image of `v` (which is equal to the image of ``self``) in ``G``.
        I = K.unit_ideal()
        for q, val in v.finite.items():
            I *= q^(val[0].valuation(q))

        return G(I)


    def is_one_mod_star(self, modulus):
        r"""
        Return whether or not ` ``self`` \equiv 1 \mod^\ast ``modulus`` ` holds

        EXAMPLES::

            sage: Q = NumberField(x-1, "one")
            sage: J = IdeleGroup(Q)
            sage: u = J(None, None, {5: (6, 2)})
            sage: u.is_one_mod_star(Modulus(Q.ideal(3)))
            False
            sage: u.is_one_mod_star(Modulus(Q.ideal(5)))
            True
            sage: u.is_one_mod_star(Modulus(Q.ideal(25)))
            False
            sage: u.is_one_mod_star(Modulus(Q.ideal(5), [0]))
            False
            sage: u.infinite[0] = RIF(1.2345)
            sage: u.is_one_mod_star(Modulus(Q.ideal(5), [0]))
            True
            sage: u.infinite[0] = RIF(-1.2345)
            sage: u.is_one_mod_star(Modulus(Q.ideal(5), [0]))
            False
            sage: u.exact = Q.one()
            sage: u.is_one_mod_star(Modulus(Q.ideal(2*3*5*7*11)))
            True
        """
        for q, e in modulus.finite_factors():
            if not q in self.finite:
                if not self._has_exact():
                    return False
                x = self.exact
            else:
                x, i = self.finite[q]
                if i < e:
                    return False
            if (x-1).valuation(q) < e:
                return False
        for i in modulus.infinite_part():
            x = self.infinite[i]
            if x is None:
                return False
            if not (x >= 0):
                return False
        return True

    def _contains_at(self, x, q):
        r"""
        Check if this idele contais the `K`-element ``x`` at the prime ``q``

        Suppose that the open subset this idele represents is `U_q^i` at ``q``.
        Then we check whether or not `x \in U_q^i` holds.

        INPUT:

        - ``x`` -- element of `K`
        - ``q`` - prime (finite or infinite) of `K`; in the infinite case, it
          should be one of the embbedings returned by ``K.places()``.

        EXAMPLES:
        
        First some examples with `\QQ`::

            sage: J = IdeleGroup(QQ)
            sage: u = J(None, None, {2: (5, 3), 3: (1/3, 0)})
            sage: u._contains_at(5, 2)
            True
            sage: u._contains_at(2, 2)
            False
            sage: u._contains_at(5/3, 3)
            True
            sage: u._contains_at(5/9, 3)
            False
            sage: u._contains_at(13, 5)
            True
            sage: u._contains_at(130, 5)
            False
            sage: u._contains_at(3.1415, QQ.places()[0])
            True

        And some examples with a non-trivial number field::

            sage: K.<a> = NumberField(x^3+x+1)
            sage: Jk = IdeleGroup(K)
            sage: v = Jk(a, [None, I], {K.prime_above(7): (a+1, 8)})
            sage: v._contains_at(a, K.places()[0])
            True
            sage: v._contains_at(a+1, K.places()[0])
            False
            sage: v._contains_at(a, K.places()[1])
            False
            sage: v._contains_at(a, K.prime_above(7))
            False
            sage: v._contains_at((a+1)*(7^8+1), K.prime_above(7))
            True

        And lastly some errors that may occur:

            sage: u._contains_at(1, 8)
            Traceback (most recent call last):
            ...
            TypeError: ``q`` should be a prime of `K`, but got 8
            sage: v._contains_at(a^2, None)
            Traceback (most recent call last):
            ...
            TypeError: ``q`` should be a prime of `K`, but got None
        """
        K = self.parent().number_field
        if q in K.places():
            i = K.places().index(q)
            if self.infinite[i] is not None:
                phi = completions(K, oo)[i][1] # phi: K --> RIF/CIF
                return phi(x) in self.infinite[i]
            elif self._has_exact():
                return x == self.exact
            else:
                # Not stored represents RR* or CC*, which ``x`` always lies in
                return True

        if (K is QQ and q in Primes()) or (K is not QQ and q in K.ideal_monoid() and q.is_prime()):
            if q in self.finite:
                y, i = self.finite[q]
                if i == ZZ(0):
                    return (x/y).valuation(q) == ZZ(0)
                else:
                    return (x/y - 1).valuation(q) >= i
            elif self._has_exact():
                return x == self.exact
            else:
                # Not stored represents Z_q^*: the elements of valuation 0.
                return x.valuation(q) == 0

        raise TypeError("``q`` should be a prime of `K`, but got {}".format(q))

    def _contains(self, x):
        """
        Check if the subset this idele represents contains the `K`-element ``x``

        INPUT:

        - ``x`` -- element of `K`, the number field to which this idele belongs

        EXAMPLES::

            sage: K.<a> = NumberField(x^7-4/7*x^3+1)
            sage: J = IdeleGroup(K)
            sage: u = J(a^5-a)
            sage: v = J(None, [completions(K, oo)[0][1](a^3), None, None, None], {K.prime_above(5): (a^3, 4)})
            sage: u._contains(a^5-a)
            True
            sage: u._contains(1)
            False
            sage: v._contains(a^3)
            True
            sage: v._contains(a)
            False
        """
        if self._has_exact() and x != self.exact:
            return False
        K = self.parent().number_field
        for phi in K.places():
            if not self._contains_at(x, phi):
                return False
        for q in self.finite:
            if not self._contains_at(x, q):
                return False
        return True

    def is_principal(self):
        """
        Return whether or not this idele is principal

        We check if there exists an "original" rational ``r`` whose image lies in
        the open subset op the idele group that this idele represents.

        ALGORITHM:

        If this idele has an exact value, that is our only candidate for ``r``.
        Else we create a generator by ... TODO

        EXAMPLES::

            sage: K.<a> = NumberField(x^2+3)
            sage: J = IdeleGroup(K)
            sage: u = J(a, None, {K.prime_above(2): (-a+1, 10)})
            sage: v = J(None, None, {K.prime_above(2): (-a+1, 4), K.prime_above(5): (-a-1, 0)})
            sage: w = J(None, None, {K.prime_above(2): (-a+1, 4), K.prime_above(5): (1/5, 0)})
            sage: u.is_principal()
            False
            sage: v.is_principal()
            True
            sage: w.is_principal()
            False
        """
        K = self.parent().number_field
        if self._has_exact():
            # The only possible ``r`` is ``self.exact``.
            return self._contains(self.exact)

        # we are not exact.
        r = K.ideal(1)
        for q, val in self.finite.items():
            x, i = val
            r *= q^(x.valuation(q))
        if not r.is_principal():
            # There does not exist an element of `K` with the right valuations.
            return False
        # The ideal ``r`` is principal, so we can take a generator.
        r = r.gens_reduced()[0]
        # This will be our candidate element of `K` that could lie in ``self``.
        # It is only uniquely determined up to a unit of `K`.
        U = K.unit_group()
        if not U.is_finite():
            raise NotImplementedError("K has infinite unit group, we can't deal with that yet")
        if len(U.gens()) > 1:
            raise NotImplementedError("K has non-cyclic unit group, we can't deal with that yet")
        u = U.gen().value()
        for e in range(U.order()):
            if self._contains(u^e * r):
                #print("DEBUG: found an original: {}".format(u^e * r))
                return True
        return False
        r"""
        old:
        if all([u_oo is None for u_oo in self.infinite]):
            # Only finite values stored.
            if K is QQ and len(self.finite) == 1:
                p = list(self.finite.keys())[0]
                x, i = self.finite[p]
                # The only candidate element is p^e, where:
                e = x.valuation(p)
                # Need to find some `z \in \ZZ_p` such that
                #     `p^e = x*(1+p^i*z)`.
                # Defining:
                y = x / p^e
                # Gives as only option:
                #     `z = (p^e-x)/(x*p^i) = (1-y)/(y*p^i)`
                # which lies in `\QQ`, but we should check that it lies in
                # `\ZZ_p \cap \QQ = \ZZ_{(p)}` (localization at `p`). 
                # That means that the valuation of `z` at `p` should be at least
                # zero. Since the valuation of `y` at `p` is zero, this boils
                # down to:
                return (1-y).valuation(p) >= i
        """

    def equivalent_modulo_principal_ideles(self, v):
        r"""
        Return whether or not `self \equiv v \mod K*`

        INPUT:

        - v -- an idele
        """
        return (self/v).is_principal()


    def increase_precision(self, primes, prec_increment=1):
        """
        Increase the precision of ``self`` at the primes given in ``primes``

        INPUT:

        - ``primes`` -- an iterable containing prime ideals of our number field
                        and/or rational prime numbers, or a prime ideal of our
                        number field
        - ``prec_increment`` -- integer (default = 1); the amount by which we
                                increase the precision at each prime in
                                ``primes``

        Let `p` be a prime ideal in ``primes``. Suppose ``self`` represents the
        open subset `x * U_p^i` at `p`. Then after calling this method, ``self``
        will represent `x * U_p^(i + prec_increment)` at `p`.

        If `p` in ``primes`` is a rational prime number, then the above is done
        for each prime `q` lying above `p`.

        .. NOTE::

            If ``self`` has an exact known value at a prime `p`, then nothing
            changes. If one sees exactness as having infinite precision, this
            just corresponds to ``oo + prec_increment == oo``.

        .. NOTE::

            Setting ``prec_increment`` to a negative value will decrease the
            precision of ``self``. If the precision drops below zero anywhere,
            we throw an exception.

        EXAMPLE::

            sage: K.<a> = NumberField(x^3-2)
            sage: J = IdeleGroup(K)
            sage: p2 = K.prime_above(2)
            sage: p5 = K.prime_above(5)
            sage: u = J(None, None, {p2: (a, 5), p5: (1/3, 1)})

        Let's increase the precision of ``p2`` and *both* prime ideals above 5
        by 3::

            sage: u.increase_precision([p2, 5], 3)
            sage: u
            Idele over Number Field in a with defining polynomial x^3 - 2 with values
                    (a,): a * U_q^8
                    (-a^2 - 1,): 1/3 * U_q^4
                    (a^2 - 2*a - 1,): 1 * U_q^3

        We can also decrease precision::

            sage: u.increase_precision(p2, -1)
            sage: u
            Idele over Number Field in a with defining polynomial x^3 - 2 with values
                    (a,): a * U_q^7
                    (-a^2 - 1,): 1/3 * U_q^4
                    (a^2 - 2*a - 1,): 1 * U_q^3

        The precision of exact values is unchanged::
        
            sage: u.exact = a+1
            sage: p3 = K.prime_above(3)
            sage: u.increase_precision([p2, p3, p5])
            sage: u
            Idele over Number Field in a with defining polynomial x^3 - 2 with values
                    (a,): a * U_q^8
                    (-a^2 - 1,): 1/3 * U_q^5
                    (a^2 - 2*a - 1,): 1 * U_q^3
            and which equals exactly a + 1 at all other primes
        """
        K = self.parent().number_field

        if primes in K.ideal_monoid() and K.ideal(primes).is_prime():
            p = primes # primes is just a single prime ideal, let's call it p
            if p in self.finite:
                x, i = self.finite[p]
            elif self._has_exact():
                # We know the value at p exactly, so we don't change
                return
            else:
                x, i = K(1), ZZ(0)
            new_prec = i + prec_increment
            if new_prec < 0:
                raise ValueError("Trying to give idele negative precision")
            self.finite[p] = (x, new_prec)
            return

        for p in primes:
            if p in K.ideal_monoid() and K.ideal(p).is_prime():
                self.increase_precision(p, prec_increment)
            elif p in Primes():
                for q in K.primes_above(p):
                    self.increase_precision(q, prec_increment)
            else:
                raise TypeError("primes should be a list of primes")



from sage.structure.unique_representation import UniqueRepresentation
from sage.groups.group import Group
class IdeleGroup(UniqueRepresentation, Group):
    Element = Idele

    def __init__(self, K):
        # TODO: implement default K=QQ, via __classcall__()
        if not is_field(K) or not K.absolute_degree() in ZZ:
            raise TypeError("K should be a number field")
        self.number_field = K
        Group.__init__(self)

    def _repr_(self):
        """
        Return a string representation of ``self``
        """
        return "Idele Group of {}".format(self.number_field)

    def _latex_(self):
        return r"\Bold{A}_{" + latex(self.number_field) + "}^*"

    def _element_constructor_(self, exact, infinite=None, finite=None):
        #print("DEBUG: IdeleGroup._element_constructor_({}, {}, {})".format(exact, infinite, finite))
        if infinite is None and finite is None:
            if exact is None:
                raise TypeError("No arguments supplied to Idele-constructor")
            if exact.parent() is Adeles(self.number_field): # conversion A_K --> J_K
                return self._from_adele(exact)
            if hasattr(exact.parent(), "_bnr"): # conversion Cl_m --> J_K
                # TODO make the check above less hacky and more robust
                # for some reason checking isinstance(exact, RayClassGroupElement) fails
                return self._from_ray_class(exact)
        return self.element_class(self, exact, infinite, finite)

    def _from_adele(self, adele):
        r"""
        Convert the adele ``adele`` to an idele

        Let U be the subset of the adele ring that ``adele`` represents. Write J
        for the idele group. Then the idele returned represents a subset that
        *contains* `U \cap J`. Note that in general this looses precision.
        TODO: example of precision loss

        EXAMPLES::
            
            sage: J = IdeleGroup(QQ)
            sage: A = Adeles(QQ)
            sage: Qhat = ProfiniteNumbers(QQ)
            sage: a = A([-1.5], Qhat(7, 24, 5))
            sage: J._from_adele(a)
            (-1.5000000000000000?, 7*(1+M_2^3), 7*(1+M_3^1), ...)
            sage: b = A([RIF(-1, 1)], Qhat(5, 25, 2))
            sage: J._from_adele(b)
            (0.?, Z_2*, Z_3*, 5*(1+M_5^1), ...)
            sage: c = A([pi.n()], 20/3)
            sage: J._from_adele(c)
            (3.1415926535897932?, 20/3, 20/3, 20/3, ...)

        If the given adele has value zero modulo on of its stored primes (i.e.
        ``4 mod 12`` has value zero modulo 2^2), then there is no idele that
        represents the adele. For example: the given adele then represents
        multiple element which have different valuation at that prime. But an
        idele always has a unique valuation at every prime. In this case, we
        throw an exception::

            sage: d = A([1], Qhat(4, 12, 3))
            sage: J._from_adele(d)
            Traceback (most recent call last):
            ...
            ValueError: adele is zero at the prime 2
        """
        K = self.number_field
        value = adele.finite.numerator.value
        modulus = adele.finite.numerator.modulus
        denominator = adele.finite.denominator

        exact = None
        if modulus.is_zero():
            exact = value / denominator
            finite = None
        else:
            x = value
            if K is not QQ:
                x = K.ideal(x)

            finite = {}
            for q, e in factor(modulus):
                v = x.valuation(q)
                if v >= e:
                    # value is zero modulo q^e, which we cannot represent as an
                    # idele
                    raise ValueError("adele is zero at the prime {}".format(q))
                finite[q] = (value, e-v)

        return self.element_class(self, exact, adele.infinite, finite)

    def _from_ray_class(self, r):
        r"""
        Convert the ray class group element ``r`` to an idele

        INPUT:

        - ``r`` -- a ray class group element

        OUPUT:
        
        Denote the ray class group to which ``r`` belongs by `G` and denote the
        modulus of `G` by `m`. So `G = I(m)/R(m)`.
        Consider the homomorphism `\phi: ` ``self`` `\to` `G` that sends a prime
        element at a finite prime `q` (not dividing `m`) to `q \mod R(m)`.
        The kernel of `\phi` is `K^* W_m` where `W_m = \prod_p U_p^{\ord_p(m)}`.

        Given ``r``, let `H` be the inverse image of ``r`` under `\phi`. We can
        try to find an idele that represents the subset `H` of ``self``. This is
        however not precisely possible. We can exactly represent `W_m`, but 
        we can not represent `K^*`. Hence what we do is the following: we find
        some `x \in H` en return an idele that represents `x \cdot W_m`.

        Although we do always return the same idele for equal inputs, the user
        should be aware that from a mathematical perspective, the output is only
        defined up to multiplication by a principal idele.
        

        EXAMPLES::

            sage: Q = NumberField(x-1, "one")
            sage: J = IdeleGroup(Q)
            sage: G = ray_class_group(Q, Modulus(Q.ideal(10), [0]))
            sage: r = G(Q.ideal(9))
            sage: factor(r.ideal())
            (Fractional ideal (3)) * (Fractional ideal (163))
            sage: J._from_ray_class(r)
            Idele over Number Field in one with defining polynomial x - 1 with values
                    infinity_0: [0.0000000000000000 .. +infinity]
                    (2,): 1 * U_q^1
                    (5,): 1 * U_q^1
                    (3,): 3 * U_q^0
                    (163,): 163 * U_q^0
            sage: s = G(Q.ideal(7))
            sage: s.ideal()
            Fractional ideal (67)
            sage: J._from_ray_class(s)
            Idele over Number Field in one with defining polynomial x - 1 with values
                    infinity_0: [0.0000000000000000 .. +infinity]
                    (2,): 1 * U_q^1
                    (5,): 1 * U_q^1
                    (67,): 67 * U_q^0
           
        :: 

            sage: K.<a> = NumberField(x^2-6)
            sage: G = ray_class_group(K, Modulus(K.ideal(10*a), [1]))
            sage: r = G([3, 0, 1])
            sage: factor(r.ideal())
            (Fractional ideal (25*a + 19)) * (Fractional ideal (28*a - 25)) * (Fractional ideal (-67*a + 109)) * (Fractional ideal (-1507*a - 5011))
            sage: Jk = IdeleGroup(K)
            sage: Jk(r)
            Idele over Number Field in a with defining polynomial x^2 - 6 with values
                    infinity_1: [0.0000000000000000 .. +infinity]
                    (-a + 2,): 1 * U_q^3
                    (a + 3,): 1 * U_q^1
                    (-a - 1,): 1 * U_q^1
                    (-a + 1,): 1 * U_q^1
                    (25*a + 19,): a + 543 * U_q^0
                    (28*a - 25,): a - 1312 * U_q^0
                    (-67*a + 109,): a - 1799 * U_q^0
                    (-1507*a - 5011,): a + 1242116 * U_q^0
        """
        K = self.number_field
        G = r.parent()  # ray class group of r
        exact = None
        infinite = [None for phi in self.number_field.places()]
        for i in G.modulus().infinite_part():
            infinite[i] = RIF(0, oo)

        finite = {}
        for q, e in G.modulus().finite_factors():
            finite[q] = (K(1), e)

        if not r.is_one():
            for q, e in factor(r.ideal()):
                if self.number_field is QQ:
                    finite[q] = (q^e, 0)
                else:
                    pi = K.uniformizer(q)
                    finite[q] = (pi^e, 0)
        
        return self.element_class(self, exact, infinite, finite)

    def cardinality(self):
        return Infinity

    def _coerce_map_from_(self, S):
        if self.number_field.has_coerce_map_from(S):
            return True
        return False


