# -*- coding: utf-8 -*-
r"""
Class Groups of Number Fields

Sage can compute class groups, ray class groups, and `S`-class groups of number
fields, and does so by wrapping the functionality from the PARI C-library. Some
of what can be computed includes the group structure, representative ideals, the
class of a given ideal, generators of the group, and products of elements. This
file also implements moduli of number fields.

AUTHORS:

- Robert Harron (2016-08-15): implemented ray class groups and moduli

EXAMPLES:

Computations with a ray class group of a quadratic field::

    sage: F = QuadraticField(40)
    sage: m = F.ideal(3).modulus([0, 1]); m
    (Fractional ideal (3)) * infinity_0 * infinity_1
    sage: R = F.ray_class_group(m); R
    Ray class group of order 8 with structure C4 x C2 of Number Field in a with defining polynomial x^2 - 40 of modulus (Fractional ideal (3)) * infinity_0 * infinity_1

Unlike for class groups and `S`-class groups, ray class group elements
do not carry around a representative ideal (for reasons of efficiency).
Nevertheless, one can be demanded. The returned ideal should be somewhat
'small.'

::

    sage: R.gens()
    (c0, c1)
    sage: R.gens_ideals()
    (Fractional ideal (430, 1/2*a + 200), Fractional ideal (-3/2*a + 2))
    sage: c = R.gens()[0]^3 * R.gens()[1]; c
    c0^3*c1
    sage: c.ideal()
    Fractional ideal (10, a)
    sage: c = R(F.ideal(2)); c
    c0^2
    sage: R.gen(0).ideal()^2
    Fractional ideal (30*a - 470)
    sage: R(R.gen(0).ideal()^2).ideal()
    Fractional ideal (2)

"""

# ****************************************************************************
#       Copyright (C) 2004, 2005, 2006, 2007 William Stein <wstein@gmail.com>
#                     2014 Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.groups.abelian_gps.abelian_group import AbelianGroup_class
from sage.groups.abelian_gps.abelian_group_element import AbelianGroupElement
from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method
from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.libs.pari import pari

def pari_real_places_to_sage(K):
    """
    Return a list converting from the ordering of real places in pari to that
    of Sage's :func:`places`.

    EXAMPLES:

    A totally real quartic field where the pari and Sage orderings are different.

    ::

        sage: x = polygen(QQ)
        sage: f = x^4 - x^3 - 3*x^2 + x + 1
        sage: F.<a> = NumberField(f(1-2*x))
        sage: F.defining_polynomial()
        16*x^4 - 24*x^3 + 8*x - 1
        sage: F.pari_polynomial()
        x^4 - 6*x^2 - 5*x - 1
        sage: pari_real_places_to_sage(F)
        (0, 3, 2, 1)

    A quintic field with three real places.

    ::

        sage: x = polygen(QQ)
        sage: f = x^5 - x^3 - 2 * x^2 + 1
        sage: F.<a> = NumberField(f(1 - x))
        sage: pari_real_places_to_sage(F)
        (2, 1, 0)
    """
    try:
        return K._pari_real_places
    except AttributeError:
        pass
    pari_conv = K._pari_absolute_structure()[1].lift()
    pari_conv = [pari_conv.polcoef(i).sage() for i in range(pari_conv.poldegree() + 1)]
    R = K.defining_polynomial().parent()
    pari_conv = R(pari_conv)
    pari_roots = [pari_conv(r.sage()) for r in K.pari_nf()[5][:K.signature()[0]]]
    pari_roots_sorted = list(pari_roots)
    pari_roots_sorted.sort()
    K._pari_real_places = tuple(pari_roots_sorted.index(r) for r in pari_roots)
    return K._pari_real_places

def pari_extended_ideal_to_sage(K, Ix):
    """
    Convert an 'extended ideal' in pari format to an ideal in Sage.

    INPUT:
    
    - ``K`` -- number field to which the ideal belongs
    - ``Ix`` -- A pair whose first entry is a pari ideal and whose second entry
      is a pari factorization matrix representing an algebraic number.

    OUTPUT:

    - The ideal that is the product of the ideal ``Ix[0]`` and the principal
      ideal generated by ``Ix[1]``.

    EXAMPLES::

        sage: F = NumberField(x^3 - 2, 'z')
        sage: I = F.ideal(27)
        sage: a = Matrix(0)
        sage: pari_extended_ideal_to_sage(F, pari([I, a])) == I
        True
        sage: factorization_matrix = pari([pari([pari([31, 0, 1]).Col(), 1]).Mat(), pari([pari([5, 1, 0]).Col(), 1]).Mat()]).Col().matconcat()
        sage: pari_extended_ideal_to_sage(F, pari([I, factorization_matrix]))
        Fractional ideal (-135*z^2 - 837*z - 4239)
        sage: pari_extended_ideal_to_sage(F, F.ray_class_group(I).gens_values()[0])
        Fractional ideal (-8*z^2 - 9*z + 1)
    """
    I = K.ideal(Ix[0])
    nf = K.pari_nf()
    a = K(nf.nfbasistoalg(nf.nffactorback(Ix[1])))
    return a * I

def ray_class_group(K, modulus, proof=True, names='c'):
    r"""
    Return the ray class group modulo ``modulus``.

    EXAMPLES:

    Ray class group of a real quadratic field of class number 1.

    ::

        sage: F = NumberField(x^2 - 5, 'a')
        sage: m = F.modulus(F.prime_above(5) * F.prime_above(29), [0, 1])
        sage: G = F.ray_class_group(m); G
        Ray class group of order 8 with structure C4 x C2 of Number Field in a with defining polynomial x^2 - 5 of modulus (Fractional ideal (-11/2*a - 5/2)) * infinity_0 * infinity_1
        sage: G.elementary_divisors()
        (2, 4)
        sage: G.gens_ideals()
        (Fractional ideal (31), Fractional ideal (12672))
        sage: G.gens_orders()
        (4, 2)

    Ray class group of an imaginary quadratic field of class number 3.

    ::

        sage: F.<a> = QuadraticField(-23)
        sage: R = F.ray_class_group(F.ideal(3/2 + a/2)); R
        Ray class group of order 6 with structure C6 of Number Field in a with defining polynomial x^2 + 23 of modulus Fractional ideal (1/2*a + 3/2)
        sage: R.gens_ideals()
        (Fractional ideal (3, 1/2*a + 1/2),)
        sage: R.modulus().finite_part().norm()
        8
        sage: F.class_group().gens_ideals()
        (Fractional ideal (2, 1/2*a - 1/2),)

    Over `\QQ`, the ray class group modulo `(m)\infty` is the unit group of `\ZZ/m\ZZ`, while
    the ray class group modulo `(m)` is the latter modulo `\{\pm1\}`.

    ::

        sage: Q = NumberField(x - 1, 'a')
        sage: Q.ray_class_group(Q.ideal(40).modulus([0])).gens_ideals()
        (Fractional ideal (17), Fractional ideal (31), Fractional ideal (21))
        sage: Zmod(40).unit_gens()
        (31, 21, 17)
        sage: Q.ray_class_group(Q.ideal(40)).gens_ideals()
        (Fractional ideal (17), Fractional ideal (21))
    """
    try:
        if modulus.parent() is K.ideal_monoid():
            modulus = modulus.modulus()
    except AttributeError:
        pass
    try:
        return K.__ray_class_group[modulus, proof, names]
    except KeyError:
        pass
    except AttributeError:
        K.__ray_class_group = {}
    Kbnf = K.pari_bnf()
    Kbnr = Kbnf.bnrinit(pari(modulus), 1)
    cycle_structure = tuple(ZZ(c) for c in Kbnr[4][1])
    #@@gens = tuple(K.ideal(hnf) for hnf in Kbnr[4][2])
    from sage.matrix.constructor import Matrix
    gens = tuple([Kbnf.idealred([hnf, pari(Matrix(0))]) for hnf in Kbnr[4][2]])#@@
    G = RayClassGroup(cycle_structure, names, modulus, gens, proof=proof, bnr=Kbnr)
    K.__ray_class_group[modulus, proof, names] = G
    return G

def _integer_n_tuple_L1_iterator(n):
    if n == 1:
        i = 1
        while True:
            yield [i]
            yield [-i]
            i += 1
    else:
        from sage.combinat.partition import OrderedPartitions
        from sage.combinat.subset import Subsets
        N = 1
        sign_options_dict = {}
        Subsets_of_n = []
        while True:
            #print "***"
            #print N
            for k in range(1, n + 1):
                Ps = OrderedPartitions(N, k)
                for P in Ps:
                    #print "--"
                    #print P
                    try:
                        Ss = Subsets_of_n[k - 1]
                    except IndexError:
                        Ss = Subsets(range(n), k)
                        Subsets_of_n.append(Ss)
                    for S in Ss:
                        #print "="
                        #print S
                        i = [0] * n
                        for j in range(k):
                            i[S[j]] = P[j]
                        yield i
                        try:
                            sign_options = sign_options_dict[S]
                        except KeyError:
                            sign_options = Subsets(S)[1:]
                            sign_options_dict[S] = sign_options
                        for signs in sign_options:
                            ii = list(i)
                            for index in signs:
                                ii[index] = - ii[index]
                            yield ii
            N += 1


class Modulus(SageObject):
    def __init__(self, finite, infinite=None, check=True):
        r"""
        Create a modulus of a number field.

        INPUT:

        - ``finite`` -- a non-zero fractional ideal in a number field.
        - ``infinite`` -- a list of indices corresponding to real places
          of the number field, sorted.
        - ``check`` (default: True) -- If ``True``, run a few checks on the input.
        """
        self._finite = finite
        if infinite is None:
            self._infinite = ()
        else:
            self._infinite = tuple(ZZ(i) for i in infinite)
        K = self._finite.number_field()
        self._number_field = K
        if check:
            #insert various checks here
            if self._finite == 0:
                raise ValueError("Finite component of a modulus must be non-zero.")
            sgn = K.signature()[0]
            for i in self._infinite:
                if i < 0 or i >= sgn:
                    raise ValueError("Infinite component of a modulus must be a list non-negative integers less than the number of real places of K")
        return

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: K.<a> = NumberField(x^2-5)
            sage: m = K.modulus(K.ideal(31), [0,1]); m
            (Fractional ideal (31)) * infinity_0 * infinity_1
        """
        if len(self._infinite) == 0:
            return str(self._finite)
        str_inf = ''
        for i in self._infinite:
            str_inf += ' * infinity_%s'%(i)
        return '(' + str(self._finite) + ')' + str_inf

    def __eq__(self, other):
        return self._number_field == other._number_field and self._finite == other._finite and self._infinite == other._infinite

    def __mul__(self, other):
        r"""
        Multiply two moduli.

        This multiplies the two finite parts and performs an exclusive or on the real places.

        EXAMPLES::

            sage: K = NumberField(x^3 - 2, 'a')
            sage: m1 = K.modulus(K.ideal(2), [0])
            sage: m2 = K.modulus(K.ideal(3), [0])
            sage: m1 * m2
            Fractional ideal (6)

        A higher degree totally real field::

            sage: K = NumberField(x^5 - x^4 - 4*x^3 + 3*x^2 + 3*x - 1, 'a')
            sage: m1 = K.modulus(K.ideal(5), [2, 3])
            sage: m2 = K.modulus(K.ideal(25), [0, 1, 3, 4])
            sage: m1 * m2
            (Fractional ideal (125)) * infinity_0 * infinity_1 * infinity_2 * infinity_4
            sage: _ == m2 * m1
            True
        """
        inf = tuple(set(self.infinite_part()).symmetric_difference(other.infinite_part()))
        return Modulus(self.finite_part() * other.finite_part(), inf, check=False)

    def lcm(self, other):
        inf = tuple(set(self.infinite_part()).union(other.infinite_part()))
        #Pe_out = []
        self_fact_P, self_fact_e = zip(*self.finite_part().factor())
        self_fact_P = list(self_fact_P)
        self_fact_e = list(self_fact_e)
        #self_facts = self.finite_part().factor()
        #of = other.finite_part()
        other_facts = other.finite_part().factor()
        mf = self._number_field.ideal_monoid().one()
        for P, e in other_facts:
            try:
                i = self_fact_P.index(P)
            except ValueError:
                #Pe_out.append([P, e])
                mf *= P**e
                continue
            #Pe_out.append([P, max(e, self_fact_e[i])])
            mf *= P**max(e, self_fact_e[i])
            del self_fact_P[i]
            del self_fact_e[i]
        for i in range(len(self_fact_P)):
            mf *= self_fact_P[i]**self_fact_e[i]
        return Modulus(mf, inf, check=False)

    def divides(self, other):
        if not set(self.infinite_part()).issubset(other.infinite_part()):
            return False
        return self.finite_part().divides(other.finite_part())

    def number_field(self):
        return self._number_field

    def finite_part(self):
        return self._finite

    def infinite_part(self):
        return self._infinite

    def finite_factors(self):
        try:
            return self._finite_factors
        except AttributeError:
            self._finite_factors = self.finite_part().factor()
            return self._finite_factors

    def equivalent_coprime_ideal_multiplier(self, I, other):
        r"""
        Given ``I`` coprime to this modulus `m`, return a number field element `\beta`
        such that `\beta I` is coprime to the modulus ``other`` and equivalent to
        ``I`` `\mathrm{mod}^\ast m`; in particular, `\beta` will be `1 \mathrm{mod}^\ast m`.

        EXAMPLES:

        An example with two prime factors difference between this modulus and ``other``.

        ::

            sage: F.<a> = QuadraticField(5)
            sage: m_small = F.modulus(3/2*a - 1/2, [0, 1])
            sage: m_big = F.modulus(2*a - 30, [0, 1])
            sage: m_small.equivalent_coprime_ideal_multiplier(F.ideal(6), m_big)
            109/54
        """
        F = self._number_field
        other_Ps = [P for P, _ in other.finite_factors() if self._finite.valuation(P) == 0]
        if len(other_Ps) == 0: #If prime factors of other is a subset of prime factors of this modulus
            return F.one()
        alpha = I.idealcoprime(other._finite)
        if self.number_is_one_mod_star(alpha):
            return alpha
        nf = F.pari_nf()
        other_Ps_facts = [nf.idealfactor(P) for P in other_Ps]
        other_Ps_fact_mat = pari(other_Ps_facts).Col().matconcat()
        self_fact_mat = nf.idealfactor(self.finite_part())
        x = pari([self_fact_mat, other_Ps_fact_mat]).Col().matconcat()
        conditions = [~alpha] * self_fact_mat.nrows() + [1] * len(other_Ps)
        gamma = F(nf.idealchinese(x, conditions))
        beta = alpha * gamma
        if self.number_is_one_mod_star(beta):
            return beta
        from sage.misc.misc_c import prod
        beta_fixed = F.modulus(self._finite * prod(other_Ps), self._infinite).fix_signs(beta) #should be able to do this more efficiently
        return beta_fixed

    def equivalent_ideal_coprime_to_other(self, I, other):
        r"""
        Given ``I`` coprime to this modulus `m`, return an ideal `J` such that `J` is coprime
        to the modulus ``other`` and equivalent to ``I`` `\mathrm{mod}^\ast m`.

        This is useful for lowering the level of a non-primitive Hecke character.

        INPUT:

        - ``I`` -- an ideal relatively prime to this modulus (not checked).
        - ``other`` -- some other modulus.

        OUTPUT:

        an ideal coprime to ``other`` and equivalent to ``I`` in the ray class
        group modulo this modulus.
        """
        return self.equivalent_coprime_ideal_multiplier(I, other) * I

    def number_is_one_mod_star(self, a):
        K = self.number_field()
        am1 = K(a - 1)
        for P, e in self.finite_factors():
            if am1.valuation(P) < e:
                return False
        inf_places = K.places()
        for i in self.infinite_part():
            if inf_places[i](a) <= 0:
                return False
        return True

    def fix_signs(self, a):
        r"""
        Given ``a`` in ``self.number_field()``, find `b` congruent to ``a`` `mod^\ast` ``self.finite_part()``
        such that `b` is positive at the infinite places dividing ``self``.
        """
        if self.is_finite() or a == 0:
            return a
        places = self.number_field().places()
        positive = []
        negative = []
        for i in self.infinite_part():
            if places[i](a) > 0:
                positive.append(i)
            else:
                negative.append(i)
        if not negative:
            return a
        t = self.get_one_mod_star_finite_with_fixed_signs(positive, negative)
        return t * a

    def get_one_mod_star_finite_with_fixed_signs(self, positive, negative):
        if len(negative) == 0:
            return self.number_field().one()
        negative = tuple(negative)
        try:
            return self._one_mod_star[negative]
        except AttributeError:
            self._one_mod_star = {}
        except KeyError:
            pass
        try:
            beta_is, Ainv = self._beta_is_Ainv
        except AttributeError:
            beta_is, Ainv = self._find_beta_is_Ainv()
        d = len(self.infinite_part())
        v = [0] * d
        for i in range(d):
            if self._infinite[i] in negative:
                v[i] = 1
        v = (GF(2)**d)(v)
        w = Ainv * v
        t = self.number_field().one()
        for i in range(d):
            if w[i] != 0:
                t *= (1 + beta_is[i])
        self._one_mod_star[negative] = t
        return t

    def _find_beta_is_Ainv(self):
        r"""
        Step 2 of Algorithm 4.2.20 of Cohen's Advanced...
        """
        from sage.matrix.special import column_matrix
        gammas = self.finite_part().basis()
        k = len(self.infinite_part())
        beta_is = []
        Acols = []
        V = GF(2)**k
        it = _integer_n_tuple_L1_iterator(k)
        while len(beta_is) < k:
            e = next(it)
            beta = sum(e[i] * gammas[i] for i in range(k))
            sbeta = V(self._signs(beta))
            Acols_new = Acols + [sbeta]
            A = column_matrix(GF(2), Acols_new)
            if A.rank() == len(Acols_new):
                Acols = Acols_new
                beta_is.append(beta)
        self._beta_is_Ainv = (beta_is, ~A)
        return self._beta_is_Ainv

    def _signs(self, b):
        if b == 0:
            raise ValueError("Non-zero input required.")
        sigmas = self._number_field.real_places()
        return [ZZ.one() if sigmas[i](b).sign() == -1 else ZZ.zero() for i in self.infinite_part()]

    def is_finite(self):
        return len(self._infinite) == 0

    def is_infinite(self):
        return self._finite.is_one()

    def __pari__(self):
        """
        Return the corresponding pari modulus.

        Note that this function performs the conversion between the ordering
        of the real places of the number field in Sage and the ordering of the
        underlying pari nf object.

        EXAMPLES:

        An example where the places in Sage and pari are in a different order.

        ::

            sage: x = polygen(QQ)
            sage: f = x^4 - x^3 - 3*x^2 + x + 1
            sage: F.<a> = NumberField(f(1-2*x))
            sage: F.modulus(1, [2, 3]).__pari__()[1]
            [0, 1, 1, 0]
        """
        inf_mod = [0] * self._number_field.signature()[0]
        conversion = pari_real_places_to_sage(self._number_field)
        for i in self._infinite:
            pari_index = conversion.index(i)
            inf_mod[pari_index] = 1
        return pari([self._finite, inf_mod])

    def __hash__(self):
        return hash((self._finite, self._infinite))

class RayClassGroupElement(AbelianGroupElement):

    def ideal(self, reduce=True):
        """
        Return an ideal representing this ray class.

        If ``reduce`` is ``True`` (by default) the returned ideal is
        reduced to 'small' (this can be slow on large inputs).

        INPUT:

        - ``reduce`` -- (default: ``True``) determine whether or not
          to output a 'small' representative.

        OUTPUT:

        An ideal representing this ray class. If ``reduce`` is ``True``,
        the ideal returned is made 'small' by the ideal's
        ``reduce_equiv`` function (and ``reduce_equiv`` is used at
        each step of computing this representative). Otherwise, the
        output is just the appropriate product of the powers of the
        generators of the ray class group.

        EXAMPLES:

        Over a real quadratic field field of class number 1::

            sage: F.<a> = NumberField(x^2 - 5)
            sage: m = F.ideal(11).modulus([0, 1])
            sage: R = F.ray_class_group(m)
            sage: c0, c1 = R.gens()
            sage: c = c0^4*c1; c.ideal()
            Fractional ideal (-a)
            sage: c.ideal(False)
            Fractional ideal (-6242265*a + 1268055)

        Over a real quadratic field of class number 2::

            sage: F = QuadraticField(40)
            sage: R = F.ray_class_group(F.prime_above(13).modulus([0, 1]))
            sage: for c in R:
            ...       if R(c.ideal()) != c:
            ...           print "Bug!"
        """
        R = self.parent()
        exps = self.exponents()
        gens = R.gens_ideals()
        L = len(exps)
        #Speed this up later using binary powering
        i = 0
        while exps[i] == 0:
            i += 1
            if i == L:
                return R.one()
        I = (gens[i]**exps[i])
        if reduce:
            I = R.ideal_reduce(I)
        i += 1
        while i < L:
            e = exps[i]
            g = gens[i]
            i += 1
            if e != 0:
                I = (I * (g**e))
                if reduce:
                    I = R.ideal_reduce(I)
        return I

class RayClassGroup(AbelianGroup_class):
    Element = RayClassGroupElement

    def __init__(self, gens_orders, names, modulus, gens, bnr, proof=True):
        r"""
        ``gens`` -- a tuple of pari extended ideals
        """
        #AbelianGroupWithValues_class.__init__(self, gens_orders, names, gens,
        #                                      values_group=modulus.number_field().ideal_monoid())
        AbelianGroup_class.__init__(self, gens_orders, names)
        self._gens = gens
        self._proof_flag = proof #TODO: use this, in _ideal_log?
        self._modulus = modulus
        self._number_field = modulus.number_field()
        self._bnr = bnr
        self._is_narrow = (len(modulus.infinite_part()) == self._number_field.signature()[0]) and modulus.finite_part().is_one()

    def _element_constructor_(self, *args, **kwds):
        try:
            L = len(args[0])
        except TypeError:
            L = -1
        if L == self.ngens():
            return self.element_class(self, args[0])
        from idele import Idele
        if isinstance(args[0], RayClassGroupElement):
            c = args[0]
            if c.parent() is self:
                return self.element_class(self, c.exponents())
            nf = self._number_field.pari_nf()
            bnr = self._bnr
            old_exps = c.exponents()
            old_gens = c.parent().gens_values()
            L = len(old_exps)
            if L ==0:
                return self.one()
            i = 0
            while old_exps[i] == 0:
                i += 1
                if i == L:
                    return self.one()
            I = nf.idealpow(old_gens[i], old_exps[i], flag=1)
            i += 1
            while i < L:
                e = old_exps[i]
                g = old_gens[i]
                i += 1
                if e != 0:
                    I = nf.idealmul(I, nf.idealpow(g, e, flag=1), flag=1)
            exps = tuple(ZZ(c) for c in bnr.bnrisprincipal(nf.idealmul(I[0], nf.nffactorback(I[1])), flag = 0))
            return self.element_class(self, exps)
        elif isinstance(args[0], Idele):
            return args[0].to_ray_class(self.modulus())
        else:
            I = self._number_field.ideal(*args, **kwds)
            if I.is_zero():
                raise TypeError("The zero ideal is not a fractional ideal")
            exps = self._ideal_log(I)
            return self.element_class(self, exps)

    def _repr_(self):
        if self._is_narrow:
            s0 = 'Narrow '
        else:
            s0 = 'Ray '
        s = 'class group of order %s '%self.order()
        if self.order() > 1:
            s += 'with structure %s '%self._group_notation(self.gens_orders())
        s += 'of %s'%(self._number_field)
        if self._is_narrow:
            return s0 + s
        s += ' of modulus %s'%(self._modulus)
        return s0 + s

    def ray_class_field(self, subgroup=None, names=None, algorithm='stark'):
        r"""
        Two different algorithms are possible: pari's :pari:`bnrstark` and
        :pari:`rnfkummer`. The first one uses the Stark conjecture and only
        deals with totally real extensions of a totally real base
        field. The second one uses Kummer theory and only deals with
        extensions of prime degree.

        INPUT:

        - algorithm -- (default: ``stark``) if the value is ``stark``,
          then pari's :pari:`bnrstark` function is tried first, and if that
          fails, :pari:`rnfkummer` will be attempted. If the value is
          ``kummer``, then pari's :pari:`rnfkummer` is tried first, with
          :pari:`bnrstark` as a backup. Using ``stark_only`` or ``kummer_only``
          will just raise an exception if the first attempt fails.

        OUTPUT:

        The class field corresponding to the given subgroup, or the
        ray class field if ``subgroup`` is ``None``, as a relative
        number field.

        EXAMPLES:

        Class fields of `\QQ(\sqrt{3})`::

            sage: F.<a> = QuadraticField(3)
            sage: m = F.ideal(7).modulus()
            sage: R = F.ray_class_group(m)
            sage: R.ray_class_field(names='b')
            Number Field in b with defining polynomial x^6 + a*x^5 - 4*x^4 - 4*a*x^3 + 2*x^2 + 2*a*x - 1 over its base field
            sage: S = R.subgroup([R.gen()^2])
            sage: R.ray_class_field(S, names='b')
            Number Field in b with defining polynomial x^2 - a*x - 1 over its base field
            sage: m = F.modulus(20)
            sage: R = F.ray_class_group(m)
            sage: S = R.subgroup([R.gens()[0]^2, R.gens()[1]])
            sage: R.ray_class_field(S, names='b')
            Number Field in b with defining polynomial x^2 + (a - 1)*x + 2*a - 4 over its base field

        An example where :pari:`bnrstark` fails, but :pari:`rnfkummer` saves the day::

            sage: F.<a> = NumberField(x^8 - 12*x^6 + 36*x^4 - 36*x^2 + 9)
            sage: m = F.ideal(2).modulus()
            sage: R = F.ray_class_group(m)
            sage: set_verbose(1)
            sage: K = R.ray_class_field(names='b'); K
            verbose 1 (...: class_group.py, ray_class_field) bnrstark failed; trying rnfkummer.
            Number Field in b with defining polynomial x^2 + (1/3*a^6 - 10/3*a^4 + 5*a^2)*x + 1/3*a^6 - 1/3*a^5 - 11/3*a^4 + 3*a^3 + 8*a^2 - 4*a - 5 over its base field
            sage: set_verbose(0)
        """
        if subgroup is not None:
            try:
                test_subgrp = (subgroup.ambient_group() is self)
            except AttributeError:
                subgroup = self.subgroup(subgroup)
                test_subgrp = (subgroup.ambient_group() is self)
            if not test_subgrp:
                raise ValueError("subgroup does not define a subgroup of this ray class group.")
            gens_coords = [h.exponents() for h in subgroup.gens()]
            if len(gens_coords) == 0:
                subgroup = None
                subgroup_mat = None
            else:
                from sage.matrix.special import column_matrix
                subgroup_mat = column_matrix(gens_coords)
        else:
            subgroup_mat = None

        from cypari2.handle_error import PariError
        from sage.misc.all import verbose

        bnr = self._bnr
        if algorithm == 'stark_only':
            if len(self._modulus.infinite_part()) > 0 or not self._number_field.is_totally_real():
                raise NotImplementedError("Stark's conjecture algorithm only implemented for totally real extensions of a totally real base field.")
            f = bnr.bnrstark(subgroup=subgroup_mat)
        elif algorithm == 'kummer_only':
            if (subgroup is None and not self.order().is_prime()) or (subgroup is not None and not self.order().divide_knowing_divisible_by(subgroup.order()).is_prime()):
                raise NotImplementedError("Kummer theory algorithm only implemented extensions of prime degree.")
            f = bnr.rnfkummer(subgp=subgroup_mat)
        elif algorithm == 'stark':
            if len(self._modulus.infinite_part()) > 0 or not self._number_field.is_totally_real():
                if (subgroup is None and not self.order().is_prime()) or (subgroup is not None and not self.order().divide_knowing_divisible_by(subgroup.order()).is_prime()):
                    raise NotImplementedError("Ray class fields only implemented for totally real extensions of totally real base fields, or for extensions of prime degree.")
                f = bnr.rnfkummer(subgp=subgroup_mat)
            else:
                try:
                    f = bnr.bnrstark(subgroup=subgroup_mat)
                except PariError:
                    if (subgroup is None and self.order().is_prime()) or (subgroup is not None and self.order().divide_knowing_divisible_by(subgroup.order()).is_prime()):
                        verbose("bnrstark failed; trying rnfkummer.")
                        f = bnr.rnfkummer(subgp=subgroup_mat)
                    else:
                        raise
        elif algorithm == 'kummer':
            if (subgroup is None and self.order().is_prime()) or (subgroup is not None and self.order().divide_knowing_divisible_by(subgroup.order()).is_prime()):
                f = bnr.rnfkummer(subgp=subgroup_mat)
            else:
                f = bnr.bnrstark(subgroup=subgroup_mat)
        else:
            raise ValueError("Value of algorithm must be one of \'stark\', \'stark_only\', \'kummer\', or \'kummer_only\'.")
        if f.type() == 't_VEC':
            raise NotImplementedError("bnrstark returned a list of polynomials. Dealing with this has not been implemented.")
        F = self._number_field
        nf = F.pari_nf()
        f = nf.rnfpolredbest(f)
        d = f.poldegree()
        cs = [F(f.polcoef(i)) for i in range(d+1)]
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        f = PolynomialRing(F, 'x')(cs)
        return F.extension(f, names=names)

    def _ideal_log(self, ideal):
        return tuple(ZZ(c) for c in self._bnr.bnrisprincipal(ideal, flag = 0))

    def ideal_reduce(self, ideal):
        from cypari2.handle_error import PariError
        ideal = pari(ideal)
        try:
            pari_ideal = self._bnr.idealmoddivisor(ideal)
        except PariError as err:
            if err.errtext().find('not coprime') != -1:
                raise ValueError('Ideal in ideal_reduce is not coprime to the modulus of this ray class group.')
            else:
                raise err
        return self._number_field.ideal(pari_ideal)

    def gens_values(self):
        return self._gens

    @cached_method
    def gens_ideals(self):
        return tuple(pari_extended_ideal_to_sage(self._number_field, v) for v in self.gens_values())

    def modulus(self):
        return self._modulus

    def number_field(self):
        return self._number_field

    def pari_bnr(self):
        return self._bnr

    def pari_gens(self):
        return self._bnr[4][2]

