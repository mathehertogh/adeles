r"""

References
--------------------

- [GS1998] -- Gee A., Stevenhagen P. (1998) Generating class fields using
  Shimura reciprocity. In: Buhler J.P. (eds) Algorithmic Number
  Theory. ANTS 1998. Lecture Notes in Computer Science, vol 1423.
  Springer, Berlin, Heidelberg. https://doi.org/10.1007/BFb0054883




-------------------------------------------------------------------
Computing Hilbert class fields using idelic Shimura reciprocity law
-------------------------------------------------------------------

Below we will demonstrate the computation of Hilbert class fields of quadratic
imaginary number fields using Shimura's reciprocity law.
We will follow the same examples given by Alice Gee and Peter Stevenhagen in
their article [GS1998].
Contrary to them, we will do the calculation using our implementation of ideles.
In order to extract an explicit action, we use ``GL2Qhat_factor()`` to factor
matrices in `GL_2(\hat{\QQ})` into a product of a matrix in `GL_2(\hat{\ZZ})`
and a matrix in `GL_2^+(\QQ)`. See [TODO] for details on this factorization.


Example 1 - `K = \QQ(\sqrt{-71})` and `f = \gamma_2`
----------------------------------------------------

The first example will compute the Hilbert class field of `K = \QQ(\sqrt{-71})`,
using the modular function `\gamma_2` of level 3 (TODO: does it have a name or
better description...? "Weber's `\gamma_2`" maybe?).
The ring of integers `\mathcal{O}` of `K` is generated by `\theta =
\frac{-1+\sqrt{-71}}{2}`, which has minimal polynomial `X^2+X+18`.

We start with loading our code and initilizing our quadratic imaginary number
field together with its idele group. ::

    sage: from shimura import *
    sage: from modular import *
    sage: R = ZZ['x']; x = R.gen()
    sage: K.<theta> = NumberField(x**2+x+18)
    sage: K.discriminant()  # Check that we have the correct field
    -71
    sage: O = K.maximal_order()
    sage: O.basis()  # Check that theta indeed generates the ring of integers
    [1, theta]
    sage: J = IdeleGroup(K)
    sage: level = 3  # the level of our modular function gamma_2

Finding the class invariant
------------------------------------------------------

Now we want to compute our class invariant. This means finding an element
`\alpha \in H_3` that is invariant under `Gal(H_3/H)`, where `H` is the
Hilbert class field of `K` and `H_3` its ray class field modulo 3.

The Artin map `(\mathcal{O}/3\mathcal{O})^* \to Gal(H_3/H)` is surjective with
kernel `\mathcal{O}^* = \{\pm 1\}`. Hence it suffices for `\alpha` to be
be invariant under `(\mathcal{O}/3\mathcal{O})^* / \{\pm 1\}`.
We compute::

    sage: Omod3 = O.quotient_ring(3, 'b')
    sage: Omod3star = K.ideal(3).idealstar(flag=2)  # flag=2 means compute generators
    sage: Omod3star.gens_values()
    (-theta + 1, theta - 1)

We see that `(\mathcal{O}/3\mathcal{O})^*` is generated by `-\theta+1` and
`\theta-1`, hence `(\mathcal{O}/3\mathcal{O})^* / \{\pm 1\}` is generated by
`\theta-1`. We compute the action of `\theta-1 \mod 3` on `\gamma_2` as follows.

First we create an idele out of `\theta-1 \mod 3`::

    sage: x = J(Omod3(theta-1)); x  # TODO: fix the ugly infinite place-printing below to CC^*
    ([-infinity .. +infinity] + [-infinity .. +infinity]*I, Z_p2*, Z_q2*, (theta - 1)*(1+M_p3^1), (theta - 1)*(1+M_q3^1), ...)
    where:
            p2 = Fractional ideal (2, theta)
            q2 = Fractional ideal (2, theta + 1)
            p3 = Fractional ideal (3, theta)
            q3 = Fractional ideal (3, theta + 1)

Next we compute the image of `x` under Shimura's connecting homomorphism, which
gives us a matrix in `GL_2(\widehat{\ZZ})`::

    sage: A = connecting_homomorphism_integral(x, level); A
    [  1 mod 3 36 mod 54]
    [  1 mod 3   2 mod 3]

.. NOTE::

    Actually, the function ``connecting_homomorphism_integral()`` above does a
    little bit more then just implementing the connecting homomorphism. It
    increases the precision of the idele ``x`` until the image matrix is
    integral and has modulus divisible by ``level``. So the returned matrix `A`
    is actually the image of a certain choice of idele represented by ``x``.
    See :func:`connecting_homomorphism_integral` for details.

As we now have `A \in GL_2(\hat{\ZZ})`, the action of `A` on `\gamma_2` only
depends on `A \mod 3`. So we project `A` to `GL_2(\ZZ/3\ZZ)`.

::

    sage: Amod3 = matrix_modulo(A, level); Amod3
    [1 0]
    [1 2]
    sage: det(Amod3)
    2

Now we factor ``Amod3`` as `\iota(d) \cdot U` with `U \in SL_2(\ZZ/3\ZZ)`,
`d = det(A)` and `\iota(d) = (1, 0; 0, d)`::

    sage: d = det(Amod3); d
    2
    sage: U = ~iota(d) * Amod3
    sage: iota(d), U
    (
    [1 0]  [1 0]
    [0 2], [2 1]
    )
    sage: iota(d)*U == A
    True

This tells us that `x` acts on `\zeta_3` as `\zeta_3 \mapsto \zeta_3^d =
\zeta_3^2`.
We determine the action of `x` on `\gamma_2`, which equals the action of `U` on
`\gamma_2`, using the following fact.
`SL_2(\ZZ)` (and hence `SL_2(\ZZ/3\ZZ)`) is generated by the "standard
generators" `S = (0, -1; 1, 0)` and `T = (1, 1; 0, 1)` and their action on
`\gamma_2` is given by

    - `\gamma_2^S = \gamma_2`;
    - `\gamma_2^T = \zeta_3^{-1}  \gamma_2`.

Hence we factor `U` as a product of `S`'s and `T`'s and find the explicit
action. This is done in the function below::

    sage: STs = ST_factor(U); STs
    (S^3*T^-1*S)^2
    sage: print_action_on_gamma_2(STs)
      gamma_2 ]--> zeta_3^2 * gamma_2

We see that `\zeta_3^x/\zeta_3 = \zeta_3` and `\gamma_2^x/\gamma_2 = \zeta_3^2`.
Hence `\zeta_3 \gamma_2` is left invariant by `x` and so we define `\alpha = 
\zeta_3 \gamma_2(\theta)`.
Because now by Shimura's reciprocity law, we have
`\alpha^x = (\zeta_3 \gamma_2)^x(\theta) = \zeta_3 \gamma_2(\theta) = \alpha`.

Computing the Hilbert class field
------------------------------------------------------------

We now know that `\alpha \in H` and in all likelihood, we will have `K(\alpha)
= H`. So we want to compute the minimal polynomial of `\alpha` over `K`.
We will do this by computing the conjugates of `\alpha` over `K`.

Recall that the Hilbert class field is the ray class field modulo 1. Write
`Cl` for the ray class group modulo 1 (which is isomorphic to the class group
of `K`). The Artin map `Cl \to Gal(H/K)` is an isomorphism and so the
conjugates of `\alpha` are given by `\alpha^x` for `x \in Cl`.

For each `x \in Cl` we will compute an approximation of `\alpha^x` in `\CC` as
follows. We replace `x` by its image in the idele group `J` of `K`.
Then we apply Shimura's connecting homomorphism to obtain a matrix `A \in
GL_2(\hat{\QQ})`. We use our ``GL2Qhat_factor()`` function to obtain matrices
`B \in GL_2(\hat{\ZZ})` and `M \in GL_2^+(\QQ)` such that `A = B \cdot M`.

The action of `B` on `\gamma_2` is determined as we did above. The action of `M`
is given by factional linear transformations. This gives us the exact value of
`\alpha^x`, which we evaluate numerically to obtain an appoximation in `\CC`.

From these approximations of `\alpha^x`, we can determine the minimal polynomial
`f_K^\alpha = \prod_{x \in Cl} (X - \alpha^x)` of `\alpha`.

We set up our embedding `\phi: K \to \CC` and compute the group `Cl`::

    sage: prec = 100  # We do our complex calculations with 100 bits of precision
    sage: CF = ComplexField(prec=prec)
    sage: zeta3 = CF(exp(2*CF(pi)*CF(I)/3))
    sage: phi = K.embeddings(CF)[1]  # the embedding where theta has positive imaginary part
    sage: Cl = ray_class_group(K, Modulus(K.ideal(1)))

Now we compute the conjugates of `\alpha`::

    sage: conjugates = []
    sage: for a in Cl:
    ....:     y = J(a)
    ....:     d, T = action_idele(y, level)  # zeta_3 * gamma_2^x ]--> zeta_3^d * gamma_2 \circ T
    ....:     tau = apply_fractional_linear_transformation(T, theta)
    ....:     conjugate_alpha = zeta3**d * gamma_2(phi(tau), prec)
    ....:     conjugates.append(conjugate_alpha)
    sage: conjugates
    [-6794.3278793864268700913699113 - 3.2311742677852643549664402034e-27*I,
     -0.036501034995017361360497864966 - 82.427712003227761738150494070*I,
     6.3732765567256884708225027190 - 3.4583745475524154096386135513*I,
     18.327164171482763936222950822 + 6.0318468311543925988327020435*I,
     18.327164171482763936222950826 - 6.0318468311543925988327020530*I,
     6.3732765567256884708225417027 + 3.4583745475524154096385923838*I,
     -0.036501034995017361360497840278 + 82.427712003227761738150494092*I]

And lastly we find the integral polynomial corresponding to these conjugates::

    sage: P = CF['y']; y = P.gen()
    sage: minimal_polynomial = prod([y - conjugate for conjugate in conjugates])
    sage: coefficients = [round(c.real()) + round(c.imag())*I for c in minimal_polynomial.coefficients()]
    sage: R = ZZ['X']; X = R.gen()
    sage: f = R(coefficients); f
    X^7 + 6745*X^6 - 327467*X^5 + 51857115*X^4 - 2319299751*X^3 + 41264582513*X^2 - 307873876442*X + 903568991567

We conclude that for the polynomial `f` above, we have `H = K[X]/(f)`.

.. NOTE::

    The polynomial `f` above is equal to the one in [GS1998], except for the
    sign at `X^3`, indicating a typo in [GS1998].

.. NOTE::

    Setting the complex precision ``prec`` too low will lead to an exception at
    the very last command ``f = R(coefficients)``; for example::

        sage: f = R(coefficients)
        Traceback (most recent call last):
        ...
        TypeError: Unable to coerce 5456*I + 903568988023 to an integer



Example 1 - `K = \QQ(\sqrt{-71})` and `f = \mathfrak{f}_2`
----------------------------------------------------------

..
    XXXXXXXXXXXXXXXXXXXXXXXXXXX

    Our field `K` is still the same. So in order to find a class invariant, we want
    to compute the action of `\theta-1 \mod 3` on `\QQ(\zeta_{48}, \mathfrak{f}_2)`.
    Because we are working with a higher level modular function this time, we
    compute the image of `\theta-1 \mod 3` in `GL_2(\hat{\ZZ})` up to a higher
    precision this time::

        sage: level = 48
        sage: x = J(Omod3(theta-1))
        sage: A = connecting_homomorphism_integral(x, level); A
        [   1 mod 48 576 mod 864]
        [  16 mod 48   17 mod 48]

    The reader can check that this `A` is consistent with the one from the previous
    example (it represents the same matrix, but with higher precision).

    We project `A` to `GL_2(\ZZ/48\ZZ)` and determine its action.

    XXXXXXXXXXXXXXXXXXXXXXXXXXXXX


Let us repeat the computation above, this time using Weber's `\mathfrak{f}_2`
modular function of level 48. ::

    sage: level = 48

Our field `K` is still the same. This time we compute
`(\mathcal{O}/48\mathcal{O})^*`::

    sage: Omod48 = O.quotient_ring(48, 'b')
    sage: Omod48star = K.ideal(48).idealstar(flag=2)  # flag=2 means compute generators
    sage: Omod48star.gens_values()
    (-12*theta + 13,
     12*theta - 23,
     6*theta + 19,
     -6*theta + 13,
     -16*theta + 1,
     16*theta + 17)

For each of these generators, we compute their action on `\zeta_{48}` and
`\mathfrak{f}_2`.
This time we use the fact that the standard generators `S` and `T` of
`SL_2(\ZZ)` act as follows on
`\QQ(\zeta_{48}, \mathfrak{f}, \mathfrak{f}_1, \mathfrak{f}_2)`:

- `S: (\mathfrak{f}, \mathfrak{f}_1, \mathfrak{f}_2) \mapsto
  (\mathfrak{f}, \mathfrak{f}_2, \mathfrak{f}_1)`;
- `T: (\mathfrak{f}, \mathfrak{f}_1, \mathfrak{f}_2) \mapsto
  (\zeta_{48}^{-1} \mathfrak{f}_1, \zeta_{48}^{-1} \mathfrak{f}, \zeta_{48}^2 \mathfrak{f}_2)`.

::

    sage: for x in Omod48star.gens_values():
    ....:     y = J(Omod48(x))
    ....:     A = connecting_homomorphism_integral(y, level)
    ....:     Amod48 = matrix_modulo(A, level)
    ....:     d = det(Amod48)  # A acts on zeta_48 by d-th powering
    ....:     U = ~iota(d) * Amod48  # U lies in SL_2(\ZZ/48\ZZ)
    ....:     STs = ST_factor(U)
    ....:     print("The action of {} on QQ(zeta_48, f2) is given by:".format(x))
    ....:     print("  zeta_48  ]--> zeta_48^{}".format(d))
    ....:     print_action_on_f2(d, STs)
    ....:
    The action of -12*theta + 13 on QQ(zeta_48, f2) is given by:
      zeta_48  ]--> zeta_48^37
      f2       ]--> zeta48^108*f2
    The action of 12*theta - 23 on QQ(zeta_48, f2) is given by:
      zeta_48  ]--> zeta_48^37
      f2       ]--> zeta48^108*f2
    The action of 6*theta + 19 on QQ(zeta_48, f2) is given by:
      zeta_48  ]--> zeta_48^31
      f2       ]--> zeta48^66*f2
    The action of -6*theta + 13 on QQ(zeta_48, f2) is given by:
      zeta_48  ]--> zeta_48^31
      f2       ]--> zeta48^18*f2
    The action of -16*theta + 1 on QQ(zeta_48, f2) is given by:
      zeta_48  ]--> zeta_48^17
      f2       ]--> zeta48^80*f2
    The action of 16*theta + 17 on QQ(zeta_48, f2) is given by:
      zeta_48  ]--> zeta_48^17
      f2       ]--> zeta48^32*f2

Let us write down these results a bit better readable:

=================   =========================   =================================
`x`                 `\zeta_{48}^x/\zeta_{48}`   `\mathfrak{f}_2^x/\mathfrak{f}_2`
=================   =========================   =================================
`-12 \theta + 13`   `\zeta_{48}^{36}`             `\zeta_{48}^{12}`
`12 \theta - 23`    `\zeta_{48}^{36}`             `\zeta_{48}^{12}`
`6 \theta + 19`     `\zeta_{48}^{30}`             `\zeta_{48}^{18}`
`-6 \theta + 13`    `\zeta_{48}^{30}`             `\zeta_{48}^{18}`
`-16 \theta + 1`    `\zeta_{48}^{16}`             `\zeta_{48}^{32}`
`16 \theta + 17`    `\zeta_{48}^{16}`             `\zeta_{48}^{32}`
=================   =========================   =================================

So this time we find as our class invariant `\alpha = \zeta_{48} \mathfrak{f}(\theta)`.

Next we compute the conjugates of this `\alpha` over `K`::

    sage: zeta48 = CF(exp(2*CF(pi)*CF(I)/48))
    sage: conjugates = []
    sage: for a in Cl:
    ....:     y = J(a)
    ....:     d, T = action_idele(y, level)  # f2^y == f2^(iota(d) * T)
    ....:     tau = apply_fractional_linear_transformation(T, theta)
    ....:     conjugate_alpha = zeta48**d * weber_f2(phi(tau), prec)
    ....:     if d % 8 in [3, 5]:
    ....:         conjugate_alpha *= -1
    ....:     conjugates.append(conjugate_alpha)
    sage: conjugates
    [0.46934985188028072404868708487 + 4.9303806576313237838233035330e-32*I,
     -1.3100509463724034874071664118 + 0.12777784412349129206137042647*I,
     -0.34481956767016108959807371471 + 1.0787151647032596770722167765*I,
     0.92019558810242421493511346058 + 0.33479097123881044615212608503*I,
     0.92019558810242662822461891407 - 0.33479097123880374817302036574*I,
     -0.34481956767014232858225163729 - 1.0787151647032821566490422829*I,
     -1.3100509463722659150544615130 - 0.12777784412463738678287924211*I]

And finally, these conjugates give us the minimal polynomial `f` of `\alpha`
over `K`::

    sage: P = CF['y']; y = P.gen()
    sage: minimal_polynomial = prod([y - conjugate for conjugate in conjugates])
    sage: coefficients = [round(c.real()) + round(c.imag())*I for c in minimal_polynomial.coefficients()]
    sage: R = ZZ['X']; X = R.gen()
    sage: f = R(coefficients); f
    X^7 + X^6 - X^5 - X^4 - X^3 + X^2 + 2*X - 1

.. NOTE::

    This polynomial also differs from the one in [GS1998]:
    the constant term has opposite sign.

.. TODO::

    In our computations of minimal polynomials of `\alpha` over `K`, we
    assume the polynomial will lie in `\ZZ[X]`. Why not `\mathcal{O}[X]`?
    Think of a clean way to do it over `\mathcal{O}[X]` in general.


Example 2 - `K = \QQ(\sqrt{-145})` and `f = \mathfrak{f}`
---------------------------------------------------------

Next up, example 2::

    sage: R = ZZ['x']; x = R.gen()
    sage: K.<theta> = NumberField(x**2+145)
    sage: K.discriminant()
    -580
    sage: O = K.maximal_order()
    sage: O.basis()
    [1, theta]
    sage: J = IdeleGroup(K)
    sage: level = 48
    sage: Omod48 = O.quotient_ring(48, 'b')
    sage: Omod48star = K.ideal(48).idealstar(flag=2)
    sage: for x in Omod48star.gens_values():
    ....:     y = J(Omod48(x))
    ....:     A = connecting_homomorphism_integral(y, level)
    ....:     Amod48 = matrix_modulo(A, level)
    ....:     d = det(Amod48) # A acts on zeta_48 by d-th powering
    ....:     U = ~iota(d) * Amod48 # U lies in SL_2(\ZZ/48\ZZ)
    ....:     STs = ST_factor(U)
    ....:     print("The action of {} on QQ(zeta_48, f) is given by:".format(x))
    ....:     print("  zeta_48 ]--> zeta_48^{}".format(d))
    ....:     print_action_on_f(STs)
    ....:
    The action of 18*theta + 1 on QQ(zeta_48, f) is given by:
      zeta_48 ]--> zeta_48^37
      f       ]--> zeta48^12*f
    The action of 16*theta + 17 on QQ(zeta_48, f) is given by:
      zeta_48 ]--> zeta_48^17
      f       ]--> f
    The action of 15*theta + 16 on QQ(zeta_48, f) is given by:
      zeta_48 ]--> zeta_48^1
      f       ]--> zeta48^-48*f
    The action of 19 on QQ(zeta_48, f) is given by:
      zeta_48 ]--> zeta_48^25
      f       ]--> f

We conclude that `\mathfrak{f}^2 / \sqrt{2}` is invariant.
We determine the minimal polynomial of `\mathfrak{f}^2(\theta) / \sqrt{2}`::


"""



from shimura import *
from modular import *


###########################
# Phase 1. Initialization #
###########################

R = ZZ['x'];
x = R.gen()
N = ZZ(3)  # level of our modular function gamma_2
K = NumberField(x**2 + x + 18, 'theta')
theta = K.gen()
J = IdeleGroup(K)

######################################
# Phase 2. Compute a class invariant #
######################################
O = K.maximal_order()
OmodN = O.quotient_ring(N, 'b')
OmodNstar = K.ideal(N).idealstar(flag=2)  # flag=2 means compute generators
print("(O/{}O)^* is generated by {}".format(N, OmodNstar.gens_values()))

for x in OmodNstar.gens_values():
    y = J(OmodN(x))
    A = connecting_homomorphism_integral(y, N)
    AmodN = matrix_modulo(A, N)
    d = det(AmodN)  # A acts on zeta_N by d-th powering
    U = ~iota(d) * AmodN  # U lies in SL_2(\ZZ/N\ZZ)
    STs = ST_factor(U)
    print("The action of {} on QQ(zeta_3, gamma_2) is given by:".format(x))
    print("  zeta_3  ]--> zeta_3^{}".format(d))
    print(STs)
    print_action_on_gamma_2(STs);

print("#################################################################")
print("# Deduce (by hand) the invarent alpha = zeta_3 * gamma_2(theta) #")
print("#################################################################")

#################################################################
# Phase 3. Compute the conjugates of our class invariant over K #
#################################################################
prec = 100  # We do our complex calculations with 100 bits of precision
CF = ComplexField(prec=prec)
zeta3 = CF(exp(2*CF(pi)*CF(I)/3))
phi = K.embeddings(CF)[1]  # the embedding where theta has positive imaginary part

conjugates = []
Cl = ray_class_group(K, Modulus(K.ideal(1))) # ray class group of modulus 1, i.e. ideal class group
for a in Cl:
    y = J(a)
    d, T = action_idele(y, N)

    tau = apply_fractional_linear_transformation(T, theta)
    conjugate_alpha = zeta3**d * gamma_2(phi(tau), prec)
    conjugates.append(conjugate_alpha)
    print("alpha^{} = {}".format(a, conjugate_alpha))

##################################################################
# Phase 4. Compute the minimal polynomial of our class invariant #
##################################################################
Cy = CF['y']; (y,) = Cy._first_ngens(1)
minimal_polynomial = prod([y - conjugate for conjugate in conjugates])
coefficients = [round(c.real()) + round(c.imag())*I for c in minimal_polynomial.coefficients()]
R = ZZ['x']; x = R.gen()
try:
    f = R(coefficients)
except TypeError:
    raise ValueError("Precision of our numerical complex computations ({} bits) is too low!".format(CF.prec()))

print("Hilbert class polynomial:")
print(f)


