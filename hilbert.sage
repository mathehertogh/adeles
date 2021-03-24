r"""
In this file we will compute the Hiblert class field of the imaginary quadratic
field K of discriminant -71 using Shimura's reciprocity law.
We do this using a class invariant obtained from the modular function `\gamma_2`
of level 3.
"""
load("matricies.sage")



def connecting_homomorphism(K, x):
    r"""
    This function implements Shimura's connecting homomorphism `\hat{K}^* \to
    GL(\hat{\QQ})`

    INPUT:

    - ``x`` -- an idele over ``K``

    OUTPUT:

    The transpose of the matrix describing the multiplication by ``x`` on the
    free `\hat{\QQ}`-module `\hat{K} = \theta \cdot \hat{\QQ} + \hat{\QQ}` with
    respect to the basis `(\theta, 1)`.
    """
    x = Adeles(K)(x).finite
    t, s = x.to_profinite_rational_vector()
    C, B, _ = K.gen().minpoly().coefficients()
    return matrix(ProfiniteNumbers(QQ), [[t-B*s, -C*s], [s, t]])


def compute_action_GL2ZmodN(A, N):
    """
    Compute the action of ``A`` on the level-``N`` modular function gamma_2
    """
    from sage.modular.arithgroup.arithgroup_perm import sl2z_word_problem
    from sage.modular.local_comp.liftings import lift_matrix_to_sl2z
    d = det(A)
    iota_d = matrix(Zmod(N), [[1, 0], [0, d]])
    B = A / iota_d
    C = SL2Z(lift_matrix_to_sl2z(B.list(), N))
    G.<S,T> = FreeGroup() # S=(0,-1;1,0), T=L=(1,1;0,1), R=(1,0;1,1)
    R = S^3 * T^-1 * S
    # We want to factor C into a product of S's and T's (as in the
    # paper of Alice Gee). sl2z_word_problem() factors it into L's and R's.
    factorization = G.one()
    for is_R, e in sl2z_word_problem(C):
        if is_R:
            factorization = factorization * R^e
        else:
            factorization = factorization * T^e
    print("The action of the matrix\n{}\nis given by:".format(A))
    print("    zeta_{} |--> zeta_{}^{}".format(N, N, d))
    # S acts trivially on gamma_2
    # T acts as multiplication by zeta^-1 on gamma_2
    # We compute the number of T's in C:
    n_Ts = sum([sign(g) for g in factorization.substitute(S=G.one()).Tietze()])
    # so we have the action:
    # gamma_2 |--C--> zeta_N^(-n_Ts) * gamma_2
    # and as zeta maps to zeta^d by iota(d), we obtain:
    # gamma_2 |--C--> zeta_N^(-n_Ts) * gamma_2 |--iota(d)--> zeta_N^(-n_Ts*d) * gamma_2
    e = -n_Ts*d
    print("    gamma_2 |--> zeta_{}^{} * gamma_2".format(N, e))
    return d, e

def compute_action(K, N, x):
    gx = connecting_homomorphism(K, x)
    gx_modN = matrix_modulo(gx, N)
    compute_action_GL2ZmodN(gx_modN, N)






###########################
# Phase 1. Initialization #
###########################
CF = ComplexField(prec=53)
D = -71 # discriminant of our imaginary quadratic field K
N = 3 # level of our modular function gamma_2
B, C = 1, 18
R.<x> = ZZ[]
K.<theta> = NumberField(x^2 + B*x + C)
sqrtD = 2*theta + 1
theta_complex = K.embeddings(CF)[1].im_gens()[0] # the image of theta in CF with positive imaginary part
complex_embedding = Hom(K, CF)([theta_complex], check=False) # K --> CF, sending theta to theta_complex

######################################
# Phase 2. Compute a class invariant #
######################################
OmodNstar = K.ideal(N).idealstar(flag=2) # flag=2 means compute generators
for x in OmodNstar.gens_values():
    compute_action(K, N, x)
print("\n#####\n# Deduce (by hand) the invarent element alpha = zeta_3 * gamma_2(theta)\n#####")

#################################################################
# Phase 3. Compute the conjugates of our class invariant over K #
#################################################################
conjugates = []
J = IdeleGroup(K)
Ak = Adeles(K)
Cl = ray_class_group(K, Modulus(K.ideal(1))) # ray class group of modulus 1, i.e. ideal class group
L = BinaryQF_reduced_representatives(D)
for bqf in L:
    a, b, c = bqf.polynomial().coefficients()

    # Use direct formula's in the paper (formula (17) at page 32 of Gee) to
    # check our own results via ideles.
    if not 3.divides(a):
        ux_accoring_to_paper = matrix(Zmod(3), [[a, (b-1)/2], [0, 1]])
    elif not 3.divides(c):
        ux_accoring_to_paper = matrix(Zmod(3), [[(-b-1)/2, -c], [1, 0]])
    else:
        ux_accoring_to_paper = matrix(Zmod(3), [[(-b-1)/2-a, (1-b)/2-c], [1, -1]])

    finite = {}
    for q, e in factor(K.ideal(a)):
        if q.divides(c):
            x_q = -b/2 + sqrtD/2 - a
        else:
            x_q = -b/2 + sqrtD/2
        finite[q] = (x_q, 1)
    x = J(None, None, finite)
    gx = connecting_homomorphism(K, x)
    beta, alpha = (-b/2+sqrtD/2).vector()
    M = matrix(QQ, [[alpha, beta], [0, a]])
    ux = gx * ~M

    n = 0
    while not is_defined_modulo(ux, N):
        d = (matrix_modulus(ux) / N).denominator()
        d *= gcd(N, matrix_denominator(ux))
        for p in d.prime_divisors():# + primes_first_n(n):
            for q in K.primes_above(p):
                if q in x.finite:
                    x_q, i_q = x.finite[q]
                    x.finite[q] = (x_q, i_q+1)
                else:
                    x.finite[q] = (a, 1)
        gx = connecting_homomorphism(K, x)
        ux = gx * ~M
        n += 1

    print("\nux at (a,b,c)=({},{},{}):\n{}".format(a, b, c, ux))
    ux_modN = matrix_modulo(ux, N)
    assert ux_modN == ux_accoring_to_paper, "ERROR!!!!!!!! According to the paper we should get\n{}".format(ux_accoring_to_paper)
    d, e = compute_action_GL2ZmodN(ux_modN, N)

    tau = (-b + sqrtD) / (2*a)
    zeta3 = CF(exp(2*CF(pi)*CF(I)/3))
    conjugate_theta = zeta3^ZZ(e) * weber_gamma_2(complex_embedding(tau))
    conjugate_alpha = zeta3^ZZ(d) * conjugate_theta
    print("conjugate of theta is: {}".format(conjugate_theta))
    conjugates.append(conjugate_alpha)

    #####
    # Third way of computing these Z/3Z-matrices, this time using our
    # ray class group to idele conversion.
    #####
    aa = K.ideal(a, (-b+sqrtD)/2)
    x = J(Cl(aa))
    gx = connecting_homomorphism(K, x)
    beta, alpha = (-b/2+sqrtD/2).vector()
    M = matrix(QQ, [[alpha, beta], [0, a]])
    ux = gx * ~M

    while not is_defined_modulo(ux, N):
        d = (matrix_modulus(ux) / N).denominator()
        d *= gcd(N, matrix_denominator(ux))
        for p in d.prime_divisors():
            for q in K.primes_above(p):
                if q in x.finite:
                    x_q, i_q = x.finite[q]
                    x.finite[q] = (x_q, i_q+1)
                else:
                    x.finite[q] = (K(1), 1)
        gx = connecting_homomorphism(K, x)
        ux = gx * ~M

    new_ux_modN = matrix_modulo(ux, N)
    print("this is ideal class {}".format(Cl(aa)))
    if new_ux_modN != ux_modN and new_ux_modN != -1*ux_modN:
        print("ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! new_ux_modN is now:\n{}".format(new_ux_modN))
    #####
    # End of the third way
    #####




##################################################################
# Phase 4. Compute the minimal polynomial of our class invariant #
##################################################################
Cy.<y> = CF[]
minimal_polynomial = prod([y - conjugate for conjugate in conjugates])
coefficients = [round(c.real()) + round(c.imag())*I for c in minimal_polynomial.coefficients()]
R.<x> = ZZ[]
try:
    f = R(coefficients)
except TypeError:
    raise ValueError("Precision of our numerical complex computations ({} bits) is too low!".format(CF.prec()))
if not NumberField(f, 'a').is_isomorphic(NumberField(hilbert_class_polynomial(-71), 'b')):
    raise ValueError("Computed hibert class field differs from the standard sage one!")
return f