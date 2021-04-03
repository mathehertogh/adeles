r"""
In this file we will compute the Hiblert class field of the imaginary quadratic
field K of discriminant -71 using Shimura's reciprocity law.
We do this using a class invariant obtained from the modular function `\gamma_2`
of level 3.
"""

from sage.all_cmdline import *

from profinite_number import ProfiniteNumbers
from adele import Adeles
from idele import IdeleGroup
from ray_class_group import Modulus, ray_class_group
from matrix import *


def weber_f(x):
    r"""
    Evaluate Weber's `f` function in ``x``
    """
    zeta48 = CF(exp(2*CF(pi)*CF(I)/48))
    return eta((x+ZZ(1))/ZZ(2)) / (zeta48 * eta(x))

def gamma_2(x):
    r"""
    Evaluate Weber's `\gamma_2` function in ``x``
    """
    fx = weber_f(x)
    return (fx**ZZ(24) - ZZ(16)) / fx**ZZ(8) 

def iota(d, N=3):
    return matrix(Zmod(N), [[1, 0], [0, d]])

def connecting_homomorphism(K, x, N=3):
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
    Qhat = ProfiniteNumbers(QQ)
    x = Adeles(K)(x).finite()
    t, s = x.to_profinite_rational_vector(enclosure=False)
    C, B, _ = K.gen().minpoly().coefficients()
    return matrix(Qhat, [[t-B*s, -C*s], [s, t]])

def apply_fractional_linear_transformation(M, tau):
    r"""
    Apply the fractional linear transformation given by ``M`` to ``tau``

    INPUT:

    - ``M`` -- a matrix in `GL_2(\QQ)` with positive determinant
    - ``tau`` -- a field element

    OUTPUT:
    
    If ``M = (a, b; c, d)``, then return the element ``(a*tau+b)/(c*tau+d)``.
    """
    a, b, c, d = M.list()
    return (a*tau + b) / (c*tau + d)

def action_SL2ZmodN(U, N=3):
    r"""
    Compute the action of the `SL_2(\ZZ/N\ZZ)`-matrix ``U``on `\gamma_2`

    INPUT:

    - ``U`` -- a matrix in `SL_2(\ZZ/N\ZZ)`, where `N` is some integer divisible
               by 3 (the level of `\gamma_3`)

    OUTPUT:

    An integer `e` \in \{0, 1, 2\}` such that `\gamma_2^U = \zeta_3^e \gamma_2`.
    Here `zeta_3` is `exp(2*\pi*i/3)`.
    """
    factorization = ST_factor(U)
    G = factorization.parent()
    # S acts trivially on gamma_2, so we are only interested in the number of
    # T's in U.
    n_Ts = sum([sign(g) for g in factorization.substitute(S=G.one()).Tietze()])
    # T acts as multiplication by zeta^-1 on gamma_2, so U acts as
    # multiplication by zeta^-n_Ts on gamma_2
    e = Zmod(N)(-n_Ts)
    """
    TODO: the above seems to be the correct e, but before in alice_gee.sage we
    did this below. So there we had e = -n_Ts*d, an extra factor d (the 
    determinant of the original matrix). Check that this is now correct and we
    should not include our "iota"-part in this by doing *d as well.
    print("The action of the matrix\n{}\nis given by:".format(A))
    print("    zeta_{} |--> zeta_{}^{}".format(N, N, d))
    # so we have the action:
    # gamma_2 |--C--> zeta_N^(-n_Ts) * gamma_2
    # and as zeta maps to zeta^d by iota(d), we obtain:
    # gamma_2 |--C--> zeta_N^(-n_Ts) * gamma_2 |--iota(d)--> zeta_N^(-n_Ts*d) * gamma_2
    e = -n_Ts*d
    print("    gamma_2 |--> zeta_{}^{} * gamma_2".format(N, e))
    """
    return ZZ(e)

def action_iota(S, N=3):
    r"""
    Compute the action of the matrix ``S`` in the image of `\iota`.

    `\iota` is the map `\hat{\ZZ}^* \to GL_2(\hat{\ZZ})` sending `d` to the
    matrix `(1, 0; 0, d)`.

    INPUT:

    - ``S`` -- a 2x2-matrix of the form `(1, 0; 0, d)` with `d \in \hat{\ZZ}^*`
               that is defined modulo 3 (the level of `\gamma_2`)

    OUTPUT:

    An integer `d \in \{0, 1, 2\}` such that `\zeta_3^S = \zeta_3^d`.
    Here `zeta_3` is `exp(2*\pi*i/3)`.
    """
    return ZZ(Zmod(N)(S[1,1].value()))

def action_GL2ZmodN(A, N=ZZ(3)):
    r"""
    Compute the action of the `GL_2(\ZZ/N\ZZ)`-matrix ``A`` on `\gamma_2`

    OUPTPUT:

    A pair of integers `(d, e)` such that the action of ``A`` is given by:
        - `\zeta_3 \mapsto \zeta_3^d`
        - `\gamma_2 \mapsto \zeta_3^e \gamma_2`
    """
    U = ~iota(det(A)) * A
    d = ZZ(det(A))
    e = action_SL2ZmodN(U)
    return d, e


def action_idele(x, N=ZZ(3)):
    r"""
    Compute the action on `\gamma_2` of the image of the idele ``x`` under
    Shimura's connecting homomorphism

    Denote the connecting homomorphism by `g: \hat{K}^* \to GL_2(\hat{\QQ})`.

    OUTPUT:

    A triple `(d, e, M)` satisfying:
        - `e` and `d` are integers
        - `M \in GL_2^+(\QQ)`
        - `\zeta_3^g(x) = \zeta_3^d
        - `\gamma_2^g(x) = \zeta_3^e \gamma_2 \circ M`
    """
    K = x.parent().number_field()

    A = connecting_homomorphism(K, x)
    B, M = GL2Zhat_GL2QQ_factor(A)

    while not (in_GL2Zhat(B) and N.divides(ZZ(matrix_modulus(B)))):
        ps = primes_missing_precision(B, N)
        x.increase_precision(ps)
        A = connecting_homomorphism(K, x)
        B, M = GL2Zhat_GL2QQ_factor(A)

    C = matrix_modulo(B, N)

    d = ZZ(det(C))
    U = ~iota(det(C)) * C
    e = action_SL2ZmodN(U)

    return d, e, M

def print_action(d, e, M, x=None):
    r"""
    Print a description of the action given by ``d`` and ``e`` as returend by
    action_idele().
    """
    if x is None:
        print("Action given by:")
    else:
        print("Action of {} given by:".format(x))
    print("  d={};  zeta_3  ]--> zeta_3^{}".format(d, d))
    print(r"  e={};  gamma_2 ]--> zeta_3^{} * gamma_2 \circ M".format(e, e))
    print("where M equals:\n{}".format(M))



###########################
# Phase 1. Initialization #
###########################
CF = ComplexField(prec=100)
zeta3 = CF(exp(2*CF(pi)*CF(I)/3))

R = ZZ['x'];
x = R.gen()
D = ZZ(-71) # discriminant of our imaginary quadratic field K
N = ZZ(3) # level of our modular function gamma_2
B, C = ZZ(1), ZZ(18) 
K = NumberField(x**2  + B*x + C, 'theta')
theta = K.gen()
O = K.maximal_order()
sqrtD = 2*theta + 1 
theta_complex = K.embeddings(CF)[1].im_gens()[0] # the image of theta in CF with positive imaginary part
complex_embedding = Hom(K, CF)([theta_complex], check=False) # K --> CF, sending theta to theta_complex
J = IdeleGroup(K)

######################################
# Phase 2. Compute a class invariant #
######################################
OmodN = O.quotient_ring(N, 'b')
OmodNstar = K.ideal(N).idealstar(flag=2) # flag=2 means compute generators
"""
for x in OmodNstar.gens_values():
    x = OmodN(x)
    d, e, M = action_idele(J(x))
    print(); print_action(d, e, M, x)
    tau = apply_fractional_linear_transformation(M, theta_complex)
    conjugate_theta = zeta3**ZZ(e) * gamma_2(tau) # TODO: e or e+d??
    conjugate_alpha = zeta3**ZZ(d) * conjugate_theta
    print("alpha^{} = {}".format(x, conjugate_alpha))

print("\n#####\n# Deduce (by hand) the invarent element alpha = zeta_3 * gamma_2(theta)\n#####")

#################################################################
# Phase 3. Compute the conjugates of our class invariant over K #
#################################################################
conjugates = []
Cl = ray_class_group(K, Modulus(K.ideal(1))) # ray class group of modulus 1, i.e. ideal class group
for a in Cl:
    print("\ncomputing alpha^{}".format(a))
    x = J(a)
    d, e, M = action_idele(x)
    print_action(d, e, M)

    tau = apply_fractional_linear_transformation(M, theta_complex)
    conjugate_theta = zeta3**ZZ(e) * gamma_2(tau) # TODO: e or e+d??
    conjugate_alpha = zeta3**ZZ(d) * conjugate_theta
    conjugates.append(conjugate_alpha)
    print("alpha^{} = {}".format(a, conjugate_alpha))

##################################################################
# Phase 4. Compute the minimal polynomial of our class invariant #
##################################################################
Cy = CF['y']; (y,) = Cy._first_ngens(1)
minimal_polynomial = prod([y - conjugate for conjugate in conjugates])
print(minimal_polynomial)
coefficients = [round(c.real()) + round(c.imag())*I for c in minimal_polynomial.coefficients()]
R = ZZ['x']; (x,) = R._first_ngens(1)
try:
    f = R(coefficients)
except TypeError:
    pass # raise ValueError("Precision of our numerical complex computations ({} bits) is too low!".format(CF.prec()))
#if not NumberField(f, 'a').is_isomorphic(NumberField(hilbert_class_polynomial(D), 'b')):
    #raise ValueError("Computed hibert class field differs from the standard sage one!")
#print(f)

"""