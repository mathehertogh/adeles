#from sage.modular.arithgroup.arithgroup_perm import eval_sl2z_word, sl2z_word_problem
#A = SL2Z([7,8,-50,-57])
#factorization = sl2z_word_problem(A)
#print(A == eval_sl2z_word(factorization))

load("ideles.sage")

def split(A):
    """
    Split the finite adele matrix ``A`` into a rational and an integral matrix

    INPUT:

    - ``A`` -- a matrix with entries profinite numbers over a number field `K`

    OUTPUT:

    A pair `(B, C)` such that:
    - `B` is a matrix with entries in the number field `K`
    - `C` is a matrix with entries profinite integers over `K`
    - ``A == B * C`` holds

    .. NOTE::

        This operation looses precision in general. Every entry of ``B *C``
        always represents at least all profinite integers that the corresponding
        entry of ``A`` represents. See the example below.

    EXAMPLES::


    """
    K = A.parent().base().base()

    # First we compute an appropriate inverse Binv of B such that Binv*A becomes
    # integral.
    a = c = lcm(A[0,0].denominator, A[0,1].denominator)
    b = d = lcm(A[1,0].denominator, A[1,1].denominator)
    a *= 2 # Make determinant of Binv non-zero: 
    Binv = MatrixSpace(K, 2).matrix([a, b, c, d])

    # Now we get our integral C:
    C = Binv*A
    # And B, the inverse of Binv:
    B = Binv.inverse()

    return B, C

# Qhat = ProfiniteNumbers(QQ)
# A = MatrixSpace(Qhat, 2).matrix([Qhat(3, 240, 2), Qhat(3, 240, 40),
#                                Qhat(1, 240, 3), Qhat(7, 240, 2)])
# B, C = split(A)



def X(b):
    return SL(2, b.parent())([1, b, 0, 1])

def Y(c):
    return SL(2, c.parent())([1, 0, c, 1])

def W(a):
    return X(a) * Y(-1/a) * X(a)

def SL2K_factor(A):
    r"""
    Factor the SL(2, K) matrix ``A`` into factors of the form `(1, b; 0, 1)` and
    `(1, 0; c, 1)`

    INPUT:

    - ``A`` -- a matrix in `SL(2, K)` with `K` a *field*

    OUTPUT:

    A list `f` of matrices of the form `(1, b; 0, 1)` and `(1, 0; c, 1)` with
    `b, c \in K`, such that ``prod(f) == A``.
    """
    print("---- in SL2K_factor()---")
    print(A)
    print("----")
    original_A = A
    K = A.parent().base()
    factorization = []
    words = []

    # Create non-zero at top-left:
    if A.matrix()[0,0] == 0:
        # Since determinant is one, A[1,0] must be non-zero.
        # So just add the bottom-row to the top-row.
        A = X(K(1)) * A
        factorization.append(X(K(-1)))
        words.append(("X", K(-1)))

    # Create zero at bottom-left:
    if A.matrix()[1,0] != 0:
        c = -A.matrix()[1,0] / A.matrix()[0,0]
        A = Y(c) * A
        factorization.append(Y(-c))
        words.append(("Y", -c))

    # Create zero at top-right:
    if A.matrix()[0,1] != 0:
        b = -A.matrix()[0,1] / A.matrix()[1,1]
        A = X(b) * A
        factorization.append(X(-b))
        words.append(("X", -b))

    a = A.matrix()[0,0]
    # Now we have A == W(a)*W(K(-1))).
    factorization.append(X(a))
    factorization.append(Y(-1/a))
    factorization.append(X(a-1))
    factorization.append(Y(K(1)))
    factorization.append(X(K(-1)))
    words.append(("X", a))
    words.append(("Y", -1/a))
    words.append(("X", a-1))
    words.append(("Y", K(1)))
    words.append(("X", K(-1)))

    return factorization, words


def Serge_Langs_method():
    K = QQ
    A = SL(2, K)([43, 30, 10, 7])
    f, w = SL2K_factor(A)
    assert prod(f) == A, "factorization in X(b)'s and Y(c)'s went wrong!"

    Qhat = ProfiniteNumbers(QQ)
    # Create a matrix with common modulus 100:
    n = 4
    N = 7^n * 100
    A = MatrixSpace(Qhat, 2).matrix([Qhat(13, 7^n * 100, 7), Qhat(14, 7^n * 100, 7),
                                     Qhat(3, 7^n * 100, 7), Qhat(7,  7^n * 100, 7)])

    x = SL(2, QQ).one()
    for p in N.prime_divisors():
        Ap = [[A[0,0][p], A[0,1][p]], [A[1,0][p], A[1,1][p]]]
        print("in loop at p={}".format(p))
        print(matrix(Ap))
        print("det(Ap) = {}".format(matrix(Ap).det()))
        Ap = SL(2, Qp(p))(Ap)
        print("...")
        f, w = SL2K_factor(Ap)
        n = sum([-2*min(0, b.valuation(p)) for Z, b in w])
        print("at p={} we get n={}".format(p, n))
        Zr_list = []
        for Z, b in w:
            e = min(0, b.valuation(p))
            b_int = b * p^(-e)
            r = QQ( ZZ(b_int) * p^e )
            Zr = X(r) if Z == "X" else Y(r)
            Zr_list.append(Zr)
        print("Zr_list:")
        for Zr in Zr_list:
            print(Zr)
        xp = prod(Zr_list)
        print(xp)
        x *= xp


def marco_skype_attempt_2021_03_01(A):
    r"""
    INPUT: a matrix `A \in SL(2; Qhat)` with all entries the same modulus
    """
    print(A); print()

    # We keep track of the transformations in the matrix ``transformation``
    transformation = matrix(QQ, [[1, 0], [0, 1]])

    # Make A integral.
    print("# Make A integral.")
    D = lcm([A[0,0].denominator, A[0,1].denominator, A[1,0].denominator, A[1,1].denominator])
    transformation *= D
    A = transformation * A
    A = MatrixSpace(ProfiniteIntegers(), 2)(A)
    print(A); print()

    # Create a zero at the bottom-left entry of A.
    print("# Create a zero at the bottom-left entry of A.")
    while A[1,0].value != 0:
        if A[0,0].value < A[1,0].value:
            trans = matrix(ZZ, [[0, 1], [1, 0]])
            A = trans * A # swap rows
            transformation = trans * transformation
            print(A); print()
            if A[1,0].value.is_zero():
                break
        scalar = A[0,0].value // A[1,0].value
        trans = matrix(ZZ, [[1, -scalar], [0, 1]])
        A = trans * A
        transformation = trans * transformation
        #A -= MSint.matrix([scalar*A[1], (0, 0)])
        print(A); print()

    # Make sure that 0 <= A[0,1].value < A[1,1].value.
    print("# Make sure that 0 <= A[0,1].value < A[1,1].value.")
    scalar = A[0,1].value // A[1,1].value # determinant of original (hence of current) A is 1, so A[1,1]!=0
    trans = matrix(ZZ, [[1, -scalar], [0, 1]])
    A = trans * A
    transformation = trans * transformation
    #A -= MSint.matrix([scalar*A[1], (0, 0)])
    print(A); print()

    # Make A[1,1] coprime to the modulus of A.
    print("# Make A[1,1] coprime to the modulus of A.")
    g = gcd(A[1,1].modulus, A[1,1].value)
    A[1,0].value //= g
    A[1,0].modulus //= g
    A[1,1].value //= g
    A[1,1].modulus //= g
    transformation = matrix(QQ, [[1, 0], [0, 1/g]]) * transformation
    # We now have: (current A) == transformation * (previous A)
    print(A); print()

    # Make A[0,1] (top-right) zero.
    print("# Make A[0,1] (top-right) zero.")
    f = A[1,1].value.inverse_mod(A[1,1].modulus)
    trans = matrix(ZZ, [[1, -f*A[0,1].value], [0, 1]])
    A = trans * A
    transformation = trans * transformation
    #A -= MSint.matrix([f*A[0,1].value*A[1], (0,0)])
    print(A); print()

    # Make A[0,0] coprime to the modulus of row 0.
    print("# Make A[0,0] coprime to the modulus of row 0.")
    g = gcd(A[0,0].modulus, A[0,0].value)
    A[0,0].value //= g
    A[0,0].modulus //= g
    A[0,1].value //= g
    A[0,1].modulus //= g
    transformation = matrix(QQ, [[1/g, 0], [0, 1]]) * transformation
    print(A); print()

    # Compute "iota"-part:
    print("# Compute 'iota'-part:")
    p = A[0,0] * A[1,1]
    t = p.value.inverse_mod(p.modulus)
    iota = MatrixSpace(ZZ, 2).matrix([1, 0, 0, t])
    A = A * iota
    print(A); print()
    return transformation


def QQhat(x, m):
    d = lcm(x.denominator(), m.denominator())
    Qhat = ProfiniteNumbers(QQ)
    return Qhat(x*d, m*d, d)

def modulus(A):
    m = [[0, 0], [0, 0]]
    all_equal = True
    for i in [0, 1]:
        for j in [0, 1]:
            m[i][j] = A[i,j].numerator.modulus / A[i,j].denominator
            if m[i][j] != m[0][0]:
                all_equal = False
    if all_equal:
        return m[0][0]
    return matrix(QQ, m)
    

def random_SL2Qhat_element():
    Qhat = ProfiniteNumbers(QQ)
    N = 20

    m = QQ.random_element(N, N)

    a = QQ.random_element(N, N)
    while a.is_zero():
        a = QQ.random_element(N, N)

    d = QQ.random_element(N, N)
    while d.is_zero():
        d = QQ.random_element(N, N)

    bc = a*d-1


    b = QQ.random_element(N, N)
    while b.is_zero() or (bc/b).height() > 2*N:
        b = QQ.random_element(N, N)

    c = bc / b

    return SL(2, Qhat)([[QQhat(a, m), QQhat(b, m)], [QQhat(c, m), QQhat(d, m)]]).matrix()



"""
Qhat = ProfiniteNumbers(QQ)
Zhat = ProfiniteIntegers(ZZ)

# We make an A with modulus 10 at every entry.
A = random_SL2Qhat_element()
print(A)
T = marco_skype_attempt_2021_03_01(A)

B = ~T
C = T*A
"""