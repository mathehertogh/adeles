

def matrix_modulus(A):
    return gcd([a.modulus() for a in A.list()])

def matrix_denominator(A):
    return lcm([a.denominator for a in A.list()])

def is_defined_modulo(A, N):
    return gcd(N, matrix_denominator(A)) == 1 and matrix_modulus(A) / N in ZZ

def matrix_modulo(A, N):
    r"""
    Project a matrix over `\hat{\ZZ}` to a matrix over `\ZZ/N\ZZ`

    INPUT:

    - ``A`` -- a matrix with entries *integral* elements of `\hat{\QQ}`
    - ``N`` -- an integer
    """
    rows_modN = []
    for row in A:
        row_modN = []
        for a in row:
            if not N.divides(a.modulus()):
                raise ValueError("not every entry of A is defined modulo {}".format(N))
            a_modN = Zmod(N)(a.value())
            row_modN.append(a_modN)
        rows_modN.append(row_modN)
    return matrix(Zmod(N), rows_modN)


def SUM_factorization(A):
    r"""
    Compute the SUM-factorization of the `GL_2(\hat{\QQ})`-matrix ``A``

    INPUT:

    - ``A`` -- a matrix in `GL_2(\hat{\QQ})`

    OUTPUT:

    A triple (S, U, M) satisfying:
    - S*U*M == A
    - S = (1, 0; 0, d) with `d \in \hat{\ZZ}^*`
    - U in `SL_2(\hat{\ZZ})`
    - M in `GL_2^+(\QQ)`

    Increasing the precision of ``A`` at a prime `p` will ultimately increase
    the precision of the returned ``S`` and ``U`` at `p`.
    """
    Zhat = ProfiniteIntegers()
    print(A); print()

    # Make all of A's entries integral.
    print("# Make all of A's entries integral.")
    denominator = lcm([a.denominator for a in A.list()])
    omega = matrix(QQ, [[denominator, 0], [0, denominator]])
    A = A * omega
    M = ~omega
    A = MatrixSpace(Zhat, 2)(A)
    #A = matrix(Zhat, [[Zhat._from_profinite_number(A[0,0]), Zhat._from_profinite_number(A[0,1])], [Zhat._from_profinite_number(A[1,0]), Zhat._from_profinite_number(A[1,1])]])
    print(A); print()

    # Create a zero at the bottom-left entry of A.
    print("# Create a zero at the bottom-left entry of A.")
    while A[1,0].value != 0:
        if A[1,1].value < A[1,0].value:
            alpha = matrix(ZZ, [[0, 1], [1, 0]])
            A = A * alpha # swap columns
            M = ~alpha * M
            print(A); print()
            if A[1,0].value.is_zero():
                break
        scalar = A[1,1].value // A[1,0].value
        alpha = matrix(ZZ, [[1, -scalar], [0, 1]])
        A = A * alpha
        M = ~alpha * M
        print(A); print()

    # Make top-left entry of A an element of `\hat{\ZZ}^*` (i.e. make its value
    # coprime to its modulus).
    print(r"# Make top-left entry of A an element of `\hat{\ZZ}^*`")
    g = gcd(A[0,0].modulus, A[0,0].value)
    beta = matrix(QQ, [[1/g, 0], [0, 1]])
    # We program the statement A = A * beta ad hoc:
    A[0,0].value //= g
    A[0,0].modulus //= g
    A[1,0].value //= g
    A[1,0].modulus //= g
    M = ~beta * M
    print(A); print()

    # Make top-right entry of A zero.
    print("# Make top-right entry of A zero.")
    f = A[0,0].value.inverse_mod(A[0,0].modulus)
    gamma = matrix(ZZ, [[1, -f*A[0,1].value], [0, 1]])
    A = A * gamma
    M = ~gamma * M
    print(A); print()

    # Make bottem-right entry of A an element of `\hat{\ZZ}^*` (i.e. make its
    # value coprime to its modulus).
    print(r"# Make bottem-right entry of A an element of `\hat{\ZZ}^*`")
    g = gcd(A[1,1].modulus, A[1,1].value)
    delta = matrix(QQ, [[1, 0], [0, 1/g]])
    # We program the statement A = A * delta ad hoc:
    A[0,1].value //= g
    A[0,1].modulus //= g
    A[1,1].value //= g
    A[1,1].modulus //= g
    M = ~delta * M
    print(A); print()

    # Compute S = iota(det(A)) and U, the left-over with determinant 1
    S = matrix(Zhat, [[1, 0], [0, det(A)]])
    det_inv = det(A).value.inverse_mod(det(A).modulus)
    S_inv = matrix(Zhat, [[1, 0], [0, det_inv]])
    U = S_inv * A

    return S, U, M
