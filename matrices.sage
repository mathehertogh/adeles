#from sage.modular.arithgroup.arithgroup_perm import eval_sl2z_word, sl2z_word_problem
#A = SL2Z([7,8,-50,-57])
#factorization = sl2z_word_problem(A)
#print(A == eval_sl2z_word(factorization))

load("adeles.sage")

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

Qhat = ProfiniteNumbers(QQ)
A = MatrixSpace(Qhat, 2).matrix([Qhat(3, 240, 2), Qhat(3, 240, 40),
	           					 Qhat(1, 240, 3), Qhat(7, 240, 2)])
B, C = split(A)

print(A)
print("=")
print(B)
print("*")
print(C)

print(A == B*C)