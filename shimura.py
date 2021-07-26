
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.arith.functions import lcm
from sage.matrix.constructor import matrix
from profinite_number import Qhat
from adele import Adeles
from matrix import modulus, denominator, factor_GLQhat#, value_matrix, denominator_matrix, 


def shimura_connecting_homomorphism(u, output_prec=None):
	"""
	Evaluate Shimura's connecting homomorphism at the idele ``u``
	"""
	if u.has_exact_finite_part():
		raise NotImplementedError()

	# First we convert u to a vector of profinite QQ-numbers.
	K = u.parent().number_field()
	a = Adeles(K)(u)
	t, s = a.finite().to_profinite_rational_vector()

	# Next we compute the entries of the transpose of the matrix
	# representing the multiplication by x map.
	theta = K.gen()
	C, B, _ = theta.minpoly().list()
	E = matrix(Qhat, [[t-B*s, -C*s], [s, t]])

	# Lastly we compute the valuations of the determinant of the exact
	# matrix as a rational number Delta.
	Delta = QQ(1)
	for P in u.stored_primes():
		Delta *= P.norm()**u.valuation(P)

	if output_prec is None:
		return E, Delta

	# We make a copy of u as to not change u itself.
	v = u.copy() # TODO check if this makes sure u does not change! TODO is does not!
	# We need to avoid that type(output_prec) is int
	output_prec = ZZ(output_prec)

	while modulus(E) / output_prec not in ZZ:
		for p in (modulus(E) / output_prec).denominator().prime_divisors():
			increment = output_prec.valuation(p) - modulus(E).valuation(p)
			v.increase_precision(p, increment)
		E, Delta = shimura_connecting_homomorphism(v)

	return v, (E, Delta)

def factored_shimura_connecting_homomorphism(u, output_prec):
	"""
	Perform Algorithm 9.2 of [Her2021]
	"""
	# Step 1.
	u_1, (E_1, Delta_1) = shimura_connecting_homomorphism(u, 1)
	# Step 2.
	P_1 = Delta_1.numerator() * denominator(E_1)
	u_2, (E_2, Delta_2) = shimura_connecting_homomorphism(u_1, P_1)
	# Step 3.
	A_0 = factor_GLQhat(E_2, Delta_2)
	# Step 4.
	P_2 = lcm(P_1, output_prec*A_0.inverse().denominator())
	u_3, (E_3, Delta_3) = shimura_connecting_homomorphism(u_2, P_2)
	# Step 5.
	A = factor_GLQhat(E_3, Delta_3)
	# Step 6.
	B = E_3 * A.inverse()
	return B, A
