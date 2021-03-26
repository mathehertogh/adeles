
DEBUG = False

def debug(msg):
	global DEBUG
	if DEBUG:
		print("    " + msg)

def random_number_field(degree=None):
	if degree is None:
		from sage.probability.probability_distribution import GeneralDiscreteDistribution
		P = [0.0, 0.05, 0.4, 0.25, 0.15, 0.05, 0.05, 0.025, 0.025]
		X = GeneralDiscreteDistribution(P)
		degree = X.get_random_element() # 1 <= degree <= 8
	R = PolynomialRing(ZZ, 'x')
	f = R.random_element(degree)
	while not f.is_irreducible():
		f = R.random_element(degree)
	K = NumberField(f.monic(), 'a')
	raise ValueError("bla")
	return K