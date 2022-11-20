
p = 37
Fp = GF(p)
Fpx.<x> = PolynomialRing(Fp)
FpXY.<X, Y> = PolynomialRing(Fp)

a = Fp(2)
b = Fp(1)
E = EllipticCurve(Fp, [a, b])
F = X^3 + a*X + b

f3 = 3*x^4 + 6*a*x^2 + 12*b*x - a^2
f4 = 2 * (x^6 + 5*a*x^4 + 20*b*x^3 - 5*a^2*x^2 - 4*a*b*x - 8*b^2 - a^3)
divpolys = {0: Fpx(0), 1: Fpx(1), 2: Fpx(1), 3: f3, 4: f4}

# f is in FpXY but univariate
def conv_FpXY_FpX(f):
	li = list(f)
	ret = 0
	for coeff, base in li:
		ret += coeff * x^base.degree()
	return ret

# f mod modf
def mod_FpXY(f, modf):
	ret = (conv_FpXY_FpX(f) % conv_FpXY_FpX(modf))(X)
	return ret

class ModInvError(Exception):
	def __init__(self, gcd):
		self.gcd = gcd
		pass
	
	def get_gcd(self):
		return self.gcd

def double(P, modf):
	x1, y1 = P[0], (P[1] // Y)
	if (x1 == 0) and (y1 == 0):
		return (FpXY(0), FpXY(0))

	# inverse check
	r_denom = conv_FpXY_FpX(2*y1*F)
	modf = conv_FpXY_FpX(modf)
	g = gcd(r_denom, modf)
	# print("gcd:", g)
	if g != 1:
		raise ModInvError((g)(X))

	r = (3*x1^2 + a) * r_denom.inverse_mod(modf)(X)
	x3 = mod_FpXY(r^2 * F - 2*x1, modf(X))
	y3 = mod_FpXY(r*(x1 - x3) - y1, modf(X))
	return (x3, (y3 * Y))

def add(P, Q, modf):
	x1, y1, x2, y2 = P[0], (P[1]//Y), Q[0], (Q[1]//Y)
	if (x1 == x2) and (y1 == y2):
		return double(P, modf)
	if (x1 == x2) and (y1 == -y2):
		return (FpXY(0), FpXY(0))
	
	# inverse check
	# print(x1, x2)
	r_denom = conv_FpXY_FpX(x1 - x2)
	modf = conv_FpXY_FpX(modf)
	g = gcd(r_denom, modf)
	# print("gcd:", g)
	if g != 1:
		raise ModInvError((g)(X))

	r = (y1 - y2) * r_denom.inverse_mod(modf)(X)
	x3 = mod_FpXY(r^2 * F - x1 - x2, modf(X))
	y3 = mod_FpXY(r*(x1 - x3) - y1, modf(X))
	return (x3, (y3 * Y))

def nP(P, n, modf):
	if n == 0:
		return (FpXY(0), FpXY(0))
	ret = P
	for i in range(n - 1):
		ret = add(ret, P, modf)
	return ret

def divpoly(m):
	if m in divpolys:
		return divpolys[m]
	
	F = 4 * (x^3 + a*x + b)
	mm = m // 2
	if not (m & 1):
		f = (divpoly(mm + 2) * divpoly(mm - 1)^2 - divpoly(mm - 2) * divpoly(mm + 1)^2) * divpoly(mm)
		divpolys[m] = f
		return f

	if (mm & 1):
		f = divpoly(mm + 2) * divpoly(mm)^3 - F^2 * divpoly(mm - 1) * divpoly(mm + 1)^3
		divpolys[m] = f
		return f

	f = F^2 * divpoly(mm + 2) * divpoly(mm)^3 - divpoly(mm - 1) * divpoly(mm + 1)^3
	divpolys[m] = f
	return f


def is_eq(f, g, modf):
	ff = mod_FpXY(f.numerator() * g.denominator(), modf)
	gg = mod_FpXY(f.denominator() * g.numerator(), modf)
	# print("ff:",ff)
	# print("gg:",gg)
	if ff == gg:
		return True
	return False

def schoof():
	M = 1
	l = 3
	A = []
	pp = floor(4 * sqrt(p))
	while M < pp:
		pl = (p % l)
		phi_l = divpoly(l)(X)
		
		while True:
			try:
				P1 = (mod_FpXY(X^p, phi_l), Y * mod_FpXY(F^((p-1)//2), phi_l))			# (x^p, y^p)
				P2 = (mod_FpXY(X^(p^2), phi_l), Y * mod_FpXY(F^((p^2-1)//2), phi_l))		# (x^(p^2), y^(p^2))
				P3 = (mod_FpXY(X, phi_l), Y)												# (x, y)
				P5 = nP(P3, pl, phi_l)		# [pl](x, y)
				P6 = add(P2, P5, phi_l)		# (x^p^2, y^p^2) + [pl](x, y)

				for n in range(l):
					# print("[n] = [{}]".format(n))

					P4 = nP(P1, n, phi_l) 		# [n](x^p, y^p)

					if is_eq(P4[0], P6[0], phi_l) and is_eq(P4[1], P6[1], phi_l):
						tl = n
						# print("tl:", tl)
			except ModInvError as e:
				phi_l = e.get_gcd()
			else:
				break
		A.append((tl, l))
		M *= l
		l = next_prime(l)
		# print("")
	return A, M

A, M = schoof()
print(A)

A = list(zip(*A))
print(A)
t = crt(list(A[0]), list(A[1]))
if t > M // 2:
	t -= M

ord = p + 1 - t
print(ord)
print(E.order())
