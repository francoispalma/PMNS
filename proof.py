from ast import literal_eval
from ops import horner_modulo, amns_montg_mult

amns = (p, n, gamma, rho, lam) = (6152896135288560374679945371974689688835168151742564408104565373600581564260451457, 5, 220855883097298041197912187592864814478435487109452369765200775161577472, 2305843009213693952, 2)
phi = 2**64
M = [9446094647596025, -89344859775265, 366378001529314, -4558175830143501, 19231042282234940]
M1 = [7045631417041842631, -6084863496536136821, 8006399337547548431, 1601279867509509692, 4355481239625866353]



if __name__ == "__main__":
	print("Starting")
	with open("log", "r") as f:
		while True:
			try:
				a = literal_eval(f.readline()[:-1])
				b = literal_eval(f.readline()[:-1])
				c = literal_eval(f.readline()[:-1])
				c_check = amns_montg_mult(a, b, *amns, phi, M, M1)
				if horner_modulo(c, gamma, p) != horner_modulo(c_check, gamma, p):
					print("False")
				for elem in c:
					if abs(elem) >= rho:
						print("More than Rho")
			except SyntaxError:
				break
	print("Finished")
