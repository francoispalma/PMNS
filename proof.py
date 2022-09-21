from math import ceil
from random import randrange
from ast import literal_eval

from convert import rho_div_convert_to_mns as conv, montgomery_convert_to_mns
from ops import list_to_poly, horner_modulo, amns_montg_mult
from pyparams import p, n, gamma, lam, phi, rho, M, M1

#amns = (p, n, gamma, rho, lam) = (6152896135288560374679945371974689688835168151742564408104565373600581564260451457, 5, 220855883097298041197912187592864814478435487109452369765200775161577472, 2305843009213693952, 2)
#phi = 2**64
#M = [9446094647596025, -89344859775265, 366378001529314, -4558175830143501, 19231042282234940]
#M1 = [7045631417041842631, -6084863496536136821, 8006399337547548431, 1601279867509509692, 4355481239625866353]

if __name__ == "__main__":
	M = M + [0] * (n - len(M))
	rho = 2**rho
	counter = 0
	visual = 0
	print("Starting")
	with open("log", "r") as f:
		while True:
			print("\b" + str(visual), end="\r")
			visual += 1
			try:
				a = literal_eval(f.readline()[:-1])
				a = [elem - phi if elem > (phi >> 1) else elem for elem in a]
				b = literal_eval(f.readline()[:-1])
				b = [elem - phi if elem > (phi >> 1) else elem for elem in b]
				c = literal_eval(f.readline()[:-1])
				c = [elem - phi if elem > (phi >> 1) else elem for elem in c]
				c_check = amns_montg_mult(a, b, p, n, gamma, rho, lam, phi, M, M1)
				if horner_modulo(c, gamma, p) != horner_modulo(c_check, gamma, p):
					counter += 1
					#print("False")
				for elem in c:
					if abs(elem) >= rho:
						print("More than Rho")
			except SyntaxError:
				break
	print("Finished")
	print("counter:", counter)
