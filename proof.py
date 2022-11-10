from ast import literal_eval

from convert import rho_div_convert_to_mns as conv, montgomery_convert_to_mns
from ops import horner_modulo, amns_montg_mult
from pyparams import p, n, gamma, lam, phi, rho, M_or_B as M, M1_or_B1 as M1

if __name__ == "__main__":
	M = M + [0] * (n - len(M))
	rho = 2**rho
	counter = 0
	visual = 0
	print("Starting")
	with open("log", "r") as f:
		while True:
			print("\btested:" + str(visual) + "\terrors: " + str(counter), end="\r")
			visual += 1
			try:
				a = literal_eval(f.readline()[:-1])
				b = literal_eval(f.readline()[:-1])
				c = literal_eval(f.readline()[:-1])
				c_check = (horner_modulo(a, gamma, p) * horner_modulo(b, gamma, p) * pow(phi, -1, p)) % p
				if horner_modulo(c, gamma, p) != c_check:
					counter += 1
				for elem in c:
					if abs(elem) >= rho:
						print("More than Rho:", abs(elem), rho)
						exit()
			except SyntaxError:
				break
	print("")
