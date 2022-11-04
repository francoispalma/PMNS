from convert import rho_div_convert_to_mns as conv, montgomery_convert_to_mns
from pyparams128 import p, n, gamma, lam, phi, rho, M_or_B, M1_or_B1
from ops import horner_modulo
from random import randrange
from ast import literal_eval

if __name__ == "__main__":
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
				a = [elem - phi if elem > (phi >> 1) else elem for elem in a]
				b = literal_eval(f.readline()[:-1])
				b = [elem - phi if elem > (phi >> 1) else elem for elem in b]
				c = literal_eval(f.readline()[:-1])
				c = [elem - phi if elem > (phi >> 1) else elem for elem in c]
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

