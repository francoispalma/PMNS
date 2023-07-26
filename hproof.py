from ast import literal_eval

from convert import rho_div_convert_to_mns as conv, montgomery_convert_to_mns
from ops import horner_modulo, amns_montg_mult
#from generated.pmns25664n5 import pmnsdict
#from generated.pmns51264n9 import pmnsdict
#primes = list(pmnsdict.keys())
#p, n, gamma, lam, rho, M, M1 = pmnsdict[primes[0]]
from hprecalcs import p, n, gamma, lam, rho, B as M, B1 as M1, phi
#phi = 2**64

if __name__ == "__main__":
	morethanrho = 0
	M = M + [0] * (n - len(M))
	rho = 2**rho
	counter = 0
	visual = 0
	print("Starting")
	with open("hlog", "r") as f:
		while True:
			print("\btested:" + str(visual) + "\terrors: " + str(counter), end="\r")
			visual += 1
			try:
				a = literal_eval(f.readline()[:-1])
				b = literal_eval(f.readline()[:-1])
				c = literal_eval(f.readline()[:-1])
				c_check = (horner_modulo(a, gamma, p) * horner_modulo(b, gamma, p) * pow(phi, -1, p)) % p
				if type(lam) == dict:
					c_check = c_check * lam["a"] % p
				if horner_modulo(c, gamma, p) != c_check:
					counter += 1
				for elem in c:
					if abs(elem) >= rho:
						morethanrho += 1
#						print("More than Rho:", abs(elem), rho)
#						exit()
			except SyntaxError:
				break
			except:
				break
	print("")
	print("more than rho: ", morethanrho)
