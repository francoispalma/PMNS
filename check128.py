from convert import rho_div_convert_to_mns as conv, montgomery_convert_to_mns
from pyparams128 import p, n, gamma, lam, phi, rho, M_or_B, M1_or_B1
from ops import list_to_poly, horner_modulo
from random import randrange
from ast import literal_eval

#RHO = 118
#rho = 2**RHO

#phinmoinsun = pow(phi, n - 1, p)
#Pi = [montgomery_convert_to_mns((rho**i) * (phi**2), p, n, gamma, rho, lam, phi, M, M1, phinmoinsun) for i in range(n)]

#a = 0x74ff560400d0105e6381e4f7cf22ba4a3d949bbe3b03e7ec1c8aebfb02a4dedf230eef099cd1ae78adf8f142cd70ed93122a5c48c5edcba658615fa2316994dce0c84e9e54c5ae9482acdc0ed6fae84eb7e83d94016d12452ad41369e33a53a676d539439488bdc8b3462c5579a432e8b579e8af9d5b2b0b8f37856fe2de7f30
#a = 0xbf6dc9f34905d4ccea18b34313d7f22412795efa0161f7ebcd5912a900ea7d255661bb894729e4fd85a477d3c575f3e97fcd1e6e2fd01d5317724f38def3c7f944162bb4ae4dcd5b1522efca1f3713a927c91f1113096ced7585edf7fef8cc9334dc56e8483a3c49f4a0fb9bb73c00b8b00e3d11435184eacbd45dd38fcbcadd
#a = 0x7025ec872fe66fc290f7a51d37cbebe888d1bf7e7cbc326cd9f54a80fac9b810ba8c44ba24ab6601f39af4248c1fa36f77be47c24f7db35db70decfc7e8d7601d3e9c2da091f85170babe8499c42dc51933310e8234eed7a07d5f8d4acaef2bf3a2c69b664c4b266bca70292ebd8c2cddede907d128efa78253c2514cb02841a
#print("HAAAAAAAAAA")
#tmp = conv(a, p, n, gamma, rho, lam, phi, M, M1, Pi)

#print([hex(elem) if elem >= 0 else hex(2**128 + elem) for elem in tmp])

#print()

#tmp = horner_modulo(tmp, gamma, p)
#print(tmp)

#print()
#print(tmp == a * phi % p)
##print([hex(elem) if elem >= 0 else hex(2**128 + elem) for elem in M1])

if __name__ == "__main__":
	rho = 2**rho
	counter = 0
	visual = 0
	print("Starting")
	with open("log", "r") as f:
		while True:
			#print("\b", end="")
			print("\b" + str(visual), end="\r")
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
#					print("a:", [hex(elem) for elem in a])
#					print("b:", [hex(elem) for elem in b])
#					print("c:", [hex(elem) for elem in c])
#					print("cc:", [hex(elem) for elem in c_check])
#					print()
#					exit()
				for elem in c:
					if abs(elem) >= rho:
						print("More than Rho")
			except SyntaxError:
				break
	print("Finished")
	print("counter:", counter)

