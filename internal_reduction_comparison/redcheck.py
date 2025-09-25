import os

from ast import literal_eval
from codegen import p, n, gamma, lam, rho

NOMONTG, NOBARRETT = False, False
if "NOMONTG" in os.environ:
	NOMONTG = True
	NOBARRETT = True

if "NOBARRETT" in os.environ:
	NOBARRETT = True

def pmns_conv_to_binary(Poly: list, Gi: list, modulo: int) -> int:
	sum_ = Poly[0]
	for i in range(1,len(Poly)):
		sum_ += Poly[i]*Gi[i-1]
	return int(sum_ % modulo)

if __name__ == "__main__":
	Gi = [pow(gamma,i,p) for i in range(1,n)]
	morethanrho = 0
	counter = 0
	visual = 0
	lastmorethanrho = 0
	lastinp = 0
	phi = 2**64
	if not NOBARRETT:
		with open("blog", "r") as f:
			lines = f.readlines()
		while lines:
			print("\btested:" + str(visual) + "\terrors: " + str(counter), end="\r")
			visual += 1
			a = literal_eval(lines.pop(0))
			b = literal_eval(lines.pop(0))
			c = literal_eval(lines.pop(0))
			c_check = pmns_conv_to_binary(a, Gi, p) * pmns_conv_to_binary(b, Gi, p) % p
			if pmns_conv_to_binary(c, Gi, p) != c_check:
				counter += 1
			toaddto = 0
			for elem in c:
				if abs(elem) >= rho:
					toaddto = 1
					lastmorethanrho = abs(elem)
					lastinp = max(lastinp, max([abs(elem) for elem in a] + [abs(elem) for elem in b]))
			morethanrho += toaddto
		print(f"Barrett-like\ttested: {visual}\terrors: {counter}\tmore than rho: {morethanrho}")
	morethanrho = 0
	counter = 0
	visual = 0
	lastmorethanrho = 0
	lastinp = 0
	phi = 2**64
	if not NOMONTG:
		with open("mlog", "r") as f:
			lines = f.readlines()
		while lines:
			print("\btested:" + str(visual) + "\terrors: " + str(counter), end="\r")
			visual += 1
			a = literal_eval(lines.pop(0))
			b = literal_eval(lines.pop(0))
			c = literal_eval(lines.pop(0))
			c_check = pmns_conv_to_binary(a, Gi, p) * pmns_conv_to_binary(b, Gi, p) * pow(phi, -1, p)  % p
			if pmns_conv_to_binary(c, Gi, p) != c_check:
				counter += 1
			toaddto = 0
			for elem in c:
				if abs(elem) >= rho:
					toaddto = 1
					lastmorethanrho = abs(elem)
					lastinp = max(lastinp, max([abs(elem) for elem in a] + [abs(elem) for elem in b]))
			morethanrho += toaddto
		print(f"Montgom-like\ttested: {visual}\terrors: {counter}\tmore than rho: {morethanrho}")
	morethanrho = 0
	counter = 0
	visual = 0
	lastmorethanrho = 0
	lastinp = 0
	phi = 2**64
	with open("plog", "r") as f:
		lines = f.readlines()
	while lines:
		print("\btested:" + str(visual) + "\terrors: " + str(counter), end="\r")
		visual += 1
		a = literal_eval(lines.pop(0))
		b = literal_eval(lines.pop(0))
		c = literal_eval(lines.pop(0))
		c_check = pmns_conv_to_binary(a, Gi, p) * pmns_conv_to_binary(b, Gi, p) * pow(phi, -2, p) % p
		if pmns_conv_to_binary(c, Gi, p) != c_check:
			counter += 1
		toaddto = 0
		for elem in c:
			if abs(elem) >= rho:
				toaddto = 1
				lastmorethanrho = abs(elem)
				lastinp = max(lastinp, max([abs(elem) for elem in a] + [abs(elem) for elem in b]))
		morethanrho += toaddto
	print(f"Plantard-like\ttested: {visual}\terrors: {counter}\tmore than rho: {morethanrho}")
#	print("more than rho: ", morethanrho)
