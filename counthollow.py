import sys, os, cypari2

from sage.all import matrix, ZZ, previous_prime, next_prime, factor, next_probable_prime, is_prime, pari
from numpy import count_nonzero


PSIZE = 256
n = 5
PHI = 2**64

if len(sys.argv) > 2:
	PSIZE = int(sys.argv[1])
	n = int(sys.argv[2])


NBZ = 0.5
NN = float(n*n)

lowgamma = int(2**((PSIZE-1)/float(n)))
while lowgamma**n > 2**(PSIZE-1):
	lowgamma -= 1

highgamma = int(2**(PSIZE/float(n)))
while highgamma**n > 2**(PSIZE):
	highgamma -= 1


w = lambda lam : lam*(n-1) + 1

delta = lambda gamma : ((((gamma)*(n-1))**2 + gamma*(6 - 4*n) + PHI*(n-1) + 1)**0.5 - (gamma*(n-1) + 1))/(2*n - 2)


lam = int(delta(lowgamma))

while 4*w(lam)*(lowgamma + lam) < PHI:
	lam += 1
while 4*w(lam)*(lowgamma + lam) > PHI:
	lam -= 1

print(lam)

gamma = lowgamma
count = 0
valid = 0
hollow = 0
hollows = []
try:
	while gamma <= highgamma:
		print(f"\b{gamma - lowgamma} count: {count}, valid: {valid}, hollow: {hollow}", end="\r")
		maxlam = int(delta(gamma))
		while 4*w(maxlam)*(gamma + maxlam) < PHI:
			maxlam += 1
		while 4*w(maxlam)*(gamma + maxlam) > PHI:
			maxlam -= 1
		gamman = gamma**n
		p = next_probable_prime(gamman - maxlam - 1)
		while True:
			print(f"\b{gamma - lowgamma} count: {count}, valid: {valid}, hollow: {hollow}", end="\r")
			count += 1
			lam = p - gamman
			if 4*w(abs(lam))*(gamma + abs(lam)) < PHI:
				try:
					if is_prime(p):
						valid += 1
						B = matrix(ZZ, [[gamma if i == j == n-1 else -gamma if i == j else 1 if j == i + 1 else lam if (i == n-1) and (j == 0) else 0 for j in range(n)] for i in range(n)])
						cz = count_nonzero(B.inverse() % PHI)
						hollow += cz / NN < NBZ
						if cz / NN < NBZ:
							hollows += [factor(gamma)]
					p = next_probable_prime(p)
				except cypari2.handle_error.PariError:
					#pari.allocatemem()
					print(f"gamm {gamma - lowgamma}")
					print()
					print(f"\b{gamma - lowgamma} count: {count}, valid: {valid}, hollow: {hollow}", end="\r")
					exit(1)
			elif p > gamman:
				break
		gamma += 1
		print(f"\b{gamma - lowgamma} count: {count}, valid: {valid}, hollow: {hollow}", end="\r")
except KeyboardInterrupt:
	print(f"count: {count}, valid: {valid}, hollow: {hollow}")
	print("hollows:")
	print(set([str(elem) for elem in hollows]))
	try:
		sys.exit(130)
	except SystemExit:
		os._exit(130)
