import sys, os, cypari2

from math import floor, ceil, log2
from sage.all import matrix, ZZ, previous_prime, next_prime, factor, next_probable_prime, is_prime, pari
from numpy import count_nonzero


PSIZE = 256
n = 5
PHI = 2**64
s = 1
#253



if len(sys.argv) > 2:
	PSIZE = int(sys.argv[1])
	n = int(sys.argv[2])
if len(sys.argv) > 3:
	s = int(sys.argv[3])

powert = 0
while not((s//(2**powert)) % 2):
	powert += 1

print("powert", powert)

print("finding p of size", PSIZE, "bits with n =", n, "such that s*p = gamma^n - lambda with s =", s)
print()

NBZ = 0.5
NN = float(n*n)

lowgamma = int(ceil(2**((log2(s) + PSIZE-1)/float(n))))
while lowgamma**n > s*2**(PSIZE-1):
	lowgamma -= 1

highgamma = int(ceil(2**((log2(s) + PSIZE)/float(n))))
while highgamma**n > s*2**(PSIZE):
	highgamma -= 1


print(highgamma - lowgamma)
print("2^" + str(floor(log2(highgamma - lowgamma))))

w = lambda lam : lam*(n-1) + 1

#w = lambda alpha, beta : max([alpha*(n-i) + i*beta for i in range(n)])


delta = lambda gamma : ((((gamma)*(n-1))**2 + gamma*(6 - 4*n) + PHI*(n-1) + 1)**0.5 - (gamma*(n-1) + 1))/(2*n - 2)


lam = int(delta(lowgamma))

#########################################################################
# w*rho**2 + n1B*PHI < rho*PHI                                          #
# (w*rho**2)/(rho-n1B/2) < PHI                                          #
# min with rho = n1B + 1                                                #
#########################################################################

# w*(n1B)² + n1B*PHI/2 < (n1B+1)*PHI
# 2*w*(n1B)² + n1B*PHI < 2*n1B*PHI + 2*PHI
# 2*w*(n1B)² < PHI*(n1B + 2)
# 2wn1B² < PHI*(n1B + 2)
# 2wn1B² < PHI*(n1B + 2)

#while 4*w(lam)*(lowgamma + lam) < PHI:
#	lam += 1
#while 4*w(lam)*(lowgamma + lam) > PHI:
#	lam -= 1

n1B = lambda lam: lowgamma + lam
while 2*w(lam)*(n1B(lam)**2) < (n1B(lam)+2)*PHI:
	lam += 1
while 2*w(lam)*(n1B(lam)**2) > (n1B(lam)+2)*PHI:
	lam -= 1

print(lam)

gamma = lowgamma
Minlab = 2**64
count = 0
valid = 0
hollow = 0
vgm = 0
hollows = []
#gamma = lowgamma + 11041120
#count = 106844272
#valid = 95803152
#hollow = 13
#vgm = 11039636
#hollows = ['2^51']
n1B = lambda lam: gamma + lam
printstr = "\b{gamma - lowgamma}, vgm: {vgm}, count: {count}, valid: {valid}, {Minlab}, {round(float(vgm)/(gamma - lowgamma + 1), 4)}, {round(float(valid)/(vgm + (vgm == 0)), 4)}"
try:
	while gamma <= highgamma:
		fcfg = False
		#print(eval(f"f'{printstr}'"), end="\r")
		maxlam = int(delta(gamma))
		while 2*w(maxlam)*(n1B(maxlam)**2) < (n1B(maxlam)+2)*PHI:
			maxlam += 1
		while 2*w(maxlam)*(n1B(maxlam)**2) < (n1B(maxlam)+2)*PHI:
			maxlam -= 1
		if maxlam < 1:
			break
		gamman = (gamma**n)//s
		p = next_probable_prime(gamman - maxlam - 1)
		while True:
			#print(eval(f"f'{printstr}'"), end="\r")
			count += 1
			lam = gamma**n - s*p
			if abs(lam) < Minlab:
				Minlab = abs(lam)
			if not(lam % (2**powert)) and not(gamma % (2**powert)) and (
					2*w(abs(lam))*(n1B(abs(lam >> powert))**2) < (n1B(abs(lam >> powert))+2)*PHI):
#			if 4*w(abs(lam))*(gamma + abs(lam)) < PHI:
				try:
					if (True or is_prime(p)):
						fcfg = True
						valid += 1
#						B = matrix(ZZ, [[gamma if i == j == n-1 else -gamma if i == j else 1 if j == i + 1 else lam if (i == n-1) and (j == 0) else 0 for j in range(n)] for i in range(n)])
#						cz = count_nonzero(B.inverse() % PHI)
#						hollow += cz / NN < NBZ
#						if cz / NN < NBZ:
#							hollows += [factor(gamma)]
				except cypari2.handle_error.PariError:
					#pari.allocatemem()
					print(f"gamm {gamma - lowgamma}")
					print()
					print(eval(f"f'{printstr}'"), end="\r")
					exit(1)
				except ZeroDivisionError:
					valid -= 1
					fcfg = False
			elif p > gamman:
				break
			p = next_probable_prime(p)
		vgm += int(fcfg)
		print(eval(f"f'{printstr}'"), end="\r")
		gamma += 1
except KeyboardInterrupt:
	print(eval(f"f'{printstr[1:]}'"))
	print("hollows:")
	print(set([str(elem) for elem in hollows]))
	try:
		sys.exit(130)
	except SystemExit:
		os._exit(130)
