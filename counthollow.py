import sys, os, cypari2

from math import floor, ceil, log2, sqrt
from sage.all import matrix, ZZ, previous_prime, next_prime, factor, next_probable_prime, is_prime, pari, is_square
from numpy import count_nonzero


PSIZE = 256
n = 5
PHI = 2**64
s = 1
alpha = 1
delta = 0
#253

WRITEOUT = False
HOLLOW = False

ARGS = sys.argv
if "write" in ARGS:
	WRITEOUT = True
	print("writing to file")
	ARGS.remove("write")
if "hollow" in ARGS:
	HOLLOW = True
	print("hollow only")
	ARGS.remove("hollow")
if "maxbeta" in ARGS:
	MAXBETA = int(ARGS.pop(ARGS.index("maxbeta") + 1))
	ARGS.remove("maxbeta")
else:
	MAXBETA = 2**64
if len(ARGS) > 2:
	PSIZE = int(ARGS[1])
	n = int(ARGS[2])
if len(ARGS) > 3:
	s = int(ARGS[3])
if len(ARGS) > 4:
	alpha = int(ARGS[4])


filename = f"generated/{'h'*int(HOLLOW)}pmns{PSIZE}{int(log2(PHI))}n{n}{'' if delta == 0 else 'delta' + str(delta)}{'' if MAXBETA == 2**64 else 'mb' + str(MAXBETA)}.py"
if WRITEOUT:
	print(filename)
	if not os.path.exists(filename):
		with open(filename, "w+") as FILE:
			FILE.write("pmnsdict = {}\n")

powert = 0
while not((s//(2**powert)) % 2):
	powert += 1

print("powert", powert)

print(f"finding p of size {PSIZE} bits with n = {n} such that {s}*p = {'' if alpha == 1 else str(alpha)+'*'}gamma^n - {'lambda' if alpha == 1 else 'beta'}")
print()

NBZ = 0.5
NN = float(n*n)

lowgamma = int(ceil(2**((log2(s) - log2(alpha) + PSIZE-1)/float(n))))
while alpha*lowgamma**n > s*2**(PSIZE-1):
	lowgamma -= 1

highgamma = int(ceil(2**((log2(s) - log2(alpha) + PSIZE)/float(n))))
while alpha*highgamma**n > s*2**(PSIZE):
	highgamma -= 1


print(highgamma - lowgamma)
print("2^" + str(floor(log2(highgamma - lowgamma))))

aleph = alpha

if alpha == 1:
	w = lambda lam : abs(lam)*(n-1) + 1
else:
	w = lambda beta : max([aleph*(n-i) + i*abs(beta) for i in range(n)])




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

gamma = lowgamma

if alpha == 1:
	n1B = lambda lam: gamma + lam
	beta = lambda gamma : ((((gamma)*(n-1))**2 + gamma*(6 - 4*n) + PHI*(n-1) + 1)**0.5 - (gamma*(n-1) + 1))/(2*n - 2)
else:
	n1B = lambda lam: alpha*gamma + 1
	beta = lambda gamma : int(PHI/((n//2)*(2*n1B(gamma)-4)) - alpha)

lam = int(beta(gamma))

while 2*w(lam)*(n1B(lam) - 2) < PHI:
	lam += 1
while 2*w(lam)*(n1B(lam) - 2) > PHI:
	lam -= 1
	if lam <= 0:
		break

print(lam)

init_lam = lam

if HOLLOW:
	gamma = ceil(gamma/2**32)*2**32
Minlab = 2**64
count = 0
valid = 0
hollow = 0
vgm = 0
mindelta, maxdelta, avdelta = 1000000,0,0
hollows = []
if alpha == 1:
	n1B = lambda lam: gamma + abs(lam)
else:
	n1B = lambda lam: aleph*gamma + 1 if aleph != 1 else gamma + abs(lam)
maxlam = 50
printstr = "\b{gamma - lowgamma}, vgm: {vgm}, count: {count}, valid: {valid}, {Minlab}, {round(float(vgm)/((gamma - lowgamma)/(1+HOLLOW*2**32) + 1), 4)}, {round(float(valid)/(vgm + (vgm == 0)), 4)}, {maxdelta}, {round(float(avdelta/(vgm + (vgm==0))),4)}"
try:
	while gamma <= highgamma:
		cdelta = 0
		fcfg = False
		#print(eval(f"f'{printstr}'"), end="\r")
		maxlam = init_lam
		while 2*w(maxlam)*(n1B(maxlam>>powert) - 2) < PHI:
			maxlam += 1
		while 2*w(maxlam)*(n1B(maxlam>>powert) - 2) > PHI:
			maxlam -= 1
			if maxlam < 0:
				print()
				print("endlam")
				break
		if maxlam < 0:
			break
		gamman = (alpha*gamma**n)//s
		p = next_probable_prime(gamman - maxlam - 1)
		while True:
			#print(eval(f"f'{printstr}'"), end="\r")
			count += 1
			lam = alpha*gamma**n - s*p
			if abs(lam) < Minlab:
				Minlab = abs(lam)
			if abs(lam) <= MAXBETA and not(lam % (2**powert)) and not(alpha*gamma % (2**powert)):
				if not(alpha % (2**powert)):
					aleph = alpha >> powert
					beth = lam >> powert
					psi = 0
				else:
					aleph = alpha
					beth = lam
					psi = powert
					while not(aleph % 2):
						aleph >>= 1
						beth >>= 1
						psi -= 1
				if 2*w(beth)*(n1B(beth >> psi) - 2) < PHI:
					currdel = floor(sqrt(PHI/(2*w(beth)*(n1B(beth >> psi) - 2)))) - 1
					if currdel > cdelta:
						cdelta = currdel
					try:
						if (True or is_prime(p)):
							fcfg = True
							valid += 1
#							if (not is_square(gamma)) and is_square(lam):
#								print()
#								print("gamma =",gamma, "lam=", lam)
#								print()
							if WRITEOUT:
								#B = matrix(ZZ, [[-alpha*gamma if i == j == (n-1) else -gamma if i == j else 1 if j == i + 1 else (lam) if (i == n-1) and (j == 0) else 0 for j in range(n)] for i in range(n)])
								B = matrix(ZZ, [[-((aleph*gamma)>>psi) if i == j == (n-1) else -gamma if i == j else 1 if j == i + 1 else (lam>>powert) if (i == n-1) and (j == 0) else 0 for j in range(n)] for i in range(n)])
								with open(filename, "a") as FILE:
									FILE.write(f"pmnsdict[{p}] = {[p, n, gamma, beth if aleph == 1 else {'a':aleph,'b':beth}, ceil(log2(n1B(beth) - 1)), list(B), list(B.inverse() % PHI)]}\n")
								if valid == 1000:
									raise KeyboardInterrupt
#							cz = count_nonzero(B.inverse() % PHI)
#							hollow += cz / NN < NBZ
#							if cz / NN < NBZ:
#								hollows += [factor(gamma)]
					except cypari2.handle_error.PariError:
						#pari.allocatemem()
						print(f"gamm {gamma - lowgamma}")
						print()
						print(eval(f"f'{printstr}'"), end="\r")
						exit(1)
					except ZeroDivisionError:
						print()
						print("oops")
						valid -= 1
						fcfg = False
			if p > gamman + maxlam:
				break
			p = next_probable_prime(p)
		vgm += int(fcfg)
		print(eval(f"f'{printstr}'"), end="\r")
		if fcfg:
			if cdelta < mindelta:
				mindelta = cdelta
			if cdelta > maxdelta:
				maxdelta = cdelta
			avdelta += cdelta
		gamma += 1
		if HOLLOW:
			gamma += 2**32 - 1
except KeyboardInterrupt:
	print(eval(f"f'{printstr[1:]}'"))
	print("hollows:")
	print(set([str(elem) for elem in hollows]))
	try:
		sys.exit(130)
	except SystemExit:
		os._exit(130)
