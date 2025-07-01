from sage.all import matrix,ZZ,xgcd,vector,ceil,prod,randrange
import os,sys
from avx512tableparams import pmnsdict


(p,n,gamma,lam,rho) = \
(113613660995108872950371506809752759690675959588751101844436533476808849634914817320631434312607403433285054417605298675248475014833114028427087425993434174182179530496401742168757932836532815342701598918429494568113346870179029590042313640659442462109916351143761861048273853218848058804709166262931087236479, 8, 28955883242582741461429713172747058474375092655532786572293269212594960972241842848949055437872923792036814087400945860615405252662305594711847454248056130253102845402302538549688787308960903140730803393040444076956884291054537928373745712154313835458896352134104675274709543803037397813420422100486108350831, 2, 1051167395993162148350165608546359597254)
s = 3
delta = 0

try:
	psize = int(os.environ["PSIZE"])
	s = int(os.environ["NBCHUNKS"])
	(p,n,gamma,lam,rho) = pmnsdict[psize][s]
except:
	try:
		psize = int(os.environ["PSIZE"])
		(p,n,gamma,lam,rho) = pmnsdict[psize]
		if psize != 1536:
			delta = 5
			if psize == 807:
				s = 3
			elif psize >= 2436:
				s = 3
			else:
				s = 2
		else:
			s = 5
	except:
		pass

NBCHUNKS = s

def findminrow(mat, lam):
	norm1 = lambda V: sum([abs(elem) for elem in V])
	compnorm = lambda X : norm1(vector(list(X)[1:]))*lam + abs(X[0])
	minnorm = p
	minM = []
	mindex = -1
	i = 0
	E = ZZ["X"](f"X^{n} - {lam}")
	for lig in mat:
		tmp = list(lig)
		while True:
			M = ZZ["X"](tmp)
			val, M1, soak = xgcd(M, E)
			try:
				val = ZZ(val)
				if val & 1:  # if val is even it won't have a gcd of 1 with phi
					norm = compnorm(tmp)
					if norm < minnorm:
						minnorm = norm
						mindex = i
						minM = tmp
			except TypeError:
				pass
			if tmp[0]%lam:
				break
			tmp = tmp[1:] + [tmp[0]//lam]
		i += 1
	return mindex, minnorm, minM

G = matrix(ZZ,[[p if (k, j) == (0, 0) else -pow(gamma, k, p) if k != 0 and j == 0 else 1 if k == j else 0 for j in range(n)] for k in range(n)]).LLL()

PHI = 2**64
ZZX = ZZ["X"]
E = ZZX(f"X^{n} - {lam}")
rho = 2*max([sum([abs(G[j][i]) for j in range(n)]) for i in range(n)])
w = abs(lam)*(n-1) + 1
phi = 2*w*rho
allpos = lambda vec: bool(prod([elem >= 0 for elem in vec]))
gha = max([max([abs(elem) for elem in lig]) for lig in G])
v = vector([2*gha]*n)
coeffs = v*G.inverse()
coeffs = vector([round(elem) for elem in coeffs])
if allpos(coeffs*G):
	M = [int(elem) for elem in coeffs * G]
else:
	G1 = G.inverse() % phi
	v = [rho//2] * n
	phi_rep = vector([int(pow(phi, n+1, p))] + [0] * (n-1))
	slashphi = lambda vec: vector([elem//phi for elem in vec])
	for _ in range(n):
		phi_rep = slashphi(phi_rep - (phi_rep*G1 % phi) * G)
	C = vector(list(ZZX(v) * ZZX(list(phi_rep)) % E))
	v2 = slashphi(C - (C*G1 % phi) * G)
	newzero = vector(v) - v2
	coeffs = newzero * G.inverse()
	gradd = lambda vec, i: vector(list(vec[:i]) + [vec[i] - (vec[i]//(abs(vec[i]) + (vec[i] == 0)))] + list(vec[i+1:]))
	norminf = lambda vec: max([abs(elem) for elem in vec])
	for i in range(n):
		if coeffs[i] == 0:
			continue
	while allpos(gradd(coeffs, i)*G) and norminf(gradd(coeffs, i)*G) < norminf(coeffs*G):
		coeffs = gradd(coeffs, i)
	M = [int(elem) for elem in coeffs * G]
while M[0] % lam == 0:
	M = M[1:] + [M[0]//lam]
norm1 = lambda V: sum([abs(elem) for elem in V])
compnorm = lambda X : norm1(vector(list(X)[1:]))*lam + abs(X[0])
minnorm = compnorm(M)
val, M1, soak = xgcd(ZZX(M), E)
M1 = (M1 * int(pow(int(val), -1, PHI) % PHI))
rho = 2*int(minnorm)
phi = 2*w*(rho-1)
philog2 = phi.bit_length()
rhover2 = philog2//NBCHUNKS
rhohi = rho//2**(rhover2*(NBCHUNKS-1))
while rhohi*2**NBCHUNKS < 2**rhover2:
	rhover2 -= 1
	rhohi = rho//2**(rhover2*(NBCHUNKS-1))
else:
	rhover2 += 1
rhover2 += (rhover2 != 49)
rhohi = rho//2**(rhover2*(NBCHUNKS-1))
philog2 = rhover2*NBCHUNKS

if __name__ == "__main__":
	with open(f"avx512params.h", "w+") as fpointer:
		write = lambda X : fpointer.write(str(X) + "\n")
		write(f"#define N {n}")
		write(f"#define NBCHUNKS {NBCHUNKS}")
		write(f"#define LAMBDA {lam}")
		write(f"#define RHOver2 {rhover2}")
		write(f"\nstatic const uint64_t RHOhi = {rho//2**(rhover2*(NBCHUNKS-1) + 1)}u;")
		write(f"\nstatic const int64_t M[{NBCHUNKS}][{2*n-1}] = {{ {{")
		for i in range(NBCHUNKS-1):
			write("\t\t" + str([(elem//2**(rhover2*i))%2**rhover2 for elem in (list(vector(M[1:])*lam) + M)])[1:-1])
			write("\t}, {")
		write(str([elem//2**(rhover2*(NBCHUNKS-1)) for elem in (list(vector(M[1:])*lam) + M)])[1:-1])
		write("} };")
		PHI = 2**rhover2
		# We switch M1 from M^-1 to -M^-1
		M1 = [(int(M1[i]) * -1) % PHI for i in range(n)]
		# We then reduce M1's coefficients
		M1 = [int(M1[i]) - PHI if M1[i] >= (PHI >> 1) else int(M1[i]) for i in range(n)]
		write("\nstatic const int64_t M1[] = {")
		write(str([elem%2**rhover2 for elem in (list(vector(M1[1:])*lam) + M1)])[1:-1])
		write("};")
