#!/usr/bin/env python3
from sage.all import matrix,ZZ,xgcd,vector,ceil

from genmppmns import findminrow

p,n,gamma,lam,rho = (90075959928965881812115109766878025624447511447676603127783616474856091308000148988333386884816606599955226522360073078552636799711181684701384803279210775865878622732460829024237522026454674316197961848642034880127947239466588020276871510535857518129559901020305397376976672385341055396383596670809791268511, 6, 86742661103592852020646276311351653153320012960072644055428490711861282748457088943269795453521326049066138959027531542035112371930820540643361740076029535781642502101159741140315633017870096551999476505620020051527465986240034247143341953599861005226419377584296588467915404211182412011109315601005117636269, 2, 6820182838396110062975394024340498657866418979727546)


NBCHUNKS = ceil(rho.bit_length()/64)


G = matrix(ZZ,[[p if (k, j) == (0, 0) else -pow(gamma, k, p) if k != 0 and j == 0 else 1 if k == j else 0 for j in range(n)] for k in range(n)]).LLL()

PHI = 2**64  # Not the real PHI, used to get M1
ZZX = ZZ["X"]
E = ZZX(f"X^{n} - {lam}")
outcouple = findminrow(G, lam)
minnorm = outcouple[1]
M = G[outcouple[0]]
M = [int(elem) for elem in M]
while not(M[0]%lam):  # We divide by X mod E as much as possible
	M = M[1:] + [M[0]//lam]
val, M1, soak = xgcd(ZZX(M), E)
M1 = (M1 * int(pow(int(val), -1, PHI) % PHI))
rho = int(minnorm - 1)
w = abs(lam)*(n-1) + 1
phi = 2*w*(rho-1)
philog2 = phi.bit_length()
rhover2 = philog2//NBCHUNKS
rhohi = rho//2**(rhover2*(NBCHUNKS-1))
#	while rhohi*2**NBCHUNKS < 2**rhover2:
#		rhover2 -= 1
#		rhohi = rho//2**(rhover2*(NBCHUNKS-1))
philog2 = rhover2*NBCHUNKS

if __name__ == "__main__":
	print(ceil(p.bit_length()))
	print(NBCHUNKS)
	print(philog2)
	G = matrix(ZZ, [list(ZZX(M)*ZZX(f"X^{i}")%E) for i in range(n)])
	with open("mpparams.h", "w+") as fpointer:
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
