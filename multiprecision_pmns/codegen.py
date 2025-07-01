#!/usr/bin/env python3
from sage.all import matrix,ZZ,xgcd,vector,ceil,prod

import os
from tableparams import pmnsdict


(p,n,gamma,lam,rho) = \
(134626326228911203302801776475443925678278255489633760896545287668303549829999313234273626117617897109776451956990190028384380667087135867956296037259427836729824065539559144053059735942997985604262768148777255323665661343041410533648843517479537384517257085145181031281905686985250491131741444027028357322183, 10, 6915070200122299431055015290546669602070576953476144957158240665058388415194596582918544415182941848060012171596172389404252354155639888469872144819773575354824721488875064345698916456320475826225085310480517765744300349926263982967354032927816587586055669373517544941504064641636165234776642356558467969924, 2, 31468319767503992344012480815574)

s = 2


if "PSIZE" in os.environ and "NBCHUNKS" in os.environ:
	PSIZE = int(os.environ["PSIZE"])
	s = int(os.environ["NBCHUNKS"])
	(p,n,gamma,lam,rho) = pmnsdict[PSIZE][s]

splitdico = {80:[2,5],60:[2,2,3],48:[2,3],40:[2,2,2],36:[2,3],32:[2,2,2],27:[3,3],24:[2,3],20:[2,2],18:[3],16:[2,2],12:[3],9:[3],6:[2]}

if n in splitdico:  # specific n splitting
	TOEPSPLIT = splitdico[n]
else:  # general guidelines (heuristic)
	TOEPSPLIT = []
	tmpn = n
	if n % 2 == 0:
		TOEPSPLIT += [2]
		tmpn //= 2
	while tmpn >= 10:
		if tmpn % 5 == 0:
			TOEPSPLIT += [5]
			tmpn //= 5
		elif tmpn % 2 == 0:
			TOEPSPLIT += [2]
			tmpn //= 2
		elif tmpn % 3 == 0:
			TOEPSPLIT += [3]
			tmpn //= 3
		else:
			break


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
shpef = findminrow(G, lam)
minnorm = shpef[1]
M = shpef[2]
val, M1, soak = xgcd(ZZX(M), E)
M1 = (M1 * int(pow(int(val), -1, PHI) % PHI))
rho = int(minnorm - 1)
w = abs(lam)*(n-1) + 1
phi = 2*w*(rho-1)
philog2 = phi.bit_length()
rhover2 = philog2//NBCHUNKS
rhohi = rho//2**(rhover2*(NBCHUNKS-1))
while rhohi*2**NBCHUNKS < 2**rhover2:
	rhover2 -= 1
	rhohi = rho//2**(rhover2*(NBCHUNKS-1))
else:
	rhover2 += 1
rhohi = rho//2**(rhover2*(NBCHUNKS-1))
philog2 = rhover2*NBCHUNKS

if __name__ == "__main__":
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
		if ceil(p.bit_length()) >= 4095 and TOEPSPLIT[0] == 2:
			M = list(vector(M[1:])*lam) + M
			Mchunks = [[int((M[j]//(2**(rhover2*i))) % 2**rhover2) for j in range(2*n-1)] for i in range(NBCHUNKS-1)]
			Mchunks += [[M[j]//(2**(rhover2*(NBCHUNKS-1))) for j in range(2*n-1)]]
			for i in range(NBCHUNKS):
				write(f"\nstatic const int64_t M{i}m0m1[{n-1}] = {{ {str(vector(Mchunks[i][n//2:n//2 + n - 1]) - vector(Mchunks[i][n:n + n - 1]))[1:-1]} }};")
				write(f"\nstatic const int64_t M{i}m0m2[{n-1}] = {{ {str(vector(Mchunks[i][n//2:n//2 + n - 1]) - vector(Mchunks[i][:n - 1]))[1:-1]} }};")
			M1 = [elem%2**rhover2 for elem in (list(vector(M1[1:])*lam) + M1)]
			M1 = [elem - 2**rhover2*(elem>2**(rhover2-1)) for elem in M1]
			write(f"\nstatic const int64_t M1_m0m1[{n-1}] = {{ {str(vector(M1[n//2:n//2 + n - 1]) - vector(M1[n:n + n - 1]))[1:-1]} }};")
			write(f"\nstatic const int64_t M1_m0m2[{n-1}] = {{ {str(vector(M1[n//2:n//2 + n - 1]) - vector(M1[:n - 1]))[1:-1]} }};")
			
		if TOEPSPLIT:
			write(f"""
static inline void mtoep{10-len(TOEPSPLIT)+1}(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	SCHOOLBOOK(N/{prod(TOEPSPLIT)})

static inline void m1toep{10-len(TOEPSPLIT)+1}(int64_t* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	M1SCHOOLBOOK(N/{prod(TOEPSPLIT)})
""")
			for i in range(len(TOEPSPLIT)-1, 0, -1):
				write(f"""
static inline void mtoep{10-i+1}(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	TOEP{TOEPSPLIT[i]}{TOEPSPLIT[i]}TOP(N/{prod(TOEPSPLIT[:i])}, mtoep{10-i})

static inline void m1toep{10-i+1}(int64_t* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	M1TOEP{TOEPSPLIT[i]}{TOEPSPLIT[i]}TOP(N/{prod(TOEPSPLIT[:i])}, m1toep{10-i})
""")
			
			write(f"""
static inline void toeplitz_vm(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	TOEP{TOEPSPLIT[0]}{TOEPSPLIT[0]}TOP(N, mtoep10)

static inline void ptoeplitz_vm(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	pTOEP{TOEPSPLIT[0]}{TOEPSPLIT[0]}TOP(N, mtoep10)

static inline void m1toep20(int64_t* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	M1TOEP{TOEPSPLIT[0]}{TOEPSPLIT[0]}TOP(N, m1toep10)
""")
		else:
			write(f"""
static inline void toeplitz_vm(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	SCHOOLBOOK(N)

static inline void ptoeplitz_vm(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	pSCHOOLBOOK(N)

static inline void m1toep20(int64_t* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	M1SCHOOLBOOK(N)
""")
		
		if ceil(p.bit_length()) >= 4095 and TOEPSPLIT[0] == 2:
			for i in range(NBCHUNKS):
				write(f"""
static inline void mtoeplitz_vm{i}(__int128* restrict rop, const int64_t* restrict vect)
{{
	__int128 t0[N/2], t1[N/2], t2[N/2];
	int64_t v0p1[N/2];
	for(int i = 0; i < N/2; i++)
	{{
		v0p1[i] = vect[i] + vect[i + N/2];
	}}
	mtoep10(t0, v0p1, M[{i}] + N/2);
	mtoep10(t1, vect, M{i}m0m1);
	mtoep10(t2, vect + N/2, M{i}m0m2);
	for(int i = 0; i < N/2; i++)
	{{
		rop[i] = t0[i] - t2[i];
		rop[i + N/2] = t0[i] - t1[i];
	}}
}}""")
				write("")
			
			write(f"""
// Utility function to extract the lower 64 bits of each vector coefficient
static inline void m1toeplitz_vm(int64_t* restrict rop, const __int128* restrict vect)
{{
	int64_t t0[N/2], t1[N/2], t2[N/2];
	int64_t v0p1[N/2], v[N];
	for(int i = 0; i < N/2; i++)
	{{
		v[i] = vect[i];
		v[i + N/2] = vect[i + N/2];
		v0p1[i] = v[i] + v[i + N/2];
	}}
	m1toep10(t0, v0p1, M1 + N/2);
	m1toep10(t1, v, M1_m0m1);
	m1toep10(t2, v + N/2, M1_m0m2);
	for(int i = 0; i < N/2; i++)
	{{
		rop[i] = t0[i] - t2[i];
		rop[i + N/2] = t0[i] - t1[i];
	}}
}}""")
		else:
			for i in range(NBCHUNKS):
				write(f"""
static inline void mtoeplitz_vm{i}(__int128* restrict rop, const int64_t* restrict vect)
{{
	toeplitz_vm(rop, vect, M[{i}]);
}}
""")
				write("")
			write(f"""
static inline void m1toeplitz_vm(int64_t* restrict rop, const __int128* restrict vect)
{{
	int64_t v[N];
	for(int i = 0; i < N; i++)
		v[i] = vect[i];
	m1toep20(rop, v, M1);
}}
""")
		write(f"""
void (* const multbym[{NBCHUNKS}]) (__int128* restrict rop, const int64_t* restrict vect) = {{ {str(["mtoeplitz_vm" + str(i) for i in range(NBCHUNKS)])[1:-1].replace("'","")} }};
""")
		
