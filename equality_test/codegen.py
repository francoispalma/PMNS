#!/usr/bin/env python3
from sage.all import matrix,ZZ,xgcd,vector,prod,ceil,log,floor

from dictofequtestpmns import pmnsdict
from dictofbasismatrices import basisdict

import os

(p,n,gamma,lam,rho,M) = \
pmnsdict[5]

if "NSIZE" in os.environ:
	NSIZE = int(os.environ["NSIZE"])
	(p,n,gamma,lam,rho,M) = pmnsdict[NSIZE]

splitdico = {189:[3,3,3],80:[2,5],60:[2,2,3],48:[2,3],40:[2,2,2],36:[2,3],32:[2,2,2],27:[3,3],24:[2,3],20:[2,2],18:[3],16:[2,2],12:[3],9:[3],6:[2]}

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

PHI = 2**64
ZZX = ZZ["X"]
E = ZZX(f"X^{n} - {lam}")
val, M1, soak = xgcd(ZZX(M), E)
M1 = (M1 * int(pow(int(val), -1, PHI) % PHI))
w = abs(lam)*(n-1) + 1


X = ZZX("X")
tmpG = [list(ZZX(M) * X**i % E) for i in range(n)]
for lig in tmpG:
	while len(lig) < n:
		lig += [0]
G = matrix(ZZ, tmpG)
beta = 1
while beta**n < p:
	beta <<= 1
trans_m = max(n*(beta-1), w*(rho-1))
Ginvnormone = G.inverse().norm(1)
trans_u = ceil(trans_m*(rho-1)*Ginvnormone)
G1 = matrix(ZZ, G.inverse() % 2**64)
G = matrix(ZZ,basisdict[n])
invG = G.inverse()
G1 = invG % 2**64

#realG = matrix(ZZ,[[p if (k, j) == (0, 0) else -pow(gamma, k, p) if k != 0 and j == 0 else 1 if k == j else 0 for j in range(n)] for k in range(n)]).LLL()
#realG = basisdict[n]

def exact_conversion_to_pmns(inp, p, n, phi, G, G1, phin):
	G = [[elem for elem in lig] for lig in G]
	G1 = [[elem for elem in lig] for lig in G1]
	vecmat = lambda vec, mat: [sum([vec[j] * mat[j][i] for j in range(len(vec))]) for i in range(len(mat))]
	modphi = lambda vec: [(elem % phi) - phi * ((elem%phi) > phi/2) for elem in vec]
	plusvec = lambda vec1, vec2: [vec1[i] + vec2[i] for i in range(len(vec1))]
	slashphi = lambda vec: [elem//phi for elem in vec]
	Gmont_like = lambda vec: slashphi(plusvec(vec, vecmat(modphi(vecmat(vec, G1)), G)))
	A = [inp * phin % p] + [0] * (n-1)
	for i in range(n):
		A = Gmont_like(A)
	return A

def doubleredmult(inp1, inp2, n, lambda_, phi, G, G1):
	G = [[elem for elem in lig] for lig in G]
	G1 = [[elem for elem in lig] for lig in G1]
	vecmat = lambda vec, mat: [sum([vec[j] * mat[j][i] for j in range(len(vec))]) for i in range(len(mat))]
	modphi = lambda vec: [(elem % phi) - phi * ((elem%phi) > phi/2) for elem in vec]
	plusvec = lambda vec1, vec2: [vec1[i] + vec2[i] for i in range(len(vec1))]
	slashphi = lambda vec: [elem//phi for elem in vec]
	Gmont_like = lambda vec: slashphi(plusvec(vec, vecmat(modphi(vecmat(vec, G1)), G)))
	matr = [inp2[i] * lambda_ for i in range(1,n)] + [inp2[i] for i in range(n)]
	rop = [0] * n
	for i in range(n):
		for j in range(n):
			rop[i] += inp1[j] * matr[n - 1 - j + i]
	return Gmont_like(Gmont_like(rop))

if __name__ == "__main__":
	with open(f"equparams.h", "w+") as fpointer:
		write = lambda X : fpointer.write(str(X) + "\n")
		write(f"#define N {n}")
		write(f"#define LAMBDA {lam}")
		write(f"#define RHO {rho}")
		write(f"#define PHIOVER2G1 {ceil(2**64/(2*(rho+1)))}")
		Theta = rho.bit_length()
		while (1<<(Theta*n) < p):
			Theta += 1
		write("#define THETA " + str(Theta) + "\n")
		tmp = w*(rho-1)**2
		write(f"#define WRHOCARRE ((((__int128){tmp>>64})<<64) | ({tmp%2**64}u))")
		write(f'#define STEECONST {ceil(rho*invG.norm(1))}')
		tmp = ceil(tmp*invG.norm(1))
		write(f"#define DTEECONST ((((__int128){tmp>>64})<<64) | ({tmp%2**64}u))")
		
		if PHI >= 2*trans_u:
			write("\n#define TRANSLATIONPOSSIBLE")
			trans_T = vector([-trans_u for i in range(n)]) * G
			trans_Thi = [elem>>64 for elem in trans_T]
			trans_Tlo = [elem%(2**64) for elem in trans_T]
			write(f"static const int64_t trans_Thi[{n}] = {{ {str(trans_Thi)[1:-1]} }};")
			write(f"static const uint64_t trans_Tlo[{n}] = {{ {str(trans_Tlo)[1:-1].replace(',','u,')}u }};\n\n")
		
		write(f"#define LENGTH_OF_P {ceil(log(p,2**64))}")
		write(f"\nstatic const uint64_t __P__[{ceil(log(p,2**64))}] = {{ {(str([(p//(2**(64*i)))%2**64 for i in range(ceil(log(p,2**64)))])[1:-1]).replace(',','u,')}u }};")
		
		Gi = [pow(gamma,i,p) for i in range(1,n)]
		write(f"static const uint64_t __Gi__[{n-1}][{ceil(log(p,2**64))}] = {{")
		for j in range(n-1):
			write(f"{{ {(str([(Gi[j]//(2**(64*i)))%2**64 for i in range(ceil(log(p,2**64)))])[1:-1]).replace(',','u,')}u }},")
		write("};\n")
#		phin = pow(PHI, n, p)
#		Pi = [0] * n
#		Tau = PHI**2 % p
##		for i in range(n-1):
##			Pi[i] = exact_conversion_to_pmns(Tau, p, n, PHI, G, G1, phin)
##			Tau = Tau * 2**Theta % p
##		Pi[n-1] = exact_conversion_to_pmns(Tau, p, n, PHI, G, G1, phin)
#		Pi[0] = exact_conversion_to_pmns(Tau, p, n, PHI, G, G1, phin)
#		# We get the other powers of Theta from Pi[0]
#		thetaphisq = exact_conversion_to_pmns(2**Theta * Tau % p, p, n, PHI, G, G1, phin)
#		for i in range(1,n):
#			Pi[i] = doubleredmult(Pi[i-1], thetaphisq, n, lam, PHI, G, G1)
#		write("\nstatic const int64_t __Pi__[N][N] = {\n")
#		for i in range(len(Pi) - 1):
#			write("\t\t{" + str(Pi[i])[1:-1] + "},\n")
#		write("\t\t{" + str(Pi[-1])[1:-1] + "}\n\t};\n\n")
		
		Gprime = [[round(elem) for elem in lig] for lig in invG * PHI]
#		val, M2, soak = xgcd(ZZX(M), E)
#		M2 = (M2 * int(pow(int(val), -1, PHI**2) % PHI**2))
		M = list(ZZX(M[1:]) * lam) + M
		write(f"\nstatic const int64_t M[{2*n-1}] = {{ {str(M)[1:-1]} }};")
		M1 = list(ZZX(list(M1)[1:]) * lam) + list(M1)
		# We switch M1 from M^-1 to -M^-1
		M1 = [(int(M1[i]) * -1) % PHI for i in range(len(M1))]
		# We then reduce M1's coefficients
		M1 = [int(M1[i]) - PHI if M1[i] >= (PHI >> 1) else int(M1[i]) for i in range(len(M1))]
		write(f"\nstatic const int64_t M1[{2*n-1}] = {{ {str(M1)[1:-1]} }};")
#		M2 = list(ZZX(list(M2)[1:]) * lam) + list(M2)
#		M2 = [(int(M2[i]) * -1) % PHI**2 for i in range(len(M2))]
#		M2 = [int(M2[i]) - PHI**2 if M2[i] >= (PHI**2 >> 1) else int(M2[i]) for i in range(len(M2))]
#		M2lo = [elem%PHI for elem in M2]
#		M2hi = [elem//PHI for elem in M2]
#		write(f"\nstatic const int64_t M2hi[{2*n-1}] = {{ {str(M2hi)[1:-1]} }};")
#		write(f"\nstatic const uint64_t M2lo[{2*n-1}] = {{ {(str(M2lo)[1:-1]).replace(',','u,')}u }};")
		#M1 = [Gprime[i][0] for i in range(n-1,0,-1)] + Gprime[0]
		#write(f"\nstatic const int64_t Gprime[{2*n-1}] = {{ {str(M1)[1:-1]} }};")
		write("\nstatic const int64_t G[N][N] = {\n")
		G = [[elem for elem in lig] for lig in G]
		for i in range(len(G) - 1):
			write("\t\t{" + str(G[i])[1:-1] + "},\n")
		write("\t\t{" + str(G[-1])[1:-1] + "}\n\t};\n\n")
		
		write("\nstatic const int64_t G1[N][N] = {\n")
		G1 = -invG%PHI
		G1 = [[int(elem) - PHI if int(elem) >= ((PHI) >> 1) else int(elem) for elem in lig] for lig in G1]
		for i in range(len(G) - 1):
			write("\t\t{" + str(G1[i])[1:-1] + "},\n")
		write("\t\t{" + str(G1[-1])[1:-1] + "}\n\t};\n\n")
		
		G2 = -invG%PHI**2
		G2 = [[int(elem) - PHI**2 if int(elem) >= ((PHI**2) >> 1) else int(elem) for elem in lig] for lig in G2]
		G2hi = [[elem//PHI for elem in lig] for lig in G2]
		G2lo = [[elem%PHI for elem in lig] for lig in G2]
		
		write("\nstatic const int64_t G2hi[N][N] = {\n")
		for i in range(len(G) - 1):
			write("\t\t{" + str(G2hi[i])[1:-1] + "},\n")
		write("\t\t{" + str(G2hi[-1])[1:-1] + "}\n\t};\n\n")
		write("\nstatic const uint64_t G2lo[N][N] = {\n")
		for i in range(len(G) - 1):
			write("\t\t{" + (str(G2lo[i])).replace(',','u,')[1:-1] + "u},\n")
		write("\t\t{" + (str(G2lo[-1])).replace(',','u,')[1:-1] + "u}\n\t};\n\n")
		
		write("\nstatic const int64_t Gprime[N][N] = {\n")
		for i in range(len(Gprime) - 1):
			write("\t\t{" + str(Gprime[i])[1:-1] + "},\n")
		write("\t\t{" + str(Gprime[-1])[1:-1] + "}\n\t};\n\n")
		if TOEPSPLIT:
			write(f"""
static inline void mtoep{10-len(TOEPSPLIT)+1}(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	SCHOOLBOOK(N/{prod(TOEPSPLIT)})

static inline void mmmtoep{10-len(TOEPSPLIT)+1}(__int128* restrict rop, const __int128* restrict vect, const int64_t* restrict matr)
	SCHOOLBOOK(N/{prod(TOEPSPLIT)})

static inline void mmtoep{10-len(TOEPSPLIT)+1}(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	SCHOOLBOOK(N/{prod(TOEPSPLIT)})

static inline void m1toep{10-len(TOEPSPLIT)+1}(int64_t* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	M1SCHOOLBOOK(N/{prod(TOEPSPLIT)})
""")
			for i in range(len(TOEPSPLIT)-1, 0, -1):
				write(f"""
static inline void mtoep{10-i+1}(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	TOEP{TOEPSPLIT[i]}{TOEPSPLIT[i]}TOP(N/{prod(TOEPSPLIT[:i])}, mtoep{10-i})

static inline void mmmtoep{10-i+1}(__int128* restrict rop, const __int128* restrict vect, const int64_t* restrict matr)
	MMTOEP{TOEPSPLIT[i]}{TOEPSPLIT[i]}TOP(N/{prod(TOEPSPLIT[:i])}, mmmtoep{10-i})

static inline void mmtoep{10-i+1}(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	MTOEP{TOEPSPLIT[i]}{TOEPSPLIT[i]}TOP(N/{prod(TOEPSPLIT[:i])}, mmtoep{10-i}, mmmtoep{10-i})

static inline void m1toep{10-i+1}(int64_t* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	M1TOEP{TOEPSPLIT[i]}{TOEPSPLIT[i]}TOP(N/{prod(TOEPSPLIT[:i])}, m1toep{10-i})
""")
			
			write(f"""
static inline void toeplitz_vm(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	TOEP{TOEPSPLIT[0]}{TOEPSPLIT[0]}TOP(N, mtoep10)

static inline void mtoeptop(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	MTOEP{TOEPSPLIT[0]}{TOEPSPLIT[0]}TOP(N, mmtoep10, mmmtoep10)

static inline void m1toep20(int64_t* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	M1TOEP{TOEPSPLIT[0]}{TOEPSPLIT[0]}TOP(N, m1toep10)
""")
		else:
			write(f"""
static inline void toeplitz_vm(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	SCHOOLBOOK(N)

static inline void mtoeptop(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	SCHOOLBOOK(N)

static inline void m1toep20(int64_t* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	M1SCHOOLBOOK(N)
""")
		
		write(f"""
static inline void mtoeplitz_vm(__int128* restrict rop, const int64_t* restrict vect)
{{
	mtoeptop(rop, vect, M);
}}
""")
		write(f"""
static inline void m1toeplitz_vm(int64_t* restrict rop, const __int128* restrict vect)
{{
	int64_t v[N];
	for(int i = 0; i < N; i++)
		v[i] = vect[i];
	m1toep20(rop, v, M1);
}}
""")
		
		if len(TOEPSPLIT) >= 1:
			write(f"""
static inline void mtoep{128-len(TOEPSPLIT)+1}(__int128* restrict rop, const __int128* restrict vect, const __int128* restrict matr)
	SCHOOLBOOK(N/{prod(TOEPSPLIT[:])})
""")
			for i in range(len(TOEPSPLIT)-1, 0, -1):
				write(f"""
static inline void mtoep{128-i+1}(__int128* restrict rop, const __int128* restrict vect, const __int128* restrict matr)
	TOEP{TOEPSPLIT[i]}{TOEPSPLIT[i]}TOP128(N/{prod(TOEPSPLIT[:i])}, mtoep{128-i})
""")
			write(f"""
void toep128(__int128* restrict rop, const __int128* restrict vect, const __int128* restrict matr)
	TOEP{TOEPSPLIT[0]}{TOEPSPLIT[0]}TOP128(N, mtoep128)
""")
			write(f"""
static inline void bigtoep{128-len(TOEPSPLIT)+1}(__int128* restrict rop, const int64_t* restrict vect, const __int128* restrict matr)
	SCHOOLBOOK(N/{prod(TOEPSPLIT[:])})
""")
			for i in range(len(TOEPSPLIT)-1, 0, -1):
				write(f"""
static inline void bigtoep{128-i+1}(__int128* restrict rop, const int64_t* restrict vect, const __int128* restrict matr)
	bigmatTOEP{TOEPSPLIT[i]}{TOEPSPLIT[i]}TOP(N/{prod(TOEPSPLIT[:i])}, bigtoep{128-i})
""")
			write(f"""
void bigmattoep(__int128* restrict rop, const int64_t* restrict vect, const __int128* restrict matr)
	bigmatTOEP{TOEPSPLIT[0]}{TOEPSPLIT[0]}TOP(N, bigtoep128)
""")
		else:
			write(f"""
void toep128(__int128* restrict rop, const __int128* restrict vect, const __int128* restrict matr)
	SCHOOLBOOK(N)
""")
			write(f"""
void bigmattoep(__int128* restrict rop, const int64_t* restrict vect, const __int128* restrict matr)
	SCHOOLBOOK(N)
""")
