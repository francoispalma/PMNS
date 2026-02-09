#!/usr/bin/env python3
from sage.all import ZZ, matrix, ceil, vector, log, prod
from dictofhequtestpmns import pmnsdict

import os

(p,n,gamma,rho,E,delta) = \
pmnsdict[5]

if "NSIZE" in os.environ:
	NSIZE = int(os.environ["NSIZE"])
	(p,n,gamma,rho,E,delta) = pmnsdict[NSIZE]

E = ZZ["X"](E)
lam = E[0]
lambda_ = lam
alpha = E[E.degree()]
w = (n-1)*abs(lam)+1
k = 1
PHI = 2**64

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

if __name__ == "__main__":
	with open(f"hequparams.h", "w+") as fpointer:
		write = lambda X : fpointer.write(str(X) + "\n")
		write(f"#define N {n}")
		write(f"#define LAMBDA {lam}")
		write(f"#define RHO {rho}")
		psi = 0
		z = k
		while not(z % 2):
			psi += 1
			z//=2
		psadi = psi
		aleph = alpha
		while psadi > 0 and not(aleph % 2):
			aleph //= 2
			psadi -= 1
		write(f"#define PARAM_K {k}")
		write(f"#define ALPHA {alpha}")
		write(f"#define GAMMA {gamma}")
		write(f"#define ALEPH {aleph}")
		write(f"#define LAMED {lambda_//2**psi}")
		write(f"#define GIMEL {gamma//2**psadi}")
		tmp = w*(rho-1)**2
		write(f"#define WRHOCARRE ((((__int128){tmp>>64})<<64) | ({tmp%2**64}u))")
		Theta = rho.bit_length()
		while (2**(Theta*n) > p):
			Theta -= 1
		while (2**(Theta*n) < p):
			Theta += 1
		write(f"#define THETALOG2 {Theta}")
		Theta = 2**Theta
		fpointer.write(f"\nstatic const int64_t lastcol[{n}] = {{")
		for i in range(n):
			tmp =int( (pow(z*p, -1, PHI)*gamma**i) % PHI)
			tmp = (tmp - PHI) if tmp > PHI//2 else tmp
			fpointer.write(str(tmp) + "u, ")
		write("};\n")
		beta = 1
		while beta**n < p:
			beta <<= 1
		trans_m = max(n*(beta-1), w*(rho-1))
#		G = [[-(alpha*gamma)//2**psi if i == j == (n-1) else -gamma if i == j else 1 if j == i + 1 else (lambda_//2**psi) if (i == n-1) and (j == 0) else 0 for j in range(n)] for i in range(n)]
		G = matrix(ZZ, [[p if i == j == 0 else 1 if i == j else -gamma if i == j + 1 else 0 for j in range(n)] for i in range(n)]).LLL()
		Ginvnormone = matrix(ZZ, G).inverse().norm(1)
		trans_u = ceil(trans_m*(rho-1)*Ginvnormone)
		if PHI >= 2*trans_u:
			write("\n#define TRANSLATIONPOSSIBLE")
		trans_T = vector([-trans_u for i in range(n)]) * matrix(ZZ, G)
		trans_Thi = [elem>>64 for elem in trans_T]
		trans_Tlo = [elem%(2**64) for elem in trans_T]
		write(f"static const int64_t trans_Thi[{n}] = {{ {str(trans_Thi)[1:-1]} }};")
		write(f"static const uint64_t trans_Tlo[{n}] = {{ {str(trans_Tlo)[1:-1].replace(',','u,')}u }};\n\n")
		
		write(f"#define LENGTH_OF_P {ceil(log(p,2**64))}")
		write(f"\nstatic const uint64_t __P__[{ceil(log(p,2**64))}] = {{ {(str([(p//(2**(64*i)))%2**64 for i in range(ceil(log(p,2**64)))])[1:-1]).replace(',','u,')}u }};")
		
		
		invG = G.inverse()
		G = [[elem for elem in lig] for lig in G]
		write("\nstatic const int64_t G[N][N] = {\n")
		for i in range(len(G) - 1):
			write("\t\t{" + str(G[i])[1:-1] + "},\n")
		write("\t\t{" + str(G[-1])[1:-1] + "}\n\t};\n\n")
		
		write("\nstatic const int64_t G1[N][N] = {\n")
		G1 = -invG%PHI
		G1 = [[int(elem) - PHI if int(elem) >= ((PHI) >> 1) else int(elem) for elem in lig] for lig in G1]
		for i in range(len(G) - 1):
			write("\t\t{" + str(G1[i])[1:-1] + "},\n")
		write("\t\t{" + str(G1[-1])[1:-1] + "}\n\t};\n\n")
		
		Gi = [pow(gamma,i,p) for i in range(1,n)]
		write(f"static const uint64_t __Gi__[{n-1}][{ceil(log(p,2**64))}] = {{")
		for j in range(n-1):
			write(f"{{ {(str([(Gi[j]//(2**(64*i)))%2**64 for i in range(ceil(log(p,2**64)))])[1:-1]).replace(',','u,')}u }},")
		write("};\n\n")
		
		if TOEPSPLIT:
			write(f"""
static inline void mtoep{10-len(TOEPSPLIT)+1}(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	SCHOOLBOOK(N/{prod(TOEPSPLIT)})
""")
			for i in range(len(TOEPSPLIT)-1, 0, -1):
				write(f"""
static inline void mtoep{10-i+1}(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	TOEP{TOEPSPLIT[i]}{TOEPSPLIT[i]}TOP(N/{prod(TOEPSPLIT[:i])}, mtoep{10-i})
""")
			
			write(f"""
static inline void toeplitz_vm(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	TOEP{TOEPSPLIT[0]}{TOEPSPLIT[0]}TOP(N, mtoep10)
""")
		else:
			write(f"""
static inline void toeplitz_vm(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)
	SCHOOLBOOK(N)
""")
