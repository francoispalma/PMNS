#!/usr/bin/env python3
from sage.all import matrix,ZZ,xgcd,vector,ceil,prod,log

from dictofredcomppmns import pmnsdict, pmnsdict128

import os

(p,n,gamma,lam,rho,M) = \
pmnsdict[232]
pmnsdict128[18]

#rho = 2*(rho+1)
#w = 2*n-1
#print(n*w*round((rho-1)**2/rho) < 2**64)
#exit()

s = 1

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


NBCHUNKS = s

PHI = 2**64
ZZX = ZZ["X"]
E = ZZX(f"X^{n} - {lam}")
val, M1, soak = xgcd(ZZX(M), E)
M1 = (M1 * int(pow(int(val), -1, PHI) % PHI))
w = abs(lam)*(n-1) + 1

if __name__ == "__main__":
	with open("redparams.h", "w+") as fpointer:
		write = lambda X : fpointer.write(str(X) + "\n")
		write(f"#define N {n}")
		write(f"#define LAMBDA {lam}")
		if rho < 2**64:
			write(f"#define RHO {rho}")
			write(f"#define RHOhi {rho}")
		else:
			write(f"#define RHO (((__int128){rho>>64})|({rho%2**64}u))")
			write(f"#define RHOhi {rho>>64}")
			write(f"#define RHOlo {rho%2**64}u")
		DIVBYRHO = rho.bit_length()
		X = ZZX("X")
		tmpG = [list(ZZX(M) * X**i % E) for i in range(n)]
		for lig in tmpG:
			while len(lig) < n:
				lig += [0]
		G = matrix(ZZ, tmpG)
		while max([abs(round(elem)) for lig in (G.inverse() * PHI * 2**(DIVBYRHO)) for elem in lig]) > 2**63:
			DIVBYRHO -= 1
		write(f"#define DIVBYRHO {DIVBYRHO-2}")
		
		write(f"#define LENGTH_OF_P {ceil(log(p,2**64))}")
		invp = pow(p,-1,2**ceil(log(p,2)))
		write(f"\nstatic const uint64_t __P__[{ceil(log(p,2**64))}] = {{ {(str([(p//(2**(64*i)))%2**64 for i in range(ceil(log(p,2**64)))])[1:-1]).replace(',','u,')}u }};")
		write(f"\nstatic const uint64_t __invP__[{ceil(log(p,2**64))}] = {{ {(str([(invp//(2**(64*i)))%2**64 for i in range(ceil(log(p,2**64)))])[1:-1]).replace(',','u,')}u }};")
		
		Gprime = [[round(elem) for elem in lig] for lig in G.inverse() * PHI * 2**(DIVBYRHO-1)]
		val, M2, soak = xgcd(ZZX(M), E)
		M2 = (M2 * int(pow(int(val), -1, PHI**2) % PHI**2))
		M = list(ZZX(M[1:]) * lam) + M
		Mhi = [elem>>64 for elem in M]
		Mlo = [elem % 2**64 for elem in M]
		if rho > 2**64:
			M = Mhi
		write(f"\nstatic const int64_t M[{2*n-1}] = {{ {str(M)[1:-1]} }};")
		write(f"\nstatic const int64_t Mhi[{2*n-1}] = {{ {str(Mhi)[1:-1]} }};")
		write(f"\nstatic const uint64_t Mlo[{2*n-1}] = {{ {(str(Mlo)[1:-1]).replace(',','u,')}u }};")
		M1 = list(ZZX(list(M1)[1:]) * lam) + list(M1)
		# We switch M1 from M^-1 to -M^-1
		M1 = [(int(M1[i]) * -1) % PHI for i in range(len(M1))]
		# We then reduce M1's coefficients
		M1 = [int(M1[i]) - PHI if M1[i] >= (PHI >> 1) else int(M1[i]) for i in range(len(M1))]
		write(f"\nstatic const int64_t M1[{2*n-1}] = {{ {str(M1)[1:-1]} }};")
		M2 = list(ZZX(list(M2)[1:]) * lam) + list(M2)
		M2 = [(int(M2[i]) * -1) % PHI**2 for i in range(len(M2))]
		M2 = [int(M2[i]) - PHI**2 if M2[i] >= (PHI**2 >> 1) else int(M2[i]) for i in range(len(M2))]
		M2lo = [elem%PHI for elem in M2]
		M2hi = [elem//PHI for elem in M2]
		write(f"\nstatic const int64_t M2hi[{2*n-1}] = {{ {str(M2hi)[1:-1]} }};")
		write(f"\nstatic const uint64_t M2lo[{2*n-1}] = {{ {(str(M2lo)[1:-1]).replace(',','u,')}u }};")
		M1 = [Gprime[i][0] for i in range(n-1,0,-1)] + Gprime[0]
		write(f"\nstatic const int64_t Gprime[{2*n-1}] = {{ {str(M1)[1:-1]} }};")
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
