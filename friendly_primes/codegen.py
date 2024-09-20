#!/usr/bin/env python3
import sys, ast

# Code generation for PHI > 2^64 not supported for now, some adjustments will
# need to be made, such as using the __int128 type.
PHI = 2**64
SQRTPHI = 2**32

def convert_to_int_tabs(num):
	L = []
	string = hex(num)[2:]
	while string:
		L += [int(string[-1:-17:-1][::-1], 16)]
		string = string[:-16]
	return L

def Gmont_like(inp, phi, G, G1, sparseflag):
	n = len(G)
	modphi = lambda vec: [(elem % phi) - phi * ((elem%phi) > phi//2) for elem in vec]
	vectimesG = lambda vec: [vec[i] * G[i][i] + vec[(i-1)%n] * G[(i-1)%n][i] for i in range(n)]
	def vectimesG1(vec):
		if sparseflag:
			ret = modphi([vec[(i+1)%n] * G1[(i+1)%n][i] + vec[(i+2)%n] * G1[(i+2)%n][i] for i in range(n)])
		else:
			ret = [0] * n
			ret[-1] = sum([vec[i] * G1[i][-1] for i in range(n)])%phi
			ret[-2] = sum([vec[i] * G1[i][-2] for i in range(n)])%phi
			for i in range(2,n):
				ret[-i-1] = (ret[-i] * -G[0][0] - vec[-i])%phi
		return modphi(ret)
	plusvec = lambda vec1, vec2: [vec1[i] + vec2[i] for i in range(len(vec1))]
	slashphi = lambda vec: [elem//phi for elem in vec]
	return slashphi(plusvec(inp, vectimesG(vectimesG1(inp))))

def exact_conversion_to_pmns(inp, p, n, phi, G, G1, phin, sparseflag):
	A = [inp * phin % p] + [0] * (n-1)
	for i in range(n):
		A = Gmont_like(A, phi, G, G1, sparseflag)
	return A

def doubleredmult(inp1, inp2, n, alpha, lambda_, phi, G, G1, sparseflag):
	matr = [inp2[i] * lambda_ for i in range(1,n)] + [inp2[i] * alpha for i in range(n)]
	rop = [0] * n
	for i in range(n):
		for j in range(n):
			rop[i] += inp1[j] * matr[n - 1 - j + i]
	return Gmont_like(Gmont_like(rop, phi, G, G1, sparseflag), phi, G, G1, sparseflag)

if __name__ == "__main__":
	if len(sys.argv) >= 2:
		if len(sys.argv) >= 7:
			try:
				sys.argv[-2] = "'" + sys.argv[-2][:-1] + "',"
				pmnstr = "".join(sys.argv[1:])
			except Exception as e:
				print("Invalid parameters")
				print(e)
				exit()
		else:
			filename = sys.argv[1]
			linenumber = 0
			if len(sys.argv) >= 3:
				try:
					linenumber = int(sys.argv[2])
				except ValueError:
					print("Invalid line number:", sys.argv[2])
					exit()
			try:
				with open(filename, "r") as fpointer:
					try:
						pmnstr = fpointer.readlines()[linenumber]
					except IndexError:
						print("No line number", linenumber, "in file", filename)
						exit()
			except FileNotFoundError as e:
				print("No such file:", filename)
		try:
			(p,n,gamma,rho,E,delta) = ast.literal_eval(pmnstr)
		except Exception:
			print("Invalid parameters:")
			print(pmnstr)
			exit()
		with open("params.h","w+") as fpointer:
			write = lambda X: fpointer.write(X+"\n")
			write("#ifndef PMNS_PARAMS_H_INCLUDED\n#define PMNS_PARAMS_H_INCLUDED\n")
			write(f"#define RHO {rho}")
			Theta = rho.bit_length()
			while (2**(Theta*n) > p):
				Theta -= 1
			while (2**(Theta*n) < p):
				Theta += 1
			write(f"#define THETALOG2 {Theta}")
			Theta = 2**Theta
			write(f"#define N {n}")
			write(f"#define PMNS_NB_FREE_ADD {delta}")
			if E[0] == "X":
				alpha = 1
			else:
				sentinel = 0
				try:
					while E[sentinel] != '*':
						sentinel += 1
					alpha = int(E[:sentinel])
				except (IndexError, ValueError):
					print("Invalid format for E:", E)
					exit()
			write(f"#define ALPHA {alpha}")
			sentinel = 0
			try:
				while E[-1-sentinel] not in [" ","+","-"]:
					sentinel += 1
				lambda_ = int(E[-1-sentinel:])
				while E[-1-sentinel] == " ":
					sentinel += 1
				if E[-1-sentinel] == "+":
					lambda_ = -lambda_
			except (IndexError, ValueError):
				print("Invalid format for E:", E)
				exit()
			write(f"#define LAMBDA {lambda_}")
			assert (alpha*gamma**n - lambda_) % p == 0
			k = (alpha*gamma**n - lambda_) // p
			write(f"#define GAMMA {gamma}")
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
			write(f"#define ALEPH {aleph}")
			write(f"#define LAMED {lambda_//2**psi}")
			write(f"#define GIMEL {gamma//2**psadi}")
			if not((gamma//2**psadi) % SQRTPHI):
				sparseflag = True
				write("\n#define SPARSEG1")
				write(f"#define GAMMALAMM1 {((alpha*gamma//2**psi)*pow(z*p, -1, PHI)) % PHI}")
				write(f"#define GAMMALAMM2 {(gamma*pow(z*p, -1, PHI)) % PHI}")
				write(f"#define ONELAMM1 {pow(z*p, -1, PHI)}")
				G1 = [[((alpha*gamma//2**psi)*pow(z*p, -1, PHI)) % PHI if (i,j) == (0, n-2) else pow(z*p, -1, PHI) if (i,j) == (0,n-1) else (gamma*pow(z*p, -1, PHI)) % PHI if (i,j) == (1, n-1) else -gamma if j == i - 2 else -1 if j == i - 1 else 0 for j in range(n)] for i in range(n)]
			else:
				sparseflag = False
				fpointer.write(f"\nstatic const uint64_t lastcol[{n}] = {{")
				for i in range(n):
					fpointer.write(str((pow(z*p, -1, PHI)*gamma**i) % PHI) + "u, ")
				write("};\n")
				lastcol = [(pow(z*p, -1, PHI)*gamma**i) % PHI for i in range(n)]
				colforelast = [elem*(alpha*gamma)//2**psi % PHI for elem in lastcol]
				colforelast[-1] -= 1
				G1 = [colforelast, lastcol]
				for i in range(n-2):
					col = [elem for elem in G1[0][1:]]
					col += [(col[-1] * gamma) % PHI]
					G1 = [col] + G1
				G1 = [[G1[j][i] for j in range(n)] for i in range(n)]
			toep3 = 0
			toep2 = 0
			if n > 7:
				write("""#define SCHOOLBOOK(X) {\\
	for(int i = 0; i < X; i++)\\
	{\\
		rop[i] = 0;\\
		for(int j = 0; j < X; j++)\\
			rop[i] += (__int128) vect[j] * matr[X - 1 - j + i];\\
	}\\
}\n""")
				check = n
				if not(check % 2):
					write("""#define TOEP22TOP(X, F) {\\
	__int128 t0[X/2], t1[X/2], t2[X/2];\\
	int64_t v0p1[X/2], m0m1[X-1], m0m2[X-1];\\
	for(int i = 0; i < X/2; i++)\\
	{\\
		v0p1[i] = vect[i] + vect[i + X/2];\\
	}\\
	for(int i = 0; i < X-1; i++)\\
	{\\
		m0m1[i] = matr[i + X/2] - matr[i + X];\\
		m0m2[i] = matr[i + X/2] - matr[i];\\
	}\\
	F (t0, v0p1, matr + X/2);\\
	F (t1, vect, m0m1);\\
	F (t2, vect + X/2, m0m2);\\
	for(int i = 0; i < X/2; i++)\\
	{\\
		rop[i] = t0[i] - t2[i];\\
		rop[i + X/2] = t0[i] - t1[i];\\
	}\\
}\n""")
					while check > 7 and not(check % 2):
						check = check//2
						toep2 += 1
				if check > 7 and not(check%3):
					write("""#define TOEP33TOP(X, F) {\\
	__int128 t0[X/3], t1[X/3], t2[X/3], t3[X/3], t4[X/3], t5[X/3];\\
	int64_t m03, v1m2[X/3], v0m2[X/3], v0m1[X/3], m034[(2*X/3) - 1], m013[(2*X/3) - 1], m012[(2*X/3) - 1], m0[(2*X/3) - 1], m1[(2*X/3) - 1], m3[(2*X/3) - 1];\\
	for(int i = 0; i < X/3; i++)\\
	{\\
		v1m2[i] = vect[X/3 + i] - vect[2*X/3 + i];\\
		v0m1[i] = vect[i] - vect[X/3 + i];\\
		v0m2[i] = vect[i] - vect[2*X/3 + i];\\
	}\\
	for(int i = 0; i < (2*X/3) - 1; i++)\\
	{\\
		m0[i] = matr[i + 2*X/3];\\
		m1[i] = matr[i + X];\\
		m3[i] = matr[i + X/3];\\
		m03 = m0[i] + m3[i];\\
		m034[i] = m03 + matr[i];\\
		m013[i] = m03 + m1[i];\\
		m012[i] = m0[i] + m1[i] + matr[i + 4*X/3];\\
	}\\
	F (t0, vect + 2*X/3, m034);\\
	F (t1, vect + X/3, m013);\\
	F (t2, vect, m012);\\
	F (t3, v1m2, m3);\\
	F (t4, v0m2, m0);\\
	F (t5, v0m1, m1);\\
	for(int i = 0; i < X/3; i++)\\
	{\\
		rop[i] = t0[i] + t3[i] + t4[i];\\
		rop[i + X/3] = t1[i] - t3[i] + t5[i];\\
		rop[i + 2*X/3] = t2[i] - t4[i] - t5[i];\\
	}\\
}\n""")
					while check > 7 and not(check % 3):
						check = check//3
						toep3 += 1
				if not(toep3 or toep2):
					write(f"static inline void toeplitz_vm(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)\n\tSCHOOLBOOK({n})")
				else:
					write(f"static inline void toeplitz{check}(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)\n\tSCHOOLBOOK({check})\n")
					while toep3:
						check *= 3
						write(f"static inline void toeplitz{check}(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)\n\tTOEP33TOP({check}, toeplitz{check//3})\n")
						toep3 -= 1
					while toep2:
						check *= 2
						write(f"static inline void toeplitz{check}(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)\n\tTOEP22TOP({check}, toeplitz{check//2})\n")
						toep2 -= 1
					write(f"static inline void toeplitz_vm(__int128* restrict rop, const int64_t* restrict vect, const int64_t* restrict matr)\n{{\n\ttoeplitz{n}(rop, vect, matr);\n}}")
			fpointer.write("\nstatic const uint64_t __P__[] = {")
			for elem in ["0x" + hex(p)[2:][::-1][i*16:(i+1)*16][::-1] for i in range(n)]:
				if elem == "0x":
					elem = "0"
				fpointer.write(elem + ", ")
			write("};")
			write(f"#define __PINVERSE__ {-pow(p, -1, PHI) % PHI}u")
			
			G = [[-(alpha*gamma)//2**psi if i == j == (n-1) else -gamma if i == j else 1 if j == i + 1 else (lambda_//2**psi) if (i == n-1) and (j == 0) else 0 for j in range(n)] for i in range(n)]
			phin = pow(PHI, n, p)
			Pi = [0] * n
			Tau = pow(alpha,-1,p) * PHI**2 % p
			Pi[0] = exact_conversion_to_pmns(Tau, p, n, PHI, G, G1, phin, sparseflag)
			# We get the other powers of Theta from Pi[0]
			thetaphisq = exact_conversion_to_pmns(Theta * Tau % p, p, n, PHI, G, G1, phin, sparseflag)
			for i in range(1,n):
				Pi[i] = doubleredmult(Pi[i-1], thetaphisq, n, alpha, lambda_, PHI, G, G1, sparseflag)
			# Since the infinity-norm of Pi[0] and thetaphisq is bounded by half
			# the one-norm of G + 1 (see proposition 1), multiplying them give us
			# something bounded by w(||G||_1 + 1)². Hence applying Gmont_like
			# twice gives us something properly bounded. We multiply by alpha^-1
			# times phi² at the same time to account for both the reduction by E
			# and the montgomery reduction.
			write("\nstatic const int64_t __Pi__[N][N] = {")
			for i in range(len(Pi) - 1):
				write("\t\t{" + str(Pi[i])[1:-1] + "},")
			write("\t\t{" + str(Pi[-1])[1:-1] + "}\n\t};")
			
			write("\n#endif")
	else:
		print("Expects file name and a line number or tuple as parameter.")
