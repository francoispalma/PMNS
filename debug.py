from ops import amns_montg_mult
from proof import amns, phi, M, M1

if __name__ == "__main__":
	A = [3175695016735605, 20859843725, -123954529873808582, 541629668316248009, -29410447444707128]
	B = [1061418265038816869, 20374760404, -477028757217305698, 161008708292031432, -62502744134330068]
	C = amns_montg_mult(A, B, *amns, phi, M, M1)
	print(C)