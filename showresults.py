import sys

from statistics import mean

if __name__ == "__main__":
	bitsizes = ["256", "512", "1024", "2048", "4096", "8192"]
	meandict = {}
	if "--latex" in sys.argv:
		print(" &", end = "")
		print(*bitsizes, sep=" &", end=" &\\\\\n\hline\n")
	else:
		print("Median clock cycles for size of p per implementation")
		print("\t\t|", end="")
		print(*bitsizes, sep="\t|")
	for bitsize in bitsizes:
		for phisize in ["", "128"]:
			for appendix in ["", "pre", "hyb"]:
				if appendix == "hyb" and phisize == "":
					pass
				else:
					suffix = bitsize + phisize + appendix
					try:
						with open("results/results" + suffix, "r") as f:
							L = f.readlines()
						LL = [eval(elem[:-1]) for elem in L]
						meandict[suffix] = mean([elem[2] for elem in LL])
					except FileNotFoundError:
						pass
	for phisize in ["", "128"]:
		for appendix in ["", "pre", "hyb"]:
			if appendix == "hyb" and phisize == "":
					pass
			else:
				print(f"pmns{'64' if phisize == '' else phisize}{' unroll' if appendix == 'pre' else ' hybrid' if appendix == 'hyb' else ''}", end=(("\t") if appendix == "" else ""))
				
				for bitsize in bitsizes:
					suffix = bitsize + phisize + appendix
					try:
						if "--latex" in sys.argv:
							print(" &", end="")
						else:
							print("\t|", end="")
						print(f"{round(meandict[suffix])}", end="")
					except KeyError:
						pass
				if "--latex" in sys.argv:
					print("\\\\\n\\hline")
				else:
					print()
