from statistics import mean

if __name__ == "__main__":
	bitsizes = ["1024", "2048", "4096", "8192"]
	meandict = {}
	for bitsize in bitsizes:
		for phisize in ["", "128"]:
			for appendix in ["", "pre", "hyb"]:
				if appendix == "hyb" and phisize == "":
					pass
				else:
					suffix = bitsize + phisize + appendix
					with open("results/results" + suffix, "r") as f:
						L = f.readlines()
					LL = [eval(elem[:-1]) for elem in L]
					meandict[suffix] = mean([elem[2] for elem in LL])
	for phisize in ["", "128"]:
		for appendix in ["pre", "", "hyb"]:
			if appendix == "hyb" and phisize == "":
					pass
			else:
				print(f"pmns{'64' if phisize == '' else phisize}{' déroulée' if appendix == 'pre' else ' hybride' if appendix == 'hyb' else ''}", end=" ")
				for bitsize in bitsizes:
					suffix = bitsize + phisize + appendix
					print(f"& {meandict[suffix]}", end=" ")
				print("\\\\\n\\hline")
