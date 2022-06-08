from commonpmns import pmnsdicts
from math import log2, floor

def countdegs(pmnsdict):
	retdict = {}
	for key in pmnsdict:
		n = pmnsdict[key][1]
		if n in retdict:
			retdict[n] += 1
		else:
			retdict[n] = 1
	return retdict

def countlambs(pmnsdict):
	retdict = {}
	for key in pmnsdict:
		lamb = pmnsdict[key][3]
		if lamb in retdict:
			retdict[lamb] += 1
		else:
			retdict[lamb] = 1
	return retdict

if __name__ == "__main__":
	print("N:")
	for key in pmnsdicts:
		dico = pmnsdicts[key]
		print(key, countdegs(dico), len(dico))
	print("#" * 50)
	print("Lambda:")
	for key in pmnsdicts:
		dico = pmnsdicts[key]
		print(key, countlambs(dico), len(dico))
	print("#" * 50)
	print("Rho:")
	for key in pmnsdicts:
		dico = pmnsdicts[key]
		lambdict = countlambs(dico)
		ndict = countdegs(dico)
		phi = 128 if "128" in str(key) else 64
		minw = 1 + (min(ndict.keys()) - 1) * min(lambdict.keys())
		maxw = 1 + (max(ndict.keys()) - 1) * max(lambdict.keys())
		maxrho = floor(phi - 1 - log2(minw))
		minrho = floor(phi - 1 - log2(maxw))
		print(key, str(minrho) + " to " + str(maxrho))
