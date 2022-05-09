from commonpmns import pmnsdicts

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
	print("Lambda:")
	for key in pmnsdicts:
		dico = pmnsdicts[key]
		print(key, countlambs(dico), len(dico))
