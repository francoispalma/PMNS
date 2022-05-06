from generated1024pmns import pmnsdict as pmnsdict1024
from generated2048pmns import pmnsdict as pmnsdict2048
from generated4096pmns import pmnsdict as pmnsdict4096
from generated1024pmns128 import pmns128dict as pmns128dict1024
from generated2048pmns128 import pmns128dict as pmns128dict2048
from generated4096pmns128 import pmns128dict as pmns128dict4096

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
	pmnsdicts = {1024:pmnsdict1024, 2048:pmnsdict2048, 4096:pmnsdict4096,
		1024128: pmns128dict1024, 2048128:pmns128dict2048, 4096128: pmns128dict4096}
	for key in pmnsdicts:
		dico = pmnsdicts[key]
		print(key, countdegs(dico), len(dico))
	for key in pmnsdicts:
		dico = pmnsdicts[key]
		print(key, countlambs(dico), len(dico))
