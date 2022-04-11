import matplotlib.pyplot as plt

if __name__ == "__main__":
	with open("results128", "r") as f:
		L = f.readlines()
	LL = [int(elem[:-1]) for elem in L]
	plt.hist(LL)
	plt.figure()
	with open("results", "r") as f:
		L = f.readlines()
	LLL = [int(elem[:-1]) for elem in L]
	plt.hist(LLL)
	plt.show()
