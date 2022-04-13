import matplotlib.pyplot as plt
from statistics import mean

if __name__ == "__main__":
	with open("results128_old", "r") as f:
		L = f.readlines()
	LL = [int(elem[:-1]) for elem in L]
	plt.hist(LL)
	with open("results128_old2", "r") as f:
		L = f.readlines()
	LL = [int(elem[:-1]) for elem in L]
	plt.hist(LL)
	#plt.figure()
	with open("results_old", "r") as f:
		L = f.readlines()
	LLL = [int(elem[:-1]) for elem in L]
	plt.hist(LLL)
	print("Avg for pmns:", mean(LLL))
	print("Avg for pmns128:", mean(LL))
	print("Ratio:", mean(LL)/mean(LLL))
	plt.title("clocks, clocks128new, clocks128")
	plt.figure()
	with open("results", "r") as f:
		L = f.readlines()
	LL = [eval(elem[:-1]) for elem in L]
	mini, maxi, meani = ([elem[0] for elem in LL], [elem[1] for elem in LL],
		[elem[2] for elem in LL])
	plt.hist(mini)
	plt.hist(maxi)
	plt.hist(meani)
	plt.title("Min/Mean/Max")
	plt.figure()
	with open("results128", "r") as f:
		L = f.readlines()
	LL = [eval(elem[:-1]) for elem in L]
	mini, maxi, meani = ([elem[0] for elem in LL], [elem[1] for elem in LL],
		[elem[2] for elem in LL])
	plt.hist(mini)
	plt.hist(maxi)
	plt.hist(meani)
	plt.title("Min/Mean/Max 128")
	plt.show()
