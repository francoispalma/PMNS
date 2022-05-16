import matplotlib.pyplot as plt
from statistics import mean

if __name__ == "__main__":
	with open("results/results1024", "r") as f:
		L = f.readlines()
	LL = [eval(elem[:-1]) for elem in L]
	mini, maxi, meani = ([elem[0] for elem in LL], [elem[1] for elem in LL],
		[elem[2] for elem in LL])
	min64 = mean(mini)
	mean64 = mean(meani)
	plt.hist(mini)
	plt.hist(maxi)
	plt.hist(meani)
	plt.title("New Min/Mean/Max")
	plt.figure()
	with open("results/results1024128", "r") as f:
		L = f.readlines()
	LL = [eval(elem[:-1]) for elem in L]
	mini, maxi, meani = ([elem[0] for elem in LL], [elem[1] for elem in LL],
		[elem[2] for elem in LL])
	min128 = mean(mini)
	mean128 = mean(meani)
	plt.hist(mini)
	plt.hist(maxi)
	plt.hist(meani)
	plt.title("Min/Mean/Max 128")
	print("Ratio1024:", mean128/mean64, "(" + str(min128/min64) + ")")
	plt.figure()
	with open("results/results2048", "r") as f:
		L = f.readlines()
	LL = [eval(elem[:-1]) for elem in L]
	mini, maxi, meani = ([elem[0] for elem in LL], [elem[1] for elem in LL],
		[elem[2] for elem in LL])
	min64 = mean(mini)
	mean64 = mean(meani)
	plt.hist(mini)
	plt.hist(maxi)
	plt.hist(meani)
	plt.title("New Min/Mean/Max 2048")
	plt.figure()
	with open("results/results2048128", "r") as f:
		L = f.readlines()
	LL = [eval(elem[:-1]) for elem in L]
	mini, maxi, meani = ([elem[0] for elem in LL], [elem[1] for elem in LL],
		[elem[2] for elem in LL])
	min128 = mean(mini)
	mean128 = mean(meani)
	plt.hist(mini)
	plt.hist(maxi)
	plt.hist(meani)
	plt.title("Min/Mean/Max 2048128")
	print("Ratio 2048:", mean128/mean64, "(" + str(min128/min64) + ")")
#	plt.figure()
#	with open("results/results4096", "r") as f:
#		L = f.readlines()
#	LL = [eval(elem[:-1]) for elem in L]
#	mini, maxi, meani = ([elem[0] for elem in LL], [elem[1] for elem in LL],
#		[elem[2] for elem in LL])
#	min64 = mean(mini)
#	mean64 = mean(meani)
#	plt.hist(mini)
#	plt.hist(maxi)
#	plt.hist(meani)
#	plt.title("Min/Mean/Max 4096")
#	plt.figure()
#	with open("results/results4096128", "r") as f:
#		L = f.readlines()
#	LL = [eval(elem[:-1]) for elem in L]
#	mini, maxi, meani = ([elem[0] for elem in LL], [elem[1] for elem in LL],
#		[elem[2] for elem in LL])
#	min128 = mean(mini)
#	mean128 = mean(meani)
#	plt.hist(mini)
#	plt.hist(maxi)
#	plt.hist(meani)
#	plt.title("Min/Mean/Max 4096128")
#	print("Ratio 4096:", mean128/mean64, "(" + str(min128/min64) + ")")
	plt.show()
