import matplotlib.pyplot as plt
from statistics import mean

if __name__ == "__main__":
#	with open("results128_old", "r") as f:
#		L = f.readlines()
#	LL = [int(elem[:-1]) for elem in L]
#	plt.hist(LL)
#	with open("results128_old2", "r") as f:
#		L = f.readlines()
#	LL = [int(elem[:-1]) for elem in L]
#	plt.hist(LL)
#	#plt.figure()
#	with open("results_old", "r") as f:
#		L = f.readlines()
#	LLL = [int(elem[:-1]) for elem in L]
#	plt.hist(LLL)
#	print("Avg for pmns:", mean(LLL))
#	print("Avg for pmns128:", mean(LL))
#	print("Ratio:", mean(LL)/mean(LLL))
#	plt.title("clocks, clocks128new, clocks128")
#	plt.figure()
	with open("results64", "r") as f:
		L = f.readlines()
	LL = [eval(elem[:-1]) for elem in L]
	mini, maxi, meani = ([elem[0] for elem in LL], [elem[1] for elem in LL],
		[elem[2] for elem in LL])
	min64 = mean(mini)
	mean64 = mean(meani)
	plt.hist(mini)
	plt.hist(maxi)
	plt.hist(meani)
	plt.title("Min/Mean/Max")
	plt.figure()
	with open("results1024", "r") as f:
		L = f.readlines()
	LL = [eval(elem[:-1]) for elem in L]
	mini, maxi, meani = ([elem[0] for elem in LL], [elem[1] for elem in LL],
		[elem[2] for elem in LL])
	nmin64 = mean(mini)
	nmean64 = mean(meani)
	plt.hist(mini)
	plt.hist(maxi)
	plt.hist(meani)
	plt.title("New Min/Mean/Max")
	plt.figure()
	with open("results128", "r") as f:
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
	print("Ratio:", mean128/mean64, "(" + str(min128/min64) + ")")
	print("New Ratio:", mean128/nmean64, "(" + str(min128/nmin64) + ")")
	plt.figure()
	with open("results204864", "r") as f:
		L = f.readlines()
	LL = [eval(elem[:-1]) for elem in L]
	mini, maxi, meani = ([elem[0] for elem in LL], [elem[1] for elem in LL],
		[elem[2] for elem in LL])
	min64 = mean(mini)
	mean64 = mean(meani)
	plt.hist(mini)
	plt.hist(maxi)
	plt.hist(meani)
	plt.title("Min/Mean/Max 2048")
	plt.figure()
	with open("results2048", "r") as f:
		L = f.readlines()
	LL = [eval(elem[:-1]) for elem in L]
	mini, maxi, meani = ([elem[0] for elem in LL], [elem[1] for elem in LL],
		[elem[2] for elem in LL])
	nmin64 = mean(mini)
	nmean64 = mean(meani)
	plt.hist(mini)
	plt.hist(maxi)
	plt.hist(meani)
	plt.title("New Min/Mean/Max 2048")
	plt.figure()
	with open("results2048128", "r") as f:
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
	print("New Ratio 2048:", mean128/nmean64, "(" + str(min128/nmin64) + ")")
#	plt.figure()
#	with open("results4096", "r") as f:
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
#	with open("results4096128", "r") as f:
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
