#!/usr/bin/env python3

from ast import literal_eval
from mpprecalcs import p, n, gamma, lam, rho, philog2

def horner_modulo(Poly: list, X: int, modulo: int) -> int:
	"""Polynomial evaluation function using Horner's algorithm, applying a
	modulo.
	
	PARAMETERS
	Poly: The polynomial in list form with the 0th degree in 0th index
	X: The value to evaluate the polynomial in
	modulo: The value to apply a modulo for during the process
	
	For our purposes we use this to convert out of the PMNS by evaluating in
	gamma modulo p.
	"""
	sum_ = 0
	for i in range(len(Poly) - 1, -1, -1):
		sum_ = int(sum_ * X % modulo)
		sum_ = sum_ + Poly[i]
	return int(sum_ % modulo)

if __name__ == "__main__":
	morethanrho = 0
	counter = 0
	visual = 0
	lastmorethanrho = 0
	lastinp = 0
	phi = 2**philog2
	print("Starting")
	with open("tmp", "r") as f:
		lines = f.readlines()
	while lines:
		print("\btested:" + str(visual) + "\terrors: " + str(counter), end="\r")
		visual += 1
		a = literal_eval(lines.pop(0))
		b = literal_eval(lines.pop(0))
		c = literal_eval(lines.pop(0))
		c_check = (horner_modulo(a, gamma, p) * horner_modulo(b, gamma, p) * pow(phi, -1, p)) % p
		if horner_modulo(c, gamma, p) != c_check:
			counter += 1
			if counter == 1:
				print(a)
				print(b)
				print(c)
		toaddto = 0
		for elem in c:
			if abs(elem) >= rho:
				toaddto = 1
				lastmorethanrho = abs(elem)
				lastinp = max(lastinp, max([abs(elem) for elem in a] + [abs(elem) for elem in b]))
		morethanrho += toaddto
	print("")
	print("more than rho: ", morethanrho)
	if morethanrho:
		print(lastmorethanrho, rho, lastinp)
