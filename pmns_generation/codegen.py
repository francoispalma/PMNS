#!/bin/python3

import sys
from ast import literal_eval

if __name__ == "__main__":
	if len(sys.argv) >= 2:
		try:
			with open(sys.argv[1]) as inpf:
				try:
					p, n, gamma, E, rho, G, G1 = literal_eval(inpf.readlines()[0])
				except Exception as e:
					print("Invalid file format")
		except FileNotFoundError as e:
			print("File not found:", sys.argv[1])
		except Exception as e:
			print("Error:", e)
	else:
		print("oops")

