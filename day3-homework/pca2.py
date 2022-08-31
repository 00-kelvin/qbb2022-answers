#!/usr/bin/env python

# import the matplotlib and numpy
import matplotlib.pyplot as plt
import numpy as np

# generate pc coordinates from text file
joined = np.genfromtxt("joined.txt", 
						dtype = None,
						encoding = None,
						names = ["indiv", "pop", "superpop", "sex", 
									"fam", "PC1", "PC2", "PC3"]) 
# generate a list of the different populations
pop_list = []
for i in joined: 
	if i["pop"] not in pop_list:
		pop_list.append(i["pop"])

# make my first plot
fig, ax = plt.subplots(figsize = (16, 12))	

# for each population...
for pop in pop_list: 

	# initialize variables, restart for each pop
	X = []
	Y = []

	# cycle through all individuals
	for i in range(len(joined)):

		# if my individual is in the current population, add their data
		# to the variable lists
		if joined[i]["pop"] == pop:
			X.append(joined[i]["PC1"])
			Y.append(joined[i]["PC2"])

	# make individual scatter plots for each pop, to go on same graph
	ax.scatter(X, Y, label = pop)
	ax.set_xlabel("PC1")
	ax.set_ylabel("PC2")
	ax.legend()

# format and save
plt.title("PCA of 1000 Genome Project Variants by Population")
plt.tight_layout()
plt.savefig("ex3_a.png")

# generate a list of the different superpopulations
superpop_list = []
for i in joined: 
	if i["superpop"] not in superpop_list:
		superpop_list.append(i["superpop"])

# start my second plot
fig, ax = plt.subplots()	

# for each superpopulation...
for superpop in superpop_list: 

	# initialize variables, restart for each superpop
	X = []
	Y = []

	# cycle through all individuals
	for i in range(len(joined)):

		# if my individual is in the current superpopulation, add their data
		# to the variable lists
		if joined[i]["superpop"] == superpop:
			X.append(joined[i]["PC1"])
			Y.append(joined[i]["PC2"])

	# make individual scatter plots for each superpop, to go on same graph
	ax.scatter(X, Y, label = superpop)
	ax.set_xlabel("PC1")
	ax.set_ylabel("PC2")
	ax.legend()

# make my second plot
plt.title("PCA of 1000 Genome Project Variants by Superpopulation")
plt.tight_layout()
plt.savefig("ex3_b.png")

# list of sexes
sex_list = ["male", "female"]

# start my third plot
fig, ax = plt.subplots()	

# for each sex...
for sex in sex_list: 

	# initialize variables, restart for each sex
	X = []
	Y = []

	# cycle through all individuals
	for i in range(len(joined)):

		# if my individual is of the current sex, add their data
		# to the variable lists
		if joined[i]["sex"] == sex:
			X.append(joined[i]["PC1"])
			Y.append(joined[i]["PC2"])

	# make individual scatter plots for each sex, to go on same graph
	ax.scatter(X, Y, label = sex)
	ax.set_xlabel("PC1")
	ax.set_ylabel("PC2")
	ax.legend()

# make my second plot
plt.title("PCA of 1000 Genome Project Variants by Sex")
plt.tight_layout()
plt.savefig("ex3_c.png")