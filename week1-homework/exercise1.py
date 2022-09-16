#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import scipy

# set seed
np.random.seed(2729)

# FOR 5X COVERAGE

# start first plot
fig, ax = plt.subplots(figsize = (8, 4))

# create array of zeros of length 1M 
coverage = np.zeros(1000000, dtype=int)

# create list of 50,000 random integers (start positions for reads) between
# 0 and 999,900
randoms = np.random.randint(999901, size = 50000)

# for each random start site...
for random in randoms:

	# increment the count at the start position and the following 99 positions
	for i in range(100): 
		pos = random + i
		coverage[pos] += 1

# list of ints from 0 to 20 for x-axis
x = np.arange(21, dtype=int)

# set x-axis ticks
plt.xticks(x)

# plot a histogram of the frequency of coverages
plt.hist(coverage, bins=x, align='left')

# make pmf of poisson dist with lambda = 5
pmf = scipy.stats.poisson.pmf(x, 5)

# plot the frequency counts from pmf by multiplying by 1,000,000
plt.plot(x, pmf * 1000000, label = "Poisson distribution, lambda = 5")

# find out how many positions have no coverage
no_coverage = np.where(coverage == 0)
print(no_coverage[0].size)

# poisson expectation for no coverage
print(pmf[0] * 1000000)

plt.title("5x coverage")
ax.set_xlabel("Coverage")
ax.set_ylabel("Frequency")
ax.legend()
plt.savefig("ex_1.2.png")

# REPEAT FOR 15x

# start second plot
fig, ax = plt.subplots(figsize = (12, 4))

# create array of zeros of length 1M 
coverage = np.zeros(1000000, dtype=int)

# create list of 150,000 random integers (start positions for reads) between
# 0 and 999,900
randoms = np.random.randint(999901, size = 150000)

# for each random start site...
for random in randoms:

	# increment the count at the start position and the following 99 positions
	for i in range(100): 
		pos = random + i
		coverage[pos] += 1

# list of ints from 0 to 40 for x-axis
x = np.arange(41, dtype=int)

# set x-axis ticks
plt.xticks(x)

#plot a histogram of the frequency of coverages
plt.hist(coverage, bins=x, align='left')

# make pmf of poisson dist with lambda = 15
pmf = scipy.stats.poisson.pmf(x, 15)

# plot the frequency counts from pmf by multiplying by 1,000,000
plt.plot(x, pmf * 1000000, label = "Poisson distribution, lambda = 15")

# find out how many positions have no coverage
no_coverage = np.where(coverage == 0)
print(no_coverage[0].size)

# poisson expectation for no coverage
print(pmf[0] * 1000000)

plt.title("15x coverage")
ax.set_xlabel("Coverage")
ax.set_ylabel("Frequency")
ax.legend()
plt.savefig("ex_1.4.png")