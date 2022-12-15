#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
# import matplotlib

# p is the starting allele frequency, N is the population size
def sim_wright_fisher(p, N):
	
	# begin list of allele frequencies at each generation
	allele_freq_list = [p]

	# as long as neither allele has fixed...
	while p != 0 and p != 1:

		# sample from a binomial distribution with n=2N, p=p
		allele_count = np.random.binomial(2*N, p)

		# calculate new allele frequency
		new_p = allele_count/(2*N)

		# add the new allele frequency to the list
		allele_freq_list.append(new_p)

		# set p to the new allele frequency for the next generation
		p = new_p

	return allele_freq_list

def allele_freq_vs_gen_number(p, N):

	# run the sim to get y values for the plot
	allele_freq = sim_wright_fisher(p, N)

	fig, ax = plt.subplots()

	# the ```plot``` function automatically sets x to the index of y! cool
	ax.plot(allele_freq)
	ax.set_xlabel("Generation number")
	ax.set_ylabel("Allele frequency")
	ax.set_title("Allele frequency vs. generation number")

	# text to add to the plot
	annotation = "Starting allele frequency: " + str(p) + "\n" + "Population size: " + str(N)

	ax.set_ylim((0,1))

	ax.annotate(annotation, (0, 0.1))

	plt.savefig("allele_freq_vs_gen_number.png", dpi=200)

def time_to_fixation_histogram(runs):
	
	times = []

	for i in range(runs):
		time_to_fix = len(sim_wright_fisher(0.5, 100))
		times.append(time_to_fix)

	fig, ax = plt.subplots()

	ax.hist(times, bins=100)

	ax.set_xlabel("Time to fixation")
	ax.set_ylabel("Number of simulations")
	ax.set_title("starting allele frequency = 0.5, population size = 100")
	plt.suptitle(
		"Time to fixation in " + str(runs) + " Wright Fisher simulations",
		weight='bold'
		)

	plt.savefig("time_to_fixation_histogram.png", dpi=200)

# takes a list of population sizes as argument
def fixation_time_vs_pop_size(N_values):

	times = []

	for N in N_values:
		time_to_fix = len(sim_wright_fisher(0.5, N))
		times.append(time_to_fix)

	fig, ax = plt.subplots()

	ax.plot(N_values, times, marker='o')

	ax.set_xlabel("Population size")
	ax.set_ylabel("Time to fixation")

	ax.set_xscale('log')

	# Made another version of the figure with log-scale y-axis too
	ax.set_yscale('log')


	# Ended up using even values so didn't need this code
	# ax.set_xticks(N_values)
	# ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
	# plt.xticks(rotation=45, ha='right')

	# Another way of annotating the N values that I didn't end up using
	# for N, t in zip(N_values, times):
	# 	ax.annotate(
	# 		"N = " + str(N),
	# 		(N,t)
	# 		)

	ax.set_title("Time to fixation vs. population size")

	plt.tight_layout()
	plt.savefig("fixation_time_vs_pop_size_log_time.png", dpi=200)

# takes a list of allele frequencies as argument
def fixation_time_vs_allele_freq(p_values, runs, N):

	data = []

	for p in p_values:
		times = []

		for i in range(runs):
			time_to_fix = len(sim_wright_fisher(p, N))
			times.append(time_to_fix)

		data.append(times)

	fig, ax = plt.subplots()

	# I can't figure out how to get the violin plots to not overlap. Oh well
	ax.violinplot(data, p_values)
	ax.set_xlabel("Starting allele frequency")
	ax.set_ylabel("Time to fixation")
	ax.set_title(
		"in " + str(runs) + " simulations with a population size of " + str(N)
		)
	plt.suptitle("Time to fixation vs. allele frequency", weight='bold')
	ax.set_xlim(-0.2, 1.2)

	ax.set_xticks(p_values)

	plt.savefig("fixation_time_vs_allele_freq.png", dpi=200)

def main():
	# allele_freq_vs_gen_number(0.3, 100)
	# time_to_fixation_histogram(1000)
	# fixation_time_vs_pop_size([100,1000,10000,100000,1000000,10000000])
	fixation_time_vs_allele_freq(np.arange(0, 1.1, 0.1), 100, 1000)


if __name__ == "__main__":
    main()