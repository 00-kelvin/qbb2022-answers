#!/usr/bin/env python

import sys
import matplotlib.pyplot as plt
import numpy as np
from vcfParser import parse_vcf 

fname = sys.argv[1]
vcf = parse_vcf(fname)

# initialize lists for read depth (dp), genotype quality (gq), 
# allele frequency (af) and predicted annotations (eff)
dp_list = []
gq_list = []
af_list = []
eff_list = []

##====================##
## POPULATE THE LISTS ##
##====================##

# loop through variants
for i in range(1, len(vcf)):
	variant = vcf[i]

	# add allele frequencies from each variant
	af_list.append(variant[7]["AF"])

	# store snpEff annotations as a list of lists
	annotations = variant[7]["ANN"].split(',')
	annotations = [item.split('|') for item in annotations]

	# add the effect annotations to the eff_list
	for item in annotations:
		eff_list.append(item[1])
		
	# look through the format column 
	for j, ID in enumerate(variant[8]):
		
		# find the index of the genotype quality ID
		if ID == 'GQ':

			# add the GQ for each sample in the variant to the GQ list
			for k in range(9, len(variant)):
				gq_list.append(variant[k][j])

		# repeat for read depth		
		elif ID == 'DP':
			for k in range(9, len(variant)):
				dp_list.append(variant[k][j])
		else:
			pass

##==================##
## FORMAT THE LISTS ##
##==================##

# split any multiple effects so they'll all be counted
eff_list = [item.split('&') for item in eff_list]

# flatten the list
eff_list = [eff for item in eff_list for eff in item]

# remove empty strings
eff_list = [eff for eff in eff_list if eff]

# create lists of the unique effect names and counts of each
eff_names, eff_counts = np.unique(np.array(eff_list), return_counts = True)

# change genotype qualities and read depts to floats and remove '.' entries
gq_list_clean = []
for gq in gq_list:
	try:
		gq = float(gq)
		gq_list_clean.append(gq)
	except:
		pass

# repeat for read depths
dp_list_clean = []
for dp in dp_list:
	try:
		dp = int(dp)
		dp_list_clean.append(dp)
	except:
		pass

fig, ax = plt.subplots(nrows = 2, ncols = 2)

# dps, dpcounts = np.unique(np.array(dp_list_clean), return_counts = True)
# print(dps,dpcounts)

ax[0,0].hist(dp_list_clean, bins = 20)
ax[0,0].set_yscale('log')
ax[0,1].hist(gq_list_clean, bins = 20)
ax[1,0].hist(af_list, bins = 10)
ax[1,1].barh(eff_names, eff_counts, log = True)

plt.show()

