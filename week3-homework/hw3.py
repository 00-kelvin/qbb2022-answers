#!/usr/bin/env python
##USAGE: hw3.py <annotated_vcf_file> 

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

##=================##
## BUILD THE PLOTS ##
##=================##

fig, ax = plt.subplots(nrows = 2, ncols = 2, figsize=(12, 10))
fig.suptitle("Variant genotype characteristics in 10 strains of S. cerevisiae",
				fontsize = 18, y = 0.95, fontweight = 'bold')

# read depths
ax[0,0].hist(dp_list_clean, bins = 35, color = 'red', ec = 'pink')
ax[0,0].set_yscale('log')
ax[0,0].set_xlabel('Read depth', fontweight = 'bold')
ax[0,0].set_ylabel('Number of variants')

# genotype qualities
ax[0,1].hist(gq_list_clean, bins = 16, color = 'orange', ec = 'yellow')
ax[0,1].set_xlabel('Genotype quality', fontweight = 'bold')
ax[0,1].set_ylabel('Number of variants')

# allele frequencies
ax[1,0].hist(af_list, bins = 10, color = 'green', ec = 'lightgreen')
ax[1,0].set_xlabel('Allele frequency', fontweight = 'bold')
ax[1,0].set_ylabel('Number of variants')

# predicted effects
plt.xticks(rotation=45, ha='right')

ax[1,1].bar(eff_names, eff_counts, ec = 'lightblue')
ax[1,1].set_yscale('log')
ax[1,1].set_xlabel('Predicted effects', fontweight = 'bold')
ax[1,1].set_ylabel('Number of variants')

plt.subplots_adjust(bottom = 0.25, wspace = 0.25, hspace = 0.25)
plt.savefig("plots.png", dpi=300)


