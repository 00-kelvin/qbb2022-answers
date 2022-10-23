#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from vcfParser import parse_vcf
import seaborn as sb

df = pd.read_csv('CB1908_IC50_gwas_results.assoc.linear', 
	delim_whitespace=True)

# only care about genotype p-values
df = df[df['TEST'] == 'ADD']

# identify the snp with the lowest p-value
top_snp = df.loc[df['P'] == df['P'].min(), 'SNP'].item()

# read in the genotypes data
gt = pd.DataFrame(parse_vcf('gwas_data/genotypes.vcf'))

# make a list out of the genotype data from the top snp...
# the gt[2] == top_snp finds the row w/ the correct snp and makes a new
# pandas df that's just that row. then .values.tolist() converts the values
# in that df into a list object. but for some reason it makes a list of lists
# with only 1 object in the super list, which is the list of all the items in
# the row... so then I select that first super list [0] and within that, select
# the 10th through final items [9:] to grab just the genotype values

gt_list = gt[gt[2] == top_snp].values.tolist()[0][9:]

# extract a list of just the IC50 values
pt = pd.read_csv('gwas_data/CB1908_IC50.txt', delim_whitespace=True)
pt_list = pt['CB1908_IC50']

# make lists of homozygous ref, alt, and hets
ref_list = []
het_list = []
alt_list = []

for i, gen in enumerate(gt_list):
	if gen == '0/0':
		ref_list.append(pt_list[i])
	elif gen == '0/1' or gen == '1/0':
		het_list.append(pt_list[i])
	elif gen == '1/1':
		alt_list.append(pt_list[i])
	else:
		pass

# d = {'Homozygous ref': ref_list, 'Heterozygous': het_list, 
# 		'Homozygous alt': alt_list}

fig, ax = plt.subplots()

ax = sb.boxplot(data = [ref_list, het_list, alt_list])
ax.set_xticklabels(['Homozygous ref', 'Heterozygous', 'Homozygous alt'])
ax.set_ylabel('IC50')
plt.title('Effect size of SNP ' + top_snp + ' on CB1908 IC50')

plt.savefig('effect_size.png')

print(top_snp)