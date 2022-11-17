#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
import numpy.lib.recfunctions as rfn
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram

############################
##                        ##
##  0a: read in the data  ##
##                        ##
############################

input_arr = np.genfromtxt("dros_gene_expression.csv", 
 							delimiter=',', names=True, 
 							dtype=None, encoding='utf-8')

# sample names only 
col_names = list(input_arr.dtype.names)[1:]
row_names = list(input_arr['t_name'])

# subset to exclude transcript names
fpkm_values = input_arr[col_names]

##############################
##                          ##
##  0b: process input data  ##
##                          ##
##############################

fpkm_values_2d = rfn.structured_to_unstructured(fpkm_values, dtype=float)

# indices of transcripts with median expression greater than 0
nonzero_ind = np.where(np.median(fpkm_values_2d, axis=1) > 0)[0]

# filter data and transcript names
fpkm_fltrd = fpkm_values_2d[nonzero_ind]
row_names_fltrd = [row_names[i] for i in nonzero_ind]

# log transform
fpkm_fltrd_log = np.log2(fpkm_fltrd + 0.1)

#####################
##                 ##
##  1: clustering  ##
##                 ##
#####################

gene_linkage_mtrx = linkage(fpkm_fltrd_log, 'complete')
gene_leaves = leaves_list(gene_linkage_mtrx)

sample_linkage_mtrx = linkage(fpkm_fltrd_log.T, 'complete')
sample_leaves = leaves_list(sample_linkage_mtrx)

# reorder the rows according to the gene clusters
fpkm_row_sorted = fpkm_fltrd_log[gene_leaves]

# reorder the columns according to the sample clusters
fpkm_sorted = fpkm_row_sorted[:,sample_leaves]

# reorder the column names as well
col_names_sorted = [col_names[i] for i in sample_leaves]

# make my figure...
fig = plt.figure(figsize=(5,8))

# using gridspec to get the dendrogram and heatmap to line up
gs = fig.add_gridspec(nrows=2, ncols=2, height_ratios=[1,4], width_ratios=[5,1])

# making the dendrogram first 
ax1 = fig.add_subplot(gs[0,0])
sample_dendr = dendrogram(sample_linkage_mtrx, ax=ax1)

# no spines >:(
for i in ax1.spines.values():
	i.set_visible(False)

# no ax ticks either
ax1.tick_params(axis='both', which='major', labelsize=10, 
					labelbottom=False, bottom=False, 
					labelleft=False, left=False)
ax1.set_xticks([])

# making the heatmap 
ax2 = fig.add_subplot(gs[1,:])
sb.heatmap(fpkm_sorted, ax=ax2)

# no spines, some ticks
for i in ax2.spines.values():
	i.set_visible(False)

ax2.tick_params(axis='both', which='major', labelsize=10, 
					labelbottom=False, bottom=False, 
					labeltop=True, top=False, left=False)
ax2.set_xticks(np.arange(len(col_names_sorted))+.5)
ax2.set_xticklabels(col_names_sorted, rotation=90)
ax2.set_yticklabels([])

# plot!
plt.tight_layout()
plt.savefig('heatmap_and_dendr.png')


#######################################
##                                   ##
##  2. differential gene expression  ##
##                                   ##
#######################################


